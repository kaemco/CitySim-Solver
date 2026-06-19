#ifndef NOISE_H
#define NOISE_H

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "tinyxml.h"
#include "util.h"

class District;
class Climate;

using namespace std;

/**
 * Steady-state outdoor noise propagation estimate for CitySim.
 *
 * Models sound propagation from point emittors to point receptors through
 * the 3D scene geometry, accounting for:
 *   1. Geometric spherical spreading  (ISO 9613-2, Adiv)
 *   2. Climate-derived atmospheric absorption (ISO 9613-1)
 *   3. Single-band ISO ground effect   (factor G 0=hard to 1=porous)
 *   4. ISO 9613-2:2024 top-edge barrier attenuation
 *   5. Approximate first-order wall reflections (image-source method)
 *
 * The propagation matrix M[e][r] - contribution of emittor e to receptor r
 * in dB - is assembled from all direct, diffracted and reflected paths.
 * Multiple contributions are combined via energy superposition.
 *
 * This is an ISO 9613-inspired screening model, not a complete implementation
 * of ISO 9613-2:2024. It uses a single representative frequency. When only
 * A-weighted source power is known, 500 Hz can be used for an attenuation
 * estimate. Octave-band source data and octave-band aggregation are not yet
 * implemented.
 *
 * Example XML:
 * <Noise id="RoadNoise" key="road" frequency="500"
 *        groundFactor="0.5"
 *        reflectionOrder="1"
 *        backgroundLevel="30">
 *     <PointEmittor id="E1" x="10" y="0" z="0.5" heightAboveGround="0.5"
 *                    powerLevel="90"/>
 *     <PointReceptor id="R1" x="100" y="20" z="1.5" heightAboveGround="1.5"/>
 * </Noise>
 *
 * XML attributes (all optional, defaults shown):
 *   frequency          [Hz]   centre frequency for all calculations (500)
 *   atmosphericAlpha   [dB/km] optional override; otherwise derived from climate
 *   groundFactor       [-]    G=0 hard (road), G=1 soft (grass) (0.0)
 *   wallAbsorptionCoeff[-]    fallback fraction absorbed per reflection when
 *                             wall composition has no acoustic coefficient (0.1)
 *   wallReflectivityCoeff[-]  fallback reflected-energy fraction. When omitted,
 *                             reflectivity is 1 - wallAbsorptionCoeff (0.9)
 *   reflectionOrder    [-]    0=direct only, 1=first-order reflections (1)
 *   backgroundLevel    [dB]   incoherent background; disabled below -300 (-300)
 *   directivityCorrection [dB] default source correction Dc (0)
 * PointEmittor and PointReceptor may define heightAboveGround [m]. When omitted,
 * the z coordinate is also used as the acoustic height above local ground.
 * PointEmittor may override directivityCorrection, e.g. +3 dB for half-space.
 * Wall compositions may define AcousticAbsorptionCoeff/AbsorptionCoeff and/or
 * AcousticReflectivityCoeff/ReflectivityCoeff. Defaults are 0.1 absorption
 * and 0.9 reflectivity.
 */
class Noise {

public:

    struct PointEmittor {
        string id, key;
        map<string,string> attributes;
        double x=0., y=0., z=0.;
        double heightAboveGround=0.;
        double directivityCorrection=0.;
        double powerLevel=0.;  // Lw, dB re 1 pW
    };

    struct PointReceptor {
        string id, key;
        map<string,string> attributes;
        double x=0., y=0., z=0.;
        double heightAboveGround=0.;
        vector<double>         level;           // total SPL (dB) per compute() call
        vector<vector<double>> levelByEmittor;  // SPL contribution per emittor
    };

private:

    // Internal polygon used for ray intersection and reflection.
    // Uses plain doubles and std::array to avoid coupling noise.h to GENPoint.
    struct SceneSurface {
        vector<array<double,3>> vertices;
        array<double,3>         normal;          // unit outward normal
        double                  absorptionCoeff; // [0,1]
        double                  reflectivityCoeff; // [0,1]
        bool                    isWall;
    };

    string id, key;
    map<string,string>  attributes;
    vector<PointEmittor>  emittors;
    vector<PointReceptor> receptors;
    vector<double> atmosphericAlpha;
    ostream logStream;

    // ── XML helpers (same pattern as Pollutant) ──────────────────────────────
    static map<string,string> readAttributes(TiXmlElement*);
    static bool   hasAttribute(const map<string,string>&, const string& name);
    static string getAttribute(const map<string,string>&, const string& name, const string& def);
    static double getAttribute(const map<string,string>&, const string& name, double def);
    static bool   getAttribute(const map<string,string>&, const string& name, bool def);
    static void   writeAttributes(ofstream&, const map<string,string>&);

    template<class T>
    static void eraseVec(vector<T>& v, unsigned int keep) {
        if (v.size() > keep)
            v.erase(v.begin(), v.end() - min(keep, (unsigned int)v.size()));
    }

    // ── Scene geometry ───────────────────────────────────────────────────────
    vector<SceneSurface> collectSurfaces(District*, double defaultWallAbsorption,
                                         double defaultWallReflectivity) const;

    // Moller-Trumbore ray-triangle intersection. t is the hit distance along D.
    static bool rayTriangle(const array<double,3>& O, const array<double,3>& D, double tMax,
                             const array<double,3>& V0, const array<double,3>& V1,
                             const array<double,3>& V2, double& t);

    // Fan-triangulate a convex polygon and test all triangles.
    static bool raySurface(const array<double,3>& O, const array<double,3>& D, double tMax,
                            const SceneSurface&, double& t);

    // Returns the index of the first surface that blocks segment P to Q, or -1.
    // The surface at index `skip` is excluded (set to -1 to include all).
    int firstObstruction(const array<double,3>& P, const array<double,3>& Q,
                          const vector<SceneSurface>&, int skip, double& tHit) const;

    // Point on the highest edge of a blocker nearest the propagation line.
    static array<double,3> diffractionEdge(const SceneSurface& blocker,
                                            const array<double,3>& source,
                                            const array<double,3>& receptor);

    // Approximate screened path over every intersected top edge. Returns an
    // empty vector if no geometrically clear top-edge path can be constructed.
    vector<array<double,3>> diffractionPath(const array<double,3>& source,
                                            const array<double,3>& receptor,
                                            const vector<SceneSurface>& surfaces) const;

    // Mirror image of point P through the plane of surface S.
    static array<double,3> imageSource(const array<double,3>& P, const SceneSurface& S);

    // ── ISO 9613-2 attenuation terms (all positive dB values = losses) ───────
    static double Adiv (double d);                            // geometric spreading
    static double Aatm (double d, double alpha);              // atmospheric absorption
    static double atmosphericAbsorptionCoefficient(double freq, double temperatureCelsius,
                                                    double relativeHumidity, double pressurePa);
    static double Agr  (double dp, double hs, double hr, double G, double freq); // ground effect
    static double Abar (double directDistance, const vector<array<double,3>>& path,
                        double lambda);                         // top-edge barrier attenuation
    static double Arefl(double reflectivityCoeff);             // reflection loss

    // ── Core computation ─────────────────────────────────────────────────────
    double computeSPL(const PointEmittor&, const PointReceptor&,
                       const vector<SceneSurface>&,
                       double freq, double alpha, double G,
                       int reflOrder, double speedOfSound) const;

public:

    Noise(TiXmlHandle hdl, ostream* pLogStr=&cout);

    // Collect scene geometry from District and compute the propagation matrix.
    // Appends one result entry per receptor (like Pollutant::compute).
    void compute(District* pDistrict, Climate* pClimate, unsigned int day, unsigned int hour);

    void   clear();
    void   eraseResults(unsigned int keepValue);
    void   writeXML(ofstream&, string tab="");
    void   writeHeaderText(fstream&);
    void   writeResultsText(fstream&, unsigned int index);
    size_t memoryUsage() const;

    string getId()         const { return id; }
    string getKey()        const { return key; }
    size_t getnEmittors()  const { return emittors.size(); }
    size_t getnReceptors() const { return receptors.size(); }
    size_t getnResults()   const { return receptors.empty() ? 0 : receptors[0].level.size(); }
};

#endif // NOISE_H
