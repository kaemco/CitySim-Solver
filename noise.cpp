#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

#include "noise.h"
#include "climate.h"
#include "district.h"  // Building, Zone, Wall, Roof, Surface iterators

// ── Inline 3-vector arithmetic ────────────────────────────────────────────────
// Using plain array<double,3> avoids pulling GENPoint templates into noise.h.

static inline double dot3(const array<double,3>& a, const array<double,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline array<double,3> cross3(const array<double,3>& a, const array<double,3>& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
static inline array<double,3> sub3(const array<double,3>& a, const array<double,3>& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
static inline array<double,3> add3(const array<double,3>& a, const array<double,3>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}
static inline array<double,3> scale3(const array<double,3>& a, double s) {
    return {a[0]*s, a[1]*s, a[2]*s};
}
static inline double norm3(const array<double,3>& a) {
    return sqrt(dot3(a, a));
}
static inline array<double,3> normalize3(const array<double,3>& a) {
    double n = norm3(a);
    return (n > 1e-12) ? scale3(a, 1.0/n) : array<double,3>{0,0,0};
}

// ── Constructor ───────────────────────────────────────────────────────────────

Noise::Noise(TiXmlHandle hdl, ostream* pLogStr) : logStream(cout.rdbuf()) {

    associate(pLogStr, logStream);

    TiXmlElement* elem = hdl.ToElement();
    if (!elem) throw string("Noise: invalid XML element.");

    attributes = readAttributes(elem);
    id  = getAttribute(attributes, "id",  string("Noise"));
    key = getAttribute(attributes, "key", getAttribute(attributes, "id", string("")));
    if (!hasAttribute(attributes, "id")) attributes["id"] = id;

    double freq           = getAttribute(attributes, "frequency",           500.0);
    double groundFactor   = getAttribute(attributes, "groundFactor",        0.0);
    double wallAbsorption = getAttribute(attributes, "wallAbsorptionCoeff", Composite::defaultAcousticAbsorptionCoeff());
    double wallReflectivity = getAttribute(attributes, "wallReflectivityCoeff", 1.0 - wallAbsorption);
    if (!hasAttribute(attributes, "wallAbsorptionCoeff") && hasAttribute(attributes, "wallReflectivityCoeff")) wallAbsorption = 1.0 - wallReflectivity;
    double reflectionOrder = getAttribute(attributes, "reflectionOrder",    1.0);
    if (!(freq > 0.)) throw string("Noise: frequency must be positive.");
    if (hasAttribute(attributes, "atmosphericAlpha") && getAttribute(attributes, "atmosphericAlpha", 0.) < 0.) throw string("Noise: atmosphericAlpha must be non-negative.");
    if (groundFactor < 0. || groundFactor > 1.) throw string("Noise: groundFactor must be between 0 and 1.");
    if (wallAbsorption < 0. || wallAbsorption > 1.) throw string("Noise: wallAbsorptionCoeff must be between 0 and 1.");
    if (wallReflectivity < 0. || wallReflectivity > 1.) throw string("Noise: wallReflectivityCoeff must be between 0 and 1.");
    if (wallAbsorption + wallReflectivity > 1.00001) throw string("Noise: wallAbsorptionCoeff plus wallReflectivityCoeff must not exceed 1.");
    if (reflectionOrder != 0. && reflectionOrder != 1.) throw string("Noise: reflectionOrder must be 0 or 1.");

    unsigned int i = 0;
    TiXmlElement* child = hdl.ChildElement("PointEmittor", i).ToElement();
    while (child) {
        PointEmittor em;
        em.attributes  = readAttributes(child);
        em.id          = getAttribute(em.attributes, "id",          toString(i));
        em.key         = getAttribute(em.attributes, "key",         string(""));
        em.x           = getAttribute(em.attributes, "x",           0.);
        em.y           = getAttribute(em.attributes, "y",           0.);
        em.z           = getAttribute(em.attributes, "z",           0.);
        em.heightAboveGround = getAttribute(em.attributes, "heightAboveGround", em.z);
        em.directivityCorrection = getAttribute(em.attributes, "directivityCorrection", getAttribute(attributes, "directivityCorrection", 0.));
        em.powerLevel  = getAttribute(em.attributes, "powerLevel",  0.);
        if (em.heightAboveGround < 0.) throw string("Noise: PointEmittor heightAboveGround must be non-negative.");
        emittors.push_back(em);
        child = hdl.ChildElement("PointEmittor", ++i).ToElement();
    }
    if (emittors.empty()) throw string("Noise: no PointEmittor defined.");

    i = 0;
    child = hdl.ChildElement("PointReceptor", i).ToElement();
    while (child) {
        PointReceptor re;
        re.attributes = readAttributes(child);
        re.id         = getAttribute(re.attributes, "id",  toString(i));
        re.key        = getAttribute(re.attributes, "key", string(""));
        re.x          = getAttribute(re.attributes, "x",  0.);
        re.y          = getAttribute(re.attributes, "y",  0.);
        re.z          = getAttribute(re.attributes, "z",  0.);
        re.heightAboveGround = getAttribute(re.attributes, "heightAboveGround", re.z);
        if (re.heightAboveGround < 0.) throw string("Noise: PointReceptor heightAboveGround must be non-negative.");
        receptors.push_back(re);
        child = hdl.ChildElement("PointReceptor", ++i).ToElement();
    }
    if (receptors.empty()) throw string("Noise: no PointReceptor defined.");

    logStream << "Noise model \"" << id << "\": "
              << emittors.size()  << " emittor(s), "
              << receptors.size() << " receptor(s) loaded." << endl << flush;
}

// ── Scene surface collection ──────────────────────────────────────────────────

vector<Noise::SceneSurface> Noise::collectSurfaces(District* pDistrict,
                                                    double defaultWallAbsorption,
                                                    double defaultWallReflectivity) const {
    vector<SceneSurface> out;

    auto addSurf = [&](Surface* s, bool isWall) {
        if (!s || s->vertexCount() < 3) return;
        SceneSurface ss;
        ss.absorptionCoeff = defaultWallAbsorption;
        ss.reflectivityCoeff = defaultWallReflectivity;
        if (s->getComposite()) {
            ss.absorptionCoeff = s->getComposite()->getAcousticAbsorptionCoeff();
            ss.reflectivityCoeff = s->getComposite()->getAcousticReflectivityCoeff();
        }
        // Convert GENPoint vertices to array<double,3>
        vector<GENPoint>* verts = s->getVertices();
        ss.vertices.reserve(verts->size());
        for (const GENPoint& gp : *verts) 
            ss.vertices.push_back({(double)gp[0], (double)gp[1], (double)gp[2]});
        // Unit normal
        GENPoint n = s->normal();
        double nlen = norm3({(double)n[0], (double)n[1], (double)n[2]});
        ss.normal = (nlen > 1e-9)
            ? array<double,3>{n[0]/nlen, n[1]/nlen, n[2]/nlen}
            : array<double,3>{0., 0., 1.};
        ss.isWall = isWall && fabs(ss.normal[2]) < 1e-3;
        out.push_back(move(ss));
    };

    // Simulated buildings: walls and roofs from all zones
    for (size_t b = 0; b < pDistrict->getnBuildings(); ++b) {
        Building* bldg = pDistrict->getBuilding(b);
        for (size_t z = 0; z < bldg->getnZones(); ++z) {
            Zone* zone = bldg->getZone(z);
            for (size_t w = 0; w < zone->getnWalls(); ++w)
                addSurf(zone->getWall(w), true);
            for (size_t r = 0; r < zone->getnRoofs(); ++r)
                addSurf(zone->getRoof(r), false);
        }
    }

    // Non-simulated buildings stored as district-level obstructing surfaces
    for (unsigned int s = 0; s < pDistrict->getnSurfaces(); ++s)
        addSurf(pDistrict->getSurface(s), true);

    return out;
}

// ── Möller-Trumbore ray-triangle intersection ─────────────────────────────────

bool Noise::rayTriangle(const array<double,3>& O, const array<double,3>& D, double tMax,
                         const array<double,3>& V0, const array<double,3>& V1,
                         const array<double,3>& V2, double& t) {
    static const double EPS = 1e-9;
    array<double,3> E1 = sub3(V1, V0);
    array<double,3> E2 = sub3(V2, V0);
    array<double,3> h  = cross3(D, E2);
    double a = dot3(E1, h);
    if (fabs(a) < EPS) return false;
    double f = 1.0 / a;
    array<double,3> s = sub3(O, V0);
    double u = f * dot3(s, h);
    if (u < 0. || u > 1.) return false;
    array<double,3> q = cross3(s, E1);
    double v = f * dot3(D, q);
    if (v < 0. || u + v > 1.) return false;
    t = f * dot3(E2, q);
    return (t > 1e-6) && (t < tMax - 1e-6);
}

bool Noise::raySurface(const array<double,3>& O, const array<double,3>& D, double tMax,
                        const SceneSurface& surf, double& t) {
    const auto& v = surf.vertices;
    const size_t n = v.size();
    if (n < 3) return false;
    // Fan-triangulate from v[0]: test (v0,v1,v2), (v0,v2,v3), …
    for (size_t i = 1; i + 1 < n; ++i) {
        double ti;
        if (rayTriangle(O, D, tMax, v[0], v[i], v[i+1], ti)) {
            t = ti;
            return true;
        }
    }
    return false;
}

int Noise::firstObstruction(const array<double,3>& P, const array<double,3>& Q,
                              const vector<SceneSurface>& surfaces,
                              int skip, double& tHit) const {
    array<double,3> D   = sub3(Q, P);
    double           len = norm3(D);
    if (len < 1e-9) return -1;
    array<double,3> Dn = scale3(D, 1.0 / len);

    int    hit  = -1;
    double tMin = len;

    for (int i = 0; i < (int)surfaces.size(); ++i) {
        if (i == skip) continue;
        double t;
        if (raySurface(P, Dn, tMin, surfaces[i], t)) {
            tMin = t;
            hit  = i;
        }
    }
    tHit = tMin;
    return hit;
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

// Returns the point on the highest edge that is closest to the horizontal
// source-receptor line. The small vertical clearance prevents the resulting
// path from numerically intersecting the edge used to construct it.
array<double,3> Noise::diffractionEdge(const SceneSurface& blocker,
                                      const array<double,3>& source,
                                      const array<double,3>& receptor) {
    if (blocker.vertices.empty()) return {0,0,0};
    double maxZ = blocker.vertices[0][2];
    for (const auto& v : blocker.vertices) maxZ = max(maxZ, v[2]);

    const double dx = receptor[0] - source[0];
    const double dy = receptor[1] - source[1];
    const double lineLength = hypot(dx, dy);
    array<double,3> closest = blocker.vertices[0];
    double bestDistance = 1e300;
    bool foundEdge = false;

    for (size_t i = 0; i < blocker.vertices.size(); ++i) {
        const auto& a = blocker.vertices[i];
        const auto& b = blocker.vertices[(i + 1) % blocker.vertices.size()];
        if (fabs(a[2] - maxZ) > 1e-3 || fabs(b[2] - maxZ) > 1e-3) continue;

        double t = 0.;
        if (lineLength > 1e-9) {
            double sideA = dx * (a[1] - source[1]) - dy * (a[0] - source[0]);
            double sideB = dx * (b[1] - source[1]) - dy * (b[0] - source[0]);
            if (sideA * sideB <= 0. && fabs(sideA - sideB) > 1e-12)
                t = sideA / (sideA - sideB);
            else
                t = fabs(sideA) <= fabs(sideB) ? 0. : 1.;
        }
        array<double,3> candidate = add3(a, scale3(sub3(b, a), t));
        double distance = lineLength > 1e-9
            ? fabs(dx * (candidate[1] - source[1]) - dy * (candidate[0] - source[0])) / lineLength
            : hypot(candidate[0] - source[0], candidate[1] - source[1]);
        if (distance < bestDistance) {
            closest = candidate;
            bestDistance = distance;
            foundEdge = true;
        }
    }

    if (!foundEdge) {
        for (const auto& v : blocker.vertices) {
            if (fabs(v[2] - maxZ) > 1e-3) continue;
            double distance = hypot(v[0] - source[0], v[1] - source[1]);
            if (distance < bestDistance) {
                closest = v;
                bestDistance = distance;
            }
        }
    }
    closest[2] += 1e-3;
    return closest;
}

vector<array<double,3>> Noise::diffractionPath(const array<double,3>& source,
                                                const array<double,3>& receptor,
                                                const vector<SceneSurface>& surfaces) const {
    struct EdgePoint {
        double position;
        array<double,3> point;
    };

    const array<double,3> direct = sub3(receptor, source);
    const double directLength2 = dot3(direct, direct);
    if (directLength2 < 1e-12) return {};

    vector<EdgePoint> edges;
    auto addEdge = [&](int surfaceIndex) {
        array<double,3> point = diffractionEdge(surfaces[surfaceIndex], source, receptor);
        double position = dot3(sub3(point, source), direct) / directLength2;
        if (position <= 1e-6 || position >= 1. - 1e-6) return false;
        for (const auto& edge : edges)
            if (norm3(sub3(edge.point, point)) < 1e-2) return false;
        edges.push_back({position, point});
        return true;
    };

    auto addObstructions = [&](const array<double,3>& from, const array<double,3>& to) {
        array<double,3> segment = sub3(to, from);
        double length = norm3(segment);
        if (length < 1e-9) return false;
        array<double,3> direction = scale3(segment, 1. / length);
        bool added = false;
        for (int i = 0; i < (int)surfaces.size(); ++i) {
            double t;
            if (raySurface(from, direction, length, surfaces[i], t))
                added = addEdge(i) || added;
        }
        return added;
    };

    addObstructions(source, receptor);
    for (size_t attempt = 0; attempt <= surfaces.size(); ++attempt) {
        sort(edges.begin(), edges.end(), [](const EdgePoint& a, const EdgePoint& b) {
            return a.position < b.position;
        });

        vector<EdgePoint> mergedEdges;
        for (const auto& edge : edges) {
            if (!mergedEdges.empty() && fabs(edge.position - mergedEdges.back().position) < 1e-9) {
                if (edge.point[2] > mergedEdges.back().point[2])
                    mergedEdges.back() = edge;
            } else {
                mergedEdges.push_back(edge);
            }
        }

        vector<EdgePoint> rubberBand;
        rubberBand.push_back({0., source});
        for (const auto& edge : mergedEdges) rubberBand.push_back(edge);
        rubberBand.push_back({1., receptor});

        vector<EdgePoint> visibleEdges;
        for (const auto& edge : rubberBand) {
            while (visibleEdges.size() >= 2) {
                const EdgePoint& a = visibleEdges[visibleEdges.size() - 2];
                const EdgePoint& b = visibleEdges.back();
                double ratio = (b.position - a.position) / (edge.position - a.position);
                double lineZ = a.point[2] + ratio * (edge.point[2] - a.point[2]);
                if (b.point[2] > lineZ + 1e-6) break;
                visibleEdges.pop_back();
            }
            visibleEdges.push_back(edge);
        }

        vector<array<double,3>> path;
        path.reserve(visibleEdges.size());
        for (const auto& edge : visibleEdges) path.push_back(edge.point);

        bool clear = true;
        bool added = false;
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            array<double,3> segment = sub3(path[i + 1], path[i]);
            double length = norm3(segment);
            if (length < 1e-9) continue;
            array<double,3> direction = scale3(segment, 1. / length);
            for (int j = 0; j < (int)surfaces.size(); ++j) {
                double t;
                if (raySurface(path[i], direction, length, surfaces[j], t)) {
                    clear = false;
                    added = addEdge(j) || added;
                }
            }
        }
        if (clear) return path;
        if (!added) break;
    }
    return {};
}

// Mirror image of P through the plane defined by surface S (unit normal + V0).
array<double,3> Noise::imageSource(const array<double,3>& P, const SceneSurface& S) {
    const auto& V0 = S.vertices[0];
    double signedDist = dot3(S.normal, sub3(P, V0));
    return sub3(P, scale3(S.normal, 2.0 * signedDist));
}

// ── ISO 9613-2 attenuation terms ─────────────────────────────────────────────

// Spherical divergence for an omnidirectional point source in a free field.
// Source directivity belongs in Dc, not in Adiv.
double Noise::Adiv(double d) {
    return 20.0 * log10(max(d, 0.1)) + 11.0;
}

// Atmospheric absorption, alpha in dB/km.
double Noise::Aatm(double d, double alpha) {
    return alpha * d / 1000.0;
}

// ISO 9613-1 atmospheric absorption coefficient in dB/km.
double Noise::atmosphericAbsorptionCoefficient(double freq, double temperatureCelsius,
                                               double relativeHumidity, double pressurePa) {
    const double referencePressure = 101325.;
    const double referenceTemperature = 293.15;
    const double triplePointTemperature = 273.16;
    const double temperature = temperatureCelsius + 273.15;
    const double pressureRatio = pressurePa / referencePressure;
    const double saturationPressureRatio = pow(10., -6.8346 * pow(triplePointTemperature / temperature, 1.261) + 4.6151);
    const double molarHumidity = 100. * max(0., min(relativeHumidity, 1.))
                               * saturationPressureRatio / pressureRatio;
    const double relaxationOxygen = (24. + 4.04e4 * molarHumidity * (0.02 + molarHumidity)
                                  / (0.391 + molarHumidity)) / pressureRatio;
    const double relaxationNitrogen = pow(temperature / referenceTemperature, -0.5)
                                    * (9. + 280. * molarHumidity
                                    * exp(-4.17 * (pow(temperature / referenceTemperature, -1. / 3.) - 1.)))
                                    / pressureRatio;
    const double classical = 1.84e-11 / pressureRatio * sqrt(temperature / referenceTemperature);
    const double molecular = pow(temperature / referenceTemperature, -2.5)
                           * (0.01275 * exp(-2239.1 / temperature)
                           / (relaxationOxygen + freq * freq / relaxationOxygen)
                           + 0.1068 * exp(-3352. / temperature)
                           / (relaxationNitrogen + freq * freq / relaxationNitrogen));
    return 1000. * 8.686 * freq * freq * (classical + molecular);
}

// ISO 9613-2:2024 general ground attenuation for one selected frequency band.
// A single global G is used for the source, middle and receptor regions.
double Noise::Agr(double dp, double hs, double hr, double G, double freq) {
    if (dp < 1e-9) return 0.;
    hs = max(hs, 0.);
    hr = max(hr, 0.);

    static const double bands[] = {63., 125., 250., 500., 1000., 2000., 4000., 8000.};
    int band = 0;
    for (int i = 1; i < 8; ++i)
        if (fabs(log(freq / bands[i])) < fabs(log(freq / bands[band]))) band = i;

    double q = dp > 30. * (hs + hr) ? 1. - 30. * (hs + hr) / dp : 0.;
    auto contribution = [&](double height) {
        if (band == 0) return -1.5;
        if (band == 1) {
            double a = 1.5 + 3. * exp(-0.12 * pow(height - 5., 2.)) * (1. - exp(-dp / 50.))
                           + 5.7 * exp(-0.09 * height * height) * (1. - exp(-2.8e-6 * dp * dp));
            return -1.5 + G * a;
        }
        if (band == 2) {
            double b = 1.5 + 8.6 * exp(-0.09 * height * height) * (1. - exp(-dp / 50.));
            return -1.5 + G * b;
        }
        if (band == 3) {
            double c = 1.5 + 14. * exp(-0.46 * height * height) * (1. - exp(-dp / 50.));
            return -1.5 + G * c;
        }
        if (band == 4) {
            double d = 1.5 + 5. * exp(-0.9 * height * height) * (1. - exp(-dp / 50.));
            return -1.5 + G * d;
        }
        return -1.5 * (1. - G);
    };

    double Aprime = contribution(hs) + contribution(hr)
                  + (band == 0 ? -3. * q : -3. * q * (1. - G));
    double Kgeo = (dp * dp + pow(hs - hr, 2.)) / (dp * dp + pow(hs + hr, 2.));
    return -10. * log10(max(1e-12, 1. + (pow(10., -Aprime / 10.) - 1.) * Kgeo));
}

// ISO 9613-2:2024 top-edge barrier attenuation Dz. Ground reflections are
// included through C2=20, so Agr is not applied separately to screened rays.
double Noise::Abar(double directDistance, const vector<array<double,3>>& path,
                   double lambda) {
    if (path.size() < 3 || directDistance <= 0. || lambda <= 0.) return 0.;

    double diffractedDistance = 0.;
    for (size_t i = 0; i + 1 < path.size(); ++i)
        diffractedDistance += norm3(sub3(path[i + 1], path[i]));
    double z = diffractedDistance - directDistance;

    double e = 0.;
    for (size_t i = 1; i + 2 < path.size(); ++i)
        e += norm3(sub3(path[i + 1], path[i]));

    const double C2 = 20.;
    double C3 = 1.;
    if (path.size() > 3 && e > 1e-9) {
        double ratio = 5. * lambda / e;
        C3 = (1. + ratio * ratio) / (1. / 3. + ratio * ratio);
    }
    double zMin = -2. * lambda / (C2 * C3);
    if (z <= zMin) return 0.;

    double dSS = norm3(sub3(path[1], path[0]));
    double dSR = norm3(sub3(path[path.size() - 1], path[path.size() - 2]));
    double Kmet = exp(-sqrt((max(dSS, dSR) + e) * min(dSS, dSR) * directDistance
                         / (2. * (z - zMin))) / 2000.);
    return 10. * log10(1. + (2. + (C2 / lambda) * C3 * z) * Kmet);
}

// Reflection loss from a wall with reflected-energy coefficient rho.
// A_refl = -10 log10(rho)  [dB, always non-negative for rho in (0,1]].
double Noise::Arefl(double reflectivityCoeff) {
    double rho = max(1e-12, min(reflectivityCoeff, 1.0));
    return -10.0 * log10(rho);
}

// ── Core SPL computation ──────────────────────────────────────────────────────

double Noise::computeSPL(const PointEmittor& em, const PointReceptor& re,
                          const vector<SceneSurface>& surfaces,
                          double freq, double alpha, double G,
                          int reflOrder, double speedOfSound) const {

    const array<double,3> O = {em.x, em.y, em.z};
    const array<double,3> R = {re.x,  re.y,  re.z};
    const double d = norm3(sub3(R, O));
    const double dp = hypot(re.x - em.x, re.y - em.y);
    if (d < 0.1) return em.powerLevel + em.directivityCorrection;  // coincident points

    const double lambda = speedOfSound / max(freq, 1.0);
    double energySum = 0.0;

    // ── 1. Direct / diffracted path ──────────────────────────────────────────
    {
        double tHit;
        int blockerIdx = firstObstruction(O, R, surfaces, -1, tHit);

        double Lp;
        if (blockerIdx < 0) {
            // Unobstructed: full direct path with ground effect
            double loss = Adiv(d)
                        + Aatm(d, alpha)
                        + Agr(dp, em.heightAboveGround, re.heightAboveGround, G, freq);
            Lp = em.powerLevel + em.directivityCorrection - loss;
        } else {
            // Obstructed: construct a clear path over all intersected top edges.
            vector<array<double,3>> path = diffractionPath(O, R, surfaces);
            if (path.empty()) {
                Lp = -300.;
            } else {
                double pathLength = 0.;
                for (size_t i = 0; i + 1 < path.size(); ++i)
                    pathLength += norm3(sub3(path[i + 1], path[i]));
                // Use direct distance for spreading (ISO 9613-2 section 7.3.2).
                double loss = Adiv(d)
                            + Aatm(pathLength, alpha)
                            + Abar(d, path, lambda)
                            + min(Agr(dp, em.heightAboveGround, re.heightAboveGround, G, freq), 0.);
                Lp = em.powerLevel + em.directivityCorrection - loss;
            }
        }
        energySum += pow(10.0, Lp / 10.0);
    }

    // ── 2. First-order reflections (image-source method) ─────────────────────
    if (reflOrder >= 1) {
        for (int i = 0; i < (int)surfaces.size(); ++i) {
            const SceneSurface& S = surfaces[i];
            if (!S.isWall) continue;                     // only vertical faces reflect
            if (S.reflectivityCoeff <= 0.0) continue;
            double sourceSide = dot3(S.normal, sub3(O, S.vertices[0]));
            double receptorSide = dot3(S.normal, sub3(R, S.vertices[0]));
            if (sourceSide * receptorSide <= 0.) continue;

            // Mirror the emittor through surface i
            array<double,3> Eimg = imageSource(O, S);

            // Ray from image source to receptor
            array<double,3> dirFull = sub3(R, Eimg);
            double dPath = norm3(dirFull);
            if (dPath < 0.1) continue;
            array<double,3> Dn = scale3(dirFull, 1.0 / dPath);

            // Does this ray hit the reflecting surface in [eps, dPath]?
            double tRefl;
            if (!raySurface(Eimg, Dn, dPath, S, tRefl)) continue;

            // Reflection point on surface i
            array<double,3> Prefl = add3(Eimg, scale3(Dn, tRefl));

            // Check O → Prefl is unobstructed (skip surface i itself)
            double t1, t2;
            if (firstObstruction(O, Prefl, surfaces, i, t1) >= 0) continue;

            // Check Prefl → R is unobstructed (skip surface i)
            if (firstObstruction(Prefl, R, surfaces, i, t2) >= 0) continue;

            // Valid first-order reflection
            double loss = Adiv(dPath)
                        + Aatm(dPath, alpha)
                        + Arefl(S.reflectivityCoeff);
            double Lp = em.powerLevel + em.directivityCorrection - loss;
            energySum += pow(10.0, Lp / 10.0);
        }
    }

    return (energySum > 0.0) ? 10.0 * log10(energySum) : -300.0;
}

// ── Public interface ──────────────────────────────────────────────────────────

void Noise::compute(District* pDistrict, Climate* pClimate, unsigned int day, unsigned int hour) {

    if (!pDistrict) throw string("Noise: no district available.");

    double wallAbsorption      = getAttribute(attributes, "wallAbsorptionCoeff",
                                             Composite::defaultAcousticAbsorptionCoeff());
    double wallReflectivity    = getAttribute(attributes, "wallReflectivityCoeff",
                                             1.0 - wallAbsorption);
    if (!hasAttribute(attributes, "wallAbsorptionCoeff")
            && hasAttribute(attributes, "wallReflectivityCoeff"))
        wallAbsorption = 1.0 - wallReflectivity;
    double freq                = getAttribute(attributes, "frequency",           500.0);
    if (!pClimate && !hasAttribute(attributes, "atmosphericAlpha"))
        throw string("Noise: climate data or atmosphericAlpha override required.");
    double temperatureCelsius  = pClimate ? pClimate->getToutCelsius(day, hour) : 20.;
    double alpha               = hasAttribute(attributes, "atmosphericAlpha")
        ? getAttribute(attributes, "atmosphericAlpha", 0.)
        : atmosphericAbsorptionCoefficient(freq, temperatureCelsius,
                                           pClimate->getRelativeHumidity(day, hour),
                                           pClimate->getPatm(day, hour));
    double G                   = getAttribute(attributes, "groundFactor",        0.0);
    int    reflOrder           = (int)getAttribute(attributes, "reflectionOrder", 1.0);
    double bg                  = getAttribute(attributes, "backgroundLevel",    -300.0);
    double speedOfSound          = 331.3 * sqrt(1. + temperatureCelsius / 273.15);

    const vector<SceneSurface> surfaces = collectSurfaces(pDistrict, wallAbsorption,
                                                          wallReflectivity);
    atmosphericAlpha.push_back(alpha);

    for (size_t ri = 0; ri < receptors.size(); ++ri) {
        PointReceptor& rec = receptors[ri];
        vector<double> byEmittor(emittors.size(), 0.);

        // Start with incoherent background (already as linear energy)
        double energyTotal = (bg > -300.) ? pow(10.0, bg / 10.0) : 0.0;

        for (size_t ei = 0; ei < emittors.size(); ++ei) {
            double lp = computeSPL(emittors[ei], rec, surfaces,
                                   freq, alpha, G, reflOrder, speedOfSound);
            byEmittor[ei] = lp;
            energyTotal  += pow(10.0, lp / 10.0);
        }

        rec.level.push_back(energyTotal > 0. ? 10.0 * log10(energyTotal) : -300.);
        rec.levelByEmittor.push_back(byEmittor);
    }
}

// ── Lifecycle ─────────────────────────────────────────────────────────────────

void Noise::clear() {
    atmosphericAlpha.clear();
    for (auto& r : receptors) {
        r.level.clear();
        r.levelByEmittor.clear();
    }
}

void Noise::eraseResults(unsigned int keep) {
    eraseVec(atmosphericAlpha, keep);
    for (auto& r : receptors) {
        eraseVec(r.level, keep);
        eraseVec(r.levelByEmittor, keep);
    }
}

// ── XML serialisation ─────────────────────────────────────────────────────────

void Noise::writeXML(ofstream& file, string tab) {
    file << tab << "<Noise";
    writeAttributes(file, attributes);
    file << ">" << endl;
    for (const auto& e : emittors) {
        file << tab << "\t<PointEmittor";
        writeAttributes(file, e.attributes);
        file << "/>" << endl;
    }
    for (const auto& r : receptors) {
        file << tab << "\t<PointReceptor";
        writeAttributes(file, r.attributes);
        file << "/>" << endl;
    }
    file << tab << "</Noise>" << endl;
}

void Noise::writeHeaderText(fstream& textFile) {
    textFile << id << "(" << key << "):atmosphericAlpha(dB/km)\t";
    for (size_t ri = 0; ri < receptors.size(); ++ri) {
        textFile << id << "(" << key << "):"
                 << receptors[ri].id << "(" << receptors[ri].key << "):level(dB)\t";
        for (size_t ei = 0; ei < emittors.size(); ++ei) {
            textFile << id << "(" << key << "):"
                     << receptors[ri].id << "(" << receptors[ri].key << "):"
                     << emittors[ei].id  << "(" << emittors[ei].key  << "):level(dB)\t";
        }
    }
}

void Noise::writeResultsText(fstream& textFile, unsigned int index) {
    textFile << fixed << setprecision(4) << atmosphericAlpha.at(index) << "\t";
    for (size_t ri = 0; ri < receptors.size(); ++ri) {
        textFile << fixed << setprecision(2) << receptors[ri].level.at(index) << "\t";
        for (size_t ei = 0; ei < emittors.size(); ++ei)
            textFile << fixed << setprecision(2)
                     << receptors[ri].levelByEmittor.at(index).at(ei) << "\t";
    }
}

size_t Noise::memoryUsage() const {
    size_t bytes = sizeof(double) * atmosphericAlpha.size();
    for (const auto& r : receptors) {
        bytes += sizeof(double) * r.level.size();
        bytes += sizeof(vector<double>) * r.levelByEmittor.size();
        for (const auto& v : r.levelByEmittor)
            bytes += sizeof(double) * v.size();
    }
    return bytes;
}

// ── XML attribute helpers ─────────────────────────────────────────────────────

map<string,string> Noise::readAttributes(TiXmlElement* elem) {
    map<string,string> m;
    for (TiXmlAttribute* a = elem->FirstAttribute(); a; a = a->Next())
        m[a->Name()] = a->Value();
    return m;
}

bool Noise::hasAttribute(const map<string,string>& attrs, const string& name) {
    return attrs.find(name) != attrs.end();
}

string Noise::getAttribute(const map<string,string>& attrs, const string& name, const string& def) {
    auto it = attrs.find(name);
    return (it != attrs.end()) ? it->second : def;
}

double Noise::getAttribute(const map<string,string>& attrs, const string& name, double def) {
    auto it = attrs.find(name);
    return (it != attrs.end()) ? to<double>(it->second) : def;
}

bool Noise::getAttribute(const map<string,string>& attrs, const string& name, bool def) {
    auto it = attrs.find(name);
    if (it == attrs.end()) return def;
    string v = it->second;
    transform(v.begin(), v.end(), v.begin(), ::toupper);
    return (v == "TRUE" || v == "1" || v == "YES");
}

void Noise::writeAttributes(ofstream& file, const map<string,string>& attrs) {
    for (const auto& kv : attrs)
        file << " " << kv.first << "=\"" << kv.second << "\"";
}
