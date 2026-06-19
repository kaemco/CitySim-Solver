#ifndef POLLUTANT_H
#define POLLUTANT_H

#include <map>
#include <string>
#include <vector>
#include <fstream>

#include "tinyxml.h"
#include "climate.h"
#include "util.h"

using namespace std;

/**
 * Hourly Gaussian plume model using CitySim coordinates and climate data.
 *
 * This is not an implementation of the full AERO 2024 regulatory methodology.
 * In particular, the wind profile, stability classes, dispersion coefficients,
 * plume rise and optional deposition attenuation are simplified inputs or
 * approximations. The XML attributes are documented where they are used in
 * pollutant.cpp.
 *
 * Example:
 * <Pollutant model="plume" windSpeed="4" windReferenceHeight="10"
 *            roughnessLength="0.1" windSpeedMin="0.5">
 *     <PointEmittor id="E1" x="0" y="0" z="2" emissionRate="0.001"/>
 *     <PointReceptor id="R1" x="10" y="0" z="1.5"/>
 * </Pollutant>
 *
 * windSpeed is optional and expressed in m/s at windReferenceHeight. When it
 * is omitted, the hourly wind speed comes from the CitySim climate file.
 */
class Pollutant {

public:
    struct PointEmittor {
        string id;
        string key;
        map<string,string> attributes;
        double x=0., y=0., z=0.;
        double emissionRate=0.; //!< pollutant mass flow rate (kg/s)
        double plumeRise=0.;    //!< user-imposed plume rise above the source point (m)
    };

    struct PointReceptor {
        string id;
        string key;
        map<string,string> attributes;
        double x=0., y=0., z=0.;
        vector<double> concentration;
        vector<vector<double> > concentrationByEmittor;
    };

private:
    string id;
    string key;
    string model;
    map<string,string> attributes;
    vector<PointEmittor> emittors;
    vector<PointReceptor> receptors;
    vector<double> windSpeed;
    vector<double> windDirection;
    vector<string> stabilityClass;
    ostream logStream;

    static map<string,string> readAttributes(TiXmlElement* elem);
    static bool hasAttribute(const map<string,string>& attrs, const string& name);
    static string getAttribute(const map<string,string>& attrs, const string& name, const string& defaultValue);
    static double getAttribute(const map<string,string>& attrs, const string& name, double defaultValue);
    static bool getAttribute(const map<string,string>& attrs, const string& name, bool defaultValue);
    static void writeAttributes(ofstream& file, const map<string,string>& attrs);
    static string upperCase(string value);
    static void eraseVector(vector<double>& values, unsigned int keepValue);
    static void eraseVector(vector<string>& values, unsigned int keepValue);

    string getStabilityClass(Climate* pClimate, unsigned int day, unsigned int hour, double u) const;
    pair<double,double> getSigma(double x, const string& stability) const;
    double getWindSpeedAtEmittorHeight(Climate* pClimate, const PointEmittor& emittor, unsigned int day, unsigned int hour) const;
    double computeConcentration(const PointEmittor& emittor, const PointReceptor& receptor, double u, double direction, const string& stability) const;

public:
    Pollutant(TiXmlHandle hdl, ostream* pLogStr=&std::cout);
    void compute(Climate* pClimate, unsigned int day, unsigned int hour);
    void clear();
    void eraseResults(unsigned int keepValue);
    void writeXML(ofstream& file, string tab="");
    void writeHeaderText(fstream& textFile);
    void writeResultsText(fstream& textFile, unsigned int index);
    size_t memoryUsage() const;

    string getModel() const { return model; }
    size_t getnEmittors() const { return emittors.size(); }
    size_t getnReceptors() const { return receptors.size(); }
    size_t getnResults() const { return windSpeed.size(); }
};

#endif
