#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>

#include "pollutant.h"

static const double TWO_PI = 2.*M_PI;

Pollutant::Pollutant(TiXmlHandle hdl, ostream* pLogStr):logStream(std::cout.rdbuf()) {

    associate(pLogStr, logStream);

    TiXmlElement* elem = hdl.ToElement();
    if (!elem) throw(string("Pollutant: invalid XML element."));

    attributes = readAttributes(elem);
    model = getAttribute(attributes, "model", string("plume"));
    id = getAttribute(attributes, "id", string("Pollutant"));
    key = getAttribute(attributes, "key", model);
    if (!hasAttribute(attributes, "model")) attributes["model"] = model;

    double roughnessLength = getAttribute(attributes, "roughnessLength", 0.1);
    double windReferenceHeight = getAttribute(attributes, "windReferenceHeight", 10.);
    double windSpeedMin = getAttribute(attributes, "windSpeedMin", 0.5);
    if (!(roughnessLength > 0.)) throw string("Error in XML pollutant configuration: roughnessLength must be strictly positive. ");
    if (!(windReferenceHeight > roughnessLength)) throw string("Error in XML pollutant configuration: windReferenceHeight must be greater than roughnessLength. ");
    if (!(windSpeedMin > 0.)) throw string("Error in XML pollutant configuration: windSpeedMin must be strictly positive. ");
    if (hasAttribute(attributes, "windSpeed") && getAttribute(attributes, "windSpeed", 0.) < 0.)
        throw string("Error in XML pollutant configuration: windSpeed must be positive or zero. ");

    unsigned int i=0;
    TiXmlElement* child = hdl.ChildElement("PointEmittor",i).ToElement();
    while (child) {
        PointEmittor emittor;
        emittor.attributes = readAttributes(child);
        emittor.id = getAttribute(emittor.attributes, "id", toString(i));
        emittor.key = getAttribute(emittor.attributes, "key", string(""));
        emittor.x = getAttribute(emittor.attributes, "x", 0.);
        emittor.y = getAttribute(emittor.attributes, "y", 0.);
        emittor.z = getAttribute(emittor.attributes, "z", 0.);
        emittor.emissionRate = getAttribute(emittor.attributes, "emissionRate", 0.);
        emittor.plumeRise = getAttribute(emittor.attributes, "plumeRise", 0.);
        emittors.push_back(emittor);
        child = hdl.ChildElement("PointEmittor",++i).ToElement();
    }

    if (emittors.empty()) throw string("Error in XML pollutant configuration: no Emittor defined. ");

    i=0;
    child = hdl.ChildElement("PointReceptor",i).ToElement();
    while (child) {
        PointReceptor receptor;
        receptor.attributes = readAttributes(child);
        receptor.id = getAttribute(receptor.attributes, "id", toString(i));
        receptor.key = getAttribute(receptor.attributes, "key", string(""));
        receptor.x = getAttribute(receptor.attributes, "x", 0.);
        receptor.y = getAttribute(receptor.attributes, "y", 0.);
        receptor.z = getAttribute(receptor.attributes, "z", 0.);
        receptors.push_back(receptor);
        child = hdl.ChildElement("PointReceptor",++i).ToElement();
    }

    if (receptors.empty()) throw string("Error in XML pollutant configuration: no Receptor defined. ");

    logStream << "Pollutant plume model: " << emittors.size() << " point emittor(s), " << receptors.size() << " point receptor(s) loaded." << endl << flush;
}

void Pollutant::compute(Climate* pClimate, unsigned int day, unsigned int hour) {

    if (!pClimate) throw(string("Pollutant: no climate data available."));

    double uSum = 0.;
    for (size_t i=0; i<emittors.size(); ++i) uSum += getWindSpeedAtEmittorHeight(pClimate, emittors.at(i), day, hour);
    double u = uSum/emittors.size(); 
    double direction = pClimate->getWindDirection(day,hour);
    string stability = getStabilityClass(pClimate, day, hour, u);

    windSpeed.push_back(u);
    windDirection.push_back(direction);
    stabilityClass.push_back(stability);

    for (size_t receptorIndex=0; receptorIndex<receptors.size(); ++receptorIndex) {
        vector<double> byEmittor(emittors.size(),0.);
        double concentration = getAttribute(attributes, "backgroundConcentration", 0.);

        for (size_t emittorIndex=0; emittorIndex<emittors.size(); ++emittorIndex) {
            double uEmittor = getWindSpeedAtEmittorHeight(pClimate, emittors.at(emittorIndex), day, hour);
            byEmittor.at(emittorIndex) = computeConcentration(emittors.at(emittorIndex), receptors.at(receptorIndex), uEmittor, direction, stability);
            concentration += byEmittor.at(emittorIndex);
        }
        receptors.at(receptorIndex).concentration.push_back(concentration);
        receptors.at(receptorIndex).concentrationByEmittor.push_back(byEmittor);
    }
}

double Pollutant::computeConcentration(const PointEmittor& emittor, const PointReceptor& receptor, double u, double direction, const string& stability) const {

    if (emittor.emissionRate <= 0.) return 0.;

    // Meteorological wind direction is the direction from which the wind blows.
    double windRad = direction*M_PI/180.;
    double windX = -sin(windRad);
    double windY = -cos(windRad);

    double dx = receptor.x - emittor.x;
    double dy = receptor.y - emittor.y;
    double x = dx*windX + dy*windY;
    double y = -dx*windY + dy*windX;

    double downwindMin = getAttribute(attributes, "downwindMin", 1.);
    if (x <= 0.) return 0.;
    x = max(x, downwindMin);

    pair<double,double> sigma = getSigma(x, stability);
    double sigmaY = sigma.first;
    double sigmaZ = sigma.second;
    if (!(sigmaY > 0.) || !(sigmaZ > 0.)) return 0.;

    // H is the effective source height. Unlike AERO, plumeRise is entered by
    // the user: it is not calculated from stack and exhaust-gas properties.
    double H = emittor.z + emittor.plumeRise;
    double z = receptor.z;
    double crossWindTerm = exp(-0.5*pow(y/sigmaY,2.));
    double verticalTerm = exp(-0.5*pow((z-H)/sigmaZ,2.));
    if (getAttribute(attributes, "groundReflection", true))
        verticalTerm += exp(-0.5*pow((z+H)/sigmaZ,2.));

    double concentration = emittor.emissionRate/(TWO_PI*u*sigmaY*sigmaZ)*crossWindTerm*verticalTerm;

    // These optional exponential attenuations are simple approximations. They
    // do not implement AERO's suspended-dust and dust-deposition calculations.
    double decayRate = getAttribute(attributes, "decayRate", 0.); // 1/s
    if (decayRate > 0.) concentration *= exp(-decayRate*x/u);

    double depositionVelocity = getAttribute(attributes, "depositionVelocity", 0.); // m/s
    if (depositionVelocity > 0.) {
        double mixingHeight = getAttribute(attributes, "mixingHeight", 1000.);
        concentration *= exp(-depositionVelocity*x/(u*max(mixingHeight,1.)));
    }

    return concentration;
}

double Pollutant::getWindSpeedAtEmittorHeight(Climate* pClimate, const PointEmittor& emittor, unsigned int day, unsigned int hour) const {

    // windSpeed, when present, overrides the hourly climate wind speed at the
    // reference height. The logarithmic profile is a CitySim model choice:
    // AERO uses a stability-dependent power law instead.
    double referenceWindSpeed = hasAttribute(attributes, "windSpeed") ? getAttribute(attributes, "windSpeed", 0.) : pClimate->getWindSpeed(day,hour);
    double roughnessLength = getAttribute(attributes, "roughnessLength", 0.1); // m
    double referenceHeight = getAttribute(attributes, "windReferenceHeight", 10.); // m
    double height = max(emittor.z + emittor.plumeRise, 1.01*roughnessLength);
    double u = referenceWindSpeed*log(height/roughnessLength)/log(referenceHeight/roughnessLength);
    return max(u,getAttribute(attributes, "windSpeedMin", 0.5));
}

string Pollutant::getStabilityClass(Climate* pClimate, unsigned int day, unsigned int hour, double u) const {

    // AUTO selects a conventional Pasquill class from the current CitySim
    // irradiance, cloud cover and wind. These are not AERO's six regulatory
    // atmospheric-equilibrium states or its prescribed meteorological cases.
    string stability = upperCase(getAttribute(attributes, "stabilityClass", string("auto")));
    if (stability!=string("AUTO")) return stability.substr(0,1);

    double Igh = pClimate->getIgh(day, hour);
    double cloud = pClimate->getCloudCoverFraction(day,hour);

    if (Igh > 10.) {
        if (Igh > 600.) {
            if (u < 2.) return "A";
            if (u < 3.) return "B";
            if (u < 5.) return "C";
            return "D";
        }
        if (Igh > 300.) {
            if (u < 2.) return "B";
            if (u < 5.) return "C";
            return "D";
        }
        if (u < 2.) return "C";
        return "D";
    }

    if (cloud > 0.5) return "D";
    if (u < 2.) return "F";
    if (u < 5.) return "E";
    return "D";
}

pair<double,double> Pollutant::getSigma(double x, const string& stability) const {

    // Generic Pasquill-Gifford coefficients for the hourly Gaussian model.
    // AERO uses the coefficients defined by the Polish reference methodology.
    if (hasAttribute(attributes, "sigmaY") && hasAttribute(attributes, "sigmaZ")) {
        return pair<double,double>(getAttribute(attributes, "sigmaY", 0.), getAttribute(attributes, "sigmaZ", 0.));
    }

    double sigmaY=0., sigmaZ=0.;
    string s = stability.empty() ? string("D") : upperCase(stability).substr(0,1);

    if (s==string("A")) {
        sigmaY = 0.22*x*pow(1.+0.0001*x,-0.5);
        sigmaZ = 0.20*x;
    }
    else if (s==string("B")) {
        sigmaY = 0.16*x*pow(1.+0.0001*x,-0.5);
        sigmaZ = 0.12*x;
    }
    else if (s==string("C")) {
        sigmaY = 0.11*x*pow(1.+0.0001*x,-0.5);
        sigmaZ = 0.08*x*pow(1.+0.0002*x,-0.5);
    }
    else if (s==string("E")) {
        sigmaY = 0.06*x*pow(1.+0.0001*x,-0.5);
        sigmaZ = 0.03*x*pow(1.+0.0003*x,-1.);
    }
    else if (s==string("F")) {
        sigmaY = 0.04*x*pow(1.+0.0001*x,-0.5);
        sigmaZ = 0.016*x*pow(1.+0.0003*x,-1.);
    }
    else {
        sigmaY = 0.08*x*pow(1.+0.0001*x,-0.5);
        sigmaZ = 0.06*x*pow(1.+0.0015*x,-0.5);
    }

    sigmaY *= getAttribute(attributes, "sigmaYFactor", 1.);
    sigmaZ *= getAttribute(attributes, "sigmaZFactor", 1.);
    return pair<double,double>(sigmaY, sigmaZ);
}

void Pollutant::clear() {
    windSpeed.clear();
    windDirection.clear();
    stabilityClass.clear();
    for (size_t i=0; i<receptors.size(); ++i) {
        receptors.at(i).concentration.clear();
        receptors.at(i).concentrationByEmittor.clear();
    }
}

void Pollutant::eraseResults(unsigned int keepValue) {
    eraseVector(windSpeed, keepValue);
    eraseVector(windDirection, keepValue);
    eraseVector(stabilityClass, keepValue);
    for (size_t i=0; i<receptors.size(); ++i) {
        eraseVector(receptors.at(i).concentration, keepValue);
        if (receptors.at(i).concentrationByEmittor.size() > keepValue)
            receptors.at(i).concentrationByEmittor.erase(receptors.at(i).concentrationByEmittor.begin(), receptors.at(i).concentrationByEmittor.end()-min(keepValue,(unsigned int)receptors.at(i).concentrationByEmittor.size()));
    }
}

void Pollutant::writeXML(ofstream& file, string tab) {

    file << tab << "<Pollutant";
    writeAttributes(file, attributes);
    file << ">" << endl;

    for (size_t i=0; i<emittors.size(); ++i) {
        file << tab << "\t<PointEmittor";
        writeAttributes(file, emittors.at(i).attributes);
        file << "/>" << endl;
    }
    for (size_t i=0; i<receptors.size(); ++i) {
        file << tab << "\t<PointReceptor";
        writeAttributes(file, receptors.at(i).attributes);
        file << "/>" << endl;
    }
    file << tab << "</Pollutant>" << endl;
}

void Pollutant::writeHeaderText(fstream& textFile) {

    textFile << id << "(" << key << "):windSpeed(m/s)\t"
             << id << "(" << key << "):windDirection(deg)\t"
             << id << "(" << key << "):stabilityClass\t";
    for (size_t receptorIndex=0; receptorIndex<receptors.size(); ++receptorIndex) {
        textFile << id << "(" << key << "):" << receptors.at(receptorIndex).id
                 << "(" << receptors.at(receptorIndex).key << "):concentration(kg/m3)\t";
        for (size_t emittorIndex=0; emittorIndex<emittors.size(); ++emittorIndex) {
            textFile << id << "(" << key << "):" << receptors.at(receptorIndex).id
                     << "(" << receptors.at(receptorIndex).key << "):" << emittors.at(emittorIndex).id
                     << "(" << emittors.at(emittorIndex).key << "):concentration(kg/m3)\t";
        }
    }
}

void Pollutant::writeResultsText(fstream& textFile, unsigned int index) {

    textFile << fixed << setprecision(4) << windSpeed.at(index) << "\t"
             << fixed << setprecision(1) << windDirection.at(index) << "\t"
             << stabilityClass.at(index) << "\t";
    for (size_t receptorIndex=0; receptorIndex<receptors.size(); ++receptorIndex) {
        textFile << scientific << setprecision(8) << receptors.at(receptorIndex).concentration.at(index) << "\t";
        for (size_t emittorIndex=0; emittorIndex<emittors.size(); ++emittorIndex) {
            textFile << scientific << setprecision(8) << receptors.at(receptorIndex).concentrationByEmittor.at(index).at(emittorIndex) << "\t";
        }
    }
}

size_t Pollutant::memoryUsage() const {

    size_t bytes = sizeof(double)*(windSpeed.size()+windDirection.size());
    bytes += sizeof(string)*stabilityClass.size();

    for (size_t receptorIndex=0; receptorIndex<receptors.size(); ++receptorIndex) {
        bytes += sizeof(double)*receptors.at(receptorIndex).concentration.size();
        bytes += sizeof(vector<double>)*receptors.at(receptorIndex).concentrationByEmittor.size();
        for (size_t resultIndex=0; resultIndex<receptors.at(receptorIndex).concentrationByEmittor.size(); ++resultIndex)
            bytes += sizeof(double)*receptors.at(receptorIndex).concentrationByEmittor.at(resultIndex).size();
    }
    return bytes;
}

map<string,string> Pollutant::readAttributes(TiXmlElement* elem) {

    map<string,string> attributes;
    for (TiXmlAttribute* attr=elem->FirstAttribute(); attr; attr=attr->Next())
        attributes[attr->Name()] = attr->Value();
    return attributes;
}

bool Pollutant::hasAttribute(const map<string,string>& attrs, const string& name) {
    return attrs.find(name)!=attrs.end();
}

string Pollutant::getAttribute(const map<string,string>& attrs, const string& name, const string& defaultValue) {
    map<string,string>::const_iterator it = attrs.find(name);
    if (it==attrs.end()) return defaultValue;
    return it->second;
}

double Pollutant::getAttribute(const map<string,string>& attrs, const string& name, double defaultValue) {
    map<string,string>::const_iterator it = attrs.find(name);
    if (it==attrs.end()) return defaultValue;
    return to<double>(it->second);
}

bool Pollutant::getAttribute(const map<string,string>& attrs, const string& name, bool defaultValue) {
    map<string,string>::const_iterator it = attrs.find(name);
    if (it==attrs.end()) return defaultValue;
    string value = upperCase(it->second);
    return (value==string("TRUE") || value==string("1") || value==string("YES"));
}

void Pollutant::writeAttributes(ofstream& file, const map<string,string>& attrs) {
    for (map<string,string>::const_iterator it=attrs.begin(); it!=attrs.end(); ++it)
        file << " " << it->first << "=\"" << it->second << "\"";
}

string Pollutant::upperCase(string value) {
    transform(value.begin(), value.end(), value.begin(), ::toupper);
    return value;
}

void Pollutant::eraseVector(vector<double>& values, unsigned int keepValue) {
    if (values.size() > keepValue)
        values.erase(values.begin(),values.end()-min(keepValue,(unsigned int)values.size()));
}

void Pollutant::eraseVector(vector<string>& values, unsigned int keepValue) {
    if (values.size() > keepValue)
        values.erase(values.begin(),values.end()-min(keepValue,(unsigned int)values.size()));
}
