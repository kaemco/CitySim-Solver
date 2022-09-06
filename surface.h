#ifndef SURFACE_H
#define SURFACE_H

#include <limits>

#include "tinyxml.h"
#include "util.h"
#include "plant.h"
#include <iomanip>

#include "DATASurfaceDelegateABC.h"
#include "GEOMPolygonInfo.h"

class Building;
class XmlScene;

using namespace std;

class Layer {

private:

  float xx, kw, Cp, rho;
  float nre, gwp, ubp; // all values per kg

public:

    Layer(float thickness, float conductivity, float specificHeat, float density, float nre=0.f, float gwp=0.f, float ubp=0.f) :
     xx(thickness), kw(conductivity), Cp(specificHeat), rho(density), nre(nre), gwp(gwp), ubp(ubp) {}

    Layer(Layer* l): xx(l->xx), kw(l->kw), Cp(l->Cp), rho(l->rho), nre(l->nre), gwp(l->gwp), ubp(l->ubp) {}

    void setxx(float meter) { xx = meter; }
    float getxx() { return xx;  }
    float getkw() { return kw;  }
    float getCp() { return Cp;  }
    float getrho(){ return rho; }
    float getNRE(){ return rho*xx*nre; /* in MJ/m2 */ }
    float getGWP(){ return rho*xx*gwp; /* in kgCO2/m2 */ }
    float getUBP(){ return rho*xx*ubp; /* in points/m2 */ }

    float getResistance() { return xx/kw; }
    float getConductance() { return kw/xx; }

    bool equals(const Layer& l) const {
        if (l.xx != xx) return false;
        if (l.kw != kw) return false;
        if (l.Cp != Cp) return false;
        if (l.rho != rho) return false;
        if (l.nre != nre) return false;
        if (l.gwp != gwp) return false;
        if (l.ubp != ubp) return false;
        return true;
    }

    friend bool operator==(const Layer& l1, const Layer& l2){ return l1.equals(l2);}
    friend bool operator!=(const Layer& l1, const Layer& l2){ return !(l1.equals(l2));}

    virtual void writeXML(ofstream& file, bool isInsulation, string tab){
        file << tab << "<Layer Thickness=\"" << xx << "\" Conductivity=\"" << kw << "\"";
        file << " Cp=\"" << Cp << "\" Density=\"" << rho << "\" NRE=\"" << nre << "\" GWP=\"" << gwp << "\" UBP=\"" << ubp << "\"";
        if(isInsulation) file << " insulation=\"true\"";
        file << "/>" << endl;
    }

    void writeGML(ofstream& file, string tab) {
        file << tab << "<energy:layer>" << endl;
        file << tab << tabs(1) << "<energy:Layer>" << endl;
        file << tab << tabs(2) << "<energy:layerComponent>" << endl;
        file << tab << tabs(3) << "<energy:LayerComponent>" << endl;
        file << tab << tabs(4) << "<energy:thickness uom=\"m\">" << xx << "</energy:thickness>" << endl;
        file << tab << tabs(4) << "<energy:material>" << endl;
        file << tab << tabs(5) << "<energy:SolidMaterial>" << endl;
        file << tab << tabs(6) << "<energy:conductivity uom=\"W/(mK)\">" << kw << "</energy:conductivity>" << endl;
        file << tab << tabs(6) << "<energy:density uom=\"kg/m3\">" << rho << "</energy:density>" << endl;
        file << tab << tabs(6) << "<energy:specificHeat uom=\"J/(kgK)\">" << Cp << "</energy:specificHeat>" << endl;
        if (gwp > 0.)
            file << tab << tabs(6) << "<energy:embodiedCarbon uom=\"kgCO2/kg\">" << gwp << "</energy:embodiedCarbon>" << endl;
        if (nre > 0.)
            file << tab << tabs(6) << "<energy:embodiedEnergy uom=\"MJ/kg\">" << nre << "</energy:embodiedEnergy>" << endl;
        file << tab << tabs(5) << "</energy:SolidMaterial>" << endl;
        file << tab << tabs(4) << "</energy:material>" << endl;
        file << tab << tabs(3) << "</energy:LayerComponent>" << endl;
        file << tab << tabs(2) << "</energy:layerComponent>" << endl;
        file << tab << tabs(1) << "</energy:Layer>" << endl;
        file << tab << "</energy:layer>" << endl;
    }

};

// To read default values in Pro more easily
class Material {

private:

    int id;
    string name;
    float conductivity, cp, density;
    float nre, gwp, ubp; // all values per kg

public:
    ostream logStream;

    Material(TiXmlHandle hdl, ostream* pLogStr=NULL);

    string getName() { return name; }
    float getConductivity() { return conductivity; }
    float getCp() { return cp;  }
    float getDensity(){ return density; }
    float getNRE(){ return nre; }
    float getGWP(){ return gwp; }
    float getUBP(){ return ubp; }

    void writeXML(ofstream& file, string tab=""){
        file << tab << "<Material name=\"" << name << "\" Conductivity=\"" << conductivity << "\" Cp=\"" << cp << "\" Density=\"" << density <<
        "\" nre=\"" << nre << "\" gwp=\"" << gwp << "\" ubp=\"" << ubp << "\"/>" << endl;
    }

};


class Composite {

 private :
    string name="";
    string category="";
    vector<Layer> vLayer;
    float Uvalue;
	int id;
    Layer* insulationLayer=NULL;

 public :

    ostream logStream;

    Composite(TiXmlHandle hdl, ostream* pLogStr=NULL);

    Composite(float Uvalue, int wtId, string name, ostream* pLogStr=NULL);

    Composite(TiXmlHandle hdl, map<string,Material*> materials, ostream* pLogStr=NULL);

    size_t getnLayers() { return vLayer.size(); }
    vector<Layer>* getLayers() { return &vLayer; }
    Layer* getLayer(unsigned int i) { return &(vLayer.at(i));}

    void addLayer(float thickness, Material* m, bool isInsulationLayer);

    int getId() { return id;}
    void setId(unsigned int newId){ id=newId;}
    string getName() { return name;}
    float getUvalue() { return Uvalue; }
    string getCategory() { return name;}
    void setCategory(string cat) { category=cat;}

    Layer* getInsulationLayer(){
        return insulationLayer;
    }

    // gets a simplified version of the composite
    void getSimplifiedNode(float& Cw, float& Kw1, float& Kw2);

    bool equals(const Composite& c) const;
    friend bool operator==(const Composite& c1, const Composite& c2){ return c1.equals(c2);}

    float getDepth() { // for a Composite obtained through getComposite (i.e. with correct insulation thickness)
        float depth = 0.f;
        for (size_t i=0; i<vLayer.size(); ++i) {
                depth += vLayer[i].getxx();
        }
        return depth;
    }

    float getDiffusivity() {
        // returns the diffusivity of the last layer (the inner one in the ground) in m^2/h
        return (vLayer.back().getkw()/(vLayer.back().getrho()*vLayer.back().getCp()))*3600.f*24.f;
    }

    float getConductance() { // W/(m≤K), without the conductive air layer(s) on the side(s), for a Composite obtained through getComposite
        float resistance = 0.f; // the resistance is in m≤K/W, sum of the resistances in series, inverse to get the conductance
        if (vLayer.empty()) { resistance = 1.f/Uvalue - 0.125f - 0.04f; }
        for (unsigned int layerIndex=0; layerIndex<vLayer.size(); ++layerIndex){
                resistance += vLayer[layerIndex].getxx()/vLayer[layerIndex].getkw();
        }
        return 1.f/resistance;
    }
    float getResistance() { // m≤K/W, without the conductive air layer(s) on the side(s), for a Composite obtained through getComposite
        float resistance = 0.f; // the resistance is in m≤K/W, sum of the resistances in series
        if (vLayer.empty()) { resistance = 1.f/Uvalue - 0.125f - 0.04f; }
        for (unsigned int layerIndex=0; layerIndex<vLayer.size(); ++layerIndex){
                resistance += vLayer[layerIndex].getxx()/vLayer[layerIndex].getkw();
        }
        return resistance;
    }
    float getCapacitance() { // J/(kg K m≤), for a Composite obtained through getComposite
        float capacitance = 0.f;
        for (unsigned int layerIndex=0; layerIndex<vLayer.size(); ++layerIndex){
                capacitance += vLayer[layerIndex].getrho()*vLayer[layerIndex].getxx()*vLayer[layerIndex].getCp();
        }
        return capacitance;
    }

    float getCapacitance(unsigned int i) { // J/(kg K m≤), for a Composite obtained through getComposite
            return vLayer.at(i).getrho()*vLayer.at(i).getxx()*vLayer.at(i).getCp();
    }
    float getNRE() { /* in MJ/m2 for a Composite obtained through getComposite (with all layers) */
        float nre = 0.f;
        for (unsigned int layerIndex=0; layerIndex<vLayer.size(); ++layerIndex){
                nre += vLayer[layerIndex].getNRE();
        }
        return nre;
    }
    float getGWP() { /* in kgCO2/m2 for a Composite obtained through getComposite (with all layers) */
        float gwp = 0.f;
        for (unsigned int layerIndex=0; layerIndex<vLayer.size(); ++layerIndex){
                gwp += vLayer[layerIndex].getGWP();
        }
        return gwp;
    }
    float getUBP() { /* in points/m2 for a Composite obtained through getComposite (with all layers) */
        float ubp = 0.f;
        for (unsigned int layerIndex=0; layerIndex<vLayer.size(); ++layerIndex){
                ubp += vLayer[layerIndex].getUBP();
        }
        return ubp;
    }

    void writeXML(ofstream& file, string tab);
    void writeGML(ofstream& file, string tab);

};

class Surface : public DATASurfaceDelegateABC {

public:
    enum SType {WALL, ROOF, FLOOR, GROUND, SURFACE};
    enum ResultsPositions {DAY_RESULTS_POSITION=0, MONTH_RESULTS_POSITION=365, YEAR_RESULTS_POSITION=365+12};
    enum DiplayableResults {SHORTWAVE, LONGWAVE, SURFACE_TEMP, PV_PRODUCTION, THERMAL_PRODUCTION, N_RESULTS};

protected:

    Building* b=nullptr;
    float area = 0.f;
    float volume = 0.f; // JK - 18.04.2015 - volume under the surface
    GENPoint norm;
    vector<GENPoint> vertices;

    //structure to hold the projected solid angles
    struct vault {
        double sky=numeric_limits<double>::quiet_NaN(), ground=numeric_limits<double>::quiet_NaN(), hemisphere=M_PI;
    } projectedSolidAngle;

    // added an id for the surface
    uint64_t id; //!< Surface ID
    string key="";

    Composite* composite=nullptr; // a pointer to the Composite

    // other values
    float shortWaveReflectance = 0.2f; //!< Surface short wave reflectance (Lambertian) \in[0,1]
    float longWaveEmissivity = 0.9f; //!< Surface long wave emissivity \in[0,1]
    float glazingRatio = 0.f;  //!< Proportion of glazing relative to the surface area \in[0,1]
    float glazingGvalue = 0.f; //!< Window transmittance (normal) for Visible
    float glazingGvalue_hemispherical = 0.f; //!< Window transmittance (hemispherical) for Visible
    float glazingUvalue = 0.f; //!< Window U-Value
    float glazingOpenableRatio = 0.f; //!< Fraction of the Window that one can open \in[0,1]

    // saves the results timestep by timestep
    vector<float> shortWaveIrradiance; //!< Hourly results for the shortwave irradiance calculation (W/m^2)
    float shortWaveIrradiance_IAM, beamIrradiance, beamAngle;
    vector<float> illuminance; //!< Hourly results for the external illuminance calculation (lux)

    // save the results timestep by timestep, internalIlluminance
    vector<float> illuminance0, illuminance1, internallyReflectedIlluminance; //!< Hourly results for the internal illuminance calculation (lux)

    // saves the results of the surface temperatuers
    vector<float> temperature; //!< Hourly results for the external surface temperature (∞C)
    vector<float> environmentalTemperature; //!< Hourly results for the environmental temperature (∞C)
    vector<float> hc; /// the convection coefficient (W/(m^2 K))

    // values to define for the PV panels
    float pvRatio = 0.f; //!< Proportion of the surface that is covered by PV panels \in[0,1]
    PhotoVoltaic *pvPanel=nullptr; //!< Pointer to a PV panel type

     // values to define for the Solar Thermal panels
    float stRatio = 0.f; //!< Proportion of the surface that is covered by Solar Thermal panels \in[0,1]
    SolarThermal *stPanel=nullptr; //!< Pointer to a Solar Thermal type

    // values to define the hybrid solar panels
    float pvtRatio = 0.f; //!< Proportion of the surface that is covered by Hybrid Solar panels \in[0,1]
    SolarHybrid *pvtPanel=nullptr;

    // values to save the pv and st production for this surface for each time step
    vector<float> pvProduction, stProduction;

    // name for E+ simulation
    string ep_id; //!< EnergyPlus id of the Surface

public:

    friend ostream& operator<<(std::ostream& os, const Surface& obj) {
        for (size_t i=0;i<obj.vertices.size();++i) {
            os << obj.vertices[i] << endl;
        }
        return os;
    }

    static const string STypeNames[5]; // ={"Wall","Roof","Floor","Ground","Surface"};

    ostream logStream;

    void writeXML(ofstream& file, string tab=""){
        // write the vertices
        for (unsigned int i=0; i< vertices.size(); ++i) {
            file << tab << "<V" << i <<" x=\"" << (vertices[i])[0] << "\" y=\"" << (vertices[i])[1] << "\" z=\"" << (vertices[i])[2] << "\"/>" << endl;
        }
        if (pvPanel) pvPanel->writeXML(file,pvRatio,tab);
        if (stPanel) stPanel->writeXML(file,stRatio,tab);
        /// TODO: write the hybrid panel
    }

    void writeGML(ofstream& file, string tab="", const vector<double>& origin={0.,0.,0.}){
        file << tab << "\t\t\t\t\t<gml:exterior>\n"
             << tab << "\t\t\t\t\t\t<gml:LinearRing>\n"
             << tab << "\t\t\t\t\t\t\t<gml:posList>\n";
        string subtab = tab + "\t\t\t\t\t\t\t\t";
        for (unsigned int i=0; i< vertices.size(); ++i){
            file << subtab << origin[0]+(vertices[i])[0] << " " << origin[1]+(vertices[i])[1] << " " << origin[2]+(vertices[i])[2] << "\n";
        }
        // repeat the first point
        file << subtab << origin[0]+(vertices[0])[0] << " " << origin[1]+(vertices[0])[1] << " " << origin[2]+(vertices[0])[2] << "\n";
        // close the tabs
        file << tab << "\t\t\t\t\t\t\t</gml:posList>\n"
             << tab << "\t\t\t\t\t\t</gml:LinearRing>\n"
             << tab << "\t\t\t\t\t</gml:exterior>" << endl;
    }

    virtual void writeGML_composedOf(ofstream& file, string tab="") {
        // opaque part
        //file << tab << "\t<energy:ThermalComponent>" << endl;
        file << tab << "<energy:area uom=\"m2\">" << (1.f-glazingRatio)*area << "</energy:area>" << endl;
        file << tab << "<energy:construction xlink:href=\"#op_" << getComposite()->getId() << "\"/>" << endl;
        //file << tab << "\t\t<energy:isGroundCoupled>false</energy:isGroundCoupled>" << endl;
        //file << tab << "\t\t<energy:isSunExposed>true</energy:isSunExposed>" << endl;
        //file << tab << "\t</energy:ThermalComponent>" << endl;
        // transparent part
        file << tab << "<energy:contains>" << endl;
        file << tab << "\t<energy:ThermalOpening>" << endl;
        file << tab << "\t\t<energy:area uom=\"m2\">" << glazingRatio*area << "</energy:area>" << endl;
        file << tab << "\t\t<energy:construction>" << endl;
        file << tab << "\t\t\t<energy:Construction>" << endl;
        file << tab << "\t\t\t\t<energy:uValue uom=\"W/K*m2\">" << glazingUvalue << "</energy:uValue>" << endl;
        file << tab << "\t\t\t\t<energy:opticalProperties>" << endl;
        file << tab << "\t\t\t\t\t<energy:OpticalProperties>" << endl;
        file << tab << "\t\t\t\t\t\t<energy:transmittance>" << endl;
        file << tab << "\t\t\t\t\t\t\t<energy:Transmittance>" << endl;
        file << tab << "\t\t\t\t\t\t\t\t<energy:fraction uom=\"ratio\">" << glazingGvalue << "</energy:fraction>" << endl;
        file << tab << "\t\t\t\t\t\t\t\t<energy:wavelengthRange>total</energy:wavelengthRange>" << endl;
        file << tab << "\t\t\t\t\t\t\t</energy:Transmittance>" << endl;
        file << tab << "\t\t\t\t\t\t</energy:transmittance>" << endl;
        file << tab << "\t\t\t\t\t</energy:OpticalProperties>" << endl;
        file << tab << "\t\t\t\t</energy:opticalProperties>" << endl;
        file << tab << "\t\t\t</energy:Construction>" << endl;
        file << tab << "\t\t</energy:construction>" << endl;
        //file << tab << "\t\t<energy:isGroundCoupled>false</energy:isGroundCoupled>" << endl;
        //file << tab << "\t\t<energy:isSunExposed>true</energy:isSunExposed>" << endl;
        file << tab << "\t</energy:ThermalOpening>" << endl;
        file << tab << "</energy:contains>" << endl;
    }

    Surface(unsigned int id=0):id(id),logStream(std::cout.rdbuf()) {}
    Surface(unsigned int id, float shortWaveReflectance, float glazingGvalue, float glazingUvalue, float glazingOpenableRatio, float glazingRatio, ostream* pLogStr=NULL)
    :id(id),shortWaveReflectance(shortWaveReflectance),glazingRatio(glazingRatio),glazingGvalue(glazingGvalue),glazingUvalue(glazingUvalue),
     glazingOpenableRatio(glazingOpenableRatio),logStream(std::cout.rdbuf()) {
        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!= NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));
        // compute once for all the hemispherical g-value
        glazingGvalue_hemispherical = computeGlazingGvalueHemispherical();
    }
    // copy constructor
    Surface(const Surface& s)
        :vertices(s.vertices),id(s.id),key(s.key),shortWaveReflectance(s.shortWaveReflectance),glazingRatio(s.glazingRatio),glazingGvalue(s.glazingGvalue),
          glazingUvalue(s.glazingUvalue),glazingOpenableRatio(s.glazingOpenableRatio),logStream(s.logStream.rdbuf()) {
        // compute once for all the hemispherical g-value
        glazingGvalue_hemispherical = computeGlazingGvalueHemispherical();
        computeNormalAndArea();
    }
    Surface(const Surface& s, string key):Surface(s) { this->key = key; }
    Surface(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr=nullptr);
    Surface(unsigned int id, TiXmlElement* pPosList, vector<double> origin, string key=""); // contructor from CityGML posList
    ~Surface() {
        //cout << "Deleting surface...";
        if (pvPanel) delete pvPanel;
        if (stPanel) delete stPanel;
        if (pvtPanel) delete pvtPanel;
        //cout << " done." << endl;
    }
    void clear() {
        shortWaveIrradiance.clear();
        illuminance.clear();
        illuminance0.clear(); illuminance1.clear(); internallyReflectedIlluminance.clear();
        temperature.clear();
        environmentalTemperature.clear();
        hc.clear();
    }

    uint64_t getId() const { return id; }
    string getKey() { return key; }
    virtual SType getType(){return SURFACE;}

    void pushVertex(float x, float y, float z) { vertices.push_back(GENPoint::Cartesian(x,y,z)); }
    void pushVertex(GENPoint *p) { vertices.push_back(*p); }
    vector<GENPoint>* getVertices(){return &vertices;}
    void popVertex() { vertices.pop_back(); }
    void popVertex(size_t i) { vertices.erase(vertices.begin()+i); }

    void reverseOrientation() {
        vector<GENPoint> newVertices;
        for(int i=vertices.size()-1; i>=0; --i){
            newVertices.push_back(vertices[i]);
        }
        vertices = newVertices;
    }

    void translate(GENPoint p) {
        vector<GENPoint> newVertices;
        for (size_t i=0; i<vertices.size(); ++i) {
            newVertices.push_back(vertices[i]+p);
        }
        vertices = newVertices;
    }

    void removeRepeatedVertices() {
        if (vertices.size() > 1) {
            // remove the double vertices
            for (size_t i=1; i<vertices.size(); ++i)
                if (vertices[vertices.size()-i]==vertices[vertices.size()-i-1])
                    vertices.erase(vertices.begin()+vertices.size()-i);
        }
        if (vertices.size() > 1) {
            // remove the last vertex if equal to the first one
            if (vertices[0]==vertices[vertices.size()-1]) vertices.erase(vertices.begin()+vertices.size()-1);
        }
        if (vertices.size() > 2) {
            // check for collinear vectors and remove the central point
            size_t i = 0;
            GEN::Point<double>  a,b,crossProdab;
            while (i<vertices.size()) {
                a = vertices[i]-vertices[(i==0)?(vertices.size()-1):(i-1)];
                b = vertices[(i==(vertices.size()-1))?(0):(i+1)]-vertices[(i==0)?(vertices.size()-1):(i-1)];
                crossProdab=cross_product(a,b);
                if (dot_product(crossProdab,crossProdab)==0.) vertices.erase(vertices.begin()+i);
                else ++i;
            }
        }
    }

    void setBuildingRef(Building* bRef){ b = bRef;}
    Building* getBuilding(){ return b;}

	unsigned int vertexCount() const { return vertices.size(); }
    GENPoint* getVertex(unsigned int i) { return &vertices[i]; }

	void sendVertices(VertexVisitor &v) const
	{
		for (std::vector<GENPoint>::const_iterator it=vertices.begin();
			it!=vertices.end();
			++it)
			v(*it);
	}

	float getArea() const { return area; }
	float getVolume() const { return volume; }
    float getRadius() const;

	void setProjectedSolidAngle(double sky, double ground, double hemisphere) { projectedSolidAngle.sky=sky; projectedSolidAngle.ground=ground; projectedSolidAngle.hemisphere=hemisphere; }
	double getProjectedSolidAngle_sky() { return projectedSolidAngle.sky; }
	double getProjectedSolidAngle_ground() { return projectedSolidAngle.ground; }
	double getProjectedSolidAngle_hemisphere() { return projectedSolidAngle.hemisphere; }
	float getSVF() const { return projectedSolidAngle.sky/projectedSolidAngle.hemisphere; }
	float getGVF() const { return projectedSolidAngle.ground/projectedSolidAngle.hemisphere; }

	GENPoint normal() const { return norm; }

	void computeNormal() {
        GEOMPolygonInfo info(vertices.begin(), vertices.end());
        norm=info.normal();
	}

	void computeArea() {
        GEOMPolygonInfo info(vertices.begin(), vertices.end());
        area=info.area();
	}

	void computeNormalAndArea() {
        GEOMPolygonInfo info(vertices.begin(), vertices.end());
        norm=info.normal();
        area=info.area();
        volume=info.volume();
	}

	GENPoint computeCentroid() {
        GEOMPolygonInfo info(vertices.begin(), vertices.end());
        return info.centroid();
	}

// added methods to get surface characteristics

    Composite* getComposite() { return composite; }
    void setComposite(Composite* c) { composite = c; }
    virtual float getShortWaveReflectance() { return shortWaveReflectance; }
    void setShortWaveReflectance(float r) { shortWaveReflectance=r; }
    virtual float getShortWaveReflectance_opaque() { return shortWaveReflectance; }

    float getLongWaveEmissivity() { return longWaveEmissivity; }
    void setLongWaveEmissivity(float emissivity) { longWaveEmissivity=emissivity; }

    float getAzimuth() { return fmod(norm.Azimuth().degrees()+360.f,360.f); } // between 0 and 360∞, the original methode returns atan2(x,y) which is 0∞ for y and 90∞ for x (clockwise)
    float getAltitude() { return norm.Altitude().degrees(); }

    //float getGlazingRatio() { return glazingRatio; }
    float getGlazingGvalue() { return glazingGvalue; }
    void setGlazingGvalue(float g) { glazingGvalue = g; }
    float getGlazingGvalue(float theta, float p=2.f, float q=4.f);
    float getGlazingGvalueHemispherical() { return glazingGvalue_hemispherical; }
    float computeGlazingGvalueHemispherical();
    void updateGlazingGvalueHemispherical() { glazingGvalue_hemispherical = computeGlazingGvalueHemispherical(); }
    float getGlazingUvalue() { return glazingUvalue; }
    void setGlazingUvalue(float u) { glazingUvalue = u; }
    float getGlazingOpenableRatio() { return glazingOpenableRatio; }
    void setGlazingOpenableRatio(float o) { glazingOpenableRatio = o; }
    float getPVRatio(){ return pvRatio; }
    void setPVRatio(float r) { pvRatio=r; }
    float getSTRatio(){ return stRatio; }
    void setSTRatio(float r) { stRatio=r; }
    bool getHasSolarThermal(){return (stPanel!=NULL);}; // Added by Max

    // method to set the PV panel
    void setPVPanel(float pvRatio, PhotoVoltaic* pv) {
        pvPanel = pv;
        this->pvRatio = pvRatio;
    }
    void setPVPanel(PhotoVoltaic* pv) {
        pvPanel = pv;
    }
    PhotoVoltaic* getPVPanel(){ return pvPanel; }

    // method to set the Solar Thermal panel
    void setSTPanel(float stRatio, SolarThermal* st) {
        stPanel = st;
        this->stRatio = stRatio;
    }
    void setSTPanel(SolarThermal* st) {
        stPanel = st;
    }
    SolarThermal* getSTPanel(){ return stPanel; }

    float getGlazingRatio() {
        if (pvRatio > 0 && stRatio > 0 )
        glazingRatio =min(1.0f - pvRatio-stRatio, glazingRatio);
        else
        if (pvRatio <= 0 && stRatio > 0 )
        glazingRatio =min(1.0f - stRatio, glazingRatio);
        else
        if (pvRatio > 0 && stRatio <=  0 )
        glazingRatio =min(1.0f - pvRatio, glazingRatio);
        return glazingRatio;
    }
    void setGlazingRatio(float g) { glazingRatio = g;}
    float getGlazingArea() { return area*glazingRatio; }

    // gets the actual surface temperature in celsius
    float getTemperature() { if (temperature.empty()) return 15.f; else return temperature.back(); } // this has to be a starting point, everything at 15°C at time step 0
    float getTemperature(unsigned int it) { if (temperature.empty()) return std::numeric_limits<float>::quiet_NaN(); else return temperature.at(it); } // returns the element it
    void setTemperature(float celsius) { temperature.push_back(celsius); }
    void eraseTemperature(unsigned int keepValue) { if (!temperature.empty()) temperature.erase(temperature.begin(),temperature.end()-min(keepValue,(unsigned int)temperature.size())); } // removes all elements but the last one
    void eraseTemperature_back() { if (!temperature.empty()) temperature.pop_back(); } // removes the last one

    float getEnvironmentalTemperature() { if (environmentalTemperature.empty()) return 15.f; else return environmentalTemperature.back(); }
    void setEnvironmentalTemperature(float celsius) { environmentalTemperature.push_back(celsius); }
    void eraseEnvironmentalTemperature(unsigned int keepValue) { environmentalTemperature.erase(environmentalTemperature.begin(),environmentalTemperature.end()-min(keepValue,(unsigned int)environmentalTemperature.size())); }
    float get_hr() { if (temperature.empty())
                        return longWaveEmissivity * 5.670373e-8 * 0.5 * pow((getEnvironmentalTemperature()+15.f+2.*273.15),3);
                    else
                        return longWaveEmissivity * 5.670373e-8 * 0.5 * pow((environmentalTemperature.back()+temperature.back()+2.*273.15),3); } // in W/(m^2 K)

    void set_hc(float value) { hc.push_back(value); }
    float get_hc() { if (hc.empty()) return 2.8f; else return hc.back(); }
    float get_hc(unsigned int it) { if (hc.empty()) return std::numeric_limits<float>::quiet_NaN(); else return hc.at(it); } // returns the element it
    void erase_hc(unsigned int keepValue) { if (!hc.empty()) hc.erase(hc.begin(),hc.end()-min(keepValue,(unsigned int)hc.size())); } // removes all elements but the last one

    void setShortWaveIrradiance(float wattPerSquareMeter) { shortWaveIrradiance.push_back(wattPerSquareMeter); }
    void addShortWaveIrradiance(float value) { if (shortWaveIrradiance.empty()) shortWaveIrradiance.push_back(value); else shortWaveIrradiance.back()+=value; }

    float getShortWaveIrradiance() { return shortWaveIrradiance.back(); } // returns the last element in W/m^2
    float getShortWaveIrradiance(unsigned int it){ return shortWaveIrradiance.at(it); } // returns the element it in W/m^2
    void eraseShortWaveIrradiance(unsigned int keepValue){ shortWaveIrradiance.erase(shortWaveIrradiance.begin(),shortWaveIrradiance.end()-min(keepValue,(unsigned int)shortWaveIrradiance.size())); } // removes all elements but the last one

    void setShortWaveIrradiance_IAM(float value) { shortWaveIrradiance_IAM = value; }

    void setBeamIrradiance(float value) { beamIrradiance = value; }
    void setBeamAngle(float value) { beamAngle = value; }
    float getBeamIrradiance() { return beamIrradiance; } // in W/m^2
    float getBeamAngle() { return beamAngle; }// in degrees

    void setLongWaveIrradiance(float wattPerSquareMeter) { if (!temperature.empty()) environmentalTemperature.push_back(pow(wattPerSquareMeter/(longWaveEmissivity * 5.670373e-8)+pow(getTemperature(environmentalTemperature.size())+273.15,4), 0.25f)-273.15); }
    float getLongWaveIrradiance()          { return longWaveEmissivity * 5.670373e-8 * (pow(environmentalTemperature.back()+273.15,4) - pow(getTemperature()+273.15,4)); }   // returns the last element in W/m^2
    float getLongWaveIrradiance(size_t it) { return longWaveEmissivity * 5.670373e-8 * (pow(environmentalTemperature.at(it)+273.15,4) - pow(getTemperature(it)+273.15,4)); }  // returns the element it in W/m^2
    float getLongWaveAbsorbed() { return longWaveEmissivity * 5.670373e-8 * pow(getEnvironmentalTemperature()+273.15,4); }  // returns the element it in W/m^2
    float getLongWaveAbsorbed(size_t it) { return longWaveEmissivity * 5.670373e-8 * pow(environmentalTemperature.at(it)+273.15,4); }  // returns the element it in W/m^2

    // sets and gets the external illuminance (in lx)
    void setIlluminance(float lux) { illuminance.push_back(lux); }
    float getIlluminance() { return illuminance.back(); } // returns the last element
    float getIlluminance(unsigned int it) { return illuminance.at(it); } // returns the element it
    void eraseIlluminance(unsigned int keepValue) { illuminance.erase(illuminance.begin(),illuminance.end()-min(keepValue,(unsigned int)illuminance.size())); } // removes all elements but the last one

    // sets and gets the internal illuminance
    void setInternalIlluminance(float *illuminance) { illuminance0.push_back(illuminance[0]); illuminance1.push_back(illuminance[1]); }
    float getInternalIlluminance0(unsigned int step) { return illuminance0.at(step); }
    float getInternalIlluminance0() { return illuminance0.back(); }
    void eraseInternalIlluminance0(unsigned int keepValue) { illuminance0.erase(illuminance0.begin(),illuminance0.end()-min(keepValue,(unsigned int)illuminance0.size())); }
    float getInternalIlluminance1(unsigned int step) { return illuminance1.at(step); }
    float getInternalIlluminance1() { return illuminance1.back(); }
    void eraseInternalIlluminance1(unsigned int keepValue) { illuminance1.erase(illuminance1.begin(),illuminance1.end()-min(keepValue,(unsigned int)illuminance1.size())); }

    // sets and gets the internally inter-reflected component
    void setInternallyReflectedIlluminance(float lux) { internallyReflectedIlluminance.push_back(lux); }
    float getInternallyReflectedIlluminance() { return internallyReflectedIlluminance.back(); }
    float getInternallyReflectedIlluminance(unsigned int step) { return internallyReflectedIlluminance[step]; }
    void eraseInternallyReflectedIlluminance(unsigned int keepValue) { internallyReflectedIlluminance.erase(internallyReflectedIlluminance.begin(),internallyReflectedIlluminance.end()-min(keepValue,(unsigned int)internallyReflectedIlluminance.size())); }

    // gets the total illuminance of the two points
    float getTotalInternalIlluminance0(unsigned int step) { return illuminance0[step] + internallyReflectedIlluminance[step]; }
    float getTotalInternalIlluminance0() { return illuminance0.back() + internallyReflectedIlluminance.back(); }
    float getTotalInternalIlluminance1(unsigned int step) { return illuminance1[step] + internallyReflectedIlluminance[step]; }
    float getTotalInternalIlluminance1() { return illuminance1.back() + internallyReflectedIlluminance.back(); }

    // gets the PV production
    float getPVElectricProduction(float Tout) {
        pvProduction.push_back(0.f);
        if (pvPanel != NULL) {
            pvProduction.back() += pvRatio*area*shortWaveIrradiance_IAM*pvPanel->getMaxPowerEfficiency(shortWaveIrradiance_IAM, Tout);
        }
        if (pvtPanel != NULL) {
            pvProduction.back() += pvtRatio*area*shortWaveIrradiance_IAM*pvtPanel->getMaxPowerEfficiency(shortWaveIrradiance_IAM, Tout);
        }
        return pvProduction.back();
    }
    float getPVElectricProduction(unsigned int it) { if (pvProduction.empty()) return std::numeric_limits<float>::quiet_NaN(); else return pvProduction.at(it); }
    void erasePVElectricProduction(unsigned int keepValue) { pvProduction.erase(pvProduction.begin(),pvProduction.end()-min(keepValue,(unsigned int)pvProduction.size())); }

    //Solar heater production
//    float getSolarThermalProduction(float tOut, float windSpeed, float heaterTemp) { // Cognet: Deleted this, so that function can be called, without saving the value of the stProduction.
    float getSolarThermalProduction(float tOut, float windSpeed, float heaterTemp, bool saveValue=true) { // Cognet: Added this.
//        stProduction.push_back(0.f); // Cognet: Deleted this.
        float prod = 0.f; // Cognet: Added this.
        if (stPanel!=NULL) {
//            stProduction.back() += stRatio*area*((stPanel->getEta0()*shortWaveIrradiance.back())-(stPanel->getA1()*(heaterTemp-tOut))-(stPanel->getA2()*(heaterTemp-tOut)*(heaterTemp-tOut))); // Cognet: Deleted this.
            prod += stRatio*area*((stPanel->getEta0()*shortWaveIrradiance.back())-(stPanel->getA1()*(heaterTemp-tOut))-(stPanel->getA2()*(heaterTemp-tOut)*(heaterTemp-tOut))); // Cognet: Added this.
        }
        if (pvtPanel!=NULL) {
//            stProduction.back() += pvtRatio*area*pvtPanel->getThermalSurfacePowerDensity(shortWaveIrradiance.back(),tOut, windSpeed, getEnvironmentalTemperature(), heaterTemp); // Cognet: Deleted this.
            prod += pvtRatio*area*pvtPanel->getThermalSurfacePowerDensity(shortWaveIrradiance.back(),tOut, windSpeed, getEnvironmentalTemperature(), heaterTemp); // Cognet: Added this.
        }
        // keeps only positive values
//        if (stProduction.back() < 0.f) stProduction.back() = 0.f; // Cognet: Deleted this.
        if (prod < 0.f) prod = 0.f; // Cognet: Added this.

        if (saveValue) { stProduction.push_back(prod); } // Cognet: Added this.

//        return stProduction.back(); // Cognet: Deleted this.
        return prod; // Cognet: Added this.
    }
    float getSolarThermalProduction(size_t it) { if (stProduction.empty()) return std::numeric_limits<float>::quiet_NaN(); else return stProduction.at(it); }
    void eraseSolarThermalProduction(unsigned int keepValue) { stProduction.erase(stProduction.begin(),stProduction.end()-min(keepValue,(unsigned int)stProduction.size())); }

    // gets the E+ id
    string getEp_id() { return ep_id; }
    void setEp_id(string value) { ep_id = value; }
};

class Wall : public Surface {
    //using Surface::writeXML;

private:

    // shading state of the glazing
    vector<float> lowerShadingState, upperShadingState;

public:

    Wall(unsigned int id, Composite* c, float shortWaveReflectance, float glazingRatio, float glazingGvalue, float glazingUvalue, float glazingOpenableRatio, ostream* logStr=NULL)
    :Surface(id,shortWaveReflectance,glazingGvalue,glazingUvalue,glazingOpenableRatio,glazingRatio,logStr){
        composite = c;
    }
    Wall(const Surface& s):Surface(s) { computeNormalAndArea(); } // Incomplete constructor for DXF reading, do not use without completing the geometry...
    Wall(int id, GENPoint p1, GENPoint p2, float elevation):Surface(id){
        // for wall to be oriented outside, floor vertices (positive orientation) must be taken backwards
        vertices.push_back(p2);
        vertices.push_back(p1);
        GENPoint p3(p1); // above p1 is p3
        p3[2]+=elevation;
        vertices.push_back(p3);
        GENPoint p4(p2); // above p2 is p4
        p4[2]+=elevation;
        vertices.push_back(p4);

        computeNormalAndArea();
    }
    Wall(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr=NULL);
    void clear() {
        lowerShadingState.clear(); upperShadingState.clear();
        Surface::clear();
    }

//    Wall(unsigned int id, WallType* wallType, float shortWaveReflectance, float glazingRatio, float glazingGvalue, float glazingUvalue, float glazingOpenableRatio, ostream* logStr=NULL)
//    :Surface(id,shortWaveReflectance,glazingRatio,glazingGvalue,glazingUvalue,glazingOpenableRatio,logStr),wallType(wallType) {}
//
//    Wall(unsigned int id, WallType* wallType, float shortWaveReflectance, float glazingRatio, float glazingGvalue, float glazingUvalue, float glazingOpenableRatio, float pvRatio, double pmp, double ac, double tref, double tcnoct, double muvoc, double vmp, ostream* logStr=NULL)
//
    // gets the shortwave reflectance of the faÁade (a mix of glazings and wall), we suppose that the glazings are reflecting nothing
    float getShortWaveReflectance() { return shortWaveReflectance*(1.f-glazingRatio); }
    float getShortWaveReflectance_opaque() { return shortWaveReflectance; }

    SType getType(){return WALL;}

    float getNRE() { return area*(1.f-glazingRatio)*getComposite()->getNRE(); }
    float getGWP() { return area*(1.f-glazingRatio)*getComposite()->getGWP(); }
    float getUBP() { return area*(1.f-glazingRatio)*getComposite()->getUBP(); }
    float getWallArea() { return (1.f-glazingRatio)*area; } // only the wall area, without the glazing

    // status of the blinds
    void setLowerShadingState(float state) { lowerShadingState.push_back(state); }
    float getLowerShadingState() { if (lowerShadingState.empty()) return 1.f; else return lowerShadingState.back(); }
    float getLowerShadingState(unsigned int step) { if (lowerShadingState.empty()) return 1.f; else return lowerShadingState.at(step); }
    void eraseLowerShadingState(unsigned int keepValue) { lowerShadingState.erase(lowerShadingState.begin(),lowerShadingState.end()-min(keepValue,(unsigned int)lowerShadingState.size())); }
    void eraseLowerShadingState_back() { lowerShadingState.pop_back(); }
    void setUpperShadingState(float state) { upperShadingState.push_back(state); }
    float getUpperShadingState() { if (upperShadingState.empty()) return 1.f; else return upperShadingState.back(); }
    float getUpperShadingState(unsigned int step) { if (upperShadingState.empty()) return 1.f; else return upperShadingState.at(step); }
    void eraseUpperShadingState(unsigned int keepValue) { upperShadingState.erase(upperShadingState.begin(),upperShadingState.end()-min(keepValue,(unsigned int)upperShadingState.size())); }
    void eraseUpperShadingState_back() { upperShadingState.pop_back(); }

    void writeXML(ofstream& file, string tab=""){
        file << tab << "<Wall id=\"" << id << "\" key=\"" << key << "\" type=\"";
        if(!composite) cout << "ERROR: no composite defined for Wall " << id << endl;
        else file << composite->getId();
        file << "\" ShortWaveReflectance=\"" << shortWaveReflectance << "\" GlazingRatio=\"" << glazingRatio <<
            "\" GlazingGValue=\"" << glazingGvalue << "\" GlazingUValue=\"" << glazingUvalue << "\" OpenableRatio=\"" << glazingOpenableRatio;
        file << "\">" << endl;
        Surface::writeXML(file, tab+"\t");
        file << tab << "</Wall>" << endl;
    }
};

class Floor : public Surface {

public:

    Floor(unsigned int id, Composite* c, ostream* logStr=NULL)
    :Surface(id,0.f,0.f,0.f,0.f,0.f,logStr){
        composite = c;
    }

    Floor(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr=NULL);
    Floor(Surface s):Surface(s){ computeNormalAndArea();} // Incomplete constructor for DXF reading, do not use without completing the geometry...
    void clear() {
        Surface::clear();
    }

    SType getType(){return FLOOR;}
    float getUAground() { return pow((double)(getComposite()->getResistance() + 1.f/3.f),(double)-1.)*area; /* hc_int = 3 from CIBSE guide */ }

    void writeXML(ofstream& file, string tab=""){
        file << tab << "<Floor id=\"" << id << "\" key=\"" << key << "\" type=\"";
        if(!composite) cout << "ERROR: no composite defined for Floor " << id << endl;
        else file << composite->getId();
        file  << "\" ShortWaveReflectance=\"" << shortWaveReflectance << "\" GlazingRatio=\"" << glazingRatio <<
             "\" GlazingGValue=\"" << glazingGvalue << "\" GlazingUValue=\"" << glazingUvalue << "\" OpenableRatio=\"" << glazingOpenableRatio << "\">" << endl;
        Surface::writeXML(file, tab+"\t");
        file << tab << "</Floor>" << endl;
    }

    void writeGML_composedOf(ofstream& file, string tab="") {
        // opaque part
        //file << tab << "\t<energy:ThermalComponent>" << endl;
        file << tab << "<energy:area uom=\"m2\">" << area << "</energy:area>" << endl;
        file << tab << "<energy:construction xlink:href=\"#op_" << getComposite()->getId() << "\"/>" << endl;
        //file << tab << "\t\t<energy:isGroundCoupled>true</energy:isGroundCoupled>" << endl;
        //file << tab << "\t\t<energy:isSunExposed>false</energy:isSunExposed>" << endl;
        //file << tab << "\t</energy:ThermalComponent>" << endl;
    }

};

class Roof : public Surface {

private:

    float shadingState = 0.f;
    MicroWindTurbine *microWindTurbine = NULL;
    GENPoint centroid = GENPoint::Origin();
    float kFactor = 0.f, X=0.f, Y=0.f;
    vector<float> waterEvapotranspiration;

public:

    Roof(unsigned int id, Composite* c, float shortWaveReflectance, float glazingRatio, float glazingGvalue, float glazingUvalue, float glazingOpenableRatio,ostream* logStr=NULL)
        :Surface(id,shortWaveReflectance,glazingGvalue,glazingUvalue,glazingOpenableRatio,glazingRatio,logStr) {
        composite = c;
    }
    Roof(Surface const & s):Surface(s){} // Incomplete constructor for DXF reading, do not use without completing the geometry...
    Roof(Floor const & f, float elevation):Surface(f){ // Incomplete constructor for DXF reading, do not use without completing the geometry...
        id = f.getId()+1;
        reverseOrientation();
        for(unsigned int i=0; i<vertices.size();++i){
            vertices[i][2] += elevation;
        }
        computeNormal();
    }
    Roof(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr=NULL);
    ~Roof();
    void clear() {
        waterEvapotranspiration.clear();
        Surface::clear();
    }

    // gets the shortwave reflectance of the faÁade (a mix of glazings and wall), we suppose that the glazings are reflecting nothing
    float getShortWaveReflectance() { return shortWaveReflectance*(1.f-glazingRatio); }
    float getShortWaveReflectance_opaque() { return shortWaveReflectance; }

    SType getType(){return ROOF;}

    float getKr() { return pow((double)getComposite()->getResistance() + 1.f/3.f,(double)-1.); /* hc_int=3 from CIBSE guide */ }
    float getKappa() { return getKr()/(getKr() + get_hc() + get_hr()); }
    float getRoofArea() { return (1.f-glazingRatio)*area; } // only the roof area, without the glazing

    float getShadingState() { return shadingState; }
    void setShadingState(float state) { shadingState = state; }

    // wind power
    float getWindElectricPower(float windSpeed, float metStationAltitude) {
        // adds the centroid of the roof (z-coordinate) and the metStationAltitude
        if (microWindTurbine) return microWindTurbine->getWindElectricPower(windSpeed, metStationAltitude, centroid[2]);
        else return 0.f;
    }

    // ET model (by GU)
    bool hasET() { return (kFactor > 0.f); }
    float get_Y() { return Y; }
    float get_X() { return X; }
    void set_YX(Climate* pClimate, unsigned int day, unsigned int hour) {
        if (hasET()) {
            Y = kFactor*pClimate->getGamma(day,hour)*(Climate::getSaturatedVapourPressure(getTemperature()) - pClimate->getSaturatedVapourPressureDerivative(getTemperature())*getTemperature() - pClimate->getVapourPressure(day,hour));
            X = kFactor*pClimate->getGamma(day,hour)*pClimate->getSaturatedVapourPressureDerivative(getTemperature()); // Qet = Ket*(Y + X*Tsurface) //Ket=0: no green Ket>0: green area, Ket =1, grass 0.12 m.
        }
        else return;
    }
    float getWaterEvapotranspiration(size_t index) { return waterEvapotranspiration.at(index); }
    void setWaterEvapotranspiration() { if (kFactor > 0.f) waterEvapotranspiration.push_back((Y + X*getTemperature())/694.5f); }
    void eraseWaterEvapotranspiration(unsigned int keepValue) { if (!waterEvapotranspiration.empty()) waterEvapotranspiration.erase(waterEvapotranspiration.begin(),waterEvapotranspiration.end()-min(keepValue,(unsigned int)waterEvapotranspiration.size())); } // removes all elements but the last one

    void writeXML(ofstream& file, string tab=""){
        file << tab << "<Roof id=\"" << id << "\" key=\"" << key << "\" type=\"";
        if(composite) file << composite->getId();
        else throw("No composite defined in Roof id=" + toString(id));
        file << "\" ShortWaveReflectance=\"" << shortWaveReflectance << "\" GlazingRatio=\"" << glazingRatio <<
             "\" GlazingGValue=\"" << glazingGvalue << "\" GlazingUValue=\"" << glazingUvalue << "\" OpenableRatio=\"" << glazingOpenableRatio << "\" ";
        file << "kFactor=\"" << kFactor << "\"";
        file << ">" << endl;
        Surface::writeXML(file, tab+"\t");
        file << tab << "</Roof>" << endl;
    }
};

class Ground : public Surface {

private:
    vector<float> layerTemperature; // contains each layer's temperature
    float C1, G1, G2; // one node equivalent values
    float kFactor = 0.f, X=0.f, Y=0.f;
    vector<float> waterEvapotranspiration;
    bool detailedSimulation = false;

public:

    Ground(unsigned int id, Composite *c, float shortWaveReflectance, ostream* logStr=NULL)
    :Surface(id,shortWaveReflectance,0.f,0.f,0.f,0.f,logStr){
        composite = c;
        initialiseModel();
    }
    Ground(Surface *s, Composite *c):Surface(*s) { composite=c; initialiseModel(); }
    Ground(Surface *s, Composite *c, string key):Surface(*s,key) { composite=c; initialiseModel(); }
    Ground(TiXmlHandle hdl, Composite* c, ostream* pLogStream=NULL);

    friend bool operator>(const Ground& lhs, const Ground& rhs) { return lhs.getId() > rhs.getId(); }

    void initialiseModel(){
        // the temperature of each layer is initialized to 15 celsius
        /// TODO: check this initial condition of 8∞C
        if (composite) layerTemperature.assign(composite->getnLayers(),8.f);
        if (composite) getComposite()->getSimplifiedNode(C1,G1,G2); ///TODO: check if no composite what happens
    }

    SType getType(){return GROUND;}

    float getLayerTemperature(unsigned int i) { return layerTemperature[i]; }
    void setLayerTemperature(unsigned int i, float value) { layerTemperature[i]=value; }

    // gets the simplified version of the ground
    float getC1() { return C1; }
    float getG1() { return G1; }
    float getG2() { return G2; }

    // gets and sets the k factor for the ET model (by GU)
    float getKfactor() { return kFactor; }
    void setKfactor(float value) { kFactor = value; }
    // ET model (by GU)
    bool hasET() { return (kFactor > 0.f); }
    float get_Y() { return Y; }
    float get_X() { return X; }
    void set_YX(Climate* pClimate, unsigned int day, unsigned int hour) {
        if (hasET()) {
            Y = kFactor*pClimate->getGamma(day,hour)*(Climate::getSaturatedVapourPressure(getTemperature()) - pClimate->getSaturatedVapourPressureDerivative(getTemperature())*getTemperature() - pClimate->getVapourPressure(day,hour));
            X = kFactor*pClimate->getGamma(day,hour)*pClimate->getSaturatedVapourPressureDerivative(getTemperature()); // Qet = Ket*(Y + X*Tsurface) //Ket=0: no green Ket>0: green area, Ket =1, grass 0.12 m.
        }
        else return;
    }
    float getWaterEvapotranspiration(size_t index) { return waterEvapotranspiration.at(index); }
    void setWaterEvapotranspiration(float value) { if (kFactor > 0.f) waterEvapotranspiration.push_back(value); }
    void eraseWaterEvapotranspiration(unsigned int keepValue) { if (!waterEvapotranspiration.empty()) waterEvapotranspiration.erase(waterEvapotranspiration.begin(),waterEvapotranspiration.end()-min(keepValue,(unsigned int)waterEvapotranspiration.size())); } // removes all elements but the last one

    // detailed simulation (multi-nodal)
    bool isDetailedSimulation() { return detailedSimulation; }
    void setDetailedSimulation(bool value) { detailedSimulation = value; }

    // gets the double conductance of layer i
    float getG(unsigned int i) { return 2.f*getComposite()->getLayer(i)->getConductance(); }

    // get the conductance between the nodes
    float getk(unsigned int i) {
        // test for range of i, nLayers+1 values
        if (i > (composite->getnLayers())) throw string("Layer index out of range.");
        else if (i == 0) return ((get_hc()+get_hr())*getG(0))/((get_hc()+get_hr())+getG(0));
        else if (i == (composite->getnLayers())) return getG(i-1);
        else return (getG(i-1)*getG(i))/(getG(i-1)+getG(i));
    }

    float getConductance(unsigned int i, unsigned int j) {
        // diagonal elements
        if (i==j) return -getk(i)-getk(i+1);
        // upper diagonal element
        else if (j==(i+1)) return getk(i+1);
        // lower diagonal element
        else if (j==(i-1)) return getk(i);
        // elsewhere nil
        else return 0.f;
    }

    float getCapacitance(unsigned int i, unsigned int j) {
        // diagonal matrix with zeroes elsewhere
        if (i==j) return getComposite()->getCapacitance(i);
        else return 0.f;
    }

    void writeXML(ofstream& file, string tab="") {
        file << tab << "<Ground id=\"" << id << "\" key=\"" << key << "\" ";
        if (composite) file << "type=\"" << composite->getId() << "\" ";
        else throw("No composite defined in Ground id=" + toString(id));
        file << "ShortWaveReflectance=\"" << shortWaveReflectance << "\" ";
        file << "kFactor=\"" << kFactor << "\" ";
        file << "detailedSimulation=\"" << ((detailedSimulation)?("true"):("false")) << "\"";
        file << ">" << endl;
        Surface::writeXML(file, tab+"\t");
        file << tab << "</Ground>" << endl;
    }
};

#endif
