#ifndef BUILDING_H
#define BUILDING_H

#include "zone.h"
#include "plant.h"
#include "util.h"
#include "surface.h"

class District;
class Zone;

// *** Building class, CitySim *** //
// *** jerome.kaempf@epfl.ch   *** //

class Building {

public:
    enum ResultsPositions {DAY_RESULTS_POSITION=0, MONTH_RESULTS_POSITION=365, YEAR_RESULTS_POSITION=365+12};
    enum DiplayableResults {HEATING_DEMAND, COOLING_DEMAND, TEMPERATURE, N_RESULTS};

private:

    // identifies the district it is in, to access global values
    District* pDistrict = nullptr;

    // the ID of the building and its key to a database
    unsigned int id=0;
    string key="";

    // indicates to simulate the building with EnergyPlus
    bool simulateEP = false;
    string fmuFile, tmpPath;

    // the vector containing the building Zones, it should be created by the constructor of building
    vector<Zone*> zones;

    // the matrices representing the building for the thermal model (constant throughout the simulation)
    vector<vector<float>> C; ///< Matrix of the capacitances, (W/K)
    vector<vector<float>> G1; ///< Matrix of the fixed conductances of the building (including the links between the zones), (W/K)

    // for the HVAC model
    bool HVACpresence=false;
    double coileff,coilHTWeff;
    bool evaporativeCooling=false;
    double TminSupply, TmaxSupply, deltaT;

    // to provide the heat needed by the HVAC model and store the tank temperature
    Tank* heatStock = nullptr;
    Tank* dhwStock = nullptr;
    Tank* coldStock = nullptr;
    vector<float> heatStockTemperature, dhwStockT, coldStockTemperature;

    // the energy conversion unit
    EnergyConversionSystem* heatingUnit = nullptr;
    EnergyConversionSystem* coolingUnit = nullptr;

    // whole building consumption on an hourly basis, emptied from time to time to prevent memory stagnation
    vector<float> electricConsumption, fuelConsumption, machinePower, solarPVProduction, solarThermalProduction;

    // save in a Compressed Row Storage sparse matrix format the links between the thermal zones
    vector<pair<float,float> > linksAn;
    vector<unsigned int> linksAj, linksAi;

    // eco-indicators
    float nre = 0.f, gwp = 0.f, ubp = 0.f;

    // blinds parameters
    float blindsLambda = 0.2f; // shape of the sigmoid curve
    float blindsIrradianceCutOff = 100.f; // 100 W/m≤ for the cut-off

    // mrt parameters
    bool mrt = false;
    float mrtEpsilon = 0.97f;

    // Cognet: Start of added code.
    // Variables for computations at each time step. (HS = hot storage, DHW = domestic hot water, CS = cold storage).
    map<string, float> imposedHeatDemand; // Stores values if one wants to impose the heat demanded to the heatingUnit. Format map["d42h4"]=3200.
    double heatingNeeds, coolingNeeds; // Sum over the building's zones' space heating/cooling needs.
    double HS_needs, DHW_needs, CS_needs; // Thermal power needed the tanks (power that will be asked to heating/cooling unit (the unit may not be powerful enough to fulfill the need) ).
    double HS_SolPp, DHW_SolPp, CS_SolPp; // Solar thermal power provided for the tanks.
    double SolTherFracLeft; // Solar thermal power is used by different applications : HS tank, DHW tank or feed-in substation. Since the temperatures change, the power available changes. When one uses a portion of the available power (only 40%), then it sets this variable to the fraction that is still available (here 60%).
    double HS_Pp, DHW_Pp, CS_Pp; // Thermal power provided for the tanks (power that the heating/cooling unit will provide).
    double VdotUsed; // Domestic hot water consumption summed over zones [m^3/s].
    float Tamb; // "Ambiant teperature", where the water tanks are located.
    // Cognet: End of added code.
    // Added by Max
    float heatingDemandUnsatisfied = 0.f;
    float coolingDemandUnsatisfied = 0.f;
    bool hasSolarThermal = false;
public:

    ostream logStream;

    // the constructor of the building, which reads in the XML file and the destructor that removes the Zones
    Building(TiXmlHandle hdl, District* pDistrict);
    Building(vector<Wall*> walls, vector<Roof*> roofs, vector<Floor*> floors, vector<Surface*> surfaces, District* pDistrict);
    ~Building();
    void clear() {
        heatStockTemperature.clear();
        dhwStockT.clear();
        coldStockTemperature.clear();
        electricConsumption.clear();
        fuelConsumption.clear();
        machinePower.clear();
        solarThermalProduction.clear();
        for (vector<Zone*>::iterator it=zones.begin();it!=zones.end();++it) (*it)->clear();
    }
    void update();

    void computeVolume();

    void addZone(Zone* z){ zones.push_back(z);}
    void removeZone(int pos){ zones[pos]=nullptr;} // to avoid deleting Zones transfered in another Building

    XmlScene* getScene();


    // returns the pointer to the creating district
    District* getDistrict() { return pDistrict; }
	unsigned int getId() { return id; }
	string getKey() { return key; }
    void setKey(string k){key = k; }

	// returns if the building is simulated with EP
	bool isEP() { return simulateEP; }
	string getFMUFile() { return fmuFile; }
	string getTMPPath() { return tmpPath; }

    // returns the links matrices
    pair<float,float> getLinksAn(unsigned int index) { return linksAn[index]; }
    unsigned int getLinksAj(unsigned int index) { return linksAj[index]; }
    unsigned int getLinksAi(unsigned int index) { if (index >= linksAi.size()) return linksAn.size(); else return linksAi[index]; }
    unsigned int getnLinksAi() { return linksAi.size(); }

    // to get the values parameters of the building
    float getNinf() {
        float volumePerHour = 0.f; // in m3/h
        for (size_t i=0;i<zones.size();++i) volumePerHour += zones.at(i)->getNinf()*zones.at(i)->getVi();
        return volumePerHour/getVolume();
    }
    void setNinf(float Ninf) { for (size_t i=0;i<zones.size();++i) zones.at(i)->setNinf(Ninf); }

    // for the thermal model
    Zone* getZone(unsigned int i) { return zones[i]; }
    Zone* getZoneFromId(unsigned int id) { return zones.at(getZoneIndexFromId(id)); }
    size_t getnZones() { return zones.size(); }
    vector<Zone*>* getZones() { return &zones; }
    unsigned int getnNodes();
    unsigned int getZonenNodes(unsigned int i);

    // gets the zone Index from the zone Id
    unsigned int getZoneIndexFromId(unsigned int id);

    // uses the zone capacitance
    double getZoneC(unsigned int zoneIndex, unsigned int nodeIndex);
    void addZoneC(unsigned int zoneIndex, float value);

    // gets the last element in the temperature vector
    double getZoneT(unsigned int i, unsigned int j);
    // sets the temperature in the zone i for node j
    void setZoneT(unsigned int i, unsigned int j, double value);

    // gets the matrices and operate on them - JK - to be put un ThermalModel class
    double getC(unsigned int i, unsigned int j) { return C[i][j]; }
    double getG1(unsigned int i, unsigned int j) { return G1[i][j]; }
    int getMatrixPosition(int zoneNumber) { int sum = 0; for (int j=0; j<zoneNumber; ++j) { sum += zones[j]->getnNodes(); } return sum; }
    int getSubMatrixPosition(int zoneNumber) { int sum = 0; for (int j=0; j<zoneNumber; ++j) { sum += zones[j]->getnNodes()-1; } return sum; }
    void loadLinkSparse(vector<vector<float>>& G) {
        // loop on the elements in the columns
        for (unsigned int i=0;i<getnLinksAi();++i) {
            // loop on the elements of the rows
            float sum = 0.f;
            for (unsigned int index=getLinksAi(i);index<getLinksAi(i+1);++index) {
                G[getMatrixPosition(i)][getMatrixPosition(getZoneIndexFromId(getLinksAj(index)))] += getLinksAn(index).first;
                sum -= getLinksAn(index).first;
            }
            G[getMatrixPosition(i)][getMatrixPosition(i)] += sum;
        }
    }

    // get the heating and cooling in Wh
    float getHeating(unsigned int step) { float watthour = 0.f;   for (size_t i=0;i<zones.size();++i)  { watthour += zones.at(i)->getHeating(step); } return watthour; }
    double getHeating(unsigned int day, unsigned int hour) { double watthour = 0.;   for (unsigned int i=0;i<zones.size();++i)  { watthour += zones[i]->getHeating(day, hour); } return watthour; }
    double getCooling(unsigned int step) { double watthour = 0.;   for (unsigned int i=0;i<zones.size();++i)  { watthour += zones[i]->getCooling(step); } return watthour; }
    double getCooling(unsigned int day, unsigned int hour) { double watthour = 0.;   for (unsigned int i=0;i<zones.size();++i)  { watthour += zones[i]->getCooling(day, hour); } return watthour; }
    // occupants number is the number of people maximum in the building
    float getOccupantsNumber() { float number = 0.f; for (size_t i=0; i<zones.size();++i) { number += zones[i]->getOccupantsNumber(); } return number; }
    // occupants count is the actual number of people
    float getOccupantsCount()  { float number = 0.f; for (size_t i=0; i<zones.size();++i) { number += zones[i]->getOccupantsCount(); } return number; }

    // method for HVAC
    bool getHVACpresence() { return HVACpresence; }
    bool getEvaporativeCooling() { return evaporativeCooling; }
    double getCoilEfficiency() { return coileff; }
    double getCoilHydroThermalWheelEfficiency() { return coilHTWeff; }
    double getTminSupply() { return TminSupply; }
    double getTmaxSupply() { return TmaxSupply; }
    double getDeltaT() { return deltaT; }

    // method for the heat and cold stocks
    Tank* getHeatStock() { return heatStock; }
    float getHeatStockTemperature() { if (heatStockTemperature.empty()) return heatStock->getTmin(); else return heatStockTemperature.back(); }
    float getHeatStockTemperature(unsigned int step) { if (heatStockTemperature.empty()) return heatStock->getTmin(); else return heatStockTemperature.at(step); }
    void setHeatStockTemperature(double celsius) { heatStockTemperature.push_back(celsius); }
    void eraseHeatStockTemperature(unsigned int keepValue) { heatStockTemperature.erase(heatStockTemperature.begin(),heatStockTemperature.end()-min(keepValue,(unsigned int)heatStockTemperature.size())); }
    void eraseHeatStockTemperature_back() { heatStockTemperature.pop_back(); }
    size_t sizeHeatStockTemperature() { return heatStockTemperature.size(); }

    Tank* getDHWHeatStock() { return dhwStock; }
    float getDHWStockT() { if (dhwStockT.empty()) return dhwStock->getTmin(); else return dhwStockT.back(); }
    float getDHWStockT(unsigned int step) { if (dhwStockT.empty()) return dhwStock->getTmin(); else return dhwStockT.at(step); }
    void setDHWStockT(double celsius) { dhwStockT.push_back(celsius); }
    void eraseDHWStockT(unsigned int keepValue) { dhwStockT.erase(dhwStockT.begin(),dhwStockT.end()-min(keepValue,(unsigned int)dhwStockT.size())); }
    void eraseDHWStockT_back() { dhwStockT.pop_back(); }

    Tank* getColdStock() { return coldStock; }
    float getColdStockTemperature() { if (coldStockTemperature.empty()) return coldStock->getTmax(); else return coldStockTemperature.back(); }
    float getColdStockTemperature(unsigned int step) { if (coldStockTemperature.empty()) return coldStock->getTmax(); else return coldStockTemperature.at(step); }
    void setColdStockTemperature(double celsius) { coldStockTemperature.push_back(celsius); }
    void eraseColdStockTemperature(unsigned int keepValue) { coldStockTemperature.erase(coldStockTemperature.begin(),coldStockTemperature.end()-min(keepValue,(unsigned int)coldStockTemperature.size())); }
    void eraseColdStockTemperature_back() { coldStockTemperature.pop_back(); }

    // method to get the energyUnit
    EnergyConversionSystem* getHeatingUnit() { return heatingUnit; }
    EnergyConversionSystem* getCoolingUnit() { return coolingUnit; }

    // Electric Consumption
    void setElectricConsumption(float joules) { electricConsumption.push_back(joules); }
    void addElectricConsumption(float joules) { electricConsumption.back() += joules; }
    float getElectricConsumption(unsigned int step) { return electricConsumption.at(step); }
    float getElectricConsumption() { return electricConsumption.back(); }
    float getElectricConsumption(unsigned int day, unsigned int hour) { return electricConsumption.at((day-1)*24 + hour -1); }
    void eraseElectricConsumption() { electricConsumption.erase(electricConsumption.begin(),electricConsumption.end()); }
    void eraseElectricConsumption_back() { electricConsumption.pop_back(); }
    float getTotalElectricConsumption() { return accumulate(electricConsumption.begin(),electricConsumption.end(),0.f); }
    // Solar PV
    void setSolarPVProduction(float joules) { solarPVProduction.push_back(joules); }
    void addSolarPVProduction(float joules) { solarPVProduction.back() += joules; }
    float getSolarPVProduction(unsigned int step) { return solarPVProduction.at(step); }
    float getSolarPVProduction() { return solarPVProduction.back(); }
    float getSolarPVProduction(unsigned int day, unsigned int hour) { return solarPVProduction.at((day-1)*24 + hour -1); }
    void eraseSolarPVProduction() { solarPVProduction.erase(solarPVProduction.begin(),solarPVProduction.end()); }
    void eraseSolarPVProduction_back() { solarPVProduction.pop_back(); }
    float getTotalSolarPVProduction() { return accumulate(solarPVProduction.begin(),solarPVProduction.end(),0.f); }
    // Solar Thermal
    void setSolarThermalProduction(float joules) { solarThermalProduction.push_back(joules); }
    float getSolarThermalProduction(unsigned int step) { return solarThermalProduction.at(step); }
    float getSolarThermalProduction() { return solarThermalProduction.back(); }
    float getSolarThermalProduction(unsigned int day, unsigned int hour) { return solarThermalProduction.at((day-1)*24 + hour -1); }
    void eraseSolarThermalProduction() { solarThermalProduction.erase(solarThermalProduction.begin(),solarThermalProduction.end()); }
    void eraseSolarThermalProduction_back() { solarThermalProduction.pop_back(); }
    float getTotalSolarThermalProduction() { return accumulate(solarThermalProduction.begin(),solarThermalProduction.end(),0.f); }
    bool getHasSolarThermal(); // Added by Max
    void computeHasSolarThermal();//Added by Max
    // Fuel Consumption
    void setFuelConsumption(float joules) { fuelConsumption.push_back(joules); }
    void addFuelConsumption(float joules) { fuelConsumption.back() += joules; }
    float getFuelConsumption(unsigned int step) { return fuelConsumption.at(step); }
    float getFuelConsumption() { return fuelConsumption.back(); }
    float getFuelConsumption(unsigned int day, unsigned int hour) { return fuelConsumption.at((day-1)*24 + hour -1); }
    void eraseFuelConsumption() { fuelConsumption.erase(fuelConsumption.begin(),fuelConsumption.end()); }
    void eraseFuelConsumption_back() { fuelConsumption.pop_back(); }
    float getTotalFuelConsumption() { return accumulate(fuelConsumption.begin(),fuelConsumption.end(), 0.f); }
    size_t sizeTotalFuelConsumption() { return fuelConsumption.size(); }

    // Machine Power
    void setMachinePower(double watts) { machinePower.push_back(watts); }
    void addMachinePower(double watts) { machinePower.back() += watts; }
    double getMachinePower(unsigned int step) { return machinePower.at(step); }
    double getMachinePower(unsigned int day, unsigned int hour);
    double getMachinePower() { return machinePower.back(); }
    void eraseMachinePower(unsigned int keepValue) { machinePower.erase(machinePower.begin(),machinePower.end()-min(keepValue,(unsigned int)machinePower.size())); } // DP: modified to keep MachinePower results (MEU webservice)
    void eraseMachinePower_back() { machinePower.pop_back(); }
    size_t sizeMachinePower() { return machinePower.size(); }

    // get the eco-indicators
    float getNRE() { return nre; }
    float getGWP() { return gwp; }
    float getUBP() { return ubp; }
    //-Floor
    float getFloorArea()   { float area = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) area += zones[zoneIndex]->getFloorArea(); return area; }
    float getFloorGWP()    { float gwp = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) gwp += zones[zoneIndex]->getFloorGWP(); return gwp; }
    float getFloorNRE()    { float nre = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) nre += zones[zoneIndex]->getFloorNRE(); return nre; }
    float getFloorUBP()    { float ubp = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) ubp += zones[zoneIndex]->getFloorUBP(); return ubp; }
    //-Roof
    float getRoofArea()    { float area = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) area += zones[zoneIndex]->getRoofArea(); return area; }
    float getRoofPVArea()  { float area = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) area += zones[zoneIndex]->getRoofPVArea(); return area; }
    float getRoofGWP()     { float gwp = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) gwp += zones[zoneIndex]->getRoofGWP(); return gwp; }
    float getRoofNRE()     { float nre = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) nre += zones[zoneIndex]->getRoofNRE(); return nre; }
    float getRoofUBP()     { float ubp = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) ubp += zones[zoneIndex]->getRoofUBP(); return ubp; }
    //-Wall
    float getWallArea()    { float area = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) area += zones[zoneIndex]->getWallArea(); return area; }
    float getWallPVArea()  { float area = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) area += zones[zoneIndex]->getWallPVArea(); return area; }
    float getWallGWP()     { float gwp = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) gwp += zones[zoneIndex]->getWallGWP(); return gwp; }
    float getWallNRE()     { float nre = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) nre += zones[zoneIndex]->getWallNRE(); return nre; }
    float getWallUBP()     { float ubp = 0.f;  for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) ubp += zones[zoneIndex]->getWallUBP(); return ubp; }
    //-Windows
    float getWindowsArea() { float area = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) area += zones[zoneIndex]->getSwi(); return area; }

    // Heated volume (m3)
    float getVolume() { float v = 0.f; for (size_t zoneIndex=0;zoneIndex<zones.size();++zoneIndex) v += zones[zoneIndex]->getVi(); return v; }

    // gets the blinds parameters
    float getBlindsLambda() { return blindsLambda; }
    void setBlindsLambda(float const& l) { blindsLambda = l;}
    float getBlindsIrradianceCutOff() { return blindsIrradianceCutOff; }
    void setBlindsIrradianceCutOff(float const& c) { blindsIrradianceCutOff = c;}

    // gets the mrt parameters
    bool isMRT() { return mrt; }
    float getMRT_Epsilon() { return mrtEpsilon; }

    void writeXML(ofstream& file, string tab="");
    void writeGML(ofstream& file, string tab="", const vector<double>& origin={0.,0.,0.});

    // method to determine the shading state of all facades and roofs
    void deterministicShadingAction(/*unsigned int day*/);

    // Cognet: Start of added code.
    // Variables for computations at each time step.
    void setHeatingNeeds(double hN) { heatingNeeds = hN; }
    double getHeatingNeeds() { return heatingNeeds; }
    void setCoolingNeeds(double cN) { coolingNeeds = cN; }
    double getCoolingNeeds() { return coolingNeeds; }

    void setHS_needs(double needs) { HS_needs = needs; }
    double getHS_needs() { return HS_needs; }
    void setDHW_needs(double needs) { DHW_needs= needs; }
    double getDHW_needs() { return DHW_needs; }
    void setCS_needs(double needs) { CS_needs= needs; }
    double getCS_needs() { return CS_needs; }

    void setHS_Pp(double pp) { HS_Pp = pp; }
    double getHS_Pp() { return HS_Pp; }
    void setDHW_Pp(double pp) { DHW_Pp = pp; }
    double getDHW_Pp() { return DHW_Pp; }
    void setCS_Pp(double pp) { CS_Pp = pp; }
    double getCS_Pp() { return CS_Pp; }

    void setHS_SolPp(double Solpp) { HS_SolPp = Solpp; }
    double getHS_SolPp() { return HS_SolPp; }
    void setDHW_SolPp(double Solpp) { DHW_SolPp = Solpp; }
    double getDHW_SolPp() { return DHW_SolPp; }
    void setCS_SolPp(double Solpp) { CS_SolPp = Solpp; }
    double getCS_SolPp() { return CS_SolPp; }

    void setSolTherFracLeft(double frac) { SolTherFracLeft = frac; }
    void multSolTherFracLeft(double d) { SolTherFracLeft *= d; }
    double getSolTherFracLeft() { return SolTherFracLeft; }

    void setVdotUsed(double vdot) { VdotUsed = vdot; }
    double getVdotUsed() { return VdotUsed; }

    void setTamb(float tamb) { Tamb = tamb; }
    float getTamb() { return Tamb; }

    bool hasImposedHeatDemand(unsigned int day, unsigned int hour, float &retValue);
    bool hasImposedHeatDemand(unsigned int day, unsigned int hour);
    // Cognet: End of added code.

    //Added by Max
    float getHeatingDemandUnsatisfied() {return heatingDemandUnsatisfied;}
    void setHeatingDemandUnsatisfied(float MissingThermalPower) {heatingDemandUnsatisfied = MissingThermalPower;}
    float getCoolingDemandUnsatisfied() {return coolingDemandUnsatisfied;}
    void setCoolingDemandUnsatisfied(float MissingThermalPower) {coolingDemandUnsatisfied = MissingThermalPower;}
};

class Tree {

private:

    unsigned int id = numeric_limits<unsigned int>::quiet_NaN();
    string name = "";
    string key = "";
    // geometrical parameters
    unsigned int layers = 1;
    float layersDistance = 1.f; // in meters
    // physical parameters
    float leafWidth = 0.1f; // in meters
    bool deciduous = false;
    string leafClass = "C3";

    vector<Surface*> leaves;
    vector<Surface*> subLeaves;
    vector<Surface*> trunc;

public:

    Tree(ostream* pLogStream=NULL):logStream(std::cout.rdbuf()) { associate(pLogStream,logStream); }
    Tree(TiXmlHandle hdl, ostream* pLogStream=NULL);
    void clear() {
        for (vector<Surface*>::iterator it=leaves.begin();it!=leaves.end();++it) (*it)->clear();
        for (vector<Surface*>::iterator it=subLeaves.begin();it!=subLeaves.end();++it) (*it)->clear();
        for (vector<Surface*>::iterator it=trunc.begin();it!=trunc.end();++it) (*it)->clear();
    }

    ostream logStream;

    vector<Surface*>* getLeaves() { return &leaves; }
    vector<Surface*>* getSubLeaves() { return &subLeaves; }
    vector<Surface*>* getTrunc() { return &trunc; }

    void setId(unsigned int id) { this->id = id; }
    unsigned int getId() { return id; }
    void setName(string name) { this->name = name; }
    string getName() { return name; }
    void setKey(string key) { this->key = key; }
    string getKey() { return key; }
    void setLeafAreaIndex(unsigned int index) { layers = index; }
    unsigned int getLeafAreaIndex() { return layers; }
    void setLeafDistance(float d) { layersDistance = d; }
    float getLeafDistance() { return layersDistance; }
    void setLeafWidth(float w) { leafWidth = w; }
    float getLeafWidth() { return leafWidth; }
    void setDeciduous(bool value) { deciduous = value; }
    bool getDeciduous() { return deciduous; }
    void setLeafClass(string value) { leafClass = value; }
    string getLeafClass() { return leafClass; }

    void addSurface(Surface* surface) {
        surface->computeNormalAndArea();
        if (surface->getArea() <= 0.f) {
            logStream << "Surface id: " << surface->getId() << " has a too small surface, removing it." << endl;
            delete surface;
        }
        else {
            // give attributes for the plant
            surface->setLongWaveEmissivity(0.95f);
            surface->setShortWaveReflectance(0.3f);
            if (surface->getAltitude() > 80.f && surface->getAltitude() <= 100.f) {
                leaves.push_back(surface);
                subLeaves.push_back(new Surface(*surface));
                subLeaves.back()->reverseOrientation();
                for (unsigned int index = 1; index < layers; ++index) {
                    subLeaves.push_back(new Surface(*surface));
                    subLeaves.back()->translate(GENPoint::Cartesian(0.f,0.f,-static_cast<float>(index)*layersDistance));
                    subLeaves.push_back(new Surface(*subLeaves.back()));
                    subLeaves.back()->reverseOrientation();
                }
            }
            else trunc.push_back(surface);
        }
    }

    size_t getnSurfaces() { return leaves.size() + subLeaves.size() + trunc.size(); }
    Surface* getSurface(size_t i) {
        if (i<leaves.size()) return leaves.at(i);
        else if (i<leaves.size()+subLeaves.size()) return subLeaves.at(i-leaves.size());
        else return trunc.at(i-leaves.size()-subLeaves.size());
    }

    void writeXML(ofstream& file, string tab="");

};

#endif
