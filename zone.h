#ifndef ZONE_H
#define ZONE_H

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>

#include "tinyxml.h"
#include "surface.h"
#include "occupants.h"

using namespace std;

class Building;

// *** Zone class, CitySim   *** //
// *** jerome.kaempf@epfl.ch *** //

class Zone {

private:

    unsigned int id;
    //Building* pBuilding; // -> protected

    // heating and cooling temperatures
    float Tmin = 20.f, Tmax = 26.f;
    unsigned int* Tprofile = nullptr;

    // Ta Foreseen
    double TaForeseen = 0.f;

    // for the HVAC
    vector<double> HVACHeat, HVACCool, HVACReheat, HVACHumidification, HVACEvaporation, HVACMassFlowRate;
    vector<double> moistureContent;

    // available from the EnergyUnit
    vector<double> HVACHeatAvailable, HVACCoolAvailable, HVACReheatAvailable, HVACHumidificationAvailable, HVACEvaporationAvailable, HVACMassFlowRateAvailable;

    // for the stochastic model
    vector<float> lumint;

    unsigned int nightVentilationBegin = 0, nightVentilationEnd = 0; //!< Night ventilation timings during the day

//  vector<double> windowAirExchangeRate;
//
//  vector<int> vPresence;
//  double mu; //parametre of mobility
//  int nbPersMax, dayOff;
//
//  bool clim, climWE, heatWE;
//  int climStartH, climStopH, climStartM, climStopM;
//  int heatStartH, heatStopH, heatStartM, heatStopM;
//
//  double  insideWallArea; //Area of each rad surface from the zone.

    string ep_id;

protected:

    // number of nodes (1,2 or 3)
    unsigned int nNodes=0;

    // pointer to the creator building
    Building* pBuilding;

    // vector of pointers to walls and roofs that belong to the zone
    vector<Wall*> walls;
    vector<Roof*> roofs;
    vector<Surface*> surfaces; // this is for shading areas
    vector<Floor*> floors; // this is for the calculation of floor areas +

    // air node temperatures
    vector<double> Ta; //!< temperature of the air node (celsius) every hour
    vector<double> TaExpl; //!< temperature of the air node (celsius) every 5 mins
    // note Tos will be stocked in the surface itself

    // variables needed by the model itself
    vector<float> Ke; //!< external conductance between the wall and the air nodes (W/K) - convective only
    float Hr; //!< external conductance for the longwave losses (W/K) - radiative only
    double Ki = 0.f; //!< internal conductance between the wall and the air nodes (W/K)
    float Kground = 0.f; //!< total conductance to the ground (W/K)
    float Kpsi = 0.f; //!< total conductance to the air due to linear and ponctual thermal losses (W/K) (by default no thermal bridge)
    float Ww = 1.0f;  // sun fraction on the walls
    float Wa = 0.0f;  // sun fraction in the air

    float Sro = 0.f; //!< Total roofs surface (m²)
    double Swa = 0.; //!< Total walls surface (m²)
    double Swi = 0.; //!< Total windows surface (m²)
    double SwiO = 0.; //!< Total openable window surface (m²)
    float Kwindow = 0.f, Kroof = 0.f;
    float Ninf = 0.1f; //!< Infiltration rate (h-1)
    float Nvent = 0.; // default no ventilation rate
    vector<float> VdotVent;
    float Vi; //!< Internal air volume (m³)
    float Cpi = 1005.f;  //!< Internal air specific heat (J/(kg*K)), initialise with air Cp inside the zone @ 300 K
    float rhoi; //!< Internal air density (kg/m³)
    float Ci; //!< capacitance of the air node (J/K)
    float Lr = 0.f, Lc= 0.f; //!< radiative and convective internal heat gains (W)

    vector<double> heating,cooling; // energy for heating and cooling to get Tmin and Tmax (in Wh)
    vector<float> Qi; // internal gains (in Wh)
    vector<double> Qs; // sensible heat that needs to be given to the room (in Wh)
    double Qs_heating = 0., Qs_cooling = 0.; // total sensible heating and cooling given to the room

    // from the explicit model due to stochastic calculation
    vector<float> windowState; // fully closed = 0.f and fully open 1.f
    vector<float> lightsState;
    bool groundFloor;

    // occupants in the Zone (& default values)
    bool occupantsStochastic = false;
    float occupantsNumber = 0.f;
    float occupantsCount = 0.f;
    float occupantsSensibleHeat = 90.f, occupantsSensibleHeatRadiantFraction = 0.6f, occupantsLatentHeat = 0.f;
    unsigned int activityType = numeric_limits<unsigned int>::signaling_NaN();
    YearProfile* occupantsYearProfile = &OccupancyProfiles::emptyYear;
    YearProfile* dhwYearProfile = &DHWProfiles::emptyYear;
    //unsigned int occupantsYearProfileID;
    Occupants *pOccupants=NULL; // Note: used only with stochastic occupants

  // à vérifier l'utilité!!
  double TminSupply = 0.f, TmaxSupply = 0.f, Deltat = 0.f;    //Temperature min max for the HVAC and loss in temperature


public:

    ostream logStream;

    Zone(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants);
    Zone(unsigned int id, Building* pBuilding, bool groundFloor):id(id),pBuilding(pBuilding),groundFloor(groundFloor),logStream(std::cout.rdbuf()){}// Incomplete constructor for DXF reading, do not use without completing the geometry...
    virtual void addSurface(Surface* s){if(s){}} // for dxf reader, does nothing except for Zone4N
    virtual ~Zone() {
        //NB: the wall and roof pointers are deleted by the GENHandle<Surface>
        //for (vector<Wall*>::iterator it=walls.begin();it!=walls.end();it++) delete *it;
        //for (vector<Roof*>::iterator it=roofs.begin();it!=roofs.end();it++) delete *it;
        if (pOccupants!=NULL) delete pOccupants;
    }
    virtual void clear() {
        HVACHeat.clear();
        HVACCool.clear();
        HVACReheat.clear();
        HVACHumidification.clear();
        HVACEvaporation.clear();
        HVACMassFlowRate.clear();
        moistureContent.clear();
        HVACHeatAvailable.clear();
        HVACCoolAvailable.clear();
        HVACReheatAvailable.clear();
        HVACHumidificationAvailable.clear();
        HVACEvaporationAvailable.clear();
        HVACMassFlowRateAvailable.clear();
        lumint.clear();
        for (vector<Wall*>::iterator it=walls.begin();it!=walls.end();++it) (*it)->clear();
        for (vector<Roof*>::iterator it=roofs.begin();it!=roofs.end();++it) (*it)->clear();
        for (vector<Surface*>::iterator it=surfaces.begin();it!=surfaces.end();++it) (*it)->clear();
        for (vector<Floor*>::iterator it=floors.begin();it!=floors.end();++it) (*it)->clear();
        Ta.clear(); TaExpl.clear();
        Ke.clear(); VdotVent.clear();
        heating.clear(); cooling.clear();
        Qi.clear(); Qs.clear();
        windowState.clear(); lightsState.clear();
    }
    virtual void update(bool constructor=false);

    vector<Surface*> getAllSurfaces();
    void computeVolume();
    void writeXML(ofstream& file, string tab);

    // gets the parent
    Building* getpBuilding() { return pBuilding; }

    // gets and sets the min and max temperatures
    float getTmin(unsigned int day=1, unsigned int hour=1);
    void setTmin(float t) { Tmin=t; }
    float getTmax(unsigned int day=1, unsigned int hour=1);
    void setTmax(float t) { Tmax=t; }
    void setTprofile(unsigned int value) { if (Tprofile) delete Tprofile; Tprofile = new unsigned int(value); }

    unsigned int getId() { return id; }
    unsigned int getnNodes() { return nNodes; }

    size_t getnWalls() { return walls.size(); }
    Wall* getWall(unsigned int i) { return walls[i]; }
    float getWallArea() { float area = 0.f; for (size_t index=0;index<walls.size();++index) area += walls[index]->getWallArea(); return area; }
    float getWallPVArea() { float area = 0.f; for (size_t index=0;index<walls.size();++index) area += walls[index]->getPVRatio()*walls[index]->getArea(); return area; }
    float getWallGWP()  { float gwp = 0.f; for (size_t index=0;index<walls.size();++index) gwp += walls[index]->getWallArea()*walls[index]->getComposite()->getGWP(); return gwp; }
    float getWallNRE()  { float nre = 0.f; for (size_t index=0;index<walls.size();++index) nre += walls[index]->getWallArea()*walls[index]->getComposite()->getNRE(); return nre; }
    float getWallUBP()  { float ubp = 0.f; for (size_t index=0;index<walls.size();++index) ubp += walls[index]->getWallArea()*walls[index]->getComposite()->getUBP(); return ubp; }
    vector<Wall*>* getWalls() { return &walls; }

    size_t getnRoofs() { return roofs.size(); }
    Roof* getRoof(unsigned int i) { return roofs[i]; }
    vector<Roof*>* getRoofs() { return &roofs; }
    float getRoofArea() { float area = 0.f; for (size_t index=0;index<roofs.size();++index) area += roofs[index]->getRoofArea(); return area; }
    float getRoofPVArea() { float area = 0.f; for (size_t index=0;index<roofs.size();++index) area += roofs[index]->getPVRatio()*roofs[index]->getArea(); return area; }
    float getRoofGWP()  { float gwp = 0.f;  for (size_t index=0;index<roofs.size();++index) gwp += roofs[index]->getRoofArea()*roofs[index]->getComposite()->getGWP(); return gwp; }
    float getRoofNRE()  { float nre = 0.f;  for (size_t index=0;index<roofs.size();++index) nre += roofs[index]->getRoofArea()*roofs[index]->getComposite()->getNRE(); return nre; }
    float getRoofUBP()  { float ubp = 0.f;  for (size_t index=0;index<roofs.size();++index) ubp += roofs[index]->getRoofArea()*roofs[index]->getComposite()->getUBP(); return ubp; }

    size_t getnSurfaces() { return surfaces.size(); }
    Surface* getSurface(unsigned int i) { return surfaces[i]; }
    vector<Surface*>* getSurfaces() { return &surfaces; }

    size_t getnFloors() { return floors.size(); }
    Floor* getFloor(unsigned int i) { return floors[i]; }
    vector<Floor*>* getFloors() { return &floors; }
    float getFloorArea() { float area = 0.f; for (size_t index=0;index<floors.size();++index) area += floors[index]->getArea(); return area; }
    float getFloorGWP()  { float gwp = 0.f;  for (size_t index=0;index<floors.size();++index) gwp += floors[index]->getArea()*floors[index]->getComposite()->getGWP(); return gwp; }
    float getFloorNRE()  { float nre = 0.f;  for (size_t index=0;index<floors.size();++index) nre += floors[index]->getArea()*floors[index]->getComposite()->getNRE(); return nre; }
    float getFloorUBP()  { float ubp = 0.f;  for (size_t index=0;index<floors.size();++index) ubp += floors[index]->getArea()*floors[index]->getComposite()->getUBP(); return ubp; }

    float getVi() { return Vi; }
    void setVi(float v) { Vi = v; }
    void setAirDensity(float value) { rhoi = value; Ci = rhoi*Cpi*Vi; /* air specific mass in kg/m3 - defines the internal air capacitance Ci */ }

    // air coupling with negligible capacitance
    void setKwindow(float wattsPerKelvin) { Kwindow = wattsPerKelvin; }
    void setKroof(float wattsPerKelvin) { Kroof = wattsPerKelvin; }
    virtual float getUA() {
        //cout << "Zone::getUA Kwindow="<< Kwindow<<" Kroof="<< Kroof<<" rhoi="<<rhoi<<" Cpi=" << Cpi <<" Vi=" << Vi <<" Ninf=" << Ninf <<" Nvent=" << Nvent << " Kpsi=" << Kpsi<<endl;
        return (Kwindow + Kroof + rhoi*Cpi*Vi*(Ninf/3600.f) + rhoi*Cpi*Vi*(Nvent/3600.f) + Kpsi); /* in W/K */ }

    // ground coupling with negligible capacitance
    float getUAground() { float value = 0.f; for (size_t index=0;index<floors.size();++index) value += floors[index]->getUAground(); return value; /* in W/K */ }

    float getSwa() { return Swa; } // the total surface area in contact with the exterior
    float getSwi() { return Swi; }
    float getSwiO() { return SwiO; }

    void setKe(double convectiveKe) { Ke.push_back(convectiveKe); }
    float getKe() { if (Ke.empty()) return 0.f; else return Ke.back(); }
    float getKe(unsigned int step) { if (Ke.empty()) return 0.f; else return Ke.at(step); }
    void eraseKe(unsigned int keepValue) { Ke.erase(Ke.begin(),Ke.end()-min(keepValue,(unsigned int)Ke.size())); }

    void setHr(float wattsPerKelvin) { Hr = wattsPerKelvin; }
    float getHr() { return Hr; }

    float getKi() { return Ki; }
    void setKi(double convectiveKi) { Ki = convectiveKi; }

    float getWw() { return Ww; }
    float getWa() { return Wa; }

    float getNinf() { return Ninf; }
    void setNinf(float Ninf) { this->Ninf = Ninf; }
    void setNvent(float value) { Nvent = value; }
    float getNvent() { return Nvent; }

    void setVdotVent(float value) { Nvent = value*3600.f/Vi; VdotVent.push_back(value); }
    float getVdotVent() { if (VdotVent.empty()) return 0.f; else return VdotVent.back(); }
    float getVdotVent(unsigned int step) { if (VdotVent.empty()) return 0.f; else return VdotVent.at(step); }
    void eraseVdotVent(unsigned int keepValue) { VdotVent.erase(VdotVent.begin(),VdotVent.end()-min(keepValue,(unsigned int)VdotVent.size())); }
    void eraseVdotVent_back() { VdotVent.pop_back(); }

    unsigned int getNightVentilationBegin() { return nightVentilationBegin; }
    void setNightVentilationBegin(unsigned int value) { nightVentilationBegin = value; }
    unsigned int getNightVentilationEnd() { return nightVentilationEnd; }
    void setNightVentilationEnd(unsigned int value) { nightVentilationEnd = value; }

    float getRadiativeInternalHeatGains() { return Lr; }
    float getConvectiveInternalHeatGains() { return Lc; }

    void setKpsi(float wattsPerKelvin) { Kpsi = wattsPerKelvin; }
    float getKpsi() { return Kpsi; }

    float getQsun1();
    float getQsun2();
    float getQsun3();

    unsigned int getnTa() { return Ta.size(); }
    double getTa(unsigned int it) { return Ta.at(it); }
    double getTa(unsigned int day, unsigned int hour);
    double getTa() { if (Ta.empty()) return 15.0; else return Ta.back(); }
    void setTa(double value) { Ta.push_back(value); }
    void setTa(double value, unsigned int step) { Ta[step] = value; }
    virtual void eraseT(unsigned int keepValue) { Ta.erase(Ta.begin(),Ta.end()-min(keepValue,(unsigned int)Ta.size())); }

    virtual void eraseT_back() { Ta.pop_back(); }

    double getTaExpl(unsigned int it) { return TaExpl.at(it); }
    double getTaExpl() { if (TaExpl.empty()) return 15.0; else return TaExpl.back(); }
    void setTaExpl(double value) { TaExpl.push_back(value); }
    virtual void eraseTExpl(unsigned int keepValue) { TaExpl.erase(TaExpl.begin(),TaExpl.end()-min(keepValue,(unsigned int)TaExpl.size())); }
    virtual void eraseTExpl_back() { TaExpl.pop_back(); }

    double getHeating(unsigned int day, unsigned int hour) { return heating.at((day-1)*24 + hour -1); }
    double getHeating(unsigned int it) { return heating.at(it); }
    double getHeating() { if (heating.empty()) return 0.0; else return heating.back(); }
    void setHeating(double watthour) { heating.push_back(watthour); }
    // erase heating only if preSimulation, otherwise keep the hourly results. Previously: else { heating.push_back(accumulate(heating.begin(),heating.end(),0.0)); heating.erase(heating.begin(),heating.end()-1) }
    void eraseHeating() { heating.erase(heating.begin(),heating.end()); }
    void eraseHeating_back() { heating.pop_back(); }
    double getTotalHeating() { return accumulate(heating.begin(),heating.end(),0.0); }

    double getCooling(unsigned int day, unsigned int hour) { return cooling.at((day-1)*24 + hour -1); }
    double getCooling(unsigned int it) { return cooling.at(it); }
    double getCooling() { if (cooling.empty()) return 0.0; else return cooling.back(); }
    void setCooling(double watthour) { cooling.push_back(watthour); }
    // erase cooling only if preSimulation, otherwise keep the hourly results. Previously: else { cooling.push_back(accumulate(cooling.begin(),cooling.end(),0.0)); cooling.erase(cooling.begin(),cooling.end()-1); }
    void eraseCooling() { cooling.erase(cooling.begin(),cooling.end()); }
    void eraseCooling_back() { cooling.pop_back(); }
    double getTotalCooling() { return accumulate(cooling.begin(),cooling.end(),0.0); }

    // Qs is the sensible heating/cooling that is provided by the HVAC & plant systems to the thermal model
    void setQs(double watthour) { Qs.push_back(watthour); }
    double getQs(unsigned int it) { return Qs.at(it); }
    double getQs() { return Qs.back(); }
    void eraseQs(bool eraseAllResults) {
        if (!eraseAllResults) {
            Qs_heating += accumulate(Qs.begin(),Qs.end(),0.0, myBinaryOperation_greater);
            Qs_cooling += accumulate(Qs.begin(),Qs.end(),0.0, myBinaryOperation_less);
        }
        Qs.erase(Qs.begin(),Qs.end());
    }
    void eraseQs_back() { Qs.pop_back(); }
    static double myBinaryOperation_greater(double a, double b) { return (b>0) ? a+b : a; }
    static double myBinaryOperation_less(double a, double b) { return (b<0) ? a+b : a; }
    double getTotalHeatingSatisfied() { return Qs_heating + accumulate(Qs.begin(),Qs.end(),0.0, myBinaryOperation_greater); }
    double getTotalCoolingSatisfied() { return Qs_cooling + accumulate(Qs.begin(),Qs.end(),0.0, myBinaryOperation_less); }

    // Qi are the internal gains of the zone
    void setQi(float watthour) { Qi.push_back(watthour); }
    float getQi(size_t index) { return Qi.at(index); }
    float getQi() { return Qi.back(); }
    void eraseQi(unsigned int keepValue) { Qi.erase(Qi.begin(),Qi.end()-min(keepValue,(unsigned int)Qi.size())); }
    void eraseQi_back() { Qi.pop_back(); }

    // methods to get the HVAC demands in the zone
    double getHVACHeat(unsigned int it) { return HVACHeat[it]; }
    double getHVACHeat() { return HVACHeat.back(); }
    void eraseHVACHeat(unsigned int keepValue) { HVACHeat.erase(HVACHeat.begin(),HVACHeat.end()-min(keepValue,(unsigned int)HVACHeat.size())); }
    void eraseHVACHeat_back() { HVACHeat.pop_back(); }
    double getHVACCool(unsigned int it) { return HVACCool[it]; }
    double getHVACCool() { return HVACCool.back(); }
    void eraseHVACCool(unsigned int keepValue) { HVACCool.erase(HVACCool.begin(),HVACCool.end()-min(keepValue,(unsigned int)HVACCool.size())); }
    void eraseHVACCool_back() { HVACCool.pop_back(); }
    double getHVACHumidification(unsigned int it) { return HVACHumidification[it]; }
    double getHVACHumidification() { return HVACHumidification.back(); }
    void eraseHVACHumidification(unsigned int keepValue) { HVACHumidification.erase(HVACHumidification.begin(),HVACHumidification.end()-min(keepValue,(unsigned int)HVACHumidification.size())); }
    void eraseHVACHumidification_back() { HVACHumidification.pop_back(); }
    double getHVACEvaporation(unsigned int it) { return HVACEvaporation[it]; }
    double getHVACEvaporation() { return HVACEvaporation.back(); }
    void eraseHVACEvaporation(unsigned int keepValue) { HVACEvaporation.erase(HVACEvaporation.begin(),HVACEvaporation.end()-min(keepValue,(unsigned int)HVACEvaporation.size())); }
    void eraseHVACEvaporation_back() { HVACEvaporation.pop_back(); }
    double getHVACReheat(unsigned int it) { return HVACReheat[it]; }
    double getHVACReheat() { return HVACReheat.back(); }
    void eraseHVACReheat(unsigned int keepValue) { HVACReheat.erase(HVACReheat.begin(),HVACReheat.end()-min(keepValue,(unsigned int)HVACReheat.size())); }
    void eraseHVACReheat_back() { HVACReheat.pop_back(); }
    double getHVACMassFlowRate(unsigned int it) { return HVACMassFlowRate[it]; }
    double getHVACMassFlowRate() { return HVACMassFlowRate.back(); }
    void eraseHVACMassFlowRate(unsigned int keepValue) { HVACMassFlowRate.erase(HVACMassFlowRate.begin(),HVACMassFlowRate.end()-min(keepValue,(unsigned int)HVACMassFlowRate.size())); }
    void eraseHVACMassFlowRate_back() { HVACMassFlowRate.pop_back(); }

    // methods to save the HVAC demands in the zone
    void setHVACHeat(double value) { HVACHeat.push_back(value); }
    void setHVACCool(double value) { HVACCool.push_back(value); }
    void setHVACReheat(double value) { HVACReheat.push_back(value); }
    void setHVACHumidification(double value) { HVACHumidification.push_back(value); }
    void setHVACEvaporation(double value) { HVACEvaporation.push_back(value); }
    void setHVACMassFlowRate(double value) { HVACMassFlowRate.push_back(value); }

    // returns the moisture content of the last time step
    unsigned int getnMoistureContent() { return moistureContent.size(); }
    double getMoistureContent() { return moistureContent.back(); }
    void setMoistureContent(double w5) { moistureContent.push_back(w5); }

    // now for the available energy
    void setHVACHeatAvailable(double value) { HVACHeatAvailable.push_back(value); }
    void setHVACCoolAvailable(double value) { HVACCoolAvailable.push_back(value); }
    void setHVACReheatAvailable(double value) { HVACReheatAvailable.push_back(value); }
    void setHVACHumidificationAvailable(double value) { HVACHumidificationAvailable.push_back(value); }
    void setHVACEvaporationAvailable(double value) { HVACEvaporationAvailable.push_back(value); }
    void setHVACMassFlowRateAvailable(double value) { HVACMassFlowRateAvailable.push_back(value); }

    double getHVACHeatAvailable(unsigned int it) { return HVACHeatAvailable[it]; }
    double getHVACHeatAvailable() { return HVACHeatAvailable.back(); }
    void eraseHVACHeatAvailable(unsigned int keepValue) { HVACHeatAvailable.erase(HVACHeatAvailable.begin(),HVACHeatAvailable.end()-min(keepValue,(unsigned int)HVACHeatAvailable.size())); }
    void eraseHVACHeatAvailable_back() { HVACHeatAvailable.pop_back(); }
    double getHVACCoolAvailable(unsigned int it) { return HVACCoolAvailable[it]; }
    double getHVACCoolAvailable() { return HVACCoolAvailable.back(); }
    void eraseHVACCoolAvailable(unsigned int keepValue) { HVACCoolAvailable.erase(HVACCoolAvailable.begin(),HVACCoolAvailable.end()-min(keepValue,(unsigned int)HVACCoolAvailable.size())); }
    void eraseHVACCoolAvailable_back() { HVACCoolAvailable.pop_back(); }
    double getHVACReheatAvailable(unsigned int it) { return HVACEvaporationAvailable[it]; }
    double getHVACReheatAvailable() { return HVACEvaporationAvailable.back(); }
    void eraseHVACReheatAvailable(unsigned int keepValue) { HVACEvaporationAvailable.erase(HVACEvaporationAvailable.begin(),HVACEvaporationAvailable.end()-min(keepValue,(unsigned int)HVACEvaporationAvailable.size())); }
    void eraseHVACReheatAvailable_back() { HVACEvaporationAvailable.pop_back(); }
    double getHVACHumidificationAvailable(unsigned int it) { return HVACHumidificationAvailable[it]; }
    double getHVACHumidificationAvailable() { return HVACHumidificationAvailable.back(); }
    void eraseHVACHumidificationAvailable(unsigned int keepValue) { HVACHumidificationAvailable.erase(HVACHumidificationAvailable.begin(),HVACHumidificationAvailable.end()-min(keepValue,(unsigned int)HVACHumidificationAvailable.size())); }
    void eraseHVACHumidificationAvailable_back() { HVACHumidificationAvailable.pop_back(); }
    double getHVACEvaporationAvailable(unsigned int it) { return HVACEvaporationAvailable[it]; }
    double getHVACEvaporationAvailable() { return HVACEvaporationAvailable.back(); }
    void eraseHVACEvaporationAvailable(unsigned int keepValue) { HVACEvaporationAvailable.erase(HVACEvaporationAvailable.begin(),HVACEvaporationAvailable.end()-min(keepValue,(unsigned int)HVACEvaporationAvailable.size())); }
    void eraseHVACEvaporationAvailable_back() { HVACEvaporationAvailable.pop_back(); }
    double getHVACMassFlowRateAvailable(unsigned int it) { return HVACMassFlowRateAvailable[it]; }
    double getHVACMassFlowRateAvailable() { return HVACMassFlowRateAvailable.back(); }
    void eraseHVACMassFlowRateAvailable(unsigned int keepValue) { HVACMassFlowRateAvailable.erase(HVACMassFlowRateAvailable.begin(),HVACMassFlowRateAvailable.end()-min(keepValue,(unsigned int)HVACMassFlowRateAvailable.size())); }
    void eraseHVACMassFlowRateAvailable_back() { HVACMassFlowRateAvailable.pop_back(); }

    // method to get the occupants
    void setOccupantsStochastic(bool value) { occupantsStochastic = value; }
    void setOccupantsNumber(float number) { occupantsNumber = number; }
    float getOccupantsNumber() { return occupantsNumber; }
    void setOccupantsSensibleHeat(float watts) { occupantsSensibleHeat = watts; }
    float getOccupantsSensibleHeat() { return occupantsSensibleHeat; }
    void setOccupantsSensibleHeatRadiantFraction(float value) { occupantsSensibleHeatRadiantFraction = value; }
    float getOccupantsSensibleHeatRadiantFraction() { return occupantsSensibleHeatRadiantFraction; }
    void setOccupantsLatentHeat(float watts) { occupantsLatentHeat = watts; }
    float getOccupantsSensibleHeatConvective() { return occupantsCount*occupantsSensibleHeat*(1.f-occupantsSensibleHeatRadiantFraction); }
    float getOccupantsSensibleHeatRadiative()  { return occupantsCount*occupantsSensibleHeat*occupantsSensibleHeatRadiantFraction; }
    float getOccupantsLatentHeat() { return occupantsLatentHeat; }
    void setActivityType(unsigned int value) { activityType = value; }
    unsigned int getActivityType() { return activityType; }
    //void setOccupantsYearProfileID(unsigned int id) { occupantsYearProfileID=id; }
    //unsigned int getOccupantsYearProfileID() { return occupantsYearProfileID; }
    void setOccupantsYearProfile(YearProfile* yd) { occupantsYearProfile=yd; }
    YearProfile* getOccupantsYearProfile() { return occupantsYearProfile; }
    void setDHWYearProfile(YearProfile* yd) { dhwYearProfile=yd; }
    YearProfile* getDHWYearProfile() { return dhwYearProfile; }

    float getOccupantsCount() { return occupantsCount; }
    void setOccupantsCountAndActivity(unsigned int day, unsigned int hour);
    float getDHWConsumption(unsigned int day, unsigned int hour);

    Occupants* getOccupants() { return pOccupants; }

    // stochastic models window state & blinds state
    float getWindowState() { if (windowState.empty()) return 0.f; else return windowState.back(); }
    float getWindowState(unsigned int step) { if (windowState.empty()) return 0.f; return windowState.at(step); }
    void setWindowState(float value) { windowState.push_back(value); }
    void eraseWindowState(unsigned int keepValue) { windowState.erase(windowState.begin(),windowState.end()-min(keepValue,(unsigned int)windowState.size())); }
    void eraseWindowState_back() { windowState.pop_back(); }
    float getLightsState() { if (lightsState.empty()) return 0.f; else return lightsState.back(); }
    float getLightsState(unsigned int step) { if (lightsState.empty()) return 0.f; return lightsState.at(step); }
    void setLightsState(float value) { lightsState.push_back(value); }
    void eraseLightsState(unsigned int keepValue) { lightsState.erase(lightsState.begin(),lightsState.end()-min(keepValue,(unsigned int)lightsState.size())); }
    bool getGroundFloor() { return groundFloor; }

    // ajouter le calcul de la consommation
    float getTotalInternalIlluminance0() {
        float illuminance = 0.f;
        for (size_t wallIndex=0; wallIndex < walls.size(); ++wallIndex)
            illuminance += walls[wallIndex]->getTotalInternalIlluminance0()*walls[wallIndex]->getLowerShadingState();
        return illuminance;
    }
    float getTotalInternalIlluminance0(unsigned int timeStep) {
        float illuminance = 0.f;
        for (size_t wallIndex=0; wallIndex < walls.size(); ++wallIndex)
            illuminance += walls[wallIndex]->getTotalInternalIlluminance0(timeStep)*walls[wallIndex]->getLowerShadingState(timeStep);
        return illuminance;
    }

    /// TODO: change the value here for reading in the XML
    float getLightsThreshold() { return 300.f; /* in lux */ }
    float getLightsPowerDensity() { return 24.f; /* in W/m² */ }

#pragma GCC diagnostic ignored "-Wunused-parameter"
    virtual double getMatrixElement(unsigned int j, unsigned int k) { return 0.0; }
    virtual double getFixedMatrixElement(unsigned int j, unsigned int k) { return 0.0; }
    virtual double getVariableMatrixElement(unsigned int j, unsigned int k) { return 0.0; }
    virtual double getC(unsigned int i) { return 0.0; }
    virtual void addC(float value) { Ci+=value; }
    virtual double getSourceTerm(unsigned int i, float Tout, float Tground, float Qs=0.f) { return 0.; }

    // returns the number of element in the vector of node i
    virtual unsigned int getnT(unsigned int i) { return 0; }
    virtual unsigned int getnTExpl(unsigned int i) { return 0; }
    // returns the last element of the temperature vector of node i
    virtual double getT(unsigned int i) { throw(string("The node type was not defined.")); }
    virtual double getTExpl(unsigned int i) { throw(string("The node type was not defined.")); }
    virtual double getTExpl(unsigned int i, unsigned int it) { throw(string("The node type was not defined.")); }
    // sets a new value to the temperature vector of node i
    virtual void setT(unsigned int i, double value) { throw(string("The node type was not defined.")); }
    virtual void setTExpl(unsigned int i, double value) { throw(string("The node type was not defined.")); }
#pragma GCC diagnostic warning "-Wunused-parameter"

    // sets the foreseen air temperature
    void setTaForeseen(double value) { TaForeseen = value; }
    double getTaForeseen() { return TaForeseen; }

    // return elements for the thermal model 1N, 2N and 3N
    virtual double getKappa1() { return 0.0; }
    virtual double getKappa2() { return 0.0; }
    virtual double getKappa3();

    virtual double getKw()  { return 0.0; }
    virtual double getG0() { return 0.0; }
    virtual double getGn() { return 0.0; }

#pragma GCC diagnostic ignored "-Wunused-parameter"
    // returns the wall temperature
    virtual double getTw(unsigned int it) { throw(string("The node type was not defined.")); }
    virtual double getTw() { throw(string("The node type was not defined.")); }

    // sets the outside surface temperature TODO: add also temperature of the ROOFS!!!
    virtual void setTos(float Tout) { throw(string("The node type was not defined.")); }
    virtual void eraseTos_back() { throw(string("The node type was not defined.")); }
#pragma GCC diagnostic warning "-Wunused-parameter"

    // gets the E+ id
    string getEp_id() { return ep_id; }
    void setEp_id(string value) { ep_id = value; }

};

class Zone2N : public Zone {

 protected:

    float hc_int=3.f;
    float Cw, Kw1, Kw2;
    float Tw=15.f, TwExpl=15.f;

 public:

    Zone2N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants);
    Zone2N(unsigned int id, Building* pBuilding, bool groundFloor):Zone(id,pBuilding,groundFloor){}// Incomplete constructor for DXF reading, do not use without completing the geometry...
    virtual void clear() {
        Tw=15.f; TwExpl=15.f;
        Zone::clear();
    }
    virtual void update(bool constructor=false);

    double getKappa1();
    double getKappa2() { return Kw2*Swa*hc_int/(Kw2+hc_int); } // if no walls, Swa=0 leading to a decoupling from the air node

    // gets the matrix elements
    double getMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return (-getUA() -getKappa2() -getKappa3() -getUAground());
        else if ( j==0 && k==1 ) return getKappa2();
        else if ( j==1 && k==1 ) return -getKappa1() -getKappa2();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getFixedMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getFixedMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return (-getKappa2() -getUAground());
        else if ( j==0 && k==1 ) return getKappa2();
        else if ( j==1 && k==1 ) return -getKappa2();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getVariableMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getVariableMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getUA()-getKappa3();
        else if ( j==0 && k==1 ) return 0.;
        else if ( j==1 && k==1 ) return -getKappa1();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getC(unsigned int node) {
        if (node == 0) return Ci;
        else if (node == 1) return Cw;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(node)));
    }

    double getSourceTerm(unsigned int i, float Tout, float Tground, float Qs=0.f) {
        // air temperature node
        if (i == 0)         return getUA()*Tout
                                   + getUAground()*Tground
                                   + getQsun3()
                                   + hc_int
                                     *(getWw()*getQsun2()+Lr)
                                     /(hc_int+Kw2)
                                   + getWa()*getQsun2()+Lc+Qs;
        // wall temperature node
        else if (i == 1) {  double term2 = Kw2*(getWw()*getQsun2()+Lr)/(Kw2+hc_int);
                            for (size_t i=0; i<walls.size(); ++i) {
                                term2 += Kw1*walls[i]->getWallArea()
                                         *(+walls[i]->get_hc()*Tout
                                           +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                                           +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature())
                                         /(Kw1+walls[i]->get_hc()+walls[i]->get_hr());
                            }
                            return term2;
        }
        // out of range
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    unsigned int getnT(unsigned int i) {
        if ( i==0 ) return Ta.size();
        else if ( i<nNodes ) return 1;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getT(unsigned int i) {
        if ( i==0 ) return getTa();
        else if ( i==1 ) return Tw;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    void setT(unsigned int i, double value) {
        if ( i==0 ) setTa(value);
        else if ( i==1 ) Tw=value;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot set element "+toString(i)));
    }

    // returns the outside surface temperature
    void setTos(float Tout);

    void eraseTos_back() {
        for (size_t i=0; i<walls.size(); ++i) walls[i]->eraseTemperature_back();
        for (size_t i=0; i<roofs.size(); ++i) roofs[i]->eraseTemperature_back();
    }

};

class Zone3N : public Zone2N {

 protected:

    float Cr, Kr1, Kr2;
    float Tr=15.f, TrExpl=15.f;

 public:

    Zone3N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants);
    Zone3N(unsigned int id, Building* pBuilding, bool groundFloor):Zone2N(id,pBuilding,groundFloor){}// Incomplete constructor for DXF reading, do not use without completing the geometry...
    virtual void clear() {
        Tr=15.f; TrExpl=15.f;
        Zone2N::clear();
    }
    virtual void update(bool constructor=false);

    float getUA() { //cout << "Zone3N::getUA Kwindow="<< Kwindow<<" Kroof="<< Kroof<<" rhoi="<<rhoi<<" Cpi=" << Cpi <<" Vi=" << Vi <<" Ninf=" << Ninf <<" Nvent=" << Nvent << " Kpsi=" << Kpsi<<endl << flush;
                    return (Kwindow + rhoi*Cpi*Vi*(Ninf/3600.f) + rhoi*Cpi*Vi*(Nvent/3600.f) + Kpsi); /* in W/K */
                    }
    double getKappa3();
    double getKappa4() { return Kr2*Sro*hc_int/(Kr2+hc_int); } // if no roofs, Sro=0 leading to a decoupling from the air node

    // gets the matrix elements
    double getMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return (-getUA() -getKappa2() -getKappa4() -getUAground());
        else if ( j<2 && k<2 ) return Zone2N::getMatrixElement(j,k);
        else if ( j==0 && k==2 ) return getKappa4();
        else if ( j==1 && k==2 ) return 0.;
        else if ( j==2 && k==2 ) return -getKappa4() -getKappa3();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getFixedMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getFixedMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return (-getKappa2() -getKappa4() -getUAground());
        else if ( j<2 && k<2 ) return Zone2N::getFixedMatrixElement(j,k);
        else if ( j==0 && k==2 ) return getKappa4();
        else if ( j==1 && k==2 ) return 0.;
        else if ( j==2 && k==2 ) return -getKappa4();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getVariableMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getVariableMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getUA();
        else if ( j<2 && k<2 ) return Zone2N::getVariableMatrixElement(j,k);
        else if ( j==0 && k==2 ) return 0.;
        else if ( j==1 && k==2 ) return 0.;
        else if ( j==2 && k==2 ) return -getKappa3();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getC(unsigned int node) {
        if (node == 0) return Ci;
        else if (node == 1) return Cw;
        else if (node == 2) return Cr;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(node)));
    }

    double getSourceTerm(unsigned int i, float Tout, float Tground, float Qs=0.f) {
        // air temperature node
        if (i == 0)         return getUA()*Tout
                                   + getUAground()*Tground
                                   + hc_int
                                     *(getWw()*getQsun2()+Lr)
                                     /(hc_int+Kw2)
                                   + getWa()*getQsun2()+Lc+Qs;
        // wall temperature node
        else if (i == 1)    return Zone2N::getSourceTerm(i,Tout,Tground,Qs);
        // roof temperature node
        else if (i == 2) {  double term3 = 0.;
                            for (size_t i=0; i<roofs.size(); ++i) {
                                term3 += Kr1*roofs[i]->getRoofArea()
                                         *(+roofs[i]->get_hc()*Tout
                                           +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())
                                           +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature()
                                           -roofs[i]->get_Y())
                                         /(Kr1+roofs[i]->get_hc()+roofs[i]->get_hr()+roofs[i]->get_X());
                            }
                            return term3;
        }
        // out of range
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    unsigned int getnT(unsigned int i) {
        if ( i==0 ) return Ta.size();
        else if ( i<nNodes ) return 1;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getT(unsigned int i) {
        if ( i==0 ) return getTa();
        else if ( i==1 ) return Tw;
        else if ( i==2 ) return Tr;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    void setT(unsigned int i, double value) {
        if ( i==0 ) setTa(value);
        else if ( i==1 ) Tw=value;
        else if ( i==2 ) Tr=value;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot set element "+toString(i)));
    }

    // returns the outside surface temperature
    void setTos(float Tout);

    void eraseTos_back() {
        for (size_t i=0; i<walls.size(); ++i) walls[i]->eraseTemperature_back();
        for (size_t i=0; i<roofs.size(); ++i) roofs[i]->eraseTemperature_back();
    }

};

class Zone3N_floor : public Zone2N {

 private:

    float Sf;
    float Cf, Kf1, Kf2;
    float Tf=15.f;//, TfExpl;

 public:

    Zone3N_floor(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants);
    virtual void clear() {
        Tf=15.f;
        Zone2N::clear();
    }
    virtual void update(bool constructor=false);

    double getKappa5() { return Kf1*Sf; }
    double getKappa6() { return Kf2*Sf*hc_int/(Kf2+hc_int); }

    // gets the matrix elements
    double getMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getUA() -getKappa2() -getKappa6();
        else if ( j<2 && k<2 ) return Zone2N::getMatrixElement(j,k);
        else if ( j==0 && k==2 ) return getKappa6();
        else if ( j==1 && k==2 ) return 0.;
        else if ( j==2 && k==2 ) return -getKappa6() -getKappa5();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getFixedMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getFixedMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getKappa2() -getKappa6();
        else if ( j<2 && k<2 ) return Zone2N::getFixedMatrixElement(j,k);
        else if ( j==0 && k==2 ) return getKappa6();
        else if ( j==1 && k==2 ) return 0.;
        else if ( j==2 && k==2 ) return -getKappa6() -getKappa5();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getVariableMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getVariableMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getUA();
        else if ( j<2 && k<2 ) return Zone2N::getVariableMatrixElement(j,k);
        else if ( k==2 ) return 0.;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getC(unsigned int i) {
        if (i<2) return Zone2N::getC(i);
        else if (i == 2) return Cf;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getSourceTerm(unsigned int i, float Tout, float Tground, float Qs=0.f) {
        // air temperature node
        if (i == 0)         return getUA()*Tout
                                   + hc_int
                                     *(Lr)
                                     /(hc_int+Kw2)
                                   + hc_int
                                     *(getWw()*getQsun2())
                                     /(hc_int+Kf2)
                                   + getWa()*getQsun2()+Lc+Qs;
        // wall temperature node
        else if (i == 1) {  double term2 = Kw2*(Lr)/(Kw2+hc_int);
                            for (size_t i=0; i<walls.size(); ++i) {
                                term2 += Kw1*walls[i]->getWallArea()
                                         *(+walls[i]->get_hc()*Tout
                                           +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                                           +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature())
                                         /(Kw1+walls[i]->get_hc()+walls[i]->get_hr());
                            }
                            return term2;
        }
        else if (i == 2)    return Kf1*Sf*Tground
                                   + Kf2*(getWw()*getQsun2())/(hc_int+Kf2);
        // out of range
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    unsigned int getnT(unsigned int i) {
        if ( i==0 ) return Ta.size();
        else if ( i<nNodes ) return 1;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getT(unsigned int i) {
        if ( i < 2 ) return Zone2N::getT(i);
        else if ( i==2 ) return Tf;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    void setT(unsigned int i, double value) {
        if ( i < 2 ) return Zone2N::setT(i,value);
        else if ( i==2 ) Tf=value;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot set element "+toString(i)));
    }

};

class Zone4N : public Zone3N {

 private:

    float Sf=0.f;
    float Cf, Kf1, Kf2;
    float Tf=15.f;//, TfExpl;

 public:

    Zone4N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants);
    // Zone extruder for DXF reading
    Zone4N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, Surface* floor, float elevation); //, vector<Surface*> surfaces
    Zone4N(Building* pBuilding, bool groundFloor);
    virtual void clear() {
        Tf=15.f;
        Zone3N::clear();
    }
    virtual void update(bool constructor=false);

    void addSurface(Surface* s);

    double getKappa5();
    double getKappa6() { return Kf2*Sf*hc_int/(Kf2+hc_int); } // if no floors, Sf=0 leading to a decoupling from the air node

    // gets the matrix elements
    double getMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getUA() -getKappa2() -getKappa4() -getKappa6();
        else if ( j<3 && k<3 ) return Zone3N::getMatrixElement(j,k);
        else if ( j==0 && k==3 ) return getKappa6();
        else if ( j==1 && k==3 ) return 0.;
        else if ( j==2 && k==3 ) return 0.;
        else if ( j==3 && k==3 ) return -getKappa6() -getKappa5();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getFixedMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getFixedMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getKappa2() -getKappa4() -getKappa6();
        else if ( j<3 && k<3 ) return Zone3N::getFixedMatrixElement(j,k);
        else if ( j==0 && k==3 ) return getKappa6();
        else if ( j==1 && k==3 ) return 0.;
        else if ( j==2 && k==3 ) return 0.;
        else if ( j==3 && k==3 ) return -getKappa6() -getKappa5();
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getVariableMatrixElement(unsigned int j, unsigned int k) {
        if ( k < j ) return getVariableMatrixElement(k,j); // symmetric matrix
        else if ( j==0 && k==0 ) return -getUA();
        else if ( j<3 && k<3 ) return Zone3N::getVariableMatrixElement(j,k);
        else if ( k==3 ) return 0.;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(j)+","+toString(k)));
    }

    double getC(unsigned int i) {
        if (i<3) return Zone3N::getC(i);
        else if (i == 3) return Cf;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getSourceTerm(unsigned int i, float Tout, float Tground, float Qs=0.f) {
        // air temperature node
        if (i == 0)         return getUA()*Tout
                                   + hc_int
                                     *(Lr)
                                     /(hc_int+Kw2)
                                   + hc_int
                                     *(getWw()*getQsun2())
                                     /(hc_int+Kf2)
                                   + getWa()*getQsun2()+Lc+Qs;
        // wall temperature node
        else if (i == 1) {  double term2 = Kw2*(Lr)/(Kw2+hc_int);
                            for (size_t i=0; i<walls.size(); ++i) {
                                term2 += Kw1*walls[i]->getWallArea()
                                         *(+walls[i]->get_hc()*Tout
                                           +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                                           +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature())
                                         /(Kw1+walls[i]->get_hc()+walls[i]->get_hr());
                            }
                            return term2;
        }
        else if (i == 2)    return Zone3N::getSourceTerm(i,Tout,Tground,Qs);
        else if (i == 3)    return Kf1*Sf*Tground
                                   + Kf2*(getWw()*getQsun2())/(hc_int+Kf2);
        // out of range
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    unsigned int getnT(unsigned int i) {
        if ( i==0 ) return Ta.size();
        else if ( i<nNodes ) return 1;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getT(unsigned int i) {
        if ( i < 3 ) return Zone3N::getT(i);
        else if ( i==3 ) return Tf;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    void setT(unsigned int i, double value) {
        if ( i < 3 ) return Zone3N::setT(i,value);
        else if ( i==3 ) Tf=value;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot set element "+toString(i)));
    }

};

class ZoneN : public Zone {

 private:

    vector<double> Tw, TwExpl;    //wall Temp. for implicit (1h) and explicit (5min)

 public:

    ZoneN(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants);
    virtual void update(bool constructor=false);

    // gets the double conductance of layer i
    double getG(unsigned int i) { return Swa*2.f*walls[0]->getComposite()->getLayer(i)->getConductance(); }
    double getG0() { return getG(0); }
    double getGn() { return getG(nNodes-2); }
    // gets the conductance between node 0 (from outside) and the outside air
    double getKappa1() { return getG0()*(getKe()+getHr())/(getG0()+(getKe()+getHr())); }
    // gets the conductance between node n-1 (from outside) and the inside air
    double getKappa2() { return getGn()*Ki/(getGn()+Ki); }
    // gets the conductance between node i-1 and i
    double getChi(unsigned int i) { return (getG(i)*getG(i+1))/(getG(i)+getG(i+1)); }

    double getMatrixElement(unsigned int j, unsigned int k) {
        // first diagonal element (air node)
        if      ( j==0 && k==0 ) return (-getUA()-getKappa2()-getUAground()-getKappa3());
        else if ( j==0 && k==(nNodes-1) ) return getKappa2();
        // second diagonal element (first wall node)
        else if ( j==1 && k==1 ) return (-getKappa1()-getChi(j-1));
        else if ( j==1 && k==2 ) return getChi(j-1);
        // other tridiagonal elements in the submatrix
        else if (j > 1 && j < (nNodes-1)) {
            if ( k==j ) return -getChi(j-2)-getChi(j-1);
            // upper diagonal element
            else if ( k==(j+1)) return getChi(j-1);
            // lower diagonal element
            else if ( k==(j-1)) return getChi(j-2);
            else return 0.f;
        }
        // last diagonal element (last wall node)
        else if ( j==(nNodes-1) && k==0 ) return getKappa2();
        else if ( j==(nNodes-1) && k==(nNodes-2) ) return getChi(j-2);
        else if ( j==(nNodes-1) && k==(nNodes-1) ) return (-getKappa2()-getChi(j-2));
        // elsewhere nil
        else return 0.f;
    }

    double getFixedMatrixElement(unsigned int j, unsigned int k) {
        // first diagonal element (air node)
        if      ( j==0 && k==0 ) return (-getKappa2()-getUAground());
        else if ( j==0 && k==(nNodes-1) ) return getKappa2();
        // second diagonal element (first wall node)
        else if ( j==1 && k==1 ) return (-getChi(j-1));
        else if ( j==1 && k==2 ) return getChi(j-1);
        // other tridiagonal elements in the submatrix
        else if (j > 1 && j < (nNodes-1)) {
            if ( k==j ) return -getChi(j-2)-getChi(j-1);
            // upper diagonal element
            else if ( k==(j+1)) return getChi(j-1);
            // lower diagonal element
            else if ( k==(j-1)) return getChi(j-2);
            else return 0.f;
        }
        // last diagonal element (last wall node)
        else if ( j==(nNodes-1) && k==0 ) return getKappa2();
        else if ( j==(nNodes-1) && k==(nNodes-2) ) return getChi(j-2);
        else if ( j==(nNodes-1) && k==(nNodes-1) ) return (-getKappa2()-getChi(j-2));
        // elsewhere nil
        else return 0.f;
    }

    double getVariableMatrixElement(unsigned int j, unsigned int k) {
        // first diagonal element (air node)
        if      ( j==0 && k==0 ) return (-getUA()-getKappa3());
        // second diagonal element (first wall node)
        else if ( j==1 && k==1 ) return (-getKappa1());
        // elsewhere nil
        else return 0.f;
    }

    double getC(unsigned int i) {
        if (i == 0) return Ci;
        else return Swa*walls[0]->getComposite()->getCapacitance(i-1);
    }

    double getSourceTerm(unsigned int i, float Tout, float Tground, float Qs=0.f) {
        // air node
        if (i == 0)             return getUA()*Tout
                                       + getUAground()*Tground
                                       + getQsun3()
                                       + getKi()
                                         *(getWw()*getQsun2()+Lr)
                                         /(getKi()+getGn())
                                       + getWa()*getQsun2()+Lc+Qs;
        // outside wall temperature node
        else if (i == 1)        return getG0()
                                       *(getKe()*Tout+getQsun1())
                                       /(getKe()+getHr()+getG0());
        // intermediate wall temperature nodes
        else if (i < nNodes-1)  return 0.;
        // inside wall temperature node
        else if (i == nNodes-1) return getGn()
                                       *(getWw()*getQsun2()+Lr)
                                       /(getGn()+getKi());
        // out of range
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    unsigned int getnT(unsigned int i) {
        if ( i==0 ) return Ta.size();
        else if ( i<nNodes ) return 1;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

    double getT(unsigned int i) {
        if ( i==0 ) return getTa();
        else if ( i<nNodes ) return Tw[i-1];
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot access element "+toString(i)));
    }

//    double getT(unsigned int i, unsigned int it) { // returns the element it of temperature vector of node i
//        if ( i==0 ) return getTa(it);
//        else if ( i==1 ) return getTw(it);
//        else if ( i==2 ) return getTw2(it);
//        else throw(string("Only three nodes, therefore cannot access element "+toString(i)));
//    }

    void setT(unsigned int i, double value) {
        if ( i==0 ) setTa(value);
        else if ( i<nNodes ) Tw[i-1]=value;
        else throw(string("Only " + toString(nNodes) + " nodes, therefore cannot set element "+toString(i)));
    }

//    // gets and sets the explicit temperature (small time step)
//
//    double getTExpl(unsigned int i) {
//        if ( i==0 ) return getTaExpl();
//        else if ( i==1 ) return getTwExpl();
//        else if ( i==2 ) return getTw2Expl();
//        else throw(string("Only three nodes, therefore cannot access element "+toString(i)));
//    }
//
//    double getTExpl(unsigned int i, unsigned int it) {
//        if ( i==0 ) return getTaExpl(it);
//        else if ( i==1 ) return getTwExpl(it);
//        else if ( i==2 ) return getTw2Expl(it);
//        else throw(string("Only three nodes, therefore cannot access element "+toString(i)));
//    }
//
//    void setTExpl(unsigned int i, double value) {
//        if ( i==0 ) setTaExpl(value);
//        else if ( i==1 ) setTwExpl(value);
//        else if ( i==2 ) setTw2Expl(value);
//        else throw(string("Only three nodes, therefore cannot set element "+toString(i)));
//    }
//

    // returns the outside surface temperature
    void setTos(float Tout);

    void eraseTos_back() {
        for (size_t i=0; i<walls.size(); ++i) walls[i]->eraseTemperature_back();
        for (size_t i=0; i<roofs.size(); ++i) roofs[i]->eraseTemperature_back();
    }
};

#endif
