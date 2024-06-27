#ifndef MODELS_H
#define MODELS_H

#include <vector>
#include <algorithm>
#include <ctime>

class Climate;
class XmlScene;
class Building;
class Zone;
class Surface;
class Wall;
class Ground;
class Tree;

//beginning of contents added by Dapeng
class DistrictEnergyCenter;
class District;
//end of contents added by Dapeng
class EnergyConversionSystem; // Cognet: added

// *** Model class, CitySim  *** //
// *** jerome.kaempf@epfl.ch *** //

class Model {

public :

    static const unsigned int dt = 3600;
    static const unsigned int dt2 = 300; //time step for stochastic model (5 min)

    static bool thermalExplicit;
    static void setThermalExplicit(bool value) { thermalExplicit = value; }

    // stability criterion
    static int ThermalWarmUpTime(Building *pBuilding);

    // the thermal model in version implicit (prevision & final calculation)
    static void ThermalStepImplicit(Building *pBuilding, Climate *pClimate, unsigned int day, unsigned int hour);
    static void ThermalStepImplicitTemperature(Building *pBuilding, Climate *pClimate, unsigned int day, unsigned int hour);
    static void ThermalStepImplicitTemperature(Ground *pGround, Climate* pClimate, unsigned int day, unsigned int hour);
    static void ThermalStepImplicitTemperature_simplified(Ground *pGround, Climate* pClimate, unsigned int day, unsigned int hour);
    static void ThermalStepTree(Tree *pTree, Climate* pClimate, unsigned int day, unsigned int hour);

    // the thermal model in version explicit (final calculation)
    static void ThermalExplicitStability(Building *pBuilding);
    static void ThermalStepExplicitTemperature(Building *pBuilding, Climate *pClimate, unsigned int day, unsigned int hour);

    // different methods for handling the matrices in the thermal model
    static int Thermal_getMatrixPosition(Building *pBuilding, int zoneNumber);
    static int Thermal_getSubMatrixPosition(Building *pBuilding, int zoneNumber);
    static double Thermal_Ke(const double& WindSpeed, const double& WindDir, const double& Tout, const std::vector<double> &azimuth, const std::vector<double> &surface);
    static double Thermal_KeClarke(const double& WindSpeed, const double& WindDir/*, const double& Tout*/, const std::vector<float> &azimuth, const std::vector<float> &surface, const double& totalSurfaceWalls);
    static float Thermal_hcClarke(const float& WindSpeed, const float& WindDir, const float &azimuth);
    static float Thermal_hcLiuHarris(const float& WindSpeed, const float& WindDir, const float& azimuth);
    static float Thermal_KeClarke(const float& WindSpeed, const float& WindDir, const std::vector<float> &azimuth, const std::vector<float> &surface);
    static float Thermal_hcCli2(Climate* pClimate, unsigned int surfaceId, unsigned int day, unsigned int hour);
    static float Thermal_KeWalls(Climate* pClimate, Zone* pZone, unsigned int day, unsigned int hour);
    static float Thermal_HrWalls(Zone* pZone);
    static float Thermal_Uprime(Surface* pSurface, float hc);
    static float Thermal_Kwindows(Climate* pClimate, Zone* pZone, unsigned int day, unsigned int hour);
    static float Thermal_Kroofs(Climate* pClimate, Zone* pZone, unsigned int day, unsigned int hour);
    static void Thermal_Kgrounds(Climate* pClimate, Ground* pGround, unsigned int day, unsigned int hour);
    static double Thermal_KiAlamdariHammond(const double& dt, const double Lvertical, const double Lhorizontal, const std::vector<double> &surface, const double& ceiling, const double& floor);
    static double Thermal_KiCibse(const double& dt, const std::vector<double> &surface, const double& ceiling, const double& floor);
    static double Thermal_KiFixed(const double& surfaceWalls);

    static void HVAC_Needs(Building *pBuilding,Climate* pClimate,unsigned int day,unsigned int hour);
    static void HVAC_Available(Building *pBuilding,Climate* pClimate,unsigned int day,unsigned int hour);
    static double HVAC_saturatedSteamPressure(double T);
    static double HVAC_moistureContent( double T, double RH, double Patm);
    static double HVAC_relativeHumidity(double T, double w, double Patm);
    static double HVAC_supplyAirTemperature(double Qs, double T);
    static double HVAC_massFlowRate(double T2, double Patm, double w2, double Qs, double &T3, double T5, double np);
    static double HVAC_totalMoistureContent(double w, double Ql, double mdot);
    static double HVAC_temperatureChange(double deltaT, double T2, double T3);
    static double HVAC_moistureControl(double Tmax, double &T2, double T3, double &deltaTHVAC, double deltaT, double T5, double &w2, double w5prime, double Patm, double &evaporation, bool evaporativeCooling);
    static double HVAC_enthalpyDryAir(double T);
    static double HVAC_enthalpyWaterVapour(double T);
    static double HVAC_humidify(/*double T2, */double w2, double T3, double w3);
    static double HVAC_reheat(double T3, double w3, double Patm);
    static double HVAC_bulbTemperature(double w, double RH, double Patm);
    static double HVAC_bulbTemperature(double cooling, double t2prime, double w2prime, double Patm);
    static void HVAC_heat(double T2, double w2, double T3, double w3, /*double Patm,*/ double &heating, double &humidification);
    static void HVAC_cool(double T2, double w2, double T3, double &w3, /*double T5, double w5,*/ double Patm, double &cooling, double& reheating, double &humidification);

    static void HVAC_Control(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour);
    static void noHVAC_Control(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour);
    static void noHVAC_Control_EnergyHub(District* pDistrict, Climate* pClimate, unsigned int day, unsigned int hour);
    // Cognet: Start of added content.
    /**
     * For a building, searches through all roofs, walls, and surfaces for thermal solar panels, adds up and returns the power produced.
     * @param pBuilding Pointer to building
     * @param pClimate Pointer to the climate
     * @param day Day of the year
     * @param hour Hour of the day
     * @param tankTemp The temperature of the heat storage tank that the solar panels are providing power to
     * @param saveValue Whether or not the the surfaces should save these as the values of solar thermal production for that time step.
     * @return Heating power provided by solar panels of building pBui to a tank at tankTemp.
     */
    static double computeSolarThermalPower(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour, double tankTemp, bool saveValue);

    /**
     * Computes the temperature of the source where energy conversion unit, if it was a heat pump would get its energy. Then adds "minus" the produced energy to the electric consumption.
     * @param pECS Pointer to EnergyConversionSystem (useful when it is a heat pump)
     * @param pClimate Pointer to the climate
     * @param day Day of the year
     * @param hour Hour of the day
     */
    static float computeHeatPumpSrcTemp(EnergyConversionSystem* pECS, Climate* pClimate, unsigned int day, unsigned int hour);

    /**
     * Sums up the electric production of PV panels (from roofs, walls and surfaces), and the wind electric prower (from roofs). Then adds "minus" the produced energy to the electric consumption. (Be careful to setElectricConsumption to zero beforehand (normally done in Model::ThermalStepImplicit).)
     * @param pBuilding = pointer to building
     * @param pClimate = pointer to the climate
     * @param day = day at which to compute the electric production
     * @param hour = hour at which to compute the electric production
     */
    static void computeAndAddPhotovoltaicAndWindElectricProduction(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour);

    /**
     * Stores in memory the DHW usage summed over all the building's zones.
     * @param pBuilding Pointer to building
     * @param day Day of the year
     * @param hour Hour of the day
     */
    static void computeAndSetDHWUsage(Building* pBuilding, unsigned int day, unsigned int hour);

    /**
     * Initial step, compute/sum over zones the heat/cooling needs, the DHW usage, the electric productions (solar and wind). Sets MachinePower and FuelConsumption to zero.
     * @param pBuilding Pointer to building
     * @param pClimate Pointer to climate
     * @param day Day of the year
     * @param hour Hour of the day
     */
    static void noHVAC_Control_Init(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour);

    /**
     * Compute for the building, the needs in heating and cooling power needs (that the heating and cooling units will be asked to supply). Distributes and sets the solar thermal power effectively used.
     * @param pBui Pointer to the building
     * @param pClim Pointer to the climate
     * @param day Day of the year
     * @param hour Hour of the day
     */
    static void noHVAC_Control_Needs(Building* pBui, Climate *pClim, unsigned int day, unsigned int hour);

    /**
     * For all buildings in district, sets the thermalPowerNeeded to the heating/cooling units. Then computes the thermal power that heating/cooling units will be able to provide. To do so it simulates the District Energy Centers.
     * @param pDis Pointer to the district
     * @param pClimate Pointer to the climate
     * @param day Day of the year
     * @param hour Hour of the day
     */
    static void noHVAC_Control_ThermalPower(District* pDis, Climate* pClimate, unsigned int day, unsigned int hour);

    /**
     * Computes and sets : the new temperatures of tanks, Machine Power (of heat/cooling units), Fuel Consumption (of heat/cooling units), Electric Consumption (of heat/cooling units). Distributes the heating/cooling to the zones of the building.
     * @param pBui Pointer to building
     * @param pClim Pointer to climate
     * @param day Day of the year
     * @param hour Hour of the day
     */
    static void noHVAC_Control_Finish(Building* pBui, Climate* pClim, unsigned int day, unsigned int hour);
    // Cognet: End of added content.

    static float deterministicShadingAction(float irradiance, float irradianceCutOff=100.f, float lambda=0.2f);
    static double deterministicWindowsNvent(double NventMax, double Tin, double Tout);

    // stochastic models
    static float ventilationFlowRate(float Ti, float Te, float V, float A);
    static void lowerShadingAction(Zone* pZone, Wall* pWall, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour);
    static void upperShadingAction(Zone* pZone, Wall* pWall, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour);

    static float randomWeibull(float scale, float shape);
    static void windowAction_Markov(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour);
    static void windowAction_Hybrid(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour);
    static void windowAction_Bernoulli(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour);
    static void windowAction_Humphreys(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour);

    static void lightAction_Lightswitch2002(Zone* pZone, unsigned int day, unsigned int hour, unsigned int fracHour);
    static void lightAction_Threshold(Zone* pZone);

    static float lightsElectricConsumption(Zone* pZone);

    static double es(double ta);
    static void computeCMIndices(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour);

};

#endif
