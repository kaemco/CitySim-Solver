#ifndef PLANT_H
#define PLANT_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>
#include <queue>
#include <stack>
#include <array>
#include <cstdint>
#include <functional>

#include "climate.h"
#include "util.h"
#include "thermal_parcel_transport.h"
#include "mass_grid_transport.h"

//beginning of contents added by Dapeng
#include "tinyxml.h"
class District;
class Building;
class DistrictEnergyCenter;
class Network;
//end of contents added by Dapeng
class NodePair;
class SubstationNodePair;
class ThermalStationNodePair;
class PipePair;
class Pipe;
class TemperatureSetpoint;
class PressureSetpoint;
class PIDController;
class PIDControllerValve;
class MassFlowSetpoint;
class MCR;
struct DECSimulationConfig;

using namespace std;

// *** Plant class, CitySim  *** //
// *** jerome.kaempf@epfl.ch *** //

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

class Polynomial {

private:
    vector<double> a; // coefficients of the polynomial regression

public:

    Polynomial(TiXmlHandle hdl);
    size_t getSize() { return a.size(); }
    double get_a(unsigned int i) { return a.at(i); }

};

class PhotoVoltaic {

protected:

    double etampref, tref, tcnoct, muvoc, vmp;
    double toutsoc = 20.;
    double gtsoc = 800.;
    Polynomial* iam_polynomial = nullptr;
    string name = "No name";

#ifdef DEBUG
    // IAM characteristics
    vector<pair<float,double>> iam;
#endif // DEBUG

    ostream logStream;

public:

    PhotoVoltaic(double etampref, double tref, double tcnoct, double muvoc, double vmp):
        etampref(etampref),tref(tref),tcnoct(tcnoct),muvoc(muvoc),vmp(vmp),logStream(std::cout.rdbuf()) {}

    PhotoVoltaic(double pmp, double ac, double tref, double tcnoct, double muvoc, double vmp):
        tref(tref),tcnoct(tcnoct),muvoc(muvoc),vmp(vmp),logStream(std::cout.rdbuf()) {etampref = pmp/(ac*1000.0);}

    PhotoVoltaic(const PhotoVoltaic& pv):
        etampref(pv.etampref),tref(pv.tref),tcnoct(pv.tcnoct),muvoc(pv.muvoc),vmp(pv.vmp),toutsoc(pv.toutsoc),gtsoc(pv.toutsoc),logStream(std::cout.rdbuf()) {}

    PhotoVoltaic(TiXmlHandle hdl, ostream* pLogStr=nullptr);
    virtual ~PhotoVoltaic();

    void writeXML(ofstream& file, float ratio, string tab);

    virtual double getMaxPowerEfficiency(double gt, double tout);
    string getName(){ return name;}

    double getIAM(float elevation); //!< Returns ths IAM with the elevation is given in radians

};

class SolarThermal {

private:

    string name = "No name";
    double eta0, a1, a2;

    ostream logStream;

public:

    SolarThermal(double eta0, double a1, double a2):eta0(eta0),a1(a1),a2(a2),logStream(std::cout.rdbuf()) {}

    SolarThermal(TiXmlHandle hdl, ostream* pLogStr=NULL);

    void writeXML(ofstream& file, float ratio, string tab);

    double getEta0() { return eta0; }
    double getA1()   { return a1; }
    double getA2()   { return a2; }
    string getName(){ return name;}

};

class SolarHybrid : public PhotoVoltaic {

private:

    // fixed parameters of the model
    const float FF=0.7595f, eta_exch=0.94f, eta0=0.554f, x_g=0.0032f, lambda_g=1.f, tau_g=0.88f, epsilon_g=0.85f, alpha_g=0.06f, x_air=0.005f, lambda_air=0.024f, x_eva=0.0005f, lambda_eva=0.007f, x_pv=3e-4, lambda_pv=150.f, alpha_pv=0.8f, tau_pv=0.95f, x_ab=5e-3, lambda_ab=25.f, epsilon_ab=0.6f, x_bs=5e-4, lambda_bs=0.007f, eff_ab=0.7289f, Cp_w=3800.f;
    // user parameters
    float Pth, massFlowRate=0.025f;
    // internal results for the iterative calculations
    float Tg=15.f+273.15f, Tpv=15.f+273.15f, Tab=15.f+273.15f, Tmw=25.f+273.15f;

public:

    SolarHybrid(TiXmlHandle hdl, ostream* pLogStr=NULL);

    float getThermalSurfacePowerDensity(float gt, float tout, float windspeed, float Tsky, float Tin);

    double getMaxPowerEfficiency(double gt, double tout) {
        return etampref*(1+(muvoc/vmp)*(Tpv-(tref+273.15f)));
    }

    //void writeXML(ofstream& file, float ratio, string tab);

};

class WindTurbine {

private:

    unsigned int id;
    float height, alpha, gamma;

    double Pr, vi, vr, vm, c;

public:

    WindTurbine(unsigned int id, float height, float alpha, float gamma):id(id),height(height),alpha(alpha),gamma(gamma) {}

    WindTurbine(double Pr, double vi, double vr, double vm, double c):Pr(Pr),vi(vi),vr(vr),vm(vm),c(c) {}

    unsigned int getId() {return id;}
    double getPr() { return Pr; }
    double getvi() { return vi; }
    double getvr() { return vr; }
    double getvm() { return vm; }
    double getc()  { return c; }

    // note: the met station is assumed to be at 10 m above ground (from Awbi, 1991, p.64) and in a rural site
    float getElectricPowerConversionEfficiency() { return alpha*pow( height/10.f, gamma)/0.85f; }

};

class MicroWindTurbine {

private:

    float cutInSpeed, ratedSpeed,cutOutSpeed, c1, c2, c3, height, testAirDensity, alpha, gamma;

public:

    MicroWindTurbine(float cutInSpeed, float ratedSpeed, float cutOutSpeed, float c1, float c2, float c3, float testAirDensity, float height, float alpha, float gamma)
    :cutInSpeed(cutInSpeed),ratedSpeed(ratedSpeed),cutOutSpeed(cutOutSpeed),c1(c1),c2(c2),c3(c3),height(height),testAirDensity(testAirDensity),alpha(alpha),gamma(gamma) {}

    float getWindElectricPower(float windSpeed, float metStationAltitude, float roofHeight) {

        float windPower = 0.f;
        // note: the met station is assumed to be at 10 m above ground (from Awbi, 1991, p.64) and in a rural site
        windSpeed = windSpeed*alpha*pow((roofHeight+height)/10.f,gamma)/0.85f;

        if (windSpeed <= cutInSpeed)
            windPower = 0.f;
        else if ((windSpeed > cutInSpeed) && (windSpeed < ratedSpeed))
            windPower = c1*pow((windSpeed-cutInSpeed), 2) + c2*(windSpeed-cutInSpeed)+ c3;
        else if (windSpeed >= ratedSpeed  && windSpeed < cutOutSpeed)
            windPower = c1*pow((ratedSpeed-cutInSpeed), 2) + c2*(ratedSpeed-cutInSpeed) + c3;
        else  if (windSpeed >= cutOutSpeed)
            windPower= 0.f;

        return windPower*Climate::getAirDensity(metStationAltitude+roofHeight+height)/testAirDensity; // standard at sea level with air density = 1.201385 kg/m3 ,  IEC 61400-12-1 more info: https://energypedia.info/wiki/Estimation_of_Wind_Energy_Production#Power_Curve_and_Air_Densit
     }

};

class Tank {

private:

    double volume, rho, Cp, phi, Tmin, Tmax, Tinlet=10., Tcritical=90.;

public:

    Tank(double tankVolume, double tankHeight, double tankDiameter, double tankThickness, double tankRho, double tankCp, double tankTmin, double tankTmax) : volume(tankVolume), rho(tankRho), Cp(tankCp), Tmin(tankTmin), Tmax(tankTmax) {

        const double Pi = 4.0*std::atan(1.0);

        double D1_heat		=	2.0*std::sqrt(volume/(Pi*tankHeight));	    //diametre interieur
        double D2_heat		=	D1_heat+2*tankThickness;		            //diametre exterieur
        double lambdA		=	0.029*1.1622;				                //coefficient de transmission de la chaleur de l'isolant (Hypothese: mousse rigide de polyurethane) en W/(m2.C)
        double alpha1_heat	=	156.0*std::pow((1/D1_heat),0.25);	        //coefficient de transmission superficielle dans l'eau du reservoir (Hypothese: la difference de temperature entre paroi et fluide est egale  1¡C,  valider)
        double alpha2_heat	=	1.25;					                    //coefficient de transmission superficielle dans l'air (Hypothese: la difference de temperature entre paroi et fluide est egale  1¡C,  valider)

        phi = Pi*tankHeight/(std::log(D2_heat/D1_heat)/(2.0*lambdA)+1.0/(D1_heat*alpha1_heat)+1.0/(D2_heat*alpha2_heat));

    }

    Tank(double tankVolume, double tankPhi, double tankRho, double tankCp, double tankTmin, double tankTmax) : volume(tankVolume), rho(tankRho), Cp(tankCp), phi(tankPhi), Tmin(tankTmin), Tmax(tankTmax) { }
    Tank(TiXmlHandle hdl) {
        volume = to<double>(hdl.ToElement()->Attribute("V"));
        phi = to<double>(hdl.ToElement()->Attribute("phi"));
        rho = to<double>(hdl.ToElement()->Attribute("rho"));
        Cp = to<double>(hdl.ToElement()->Attribute("Cp"));
        Tmin = to<double>(hdl.ToElement()->Attribute("Tmin"));
        Tmax = to<double>(hdl.ToElement()->Attribute("Tmax"));
        if (hdl.ToElement()->Attribute("Tinlet")) Tinlet = to<double>(hdl.ToElement()->Attribute("Tinlet"));
        if (hdl.ToElement()->Attribute("Tcritical")) Tcritical = to<double>(hdl.ToElement()->Attribute("Tcritical")); // Cognet: By default, this is initialized to 90.
    }
    virtual ~Tank() {}

    double getVolume() { return volume; }
    double getRho() { return rho; }
    double getCp() { return Cp; }
    double getTmin() { return Tmin; }
    double getTmax() { return Tmax; }
    double getTinlet() { return Tinlet; }
    double getTcritical() { return Tcritical; }

    double getPhi() { return phi; }

    void setTmin(double tankTmin) { Tmin = tankTmin; }
    void setTmax(double tankTmax) { Tmax = tankTmax; }

    // models

    virtual double temperature(double t, double VdotUsed, double Pp2, double Pup2, double T0, double Tinlet, double Tamb);
    virtual double power(double t, double Tf, double VdotUsed, double Pup2, double T0, double Tinlet, double Tamb);
    virtual double time(double Tf, double VdotUsed,double Pp2, double Pup2, double T0, double Tinlet, double Tamb);
    virtual double domesticHotWater(double t, double Tf, double VdotUsedUp, double Pp2, double Pup2, double T0, double Tinlet, double Tamb);
    virtual double maxSolPowerToNotExceedTcrit(double t, double VdotUsed, double power, double solPower, double T0, double Tinlet, double Tamb); // Cognet: Added this.

    void writeXML(ofstream& file, string tag, string tab=""){
        file << tab << "<" << tag <<" V=\""<< volume << "\" phi=\"" << phi << "\" rho=\"" << rho << "\" Cp=\"" << Cp << "\" Tmin=\"" << Tmin << "\" Tmax=\"" << Tmax << "\"/>" << endl;
    }
};

class TankPCM : public Tank {

private:

    double mass, a, b, Tm, Cs;

public:

    TankPCM(double tankVolume, double tankPhi, double tankRho, double tankCp, double tankTmin, double tankTmax, double massPCM, double aPCM, double bPCM, double TmPCM, double CsPCM) : Tank(tankVolume, tankPhi, tankRho, tankCp, tankTmin, tankTmax), mass(massPCM), a(aPCM), b(bPCM), Tm(TmPCM), Cs(CsPCM) { }

    double Cp(double T) { if (b != 0.0) return Cs + a*std::exp(-(1./2.)*std::pow((T-Tm)/b, 2.0));
        else return 0.0; }

    double temperature(double t, double VdotUsed, double Pp2, double Pup2, double T0, double Tinlet, double Tamb);
    double power(double t, double Tf, double VdotUsed, double Pup2, double T0, double Tinlet, double Tamb);

};

class Equations {

public:

    // PV
    static double photoVoltaicMaxPowerEfficiency(PhotoVoltaic *panel, double gt, double ta); // Duffie and Beckman

    // Solar Heaters
    static double solarHeaterEfficiency(SolarThermal *panel, double gt, double xsi);

    // Wind Turbines
    static double windTurbinePower(WindTurbine *turbine, double v);

    // Wind Speed Model
    static double windSpeedRatio(int type, double height, int typeRef, double heightRef);

};

class EnergyConversionSystem {

protected:

    // run period of the system
    unsigned int beginDay, endDay;

    double boilerCO2EmissionCoefficient;
    double electricCO2EmissionCoefficient;
    double coGenCO2EmissionCoefficient;

    double thermalPowerNeeded; // Cognet: Added this. Now first setThermalPowerNeeded(), then use getThermalPower() to see thermal power that can be provided.

    double thermalPowerNeededHS; // Added by Max
    double thermalPowerNeededDHW; // Added by Max

public:
    ostream logStream;

    EnergyConversionSystem(unsigned int beginDay=1, unsigned int endDay=365, ostream* pLogStr = &std::cout):beginDay(beginDay),endDay(endDay),logStream(std::cout.rdbuf()) {
        // CO2 coefficients
        logStream.rdbuf(pLogStr->rdbuf());

        coGenCO2EmissionCoefficient = 0.238/3.6e6;
        boilerCO2EmissionCoefficient = 0.238/3.6e6; //238g/kWh gas naturel = 0.238 kg/kWh = 0.238e-3 kg/Wh .. in kg/J!
        electricCO2EmissionCoefficient = 0.129/3.6e6; // 129g/Kwh swiss electricity mix = 0.129 kg/Wh = 0.129e-3 kg/Wh .. in kg/J
    }
    virtual ~EnergyConversionSystem() { /*cerr << "Destructor of EnergyConversionSystem" << endl;*/ }

    unsigned int getBeginDay(){ return beginDay;}
    unsigned int getEndDay(){ return endDay;}

    // check if the system is working
    bool isWorking(unsigned int day) {
        // check the given period
        if (beginDay < endDay) { // standard order
            if ( day >= beginDay && day <= endDay) return true;
        }
        else { // reversed order, for example given from september until april
            if ( day >= beginDay || day <= endDay ) return true;
        }
        return false;
    }

    // version domestique
    virtual void   getMaxThermalPower(double thermalPower1Needed, double thermalPower2Needed, double &thermalPower1Available, double &thermalPower2Available, double sourceTemp=0.0) { return; }
    virtual double getElectricProduction(double thermalPower1, double thermalPower2, double sourceTemp) { return 0.0; }
    virtual double getCO2Production(double time, double thermalPower1, double thermalPower2, double sourceTemp) { return 0.0; }
    virtual double getFuelConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) { return 0.0; }
    virtual double getElectricConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) { return 0.0; }
    virtual float getThermalPowerMax(double sourceTemp){return 0.f;} // Added by Max. Useful for the MCR

    // version simplifi�e
    virtual string getLabel() { return "EnergyConversionSystem"; } // Cognet: Corrected spelling.

    //    virtual double getThermalPower(double thermalPowerNeeded, double sourceTemperature) { return 0.0; } // Cognet: Deleted this.
    virtual double getThermalPower(double sourceTemperature) { return 0.0; } // Cognet: Added this. Now first setThermalPowerNeeded(), then use getThermalPower() to see thermal power that can be provided.
    virtual double getThermalPowerHS(double sourceTemperature) { return 0.0; } // Added by Max. Useful for the 2stages HP.
    virtual double getThermalPowerDHW(double sourceTemperature) { return 0.0; } // Added by Max. Useful for the 2stages HP.

    virtual double getFuelConsumption(double time, double thermalPower, double sourceTemp) { return 0.0; }
    virtual double getElectricConsumption(double time, double thermalPower, double sourceTemp) { return 0.0; }

    // version distribu�e par le r�seau
    virtual void   getMaxThermalPower(vector<double> thermalPowerNeeded, vector<double> &thermalPowerAvailable, double sourceTemp) { return; }
    virtual double getElectricProduction(vector<double> thermalPower, double sourceTemp) { return 0.0; }
    virtual double getCO2Production(double time, vector<double> thermalPower, double sourceTemp) { return 0.0; }

    static double epsilonC(double sourceTemp, double outputTemp) { return ((273.15+outputTemp)/(outputTemp - sourceTemp)); }

    // special aux heat pumps
    virtual float getSourceTemperature(Climate* pClimate, unsigned int day, unsigned int hour) { return numeric_limits<float>::quiet_NaN(); }

    virtual void setThermalPowerNeeded(double tPN) { thermalPowerNeeded=tPN; } // Cognet: Added this.
    float getThermalPowerNeeded() const { return thermalPowerNeeded; }

    virtual void setThermalPowerNeededHS(double tPN) { thermalPowerNeededHS=tPN; } // Added by Max.
    virtual void setThermalPowerNeededDHW(double tPN) { thermalPowerNeededDHW=tPN; } // Added by Max.

    virtual void writeXML(ofstream& file, string tab)=0;
    virtual void writeGML(ofstream& file, string tab) {}
};

class Boiler : public EnergyConversionSystem {

private:

    double boilerThermalPower, boilerThermalEfficiency;
    string name="";

public:

    Boiler(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);

    Boiler(double boilerThermalPower, double boilerThermalEfficiency,unsigned int beginDay=1,unsigned int endDay=365):EnergyConversionSystem(beginDay,endDay),boilerThermalPower(boilerThermalPower),boilerThermalEfficiency(boilerThermalEfficiency) {}
    ~Boiler() { /*cerr << "Destructor of Boiler..." << endl;*/ }

    void writeXML(ofstream& file, string tab);
    void writeGML(ofstream& file, string tab);

    friend ostream& operator<< (ostream& s, Boiler& unit) {

        s << "\n\nBoiler:\nThermal Power: " << unit.boilerThermalPower << " W(th)\nEfficiency: " << unit.boilerThermalEfficiency << endl;
        return s;

    }

    string getLabel() { return "Boiler"; }

//    double getThermalPower(double thermalPowerNeeded, double sourceTemp) { // Cognet: Deleted this.
    double getThermalPower(double sourceTemp) override { // Cognet: Added this.

        if ( thermalPowerNeeded<=0 ) return 0.; // Cognet: Added this.
        if ( thermalPowerNeeded <= boilerThermalPower ) return thermalPowerNeeded;
        else return boilerThermalPower;

    }

    double getFuelConsumption(double time, double thermalPower, double sourceTemp) {

        if (thermalPower <= boilerThermalPower) return (time*(thermalPower)/boilerThermalEfficiency);
        else return (time*boilerThermalPower/boilerThermalEfficiency);

    }

    float getThermalPowerMax(double sourceTemp)override{return boilerThermalPower;} //Added by Max
};

class CoGeneration: virtual public EnergyConversionSystem {

protected:

    double coGenThermalPower, coGenElectricalEfficiency, coGenThermalEfficiency, coGenMinPartLoadCoefficient;

public:

    CoGeneration(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    CoGeneration(double coGenThermalPower,double coGenElectricalEfficiency,double coGenThermalEfficiency,double coGenMinPartLoadCoefficient,unsigned int beginDay=1,unsigned int endDay=365);

    ~CoGeneration() { /*cerr << "Destructor of CoGeneration..." << endl;*/ }

    void writeXML(ofstream& file, string tab);

    friend ostream& operator<< (ostream& s, CoGeneration& unit){
        s << "\nCoGeneration:\nThermal Power: \t" << unit.coGenThermalPower << " W(th)\nElectrical efficiency: " << unit.coGenElectricalEfficiency << "\nThermal efficiency: " << unit.coGenThermalEfficiency << "\nMin. part-load coefficient: " << unit.coGenMinPartLoadCoefficient << endl;
        return s;
    }

    string getLabel() { return "CoGeneration"; }

    virtual double getThermalPower(double sourceTemp) override; // Cognet: Added this.
    virtual double getFuelConsumption(double time, double thermalPower, double sourceTemp);
    virtual double getElectricConsumption(double time, double thermalPower, double sourceTemp);
    virtual float getThermalPowerMax(double sourceTemp) override {return coGenThermalPower;}
};

class HeatPump : virtual public EnergyConversionSystem {

protected:

    double heatPumpElectricPower;
    double targetTemp;
    double etaTech;

    // for the ground heat pump
    bool ground;
    float z0,z1,alpha;

    // for the fixed source temperature heat pump
    float Tsource=numeric_limits<float>::quiet_NaN();

public:

    HeatPump(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    HeatPump(double heatPumpElectricPower,double heatPumpCOP,double heatPumpSrcTemp,double heatPumpOutputTemp,unsigned int beginDay,unsigned int endDay);
    HeatPump(double heatPumpElectricPower,double heatPumpEtaTech,double targetTemp,unsigned int beginDay,unsigned int endDay);

    ~HeatPump() { /*cerr << "Destructor of HeatPump..." << endl;*/ }

    void writeXML(ofstream& file, string tab);
    void writeGML(ofstream& file, string tab);

    friend ostream& operator<< (ostream& s, HeatPump& unit){
        s << "\nHeat Pump:\nElectrical Power: \t" << unit.heatPumpElectricPower << " W(el)\nEta tech: " << unit.etaTech << endl;
        return s;
    }

    void setGround(float z0, float z1, float alpha) { this->z0=z0; this->z1=z1; this->alpha=alpha; ground = true; }
    float getSourceTemperature(Climate* pClimate, unsigned int day, unsigned int hour) {
        if (isnan(Tsource)) { // in this case, we compute the temperature based on the Climate
            if (ground) return pClimate->getTgroundCelsius(day,hour,z0,alpha,z1); // ground temperature
            else return pClimate->getToutCelsius(day,hour); // air temperature
        }
        else return Tsource;
    }

    string getLabel() { return "HeatPump"; }

    double getHeatProduced(double work, double sourceTemp, double outputTemp);
    double getWorkNeeded(double thermalPower, double sourceTemp, double outputTemp);
    double getWorkNeededEvap(double PowerEvap, double sourceTemp, double outputTemp); // Added by Max. This function allows to compute the electricity power knowing the power at the evaporator.
    double getTargetTemp(){return targetTemp;}; // Added by Max.
    double getHeatPumpElectricPower(){return heatPumpElectricPower;}; // Added by Max.
    double getEtaTech(){return etaTech;}; // Added by Max.

//    virtual double getThermalPower(double thermalPowerNeeded, double sourceTemp); // Cognet: Deleted this.
    virtual double getThermalPower(double sourceTemp) override; // Cognet: Added this.
    virtual double getElectricConsumption(double time, double thermalPower, double sourceTemp);
    virtual float getThermalPowerMax(double sourceTemp) override{return getHeatProduced(heatPumpElectricPower,sourceTemp,targetTemp);}
};

class CoGenerationHeatPump : public CoGeneration, public HeatPump {

public:

    CoGenerationHeatPump(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);

    CoGenerationHeatPump(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient,
                         double heatPumpCOP, double heatPumpSrcTemp, double heatPumpOutputTemp,unsigned int beginDay=1,unsigned int endDay=365);
    CoGenerationHeatPump(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient,
                         double heatPumpEtaTech, double targetTemp,unsigned int beginDay=1,unsigned int endDay=365);

    ~CoGenerationHeatPump() { /*cerr << "Destructor of CoGenerationHeatPump..." << endl;*/ }

    void writeXML(ofstream& file, string tab);

    friend ostream& operator<< (ostream& s, CoGenerationHeatPump& unit) {

        s << "\nCoGeneration:\nThermal Power: \t" << unit.coGenThermalPower << " W(th)\nElectrical efficiency: " << unit.coGenElectricalEfficiency << "\nThermal efficiency: " << unit.coGenThermalEfficiency << "\nMin. part-load coefficient: " << unit.coGenMinPartLoadCoefficient;
        s << "\nHeat pump:\nEtatech: " << unit.etaTech << "\nTarget Temperature: " << unit.targetTemp << endl;
        return s;

    }


    string getLabel() { return "CoGenerationHeatPump"; }

//    double getThermalPower(double thermalPowerNeeded, double sourceTemp); // Cognet: Deleted this.
    double getThermalPower(double sourceTemp) override; // Cognet: Added this.
    double getFuelConsumption(double time, double thermalPower, double sourceTemp);
    double getElectricConsumption(double time, double thermalPower, double sourceTemp) { return 0.0; }
    float getThermalPowerMax(double sourceTemp) override {return coGenThermalPower;}
};

// Added by Max
class EfficiencyPump{
public:
    static EfficiencyPump* createNewEfficiencyPump(TiXmlHandle hdl);
    EfficiencyPump() {}
    virtual ~EfficiencyPump() {}
    virtual float computeEfficiency(float const& massFlow) = 0;
};

class ConstantEfficiencyPump: public EfficiencyPump {
private:
    float efficiencyPump; // must be in between 0 and 1.
public:
    ConstantEfficiencyPump(TiXmlHandle hdl);
    ConstantEfficiencyPump(float efficiencyPump):efficiencyPump(efficiencyPump) {}
    virtual float computeEfficiency(float const& massFlow) override { return efficiencyPump; }
};

class AffineEfficiencyPump : public EfficiencyPump {
protected:
    // The setpoint function of the pressure difference is a function of mass flow, affine by parts.
    vector<float> massFlows; // [kg/s]
    vector<float> efficiencyPumps; // [Pa]
public:
    AffineEfficiencyPump(TiXmlHandle hdl);
    virtual float computeEfficiency(float const& massFlow) override;
}; // end of Added by Max

class Pump {
private:
    // Stays constant, define the Pump.
    EfficiencyPump* efficiencyPump = nullptr;
    float n0; // Nominal rotational speed [rotations/min].
    float a0; // Coefficient of polynomial for flow/pressure curve [Pa] (eg. 1247180.f).
    float a1; // Coefficient of polynomial for flow/pressure curve [Pa*s/kg] (eg. -1640.236f).
    float a2; // Coefficient of polynomial for flow/pressure curve [Pa*(s/kg)^2] (eg. -0.00016031f).

    // Evolves at each time step.
    float n; // Current rotational speed [rotations/min].

    float computeNMin() { return n0*0.001f; } // One one thousandth of max.

public:
    Pump(TiXmlHandle hdl);
    ~Pump() { if (efficiencyPump) delete efficiencyPump; }
    float getn() { return n; }

    /**
     * @brief computeElectricAndThermalPower
     * @param pressureDiff Pressure loss through the pump (negative in normal conditions, since the pump increases pressure) [Pa].
     * @param massFlow Mass flow through the pump (positive in normal conditions) [kg/s].
     * @param rho Mass density of fluid going through pump [kg/m^3].
     * @param electricPow Return variable for electric power (positive when pump consumes electricity) [W].
     * @param thermalPow Return variable for thermal power (positive when fluid temperature increases in the direction of the flow) [W].
     */
    void computeElectricAndThermalPower(float const& pressureDiff, float const& massFlow, float const& rho, float& electricPow, float& thermalPow);
    float cpdT(float const& pressureDiff, float const& rho, float const& massFlow); //Added by Max
    /**
     * @brief computePressureDiff
     * @param massFlow Mass flow through the pump (positive in normal conditions) [kg/s].
     * @param deltaP Pressure loss through the pump (negative in normal conditions, since the pump increases pressure) [Pa].
     * @param dDeltaP_dm Derivative of pressure loss through the pump [Pa/(kg/s)].
     */
    void computePressureDiff(float const& massFlow, float& deltaP, float& dDeltaP_dm);

    /**
     * @brief computeIdealN Solves the rpm needed to get the pressureDiff, and the mass flow. If the mass flow is negative, there is no solution, simply try increasing to (n+n0)/2.
     * @param massFlow Mass flow through the pump (positive in normal conditions) [kg/s].
     * @param pressureDiff Pressure loss through the pump (negative in normal conditions, since the pump increases pressure) [Pa].
     * @return The value of n needed to get the pressureDiff and massFlow  [rotation/min].
     */
    float computeIdealN(float const& massFlow, float const& pressureDiff);

    void updateRpms(float& sumDeltaRpm, float& sumRpm, float const& learningRate, float const& targetMassFlow, float const& targetPressureDiff);
    void setNToMax(float& sumDeltaRpm, float& sumRpm, float const& learningRate);

    bool nIsMin() { return (n==computeNMin()); }
    bool nIsAlmostMax() { return n>0.99*n0; }

};

class Valve {
private:
    static float deltaP0_invRho0_36002; // deltaP0=1e5, invRho0=1/1e3, 3600^2=12960000, multiply them. Computed only once.

    // Stays constant, define the Valve.
    float kvMax; // Maximal valve flow coefficient (max opening of valve) [m^3/h].

    // Evolves at each time step.
    float kv; // Valve flow coefficient (opening of valve) [m^3/h].

    float computeKvMin() { return kvMax*0.0001f; } // A fraction of the max. Modified by Max. The range of operation is reduced since for high mass flows it might give huge pressure differences from 1e-5 to 1e-4.

public:
    Valve(float const& kvMax_) : kvMax(kvMax_), kv(kvMax_*0.5f) { }
    virtual ~Valve() {}
    void computePressureDiffAndDerivative(float const& m, float const& rho, float& deltaP, float& dDeltaP_dm);
    float computeIdealKv(float const& rho, float const& massFlow, float const& pressureDiff);
    void updateKv(float const& rho, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& targetMassFlow, float const& targetPressureDiff, PIDControllerValve& pid, float& Targetkv, bool& ImposedValve);
    /**
     * @brief computeTemperatureIncrease
     * @param pressureDiff Pressure loss through the valve (positive in normal conditions, since the valve decreases pressure) [Pa].
     * @param cp Fluid heat capacity [J/(kg*K)].
     * @param rho Fluid mass density [kg/(m^3)].
     * @return Temperature increase (positive in normal conditions, since pressure loss due to friction increases temperature) [K]
     */
    float computeTemperatureIncrease(float const& pressureDiff, float const& cp, float const& rho) { return pressureDiff/(cp*rho); }
    void setKvToMax(float& sumDeltaKv, float& sumKv) { sumDeltaKv += abs(kv-kvMax); sumKv += kvMax; kv = kvMax; }
    void setKvToMin(float& sumDeltaKv, float& sumKv, float const& learningRate);
    void keepKvConst(float& sumDeltaKv, float& sumKv) { sumDeltaKv += 0.f; sumKv += kv; }
    virtual int nbEdges() { return 1; } // It's only a valve.
    float getkv(){return kv;}

    bool kvIsMin() { return kv==computeKvMin(); }
    bool kvIsAlmostMin() { return kv<1.5*computeKvMin(); }
};


class PIDController {
private:
    float integralErr;
    float prevErr;

protected:
    float computeControlVariable(float const& desiredSetpoint, float const& processVariable, float const& kp, float const& ki, float const& kd);

public:
    PIDController() : integralErr(0.f), prevErr(0.f) { }
    virtual ~PIDController() { }
};

class PIDControllerValve : public PIDController {
public:
    PIDControllerValve() : PIDController() { }
    ~PIDControllerValve() { }
    float computeControlVariable(float const& desiredSetpoint, float const& processVariable) { return PIDController::computeControlVariable(desiredSetpoint, processVariable, 0.5f, 0.3f, 0.f); }
};

class PIDControllerPump : public PIDController {
public:
    PIDControllerPump() : PIDController() { }
    ~PIDControllerPump() { }
    float computeControlVariable(float const& desiredSetpoint, float const& processVariable) { return PIDController::computeControlVariable(desiredSetpoint, processVariable, 0.3f, 0.2f, 0.f); }
};

class MemoryManager {
private:
    deque<float> q;
    size_t r;
public:
    MemoryManager(size_t r) : q(), r(r) { }
    ~MemoryManager() { }
    size_t nbMemorized() { return q.size(); }
    void addElement(float e) { q.push_front(e); if(q.size()>=r) { q.pop_back(); } }
    float findMin() { float min = q.front(); for (auto const& el : q) { if (el<min) { min = el; } } return min; }
    float findMedian() { vector<float> v (q.size()); for(size_t i=0; i<q.size(); i++) { v[i]=q[i]; } nth_element(v.begin(), v.begin()+v.size()/2, v.end()); return v[v.size()/2]; }
};

class Carla {
private:
    vector<float> position_x;
    vector<float> probaDensity_fx;
    float actionChoice_r;

    float xmin() { return position_x.front(); }
    float xmax() { return position_x.back(); }
    float integrate(vector<float> const& f, float const& xmin, float const& xmax);
    float selectAction();
    float solveQuadratic(float const& xi, float const& xip1, float const& fi, float const& fip1, float const& area);

public:
    Carla(size_t const& n, float const& xmin, float const& xmax);
    ~Carla() { }
    float step(float const& beta);

};

class PIDControllerCarla : public PIDController {
private:
    Carla carlaKp;
    Carla carlaKi;
//    Carla carlaKd;
    MemoryManager mm;

public:
    PIDControllerCarla(size_t const& n, size_t const& r, float const& xmin, float const& xmax);
    ~PIDControllerCarla() { }
    float evaluatePerformanceBeta(float const& costJ);
    float computeControlVariable(float const& desiredSetpoint, float const& processVariable);
};


//beginning of contents added by Dapeng
class Substation : virtual public EnergyConversionSystem {
protected:

    // Stays constant, define the substation.
    float designThermalPower; // Must be positive. [W]
    float designTempDifference; // Must be positive. [degree C]
    float designEpsilon; // Must be positive and smaller or equal to one. []

    DistrictEnergyCenter* pDEC; // The DEC that the substation is linked to.
    SubstationNodePair* pNode; // The node where the substation links to in the DEC network.
    Building* pBuilding; // The building that the substation is linked to.


    // Evolve at each time step.
    float desiredMassFlow; // on primary side [kg/s]
    float secondarySideInputTemp;
    float thermalPowerExchanged; // Heat exchanged at heat exchanger. Positive means primary side gives heat to secondary, if everything is ok, when the substation consumes, it should reach thermalPowerNeeded [W].
    float primarySideOutputTemp; // Temperature before mixing with other fluxes, downtream from the substation heat exchanger [degree C].

    // Heat exchanged over the thermal intervals of one hydraulic step [J] and their duration [s].
    double exchangedEnergyOverStep = 0.;
    double exchangeDurationOverStep = 0.;

    PIDControllerValve pid; // To control the valve opening, to get the desired mass flow.
//    PIDControllerCarla pid;

    Valve* valve;

    // Done in order to have the case where the valve is controlled manually
    float Targetkv; // Added by Max
    bool ImposedValve; // Added by Max

    // Control directly the mass flow instead of going through the valve
    bool ImposedMassFlow;
    MassFlowSetpoint* massFlowSetpoint;

    void setPrimarySideOutputTemp(float outputTemp){primarySideOutputTemp=outputTemp;}

public:
    static Substation* createNewSubstation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    Substation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    ~Substation() override { delete valve; /*logStream << "Destructor of Substation" << endl;*/ }

    void writeXML(ofstream& file, string tab){ file << tab << "Substation saving not supported yet" << endl; }
    string getLabel() { return "Substation"; }
    float getPrimarySideOutputTemp() { return primarySideOutputTemp; }
    float getDesiredMassFlow() { return desiredMassFlow; }

    virtual void setThermalPowerNeeded(double tPN) override;

    /**
     * Returns value of thermal power provided to the user, must compute it beforehand using computeHeatExchanged().
     * @param thermalPowerNeeded (for other EnergyConversionSystems class, eventually change this)
     * @return Thermal power that the substation can provide to the user.
     */
    virtual double getThermalPower(double sourceTemp) override;
    virtual double getFuelConsumption(double time, double thermalPower, double sourceTemp) override;
    virtual double getElectricConsumption(double time, double thermalPower, double sourceTemp) override;

    /**
     * @brief mc Formula to compute primary network side mass flow rate.
     * @param thermalPowerNeeded (positive or negative) [W]
     * @param designThermalPower (always positive) [W]
     * @param m_c_n Design mass flow [kg/s]
     * @return Primary network side mass flow rate [kg/s]
     */
    virtual double mc(double thermalPowerNeeded, double designThermalPower, double m_c_n);

    Building* getBuilding() { return pBuilding; }
    void setDEC(DistrictEnergyCenter* dec) { pDEC = dec; }
    void setNode(SubstationNodePair* node) { pNode = node; }

    /**
     * Computes the heat exchanged for a given primary side mass flow and input temperature. So computes thermalPowerExchanged and primarySideOutputTemp.
     */
    virtual void computeHeatExchanged(float const& primarySideCp, float const& primarySideRho, float const& primarySideInputTemp, float const& primarySideMassFlow, float const& primarySidePressureDiff, Climate* pClim, unsigned int day, unsigned int hour);
    virtual float maxKv(float const& rho);
    virtual void updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& massFlowSupplyToReturn);

    float relativeErrorMassFlow(float const& massFlow) { return (massFlow-desiredMassFlow)/desiredMassFlow; }
    virtual void errorMassFlow(float const& massFlow, float& relErr, float& absErr, float& sumErr);
    virtual float computeDesignMassFlow();

    /**
     * Gives the number of graph edges that the substation is made up of. Eg 1 edge if the substation is only a valve, 3 edges if the substation is a valve, a pipe and a differential pressure regulator.
     */
    virtual int nbEdges() { return 1; } // A plain substation has only a valve.
    virtual bool hasRegEle() { return false; } // A plain substation has only a valve.
    virtual float computePressureDiffControlElement(float pressureDiffSubstationNodePair) { return pressureDiffSubstationNodePair; } // All pressure loss is in the control element (valve).
    virtual void computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm);
    float computeOutputTemp(float const& inputTemp, float const& thermalPowerExchanged, float const& massFlow, float const& pressureDiff, float const& rho, float const& cp);

    virtual void updateDesiredMassFlow(float const& cp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour); // Only prosumers need this information, not substations that only consume.

    /**
     * The dynamic thermal models evaluate the transfer once per thermal interval. These average it
     * over the hydraulic step, so the building and the hourly records get the heat exchanged over
     * the step rather than the value of its last interval.
     */
    void resetTransportAverage() { exchangedEnergyOverStep = 0.; exchangeDurationOverStep = 0.; }
    void accumulateTransportAverage(double dt) { exchangedEnergyOverStep += thermalPowerExchanged * dt; exchangeDurationOverStep += dt; }
    void applyTransportAverage() { if (exchangeDurationOverStep > 0.) { thermalPowerExchanged = static_cast<float>(exchangedEnergyOverStep / exchangeDurationOverStep); } }

    virtual void recordTimeStep() { }
    virtual void eraseRecords(unsigned int keepValue) { }
    virtual void eraseRecords_back() { }

    virtual void writeTHHeaderText(fstream& textFile, string prefix) { }
    virtual void writeTHResultsText(fstream& textFile, unsigned int i) { }

};
//end of contents added by Dapeng



class RegulatingElement {
private:
    float targetRegulatedPathPressureDiff; // [Pa]
    float pressureDiff; // [Pa]
public:
    RegulatingElement(TiXmlHandle hdl);
    ~RegulatingElement() { }
    float getTargetRegulatedPathPressureDiff() { return targetRegulatedPathPressureDiff; }
    float getPressureDiff() { return pressureDiff; }
    void setPressureDiff(float const& p) { pressureDiff = p; }
    void subtractPressure(float const& press) { pressureDiff -= press; }
};



class RegulatedPressureSubstation : virtual public Substation {
private:
    RegulatingElement* regEle;
public:
    RegulatedPressureSubstation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    ~RegulatedPressureSubstation() { delete regEle; }
    virtual float maxKv(float const& rho) override { return 3600.f*computeDesignMassFlow()*1.1/(rho*sqrt(regEle->getTargetRegulatedPathPressureDiff()*0.00001f)); } // At the target pressure, the max valve opening can let through the nominal mass flow times 1.1 (this factor is added in case of imprecisions).
    virtual int nbEdges() override { return 2; } // There is a valve and a differential pressure regulator.
    virtual bool hasRegEle() override { return true; } // There is a valve and a differential pressure regulator.
    virtual float computePressureDiffControlElement(float pressureDiffSubstationNodePair) override { return pressureDiffSubstationNodePair - regEle->getPressureDiff(); } // Total pressure loss = valve + regulating element pressure losses. (the control element is the valve)
    virtual void computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm) override;
};




class ProsumerSubstation : public Substation {
private:
    TemperatureSetpoint* temperatureSetpoint; // Determines the target temperature.
    PressureSetpoint* pressureSetpoint; // Determines the target pressure.

    Pump* pump;
    Valve* pumpFlowControlValve;
    bool producerModeOn;
    float pumpElectricPower;

//    bool turnOfPumpFlowControlValve;

    void deleteDynAllocated();
    bool isHeatSource();
    float computeSolarThermalPower(float const& targetSupplyTemp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour);
    float minDesiredPumpFlow();

public:
    ProsumerSubstation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    ~ProsumerSubstation() { deleteDynAllocated(); }

    virtual void setThermalPowerNeeded(double tPN) override;
    virtual double getThermalPower(double sourceTemp) override;
    virtual double getElectricConsumption(double time, double thermalPower, double sourceTemp) override;

    virtual void updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& massFlowSupplyToReturn) override;

    virtual void computeHeatExchanged(float const& primarySideCp, float const& primarySideRho, float const& primarySideInputTemp, float const& primarySideMassFlow, float const& primarySidePressureDiff, Climate* pClim, unsigned int day, unsigned int hour) override;
    virtual void computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm) override;
    virtual void updateDesiredMassFlow(float const& cp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour) override;
    void setProsumerSolarThermal();

    virtual void errorMassFlow(float const& massFlow, float& relErr, float& absErr, float& sumErr) override;

    bool getProducerMode(){return producerModeOn;}
    void setProducerMode(bool power) {producerModeOn = power;}

    virtual void recordTimeStep() override { Substation::recordTimeStep(); setProsumerSolarThermal(); } // setProsumerSolarThermal
    virtual void eraseRecords(unsigned int keepValue) override { Substation::eraseRecords(keepValue); }
    virtual void eraseRecords_back() override { Substation::eraseRecords_back(); }

};


// Added by Max
class SubstationHeatPump: public ProsumerSubstation {
private:
    HeatPump* heatPump;
public:
    SubstationHeatPump(TiXmlHandle hdl, Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    ~SubstationHeatPump() override { /*logStream << "Destructor of SubstationHP" << endl;*/ }
    void writeXML(ofstream& file, string tab) override{ file << tab << "SubstationHeatPump saving not supported yet" << endl; }
    string getLabel() override{ return "SubstationHeatPump"; }
    virtual double getThermalPower(double sourceTemp) override;
    virtual double getElectricConsumption(double time, double thermalPower, double sourceTemp) override;
    HeatPump* getHeatPump(){return heatPump;};
    //double mc(double thermalPowerNeeded, double sourceTemp);
    float computeDesignMassFlow() override;
    void setThermalPowerNeeded(double tPN) override;
    void updateDesiredMassFlow(float const& cp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour) override;
    virtual void computeHeatExchanged(float const& primarySideCp, float const& primarySideRho, float const& primarySideInputTemp, float const& primarySideMassFlow, float const& primarySidePressureDiff, Climate* pClim, unsigned int day, unsigned int hour) override;
    double logMeanTemperatureDifference(float Tin,float Tout);
};
//Added by Max
class SubstationHeatPump2stages: public SubstationHeatPump{
    double targetTemp2; // The target temperature for heating the surface.
    double heatPumpElectricPower2;
    double etaTech2;
    // We suppose etatech is identical for both stages.
public:
    SubstationHeatPump2stages(TiXmlHandle hdl, Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr = &std::cout);
    void writeXML(ofstream& file, string tab) override{ file << tab << "SubstationHP2stages saving not supported yet" << endl; }
    string getLabel() override{ return "SubstationHP2stages"; }
    double getThermalPower(double sourceTemp) override;
    double getElectricConsumption(double time, double thermalPower, double sourceTemp) override;
    double getHeatProduced(double workDHW, double workHS, double sourceTemp, double outputTemp);
    vector<double> getWorkNeeded(double thermalPower, double sourceTemp, double outputTemp);
    virtual void computeHeatExchanged(float const& primarySideCp, float const& primarySideRho, float const& primarySideInputTemp, float const& primarySideMassFlow, float const& primarySidePressureDiff, Climate* pClim, unsigned int day, unsigned int hour) override;
    virtual void updateDesiredMassFlow(float const& cp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour) override;
};

class TemperatureSetpoint {
public:
    static TemperatureSetpoint* createNewTemperatureSetpoint(TiXmlHandle hdl);
    TemperatureSetpoint() { }
    virtual ~TemperatureSetpoint() { }
    virtual float computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) = 0;
    virtual float getInitTemp() = 0;
};

class ConstantTemperatureSetpoint: public TemperatureSetpoint {
private:
    float targetSupplyTemp;
public:
    ConstantTemperatureSetpoint(TiXmlHandle hdl);
    ~ConstantTemperatureSetpoint() override { }
    virtual float computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) override { return targetSupplyTemp; }
    virtual float getInitTemp() override { return targetSupplyTemp; }
};

class AffineTemperatureSetpoint : public TemperatureSetpoint {
protected:
    float lowExtTemp, highExtTemp, lowExtTempSupplyTemp, highExtTempSupplyTemp;
public:
    AffineTemperatureSetpoint(TiXmlHandle hdl);
    ~AffineTemperatureSetpoint() override { }
    virtual float computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) override;
    virtual float getInitTemp() override { return highExtTempSupplyTemp; }
    float avgExtTempLast24Hours(Climate* pClimate, unsigned int day, unsigned int hour);
};

class AffineWinterConstantSummerSetpoint : public AffineTemperatureSetpoint {
private:
    bool notYetActivated;
    bool summerModeOn;
    float startSummerTempThreshold, endSummerTempThreshold;
public:
    AffineWinterConstantSummerSetpoint(TiXmlHandle hdl);
    ~AffineWinterConstantSummerSetpoint() override { }
    virtual float computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) override;
};

class ImposedValuesOrConstantSetpoint : public TemperatureSetpoint {
private:
    map<string, float> imposedValues; // Stores values if one wants to impose the supply temperature. Format map["d42h4"]=75.
    float constantTempIfNoImposed;

    bool hasImposedValue(unsigned int day, unsigned int hour, float& retValue);
public:
    ImposedValuesOrConstantSetpoint(TiXmlHandle hdl);
    ~ImposedValuesOrConstantSetpoint() override { }
    virtual float computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour);
    virtual float getInitTemp() { return constantTempIfNoImposed; }
};





class PressureSetpoint {
public:
    static PressureSetpoint* createNewPressureSetpoint(TiXmlHandle hdl);
    PressureSetpoint() { }
    virtual ~PressureSetpoint() { }
    virtual float computeTargetPressureDiff(float const& massFlow) = 0; // By convention they are negative when the pump increases the pressure (normal case).
};

class ConstantPressureSetpoint: public PressureSetpoint {
private:
    float targetPressureDiff; // Usually negative for pumps.
public:
    ConstantPressureSetpoint(TiXmlHandle hdl);
    ~ConstantPressureSetpoint() override { }
    virtual float computeTargetPressureDiff(float const& massFlow) override { return targetPressureDiff; }
};

class AffinePressureSetpoint : public PressureSetpoint {
protected:
    // The setpoint function of the pressure difference is a function of mass flow, affine by parts.
    vector<float> massFlows; // [kg/s]
    vector<float> pressureDiffs; // [Pa]
public:
    AffinePressureSetpoint(TiXmlHandle hdl);
    ~AffinePressureSetpoint() override { }
    virtual float computeTargetPressureDiff(float const& massFlow) override;
};

//Added by Max
class MassFlowSetpoint {
public:
    static MassFlowSetpoint* createNewMassFlowSetpoint(TiXmlHandle hdl);
    MassFlowSetpoint() { }
    virtual ~MassFlowSetpoint() { }
    virtual double computeTargetMassFlow() = 0;
};

class ConstantMassFlowSetPoint: public MassFlowSetpoint {
private:
    double massFlow; // must be positive.
public:
    ConstantMassFlowSetPoint(TiXmlHandle hdl);
    ~ConstantMassFlowSetPoint() override { }
    virtual double computeTargetMassFlow() override { return massFlow; }
};


class Storage {
public:
    static Storage* createNewStorage(TiXmlHandle hdl);
    Storage() { }
    virtual ~Storage() { }
    virtual float computeOutputTemperature(bool storageHeatsUp, float const& inputTemp, float const& m, float const& cp) = 0;
    virtual void confirmStoredHeat(float const& tempDiff, float const& m, float const& cp) = 0;
    virtual void recordTimeStep() = 0;
    virtual void eraseRecords(unsigned int keepValue) = 0;
    virtual void eraseRecords_back() = 0;
    virtual void writeTHHeaderText(fstream& textFile, string prefix) = 0;
    virtual void writeTHResultsText(fstream& textFile, unsigned int i) = 0;
};



class SimpleStorage : public Storage {
private:
    float temperature; // [degree C]
    float heatCapacity; // [J/K]

    // Stored information about simulation.
    vector<float> temperatureRecord;

    void recordTemperature() { temperatureRecord.push_back(temperature); }
    float getTemperature(unsigned int step) {return temperatureRecord.at(step); }
    void eraseTemperature(unsigned int keepValue) { temperatureRecord.erase(temperatureRecord.begin(),temperatureRecord.end()-min(keepValue,(unsigned int)temperatureRecord.size())); }
    void eraseTemperature_back() { temperatureRecord.pop_back(); }

public:
    SimpleStorage(TiXmlHandle hdl);
    ~SimpleStorage() override { }
    virtual float computeOutputTemperature(bool storageHeatsUp, float const& inputTemp, float const& m, float const& cp) override;
    virtual void confirmStoredHeat(float const& tempDiff, float const& m, float const& cp) override;

    virtual void recordTimeStep() override { recordTemperature(); }
    virtual void eraseRecords(unsigned int keepValue) override { eraseTemperature(keepValue); }
    virtual void eraseRecords_back() override { eraseTemperature_back(); }
    virtual void writeTHHeaderText(fstream& textFile, string prefix) override;
    virtual void writeTHResultsText(fstream& textFile, unsigned int i) override;

};




class ThermalStation {
protected:
    Pump* pump=nullptr;

    vector<EnergyConversionSystem*> ecs;
    vector<float> thermalPowerProvided;
    ThermalStationNodePair* node; // Node where in connects to the Network.

    // Changes at every time step
    float pumpPower; // Power consumed by pump (usually electric) [W].
    float outputTemperature; // Temperature that exits / is outputed by the thermal station [degree C].
    TemperatureSetpoint* temperatureSetpoint=nullptr; // Determines the target temperature.
    PressureSetpoint* pressureSetpoint=nullptr; // Determines the target pressure.
//    PIDControllerPump pid;
//    PIDControllerCarla pid;

    // Stored information about simulation.
    vector<float> pumpPowerRecord;
    vector<vector<float>> ecs_thermalPowerProvided, ecs_electricConsumption, ecs_fuelConsumption;

    //Added by Max
    unsigned int linkedNodeId;
    float thermalPowerNeeded=0.f;

    // Energies over the thermal intervals of one hydraulic step [J] and their duration [s].
    vector<double> providedEnergyOverStep;
    double pumpEnergyOverStep = 0., neededEnergyOverStep = 0., durationOverStep = 0.;

    void deleteDynAllocated();

    // Pump Power
    void recordPumpPower() { pumpPowerRecord.push_back(pumpPower); }
    float getPumpPower(unsigned int step) {return pumpPowerRecord.at(step); }
    void erasePumpPower(unsigned int keepValue) { pumpPowerRecord.erase(pumpPowerRecord.begin(),pumpPowerRecord.end()-min(keepValue,(unsigned int)pumpPowerRecord.size())); }
    void erasePumpPower_back() { pumpPowerRecord.pop_back(); }

    // Electric Consumption in J
    float getElectricConsumption_back() {
        float value = pumpPowerRecord.back()*static_cast<double>(Model::dt); // to convert the W in J
        for (size_t i=0;i<ecs.size();++i)
            value += ecs_electricConsumption[i].back();
        return value;
    }
    float getElectricConsumption(unsigned int step) {
        float value = pumpPowerRecord.at(step)*static_cast<double>(Model::dt); // to convert the W in J
        for (size_t i=0;i<ecs.size();++i)
            value += ecs_electricConsumption[i].at(step);
        return value;
    }
    void eraseElectricConsumption(unsigned int keepValue) {
        for (size_t i=0;i<ecs.size();++i)
            ecs_electricConsumption[i].erase(ecs_electricConsumption[i].begin(),ecs_electricConsumption[i].end()-min(keepValue,(unsigned int)ecs_electricConsumption[i].size()));
    }
    void eraseElectricConsumption_back() {
        for (size_t i=0;i<ecs.size();++i)
            ecs_electricConsumption[i].pop_back();
    }

    // Fuel Consumption
    float getFuelConsumption_back() {
        float value = 0.f;
        for (size_t i=0;i<ecs.size();++i)
            value += ecs_fuelConsumption[i].back();
        return value;
    }
    float getFuelConsumption(unsigned int step) {
        float value = 0.f;
        for (size_t i=0;i<ecs.size();++i)
            value += ecs_fuelConsumption[i].at(step);
        return value;
    }
    void eraseFuelConsumption(unsigned int keepValue) {
        for (size_t i=0;i<ecs.size();++i)
            ecs_fuelConsumption[i].erase(ecs_fuelConsumption[i].begin(),ecs_fuelConsumption[i].end()-min(keepValue,(unsigned int)ecs_fuelConsumption[i].size()));
    }
    void eraseFuelConsumption_back() {
        for (size_t i=0;i<ecs.size();++i)
            ecs_fuelConsumption[i].pop_back();
    }

    // Machine Power
    float getMachinePower_back() {
        float value = 0.f;
        for (size_t i=0;i<ecs.size();++i)
            value += ecs_thermalPowerProvided[i].back();
        return value;
    }
    float getMachinePower(unsigned int step) {
        float value = 0.f;
        for (size_t i=0;i<ecs.size();++i)
            value += ecs_thermalPowerProvided[i].at(step);
        return value;
    }
    void eraseMachinePower(unsigned int keepValue) {
        for (size_t i=0;i<ecs.size();++i)
            ecs_thermalPowerProvided[i].erase(ecs_thermalPowerProvided[i].begin(),ecs_thermalPowerProvided[i].end()-min(keepValue,(unsigned int)ecs_thermalPowerProvided[i].size()));
    }
    void eraseMachinePower_back() {
        for (size_t i=0;i<ecs.size();++i)
            ecs_thermalPowerProvided[i].pop_back();
    }
    void setMachinePower_FuelConsumption_ElectricConsumption(Climate* pClim, unsigned int day, unsigned int hour);


public:
    ostream logStream; // Cognet: Added, to do like other classes.

    static ThermalStation* createNewThermalStation(TiXmlHandle hdl, Network* net, ostream* pLogStr);

    ThermalStation(TiXmlHandle hdl, Network* net, ostream* pLogStr);
    virtual ~ThermalStation() { deleteDynAllocated(); }

    Pump* getPump() { return pump; }

    EnergyConversionSystem* getEcs() { return ecs[0]; } /// NOTE: temporary set to ecs[0]
    float getThermalPowerProvided() { return accumulate(thermalPowerProvided.begin(), thermalPowerProvided.end(), 0.f); }
    float getOutputTemperature() { return outputTemperature; }

    virtual void computeThermal(float pressureDiff, float massFlow, float rho, float cp, float inputTemp, Climate *pClimate, unsigned int day, unsigned int hour);

    virtual void computePressureDiff(float const& massFlow, float const& rho, float& deltaP, float& dDeltaP_dm);
    float computeTargetPressureDiff(float const& massFlow) { return pressureSetpoint->computeTargetPressureDiff(massFlow); }
    virtual void updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate);
    //Added by Max
    virtual void updateStage(float const& totalThermalPowerNeeded, float const& totalThermalPowerMax, float prevCurrentStage, bool prevSlaveShutDown, float rho, float cp, float sourceTemp, Climate* pClimate, unsigned int day, unsigned int hour){throw string("CitySim is trying to update the stage from a non slave ThermalStation.");}
    unsigned int getLinkedNodeId(){return linkedNodeId;}
    virtual void updateDesiredMassFlow(float cp, Climate* pClimate, unsigned int day, unsigned int hour){};
    virtual vector<float> getRampingStages(){throw string("CitySim is trying to get a rampingStages from a non slave ThermalStation.");return {-1.f};}
    virtual float getCurrentStage(){throw string("CitySim is trying to get a current stage from a non slave ThermalStation."); return 0.f;}
    virtual void setCurrentStage(float stage){throw string("CitySim is trying to set a current stage to a non slave ThermalStation.");}
    virtual bool getMaster(){return true;} // Added by Max (by convention, the thermalStation is a master. Test in MCR that forces to have thermalStation master or slave if they are several).
    float getThermalPowerNeeded() const { return thermalPowerNeeded; }

    /** Same role as Substation::applyTransportAverage, for the powers recorded by the station. */
    void resetTransportAverage() {
        providedEnergyOverStep.assign(thermalPowerProvided.size(), 0.);
        pumpEnergyOverStep = 0.; neededEnergyOverStep = 0.; durationOverStep = 0.;
    }
    void accumulateTransportAverage(double dt) {
        if (providedEnergyOverStep.size() != thermalPowerProvided.size()) { providedEnergyOverStep.resize(thermalPowerProvided.size(), 0.); }
        for (size_t i = 0; i < thermalPowerProvided.size(); ++i) { providedEnergyOverStep[i] += thermalPowerProvided[i] * dt; }
        pumpEnergyOverStep += pumpPower * dt;
        neededEnergyOverStep += thermalPowerNeeded * dt;
        durationOverStep += dt;
    }
    void applyTransportAverage() {
        if (durationOverStep <= 0.) { return; }
        for (size_t i = 0; i < thermalPowerProvided.size() && i < providedEnergyOverStep.size(); ++i) {
            thermalPowerProvided[i] = static_cast<float>(providedEnergyOverStep[i] / durationOverStep);
        }
        pumpPower = static_cast<float>(pumpEnergyOverStep / durationOverStep);
        thermalPowerNeeded = static_cast<float>(neededEnergyOverStep / durationOverStep);
    }
    virtual unsigned int getTimeClock(){throw string("Error in the code: A timeClock is trying to be accessed in a ThermalStation of different type than 'slave'.");}
    virtual void setTimeClock(unsigned int count){throw string("Error in the code: A timeClock is trying to be fixed in a ThermalStation of different type than 'slave'.");}
    virtual unsigned int getLatency(){throw string("Error in the code: A Latency is trying to be accessed in a ThermalStation of different type than 'slave'.");}
    float getTotalPumpPower() {return std::accumulate(pumpPowerRecord.begin(),pumpPowerRecord.end(),0.0);}

    virtual void errorPressureDiff(float const& massFlow, float const& pressureDiff, float& relErr, float& absErr,float& sumErrP, float& sumDesiredP);

    virtual void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& massFlow, float const& cp) { setMachinePower_FuelConsumption_ElectricConsumption(pClim, day, hour); recordPumpPower(); } // setMachinePower_FuelConsumption_ElectricConsumption setPumpPowers  confirmStorage
    virtual void eraseRecords(unsigned int keepValue) { eraseMachinePower(keepValue); eraseFuelConsumption(keepValue); eraseElectricConsumption(keepValue); erasePumpPower(keepValue); }
    virtual void eraseRecords_back() { eraseMachinePower_back(); eraseFuelConsumption_back(); eraseElectricConsumption_back(); erasePumpPower_back(); }


    virtual void updateOperationMode(float const& sumSubstationDesiredMassFlow, size_t const& nbThermalStations) { }
    void computeThermalPowerProvided(Climate* pClimate, unsigned int day, unsigned int hour);

    virtual void writeTHHeaderText(fstream& textFile, string prefix);
    virtual void writeTHResultsText(fstream& textFile, unsigned int i);
};

// Added by Max
class ThermalStationMaster: public ThermalStation{ // One thermal station to rule them all.
public:
    ThermalStationMaster(TiXmlHandle hdl, Network* net, ostream* pLogStr):ThermalStation(hdl,net,pLogStr){}
    bool getMaster(){return true;}
};

//Added by Max
class ThermalStationSlave: public ThermalStation{
private:
    float desiredMassFlow; // It is useful for the slave.

    Valve* valve;
    PIDControllerValve pid;
    float Targetkv;
    bool ImposedValve;

    vector<float> rampingStages; // [%]
    float currentStage;

    unsigned int latency; // [hour]
    unsigned int timeClock; //hour

    void deleteDynAllocated() { if (valve) delete valve; }

public:
    ThermalStationSlave(TiXmlHandle hdl, Network* net, ostream* pLogStr);
    ~ThermalStationSlave() override { deleteDynAllocated(); }

    void computePressureDiff(float const& massFlow, float const& rho, float& deltaP, float& dDeltaP_dm);
    void setDesiredMassFlow(float rho, float cp, double sourceTemp, Climate* pClimate, unsigned int day, unsigned int hour);
    void updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate);
    bool getMaster(){return false;}
    vector<float> getRampingStages(){return rampingStages;}
    float getCurrentStage(){return currentStage;}
    void setCurrentStage(float stage){currentStage=stage;}
    void errorPressureDiff(float const& massFlow, float const& pressureDiff, float& relErr, float& absErr, float& sumErrP, float& sumDesiredP);
    unsigned int getLatency(){return latency;}
    unsigned int getTimeClock(){return timeClock;}
    void setTimeClock(unsigned int count){timeClock=count;}
    void updateStage(float const& totalThermalPowerNeeded, float const& totalThermalPowerMax, float prevCurrentStage, bool prevSlaveShutDown, float rho, float cp, float sourceTemp, Climate* pClimate, unsigned int day, unsigned int hour);
};

class SeasonalStorageHeatingThermalStation : public ThermalStation { // Thermal station that only heats, and that has thermal storage. Goes with prosumer thermal stations.
private:
    Storage* storage;
    Valve* valve;
    bool storageModeOn;
    float tempDiffAroundStorage; // Temperature leaving storage - temperature going to storage [K].
    float flowToStore;

    PIDControllerValve pid;
    float Targetkv;
    bool ImposedValve;

    void deleteDynAllocated() { if (storage!=nullptr) { delete storage; } if (valve!=nullptr) { delete valve; } }

public:
    SeasonalStorageHeatingThermalStation(TiXmlHandle hdl, Network* net, ostream* pLogStr);
    ~SeasonalStorageHeatingThermalStation() override { deleteDynAllocated(); }

    void computePressureDiff(float const& massFlow, float const& rho, float& deltaP, float& dDeltaP_dm) override;

    void updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate) override;

    void computeThermal(float pressureDiff, float massFlow, float rho, float cp, float inputTemp, Climate *pClimate, unsigned int day, unsigned int hour);

    void updateOperationMode(float const& sumSubstationDesiredMassFlow, size_t const& nbThermalStations) override;

    void confirmStorage(float const& massFlow, float const& cp) { storage->confirmStoredHeat(tempDiffAroundStorage, abs(massFlow), cp); }

    void errorPressureDiff(float const& massFlow, float const& pressureDiff, float& relErr, float& absErr, float& sumErrP, float& sumDesiredP) override;

    void writeTHHeaderText(fstream& textFile, string prefix) override;
    void writeTHResultsText(fstream& textFile, unsigned int i) override;


    void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& massFlow, float const& cp) override { ThermalStation::recordTimeStep(pClim, day, hour, massFlow, cp);  confirmStorage(massFlow, cp); storage->recordTimeStep(); }
    void eraseRecords(unsigned int keepValue) override { ThermalStation::eraseRecords(keepValue); storage->eraseRecords(keepValue); }
    void eraseRecords_back() override { ThermalStation::eraseRecords_back(); storage->eraseRecords_back(); }
};

class MCR{
protected:
    vector<ThermalStation*> thermalstations; // vector of external references to ThermalStation objects created elsewhere

public:
    static MCR* createNewMCR(TiXmlHandle hdl, vector<ThermalStation*> thermalStations);
    MCR(TiXmlHandle hdl, vector<ThermalStation*> thermalStations);
    virtual ~MCR() {}
    virtual void MasterControlSlave(float rho, float cp, Climate* pClimate, unsigned int day, unsigned int hour, vector<float> prevCurrentStage){};
    virtual void MasterControlSlave2(float rho, float cp, Climate* pClimate, unsigned int day, unsigned int hour, size_t index){};
    void getSourceTemp(double& sourceTemp, Climate* pClimate, EnergyConversionSystem* ecs, unsigned int day, unsigned int hour);
    vector<ThermalStation*> getThermalStations(){return thermalstations;}
    virtual void updateThermalStationSlaveClock(vector<float> PrevCurrentStages){};

};

class SingleMCR: public MCR{
public:
    SingleMCR(TiXmlHandle hdl, vector<ThermalStation*> thermalStations);
};

class SimpleMCR:public MCR{
public:
    SimpleMCR(TiXmlHandle hdl, vector<ThermalStation*> thermalStations);
    virtual void MasterControlSlave(float rho, float cp, Climate* pClimate, unsigned int day, unsigned int hour, vector<float> prevCurrentStage);
    virtual void MasterControlSlave2(float rho, float cp, Climate* pClimate, unsigned int day, unsigned int hour, size_t index);
    void updateStageSlave(float const& desiredStage, size_t index, float prevCurrentStage);
    void updateThermalStationSlaveClock(vector<float> PrevCurrentStages);

};


struct WaterVolume {
    float temperature;
    float volume;
    //float position;

    WaterVolume() : temperature(0.f), volume(0.f) {}
    WaterVolume(float temperature_, float volume_) :
        temperature(temperature_),
        volume(volume_) {}
};

/**
 * @brief Ring buffer for WaterVolume elements used in the Lagrangian plug-flow model.
 *
 * Replaces std::deque<WaterVolume> with a pre-allocated circular array that provides:
 *   - O(1) push/pop at both ends without heap allocation during simulation
 *   - Cache-friendly, contiguous memory iteration via the two-segment interface
 *     (seg1Begin/seg1Len/seg2Len), which enables compiler auto-vectorisation (SIMD)
 *
 * Physical layout convention (same as the old deque):
 *   front() == tailNode side of the pipe  (upstream for massFlow > 0)
 *   back()  == headNode side of the pipe  (upstream for massFlow < 0)
 * This ordering is preserved across flow reversals, so no reordering is needed.
 *
 * The backing vector grows by doubling when full (amortised O(1) push). After a
 * grow() the elements are compacted so headIdx == 0, restoring the single-segment layout.
 */
struct WaterVolumeBuffer {
    vector<WaterVolume> data;
    unsigned int headIdx = 0;
    unsigned int sz      = 0;

    WaterVolumeBuffer() { data.resize(16); }

    unsigned int size()     const { return sz; }
    bool         empty()    const { return sz == 0; }
    unsigned int capacity() const { return static_cast<unsigned int>(data.size()); }

    WaterVolume&       front()       { return data[headIdx]; }
    const WaterVolume& front() const { return data[headIdx]; }
    WaterVolume&       back()        { return data[(headIdx + sz - 1u) % capacity()]; }
    const WaterVolume& back()  const { return data[(headIdx + sz - 1u) % capacity()]; }

    WaterVolume& operator[](unsigned int i) { return data[(headIdx + i) % capacity()]; }
    const WaterVolume& operator[](unsigned int i) const { return data[(headIdx + i) % capacity()]; }

    void push_front(const WaterVolume& wv) {
        if (sz == capacity()) grow();
        headIdx = (headIdx == 0u ? capacity() : headIdx) - 1u;
        data[headIdx] = wv;
        ++sz;
    }
    void push_back(const WaterVolume& wv) {
        if (sz == capacity()) grow();
        data[(headIdx + sz) % capacity()] = wv;
        ++sz;
    }
    void pop_front() { headIdx = (headIdx + 1u) % capacity(); --sz; }
    void pop_back()  { --sz; }
    void clear() { headIdx = 0u; sz = 0u; }

    // Two-segment interface for contiguous, SIMD-friendly iteration.
    // Segment A: data[seg1Begin() .. seg1Begin()+seg1Len()-1]
    // Segment B: data[0 .. seg2Len()-1]  (non-empty only when the buffer wraps)
    unsigned int seg1Begin() const { return headIdx; }
    unsigned int seg1Len()   const { return min(sz, capacity() - headIdx); }
    unsigned int seg2Len()   const { return sz - min(sz, capacity() - headIdx); }

private:
    void grow() {
        unsigned int newCap = capacity() * 2u;
        vector<WaterVolume> tmp(newCap);
        for (unsigned int i = 0; i < sz; i++) tmp[i] = data[(headIdx + i) % capacity()];
        data    = std::move(tmp);
        headIdx = 0;
    }
};


/**
 * @brief The Pipe class represents a pipe, which is part of a PipePair a DEC.
 */
class Pipe {
    friend class PipePair;

    private:

        // Dubon add
        float innerRadius; // autant descendre tout ici
        float innerVolume;
        float length;

        float wallCp = Model::defaultPipeWallCp;
        float wallRho = Model::defaultPipeWallRho;
        float wallThickness = 0.f;
        float wallkValue = Model::defaultPipeWallThermalConductivity;
        float insulationCp = Model::defaultPipeInsulationCp;
        float insulationRho = Model::defaultPipeInsulationRho;
        float insulationThickness = 0.f;
        float insulationkValue = Model::defaultPipeInsulationThermalConductivity;

        float casingThickness = 0.f;
        float casingkValue = Model::defaultPipeCasingThermalConductivity;
        float casingCp = Model::defaultPipeCasingCp;
        float casingRho = Model::defaultPipeCasingRho;

        float airConvectiveCoefficient = Model::defaultAirConvectiveCoefficient; // [W/m2/K] External heat transfer of air

        float buriedDepth; // [m] Pipe depth to center line of the pipe

        float thermalResistance; // Thermal resistance of pipe, from pipe wall to outside (soil or air). [K*m/W]
        float roughness;

        bool isAirPipe;

        // Evolve at each time step.
        float massFlow; // [kg/s] considered position if the flow goes from connectedNode[0] (tail) to connectedNode[1] (head).
        float pressureDiff; // [Pa] considered position if pressure drops from connectedNode[0] (tail) to connectedNode[1] (head).
        float thermalLoss; // [W]
        float downstreamTemp; // [degree C] the temperature just before the node where the flow outputs (the node temperature is a wheighted sum of incoming flows).
        float pressureLossHeat; // Heat from the pressure loss (always positive) per unit length of pipe [W/m]
        bool lastFlowWasPositive = true;

        unsigned int thermalSegmentCount = 1u;
        float thermalSegmentLength = 0.f;
        float thermalSegmentVolume = 0.f;
        vector<float> pipeWallTemperatureSegments; // [degree C] stationary pipe wall temperatures along the pipe.
        vector<float> casingTemperatureSegments; // [degree C] stationary casing temperatures along the pipe.
        vector<float> newPipeWallTemperatureSegments;
        vector<float> newCasingTemperatureSegments;
        vector<float> segmentSumGWaterWall;
        vector<float> segmentSumGWaterWallTimesWaterResponse;
        vector<float> segmentRhsPipeWallFromWater;
        vector<WaterVolume> waterVolumeScratch;

        // Stored information about simulation.
        vector<float> massFlowRecord; // Recorded values of the mass flow during the whole simulation.
        vector<float> pressureDiffRecord; // Recorded values of the pressure differences during the whole simulation.
        vector<float> thermalLossRecord; // Recorded values of the thermal losses during the whole simulation.

        // Dubon add — replaced deque with cache-friendly ring buffer (WaterVolumeBuffer)
        WaterVolumeBuffer waterVolumes; // Water volumes for the Lagrangian method
        // massGrid model: fixed equal-mass cells, built once from the first hydraulic solution.
        // Replaces waterVolumes for that model; its cells are also the wall segments.
        MassGridColumn massGrid;
        unsigned int massGridCellsPerSegment = 1u; // water cells per wall segment
        float massGridCellVolume = 0.f; // [m3]

        // Methods
        float computeSoilResistance(float outerPipeRadius, float soilThermalConductivity);
        float computeCylindricalResistance(float innerRadius, float outerRadius, float kValue);
        void applyHeatTransportToWaterVolumes(DistrictEnergyCenter* pDEC, float const& exteriorTemp, unsigned int dt);
        void applyThreeNodeHeatTransportToWaterVolumes(DistrictEnergyCenter* pDEC, float const& exteriorTemp, unsigned int dt);
        void advectWaterVolumes(float const& rhoWater, float const& tempIn, unsigned int dt);
        void initMassGrid(float rhoWater, unsigned int referenceDt);
        void applyMassGridHeatTransport(DistrictEnergyCenter* pDEC, float const& exteriorTemp, unsigned int dt);
        void advectMassGrid(float rhoWater, float tempIn, unsigned int dt);
        void advectPointerMerged(float rhoWater, float tempIn, unsigned int dt, float temperatureThreshold);

    public:
        Pipe(TiXmlHandle hdl, PipePair* parent, float const& soilThermalConductivity, float const& downstreamTemp, float const& airThermalConductivity);

        float getMassFlow() { return massFlow; }
        void setMassFlow(float mf) { massFlow = mf; }
        float getPressureDiff() { return pressureDiff; }
        void setPressureDiff(float dp) { pressureDiff = dp; }
        float getThermalLoss() { return thermalLoss; }
        void setThermalLoss(float tl) { thermalLoss = tl; }
        float getDownstreamTemp() { return downstreamTemp; }
        void setDownstreamTemp(float temp) { downstreamTemp = temp; }

        float getBuriedDepth() {return buriedDepth;} // Added by Max
        // int getNbWaterVolumes() { return waterVolumes.size();}

        // Dubon add
        float getPipeProperty(string propertyName);
        bool getIsAirPipe() { return isAirPipe; }
        void setRoughness(float value) { roughness = value; }

        /**
         * @brief pipeWall Computes the thermal resistance of insulation materials [K*m/W]
         * @param pipeRadius inner radius without the insulation
         * @return Insulation material resistance [K*m/W]
         */
        //float computeInsulationResistance(float pipeRadius);

        void recordThermalLoss() { thermalLossRecord.push_back(thermalLoss); }
        void eraseThermalLoss(unsigned int keepValue) { thermalLossRecord.erase(thermalLossRecord.begin(),thermalLossRecord.end()-min(keepValue,(unsigned int)thermalLossRecord.size())); }
        void eraseThermalLoss_back() { thermalLossRecord.pop_back(); }
        float getThermalLossRecord(unsigned int step)  { return thermalLossRecord.at(step); }

        void recordMassFlow() { massFlowRecord.push_back(massFlow); }
        void eraseMassFlow(unsigned int keepValue) { massFlowRecord.erase(massFlowRecord.begin(),massFlowRecord.end()-min(keepValue,(unsigned int)massFlowRecord.size())); }
        void eraseMassFlow_back() { massFlowRecord.pop_back(); }
        float getMassFlowRecord(unsigned int step)  { return massFlowRecord.at(step); }

        void recordPressureDiff() { pressureDiffRecord.push_back(pressureDiff); }
        void erasePressureDiff(unsigned int keepValue) { pressureDiffRecord.erase(pressureDiffRecord.begin(),pressureDiffRecord.end()-min(keepValue,(unsigned int)pressureDiffRecord.size())); }
        void erasePressureDiff_back() { pressureDiffRecord.pop_back(); }
        float getPressureDiffRecord(unsigned int step)  { return pressureDiffRecord.at(step); }


        void computeThermalLoss(float inputTemp, float twinNearInputTemp, float twinNearOutputTemp, float soilTemp, float cp, float length, float interPipeThermalResistance);
        void hydraulicConverged(float const& massFlow_, float const& deltaP, float const& rho, float const& length, float const& altitudePressureLoss);

        void injectAndAdvectWaterVolumes(DistrictEnergyCenter* pDEC, float const& exteriorTemp, float const& tempIn, unsigned int dt);

        unsigned int getMassGridCellCount() const { return massGrid.cells(); }
        uint64_t getPointerMergeBuildCount() const { return massGrid.mergedViewBuildCount(); }
        uint64_t getPointerMergeGroupCountSum() const { return massGrid.mergedGroupCountSum(); }
        unsigned int getPointerMergeMaxGroupCount() const { return massGrid.maxMergeGroupCount(); }

        float getEquivalentThermalResistanceToOuterCasing();
};



/**
 * @brief The Pipe class represents a pair of parallel supply and return pipes between two network nodes of a DEC.
 */
class PipePair {
    private:
        // Dubon add
        Network* parentNetwork;

        static vector<unsigned int> ids; // = {}; // To make sure all ids are unique.

        // Stay constant, define the pipe pair.
        unsigned int id;
        Pipe* supplyPipe;
        Pipe* returnPipe;
        float length; // [m]
        float innerRadius; // [m]
        float interPipeDistance; // [m]
        array<NodePair*,2> connectedNodes;
        array<string,2> singulars; //Added by Max

        float interPipeThermalResistance = 0.f; // Resistance of both pipe insulations and soil [K*m/W]

        // Evolve at each time step.
        bool alreadyTraversed;

        unsigned int networkIndex = 0u; // Position in Network::pipePairs, for the transport schedule.

        float computeInterPipeSoilThermalResistance();
        float computeInterPipeMediumThermalResistance(float supplyOuterRadius, float returnOuterRadius);
        float computeCasingCenterToCasingCenterResistance();
        void applyThreeNodeFullHeatTransportToWaterVolumes(DistrictEnergyCenter* pDEC, float const& exteriorTempS, float const& exteriorTempR, unsigned int dt);

        void recordMassFlows() { supplyPipe->recordMassFlow(); returnPipe->recordMassFlow(); }
        void eraseMassFlows(unsigned int keepValue) { supplyPipe->eraseMassFlow(keepValue); returnPipe->eraseMassFlow(keepValue); }
        void eraseMassFlows_back() { supplyPipe->eraseMassFlow_back(); returnPipe->eraseMassFlow_back(); }

        void recordPressureDiffs() { supplyPipe->recordPressureDiff(); returnPipe->recordPressureDiff(); }
        void erasePressureDiffs(unsigned int keepValue) { supplyPipe->erasePressureDiff(keepValue); returnPipe->erasePressureDiff(keepValue); }
        void erasePressureDiffs_back() { supplyPipe->erasePressureDiff_back(); returnPipe->erasePressureDiff_back(); }

        void recordThermalLosses() { supplyPipe->recordThermalLoss(); returnPipe->recordThermalLoss(); }
        void eraseThermalLosses(unsigned int keepValue) { supplyPipe->eraseThermalLoss(keepValue); returnPipe->eraseThermalLoss(keepValue); }
        void eraseThermalLosses_back() { supplyPipe->eraseThermalLoss_back(); returnPipe->eraseThermalLoss_back(); }

public:
        PipePair(TiXmlHandle hdl, Network* net);
        ~PipePair();

        unsigned int getPipePairsIdx(); // Added by Max (give the index of the pipepair in the vector pipepairs)
        unsigned int getNetworkIndex() const { return networkIndex; } // O(1) equivalent, cached at construction.
        void setNetworkIndex(unsigned int idx) { networkIndex = idx; }
        unsigned int getId() { return id; }
        Pipe* getSupplyPipe() { return supplyPipe; }
        Pipe* getReturnPipe() { return returnPipe; }
        Pipe* getPipe(bool isSupply) { if(isSupply) { return supplyPipe; } else { return returnPipe; } }
        float getLength() { return length; }
        float getInnerRadius() { return innerRadius; }
        Network* getParentNetwork() { return parentNetwork; }
        NodePair* getTailNode() { return connectedNodes[0]; }
        NodePair* getHeadNode() { return connectedNodes[1]; }
        NodePair* getOtherNode(NodePair* np);
        bool getAlreadyTraversed() { return alreadyTraversed; }
        void setAlreadyTraversed(bool b) { alreadyTraversed = b; }
        array<string,2> getSingularities() {return singulars;} // Added by Max
        void setRoughnessFromNetwork(float roughness);

        void recordTimeStep() { recordThermalLosses();  recordMassFlows();  recordPressureDiffs(); }
        void eraseRecords(unsigned int keepValue) { eraseThermalLosses(keepValue); eraseMassFlows(keepValue); erasePressureDiffs(keepValue); }
        void eraseRecords_back() { eraseThermalLosses_back(); eraseMassFlows_back(); erasePressureDiffs_back(); }

        void computeThermalLoss(float soilTemp, float cp, bool isSupply);

        bool massFlowsToMe(NodePair* n, bool isSupply);
        float computeAltitudePressureLoss(float const& rho); // Positive if head is higher than tail.
        void hydraulicConverged(float const& supplyMassFlow, float const& supplyDeltaP, float const& returnMassFlow, float const& returnDeltaP, float const& rho);

        void injectAndAdvectWaterVolumes(Climate* pClimate, unsigned int dtSeconds, unsigned int day, unsigned int hour);

        /**
         * @brief Applies only the wall/soil (and, for the coupled equations, the inter-pipe)
         * heat transport to the parcels already in both pipes, without advecting them.
         * Purely local to the pipe pair, so it is independent of the node update order.
         */
        void applyHeatTransportStep(Climate* pClimate, unsigned int dt, unsigned int day, unsigned int hour);

        /**
         * @brief Advects one side of the pair over dt, injecting the CURRENT temperature of its
         * inlet node. Must be called only once that inlet node has been finalised for this
         * interval, which is what the direction-aware sweep in Network guarantees.
         */
        void advectSide(bool isSupply, unsigned int dt);

        /**
         * @brief Heat transport + advection of one side, for the models that sub-step each pipe
         * side independently. Only valid when the heat transport equation is not one of the
         * supply/return coupled ones (those must use applyHeatTransportStep on the pair).
         */
        void applyHeatTransportAndAdvectSide(bool isSupply, unsigned int dt, float exteriorTemp);
};


/**
 * @brief The NodePair class represents a node in the pipe network of a DEC, PipePairs connect these network nodes. It contains both the node linking the supply pipes and the node linking the return pipes.
 */
class NodePair {
    private:
        static vector<unsigned int> ids; // = {}; // To make sure all ids are unique.

        // Stay constant, define the node.
        unsigned int id; // As written in the xml.
        vector<PipePair*> connectedPipes; // PipePairs that are connected to this NodePair.
        vector<PipePair*> connectedPipesOut; //PipePairs that are connected to this NodePair and that are outgoing. Added by Max
        float z;

        unsigned int transportIndex = 0u; // Position in Network::allNodePairs, for the transport schedule.

        // Evolve at each time step.
        float supplyTemperature, returnTemperature; // [degree C]

        // Stored information about simulation.
        vector<float> supplyTemperatureRecord, returnTemperatureRecord; // [degree C]

        void recordTemperatures() { supplyTemperatureRecord.push_back(supplyTemperature);  returnTemperatureRecord.push_back(returnTemperature); }
        bool temperatureRecordsAreEmpty() { return (supplyTemperatureRecord.empty() and returnTemperatureRecord.empty()); } // They should both have the same size.
        void eraseTemperatures(unsigned int keepValue) { supplyTemperatureRecord.erase(supplyTemperatureRecord.begin(),supplyTemperatureRecord.end()-min(keepValue,(unsigned int)supplyTemperatureRecord.size()));   returnTemperatureRecord.erase(returnTemperatureRecord.begin(),returnTemperatureRecord.end()-min(keepValue,(unsigned int)returnTemperatureRecord.size())); }
        void eraseTemperatures_back() { supplyTemperatureRecord.pop_back();   returnTemperatureRecord.pop_back(); }
        float getSupplyTemperatureRecord(unsigned int step) {return supplyTemperatureRecord.at(step); }
        float getReturnTemperatureRecord(unsigned int step) {return returnTemperatureRecord.at(step); }

    protected:
        bool areAllUpstreamsComputed(bool isSupply);
        void propagateDownstream(Climate* pClim, float cp, bool isSupply, unsigned int hour, unsigned int day); // Added by Max


    public:
        NodePair(TiXmlHandle hdl, float initTemp);
        virtual ~NodePair();

        unsigned int getId() { return id; }
        vector<PipePair*> getPipes() { return connectedPipes; }
        size_t getNbPipes() { return connectedPipes.size(); }
        vector<PipePair*> getPipesOut(){return connectedPipesOut;}// Added by Max
        float getSupplyTemperature() { return supplyTemperature; }
        float getReturnTemperature() { return returnTemperature; }
        float getTemperature(bool isSupply) { if(isSupply) { return supplyTemperature; } else { return returnTemperature; } }
        void setSupplyTemperature(float temp) { supplyTemperature = temp; }
        void setReturnTemperature(float temp) { returnTemperature = temp; }
        void setTemperature(float temp, bool isSupply) { if(isSupply) { supplyTemperature = temp; } else { returnTemperature = temp; } }
        float getZ() { return z; }

        void addConnectedPipe(PipePair *pipe);
        void addConnectedPipeOut(PipePair* pipe); // Added by Max for the singularities.

        void computeAndSetTemperature(float initSumTM, float initSumM, bool isSupply); // Moved here by Dubon

        // --- Direction-aware transport sweep (Network::advanceTransportInterval) ---
        // Index of this node in Network::allNodePairs, assigned once the network is built.
        unsigned int getTransportIndex() const { return transportIndex; }
        void setTransportIndex(unsigned int idx) { transportIndex = idx; }

        /**
         * @brief Finalises this node's temperature on the given side, knowing that every pipe
         * flowing into it on that side has already been advected over the current interval.
         * A node that also crosses supply and return computes its own transfer here, as soon as
         * its inlet side is known, and stores the enthalpy it injects into the other side.
         */
        virtual void transportUpdate(bool isSupply, float rho, float cp, Climate* pClim, unsigned int day, unsigned int hour);

        /**
         * @brief Seeds the cross-injection terms from the temperatures left by the previous
         * interval, so that the first sweep of an interval does not start from a zero source.
         */
        virtual void primeCrossInjection() { }

        float changeRadiusPressureLoss(float m, float innerRadiusIn, float innerRadiusOut);
        void TPressureLoss(vector<float> m, vector<float>& deltaP, bool supply);
        void TPressureLossSplit(vector<float> m, bool supply, vector<float>& deltaP, vector<float>& innerRadiusIn, vector<float>& innerRadiusOut); // Added by Max
        void TPressureLossJunction(vector<float> m, bool supply, vector<float>& deltaP, vector<float>& innerRadius, vector<float>& innerRadiusOut); // Added by Max
        virtual void computePressureDiffAndDerivative(vector<float> m, float const& rho, vector<float>& deltaP);
        virtual int nbEdges() { return connectedPipesOut.size(); }

        virtual void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp) { recordTemperatures();}
        virtual void eraseRecords(unsigned int keepValue) { eraseTemperatures(keepValue); }
        virtual void eraseRecords_back() { eraseTemperatures_back(); }
        virtual string getPrefix(unsigned int decId) { return "DEC" + to_string(decId) + ":NodePair" + to_string(getId()); }
        virtual void writeTHHeaderText(fstream& textFile, unsigned int decId);
        virtual void writeTHResultsText(fstream& textFile, unsigned int i);


        virtual bool propagateNetwork(Climate* pClim, float cp, bool isSupply, unsigned int day, unsigned int hour);// Added by Max

        virtual void checkNbPipes();
};



/**
 * @brief The NodePairConnectingSupplyReturn class represents a NodePair of the DEC network that connects the supply network to the return network (such as a thermal station or a substation), so has a mass flow.
 */
class NodePairConnectingSupplyReturn : public NodePair {
    protected:
        float massFlow; // Positive for substations when flows from supply to return, positive for thermal stations when flows from return to supply [kg/s]
        vector<float> massFlowRecord; // [kg/s]

        float pressureDiff; // Positive when the pressure drops in the same direction as the mass flow (should be positive for substations, negative for thermal station with pumps) [Pa]
        vector<float> pressureDiffRecord; // [Pa]

        bool alreadyTraversed;

        // Enthalpy this element injects into the side it feeds, produced when its inlet side is
        // finalised during the direction-aware sweep. [kg/s * degree C] and [kg/s].
        float crossInjectionSumTM = 0.f;
        float crossInjectionSumM = 0.f;

        void recordMassFlow() { massFlowRecord.push_back(massFlow); }
        void eraseMassFlow(unsigned int keepValue) { massFlowRecord.erase(massFlowRecord.begin(),massFlowRecord.end()-min(keepValue,(unsigned int)massFlowRecord.size())); }
        void eraseMassFlow_back() { massFlowRecord.pop_back(); }
        float getMassFlowRecord(unsigned int step)  { return massFlowRecord.at(step); }

        void recordPressureDiff() { pressureDiffRecord.push_back(pressureDiff); }
        void erasePressureDiff(unsigned int keepValue) { pressureDiffRecord.erase(pressureDiffRecord.begin(),pressureDiffRecord.end()-min(keepValue,(unsigned int)pressureDiffRecord.size())); }
        void erasePressureDiff_back() { pressureDiffRecord.pop_back(); }
        float getPressureDiffRecord(unsigned int step) {return pressureDiffRecord.at(step); }

    public:
        NodePairConnectingSupplyReturn(TiXmlHandle hdl , float const&  initTemp, float const& pressureDiff) : NodePair(hdl , initTemp), massFlow(0.01f), massFlowRecord(0), pressureDiff(pressureDiff), pressureDiffRecord(0) {  }
        virtual ~NodePairConnectingSupplyReturn() override { }

        void setAlreadyTraversed(bool b) { alreadyTraversed = b; }

        /**
         * @brief The side this element feeds, deduced from the sign of its mass flow. Substations
         * and valves are positive from supply to return, thermal stations from return to supply.
         */
        virtual bool outletSideIsSupply() const { return massFlow < 0.f; }
        virtual void primeCrossInjection() override {
            const float m = fabs(massFlow);
            crossInjectionSumTM = m * getTemperature(outletSideIsSupply());
            crossInjectionSumM = m;
        }

        float getMassFlow() { return massFlow; }
        void setMassFlow(float mf) { massFlow = mf; }
        float getPressureDiff() { return pressureDiff; }
        void setPressureDiff(float dp) { pressureDiff = dp; }

        virtual void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp) override { NodePair::recordTimeStep(pClim, day, hour, cp);  recordMassFlow(); recordPressureDiff(); }
        virtual void eraseRecords(unsigned int keepValue) override { NodePair::eraseRecords(keepValue);  eraseMassFlow(keepValue); erasePressureDiff(keepValue); }
        virtual void eraseRecords_back() override { NodePair::eraseRecords_back();  eraseMassFlow_back(); erasePressureDiff_back(); }
        virtual void writeTHHeaderText(fstream& textFile, unsigned int decId) override;
        virtual void writeTHResultsText(fstream& textFile, unsigned int i) override;
};


/**
 * @brief The SubstationNodePair class represents a NodePair of the DEC network to which a substation is connected.
 */
class SubstationNodePair : public NodePairConnectingSupplyReturn {
    private:
        Substation* substation;

    public:
        SubstationNodePair(TiXmlHandle hdl , float initTemp) : NodePairConnectingSupplyReturn(hdl, initTemp, 100000.f), substation(nullptr) {  }
        virtual ~SubstationNodePair() override { }

        void setSubstation(Substation* sub);
        bool substationNotSet() { return substation==nullptr; }
        Substation* getSubstation() { return substation; }

        int nbEdges() override { return substation->nbEdges(); }
        bool hasRegEle() { return substation->hasRegEle(); }

        float getDesignMassFlow() { return substation->computeDesignMassFlow(); }
        float getDesiredMassFlow() { return substation->getDesiredMassFlow(); }
        float getPrimarySideOutputTemperature() { return substation->getPrimarySideOutputTemp(); }

        void updateControlVariable(float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& massFlowSupplyToReturn) { substation->updateControlVariable(massFlow, pressureDiff, rho, sumDeltaRpm, sumRpm, sumDeltaKv, sumKv, learningRate, massFlowSupplyToReturn); }

        void errorMassFlow(float& relErr, float& absErr,float& sumErr) { substation->errorMassFlow(massFlow, relErr, absErr, sumErr); }
        void computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm);
        void hydraulicConverged(vector<float>::const_iterator m, vector<float>::const_iterator deltaP); //override;
        void computeThermal(float const& cp, float const& rho, Climate* pClim, unsigned int day, unsigned int hour, bool supplyToReturn);

        virtual bool propagateNetwork(Climate* pClim, float cp, bool isSupply, unsigned int day, unsigned hour) override; // Added by Max
        // A substation draws from the supply when it consumes (massFlow>0) and feeds it back when
        // it behaves as a producer (massFlow<=0), matching the branches of computeThermal.
        virtual bool outletSideIsSupply() const override { return massFlow <= 0.f; }
        virtual void transportUpdate(bool isSupply, float rho, float cp, Climate* pClim, unsigned int day, unsigned int hour) override;

        virtual void checkNbPipes() override;

        void updateDesiredMassFlow(float const& cp, Climate* pClim, unsigned int day, unsigned int hour) { substation->updateDesiredMassFlow(cp, getReturnTemperature(), pClim, day, hour); }
        bool hasSubstation(Substation* sub) { return sub==substation; }

        virtual void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp) override { NodePairConnectingSupplyReturn::recordTimeStep(pClim, day, hour, cp); substation->recordTimeStep();  }  // recordTemperatures  recordMassFlow  recordPressureDiff  setProsumerSolarThermal
        virtual void eraseRecords(unsigned int keepValue) override { NodePairConnectingSupplyReturn::eraseRecords(keepValue); substation->eraseRecords(keepValue);  }
        virtual void eraseRecords_back() override { NodePairConnectingSupplyReturn::eraseRecords_back(); substation->eraseRecords_back();  }
        virtual string getPrefix(unsigned int decId) override { return "DEC" + to_string(decId) + ":SubstationNodePair" + to_string(getId()); }
        virtual void writeTHHeaderText(fstream& textFile, unsigned int decId) override;
        virtual void writeTHResultsText(fstream& textFile, unsigned int i) override;
};


/**
 * @brief The ThermalStationNodePair class represents a NodePair of the DEC network to which a thermal station is connected.
 */
class ThermalStationNodePair : public NodePairConnectingSupplyReturn {
    private:
        ThermalStation* thermalStation;
    public:
        ThermalStationNodePair(TiXmlHandle hdl , float initTemp) : NodePairConnectingSupplyReturn(hdl , initTemp, -100000.f), thermalStation(nullptr) { }
        virtual ~ThermalStationNodePair() override { }

        void setThermalStation(ThermalStation* st);
        ThermalStation* getThermalStation() { return thermalStation; }
        bool thermalStationNotSet() { return thermalStation==nullptr; }

        void updateControlVariable(float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate) { thermalStation->updateControlVariable(massFlow, pressureDiff, rho, sumDeltaRpm, sumRpm, sumDeltaKv, sumKv, learningRate); }
        void updateDesiredMassFlow(float const& cp, Climate* pClimate, unsigned int day, unsigned hour) { thermalStation->updateDesiredMassFlow(cp, pClimate, day ,hour); } // Added by Max
        void errorPressureDiff(float& relErr, float& absErr, float& sumErrP, float& sumDesiredP) { thermalStation->errorPressureDiff(massFlow, pressureDiff, relErr, absErr, sumErrP, sumDesiredP); }
        void computeThermal(float rho, float cp, Climate* pClim, unsigned int day, unsigned int hour, bool supplyToReturn);

        virtual bool propagateNetwork(Climate* pClim, float cp, bool isSupply, unsigned int day, unsigned int hour) override; //Added by Max
        // A thermal station injects into the supply when producing (massFlow>=0) and behaves as a
        // consumer or slack otherwise, matching the branches of Network::updateNodesTemperatures.
        virtual bool outletSideIsSupply() const override { return massFlow >= 0.f; }
        virtual void transportUpdate(bool isSupply, float rho, float cp, Climate* pClim, unsigned int day, unsigned int hour) override;

        virtual void checkNbPipes() override;

        virtual void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp) override { NodePairConnectingSupplyReturn::recordTimeStep(pClim, day, hour, cp); thermalStation->recordTimeStep(pClim, day, hour, massFlow, cp);  }
        virtual void eraseRecords(unsigned int keepValue) override { NodePairConnectingSupplyReturn::eraseRecords(keepValue); thermalStation->eraseRecords(keepValue); }
        virtual void eraseRecords_back() override { NodePairConnectingSupplyReturn::eraseRecords_back(); thermalStation->eraseRecords_back(); }
        virtual string getPrefix(unsigned int decId) override { return "DEC" + to_string(decId) + ":ThermalStationNodePair" + to_string(getId()); }
        virtual void writeTHHeaderText(fstream& textFile, unsigned int decId) override;
        virtual void writeTHResultsText(fstream& textFile, unsigned int i) override;
};

// Added by Max
class ValveNodePair: public NodePairConnectingSupplyReturn{
    private:
        Valve* valve;
        float kvMax;
        float outputTemperature;

        PIDControllerValve pid; // To control the valve opening, to get the desired mass flow.
        // Done in order to have the case where the valve is controlled manually
        float Targetkv; // Added by Max
        bool ImposedValve; // Added by Max
    public:
        ValveNodePair(TiXmlHandle hdl, float initTemp);
        virtual ~ValveNodePair() override { }

        Valve* getValve() { return valve; }

        int nbEdges() override{ return valve->nbEdges(); }

        virtual float maxKv();
        void updateControlVariable(float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate);
        virtual float computePressureDiffControlElement(float pressureDiffSubstationNodePair) { return pressureDiffSubstationNodePair; } // All pressure loss is in the control element (valve).
        virtual void computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm);
        void hydraulicConverged(vector<float>::const_iterator m, vector<float>::const_iterator deltaP); //override;
        void computeThermal(float const& cp, float const& rho, Climate* pClim, unsigned int day, unsigned int hour, bool supplyToReturn);
        float getOutputTemp(){return outputTemperature;}

        virtual bool propagateNetwork(Climate* pClim, float cp, bool isSupply, unsigned int day, unsigned int hour) override; // Added by Max
        // A valve sends supply to return for a positive mass flow, and the other way round below.
        virtual bool outletSideIsSupply() const override { return massFlow < 0.f; }
        virtual void transportUpdate(bool isSupply, float rho, float cp, Climate* pClim, unsigned int day, unsigned int hour) override;

        virtual void checkNbPipes() override;

        virtual void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp) override { NodePair::recordTimeStep(pClim, day, hour, cp); }
        virtual void eraseRecords(unsigned int keepValue) override { NodePair::eraseRecords(keepValue); }
        virtual void eraseRecords_back() override { NodePair::eraseRecords_back(); }
        virtual string getPrefix(unsigned int decId) override { return "DEC" + to_string(decId) + ":ValveNodePair" + to_string(getId()); }
        virtual void writeTHHeaderText(fstream& textFile, unsigned int decId) override;
        virtual void writeTHResultsText(fstream& textFile, unsigned int i) override;
};

class Network {
private:
    DistrictEnergyCenter *pDEC;

    vector<NodePair*> separateNodePairs;
    vector<SubstationNodePair*> substationNodePairs;
    vector<ThermalStationNodePair*> thermalStationNodePairs;
    vector<ValveNodePair*> valveNodePairs;

    vector<PipePair*> pipePairs;

    vector<RegulatingElement*> regulatingElements;

    // --- Direction-aware thermal transport schedule ---
    // The dynamic models used to advect every pipe from the node temperatures of the previous
    // sub-step and only then update all the nodes. That Jacobi update adds one sub-step of
    // artificial delay per node crossed, so a thermal front arrived late by (number of nodes on
    // the path) x dt and the result depended on how a pipe run was split into segments. The
    // schedule below orders the updates along the flow direction instead: a pipe is advected only
    // once its inlet node is final for the interval, which removes that accumulated lag.
    struct TransportLevel {
        vector<NodePair*> nodes;  // Finalised at this depth: all their inflowing pipes are done.
        vector<PipePair*> pipes;  // Advected at this depth: their inlet node has just been set.
    };
    struct TransportSchedule {
        vector<TransportLevel> levels;
        vector<PipePair*> stagnantPipes; // No mass flow on this side, so outside the flow graph.
        unsigned int nbBackEdges = 0u;   // Loop edges still relaxed with the previous sub-step.
        unsigned int nbActivePipes = 0u;
        unsigned int sweepThreads = 1u;  // 1 when the levels are too thin to pay for a team.
    };
    // Pipes a level must hold per thread before threading the sweep pays for its barriers.
    static constexpr unsigned int MIN_PIPES_PER_LEVEL_THREAD = 24u;
    vector<NodePair*> allNodePairs; // Every node of the network, in a stable order.
    TransportSchedule supplyTransport, returnTransport;
    bool transportSchedulesLogged = false;

    // Optional sub-hourly temperature trace (see DECSimulationConfig::subHourlyDiagnosticFile).
    fstream subHourlyDiagnosticStream;
    vector<NodePair*> subHourlyDiagnosticNodes;
    bool subHourlyDiagnosticReady = false;
    bool subHourlyDiagnosticEnabled = false;
    unsigned int subHourlyDiagnosticCounter = 0u;

    float soilkValue;
    float airHValue = Model::defaultAirConvectiveCoefficient; // Thermal convective coefficient
    float roughness = Model::defaultPipeRoughness; // absolute rougness in (m)

    // Parallelism parameters — computed once in initDynamicParallelism(), read-only during simulation.
    // If not set, we are in Debug mode
    unsigned int numDynamicThreads  = 1;  // Effective OMP thread count for the pipe advection loop.
    unsigned int dynamicPipeChunkSize = 1; // Minimum chunk size for schedule(guided, N).
    unsigned int maxThreadsFromXml  = 0;  // 0 = auto (no cap); set via XML attribute 'maxThreads'.
    // Minimum number of pipes required per thread before spawning an additional thread.
    // Below this ratio, thread synchronisation overhead exceeds the parallelism benefit.
    // This value may be tuned if necessary
    static constexpr unsigned int MIN_PIPES_PER_THREAD = 1u;
    // fixedDV maxDtMode: pipes carrying less than this fraction of the largest pipe mass flow do not
    // set the synchronisation interval. A near-stagnant dead end needs hours to push the fixed
    // volume, which would otherwise pin the interval of the whole network to its cap.
    static constexpr float FIXED_DV_SIGNIFICANT_FLOW_FRACTION = 1.e-3f;

    vector<float> loopMassFlows; // Mass flow through the loops. Values stored after convergence, to be used as initial condition in the next iteration [kg/s]
    vector<vector<float>> loopMatrix; // nb of rows=nb of loops, nb of columns=nb of edges. If [i,j]==1, then edge j is in loop i. If [i,j]==-1, then edge j is in loop i, but in the opposite direction. If [i,j]==0, then edge j is not in loop i.
    vector<vector<float>> regulatedPathMatrix; // nb of rows=nb of regulated paths, nb of columns=nb of edges. If [i,j]==1, then edge j is in regulated path i. If [i,j]==-1, then edge j is in regulated path i, but in the opposite direction. If [i,j]==0, then edge j is not in regulated path i.
    vector<vector<float>> jacobianRegulatedPathColumns; // nb of rows=nb of loops+nb of regulated paths, nb of columns=nb of regulated paths.

    // Private methods.
    void deleteNodesAndPipes();
    void convergeThermal(float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour);
    void computeResiduals(vector<float>& edgeMassFlows, vector<double>& residuals, vector<float>& deltaP, vector<float>& dDeltaP_dm, float const& rho);
    double computeResidualNorm(vector<double> const& residuals);
    double computeMuPrime(vector<double> const& residuals_n, vector<double> const& residuals_np1, double const& lambda_n, double const& residualNorm_n);
    bool hasConverged(vector<double> const& residuals);
    void convergeHydraulic(float rho);
    void initDynamicParallelism(); // Computes numDynamicThreads and dynamicPipeChunkSize from network size and hardware.
    void solveSteadyStateHydraulic(MCR* mcr, float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour);
    void updateNodesTemperatures(float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour);

    /**
     * @brief Builds the flow-ordered update schedule of one side from the mass flows of the
     * current hydraulic solution. Nodes are levelled by longest path from the sources (Kahn), so
     * every pipe within a level is independent and can be advected in parallel. Meshed parts that
     * cannot be levelled are broken at the edge with the fewest unresolved predecessors; those
     * back edges keep the previous sub-step value, as the whole network did before.
     */
    void buildTransportSchedule(bool isSupply, TransportSchedule& schedule);
    void buildTransportSchedules();
    void sweepTransportDirection(bool isSupply, float rho, float cp, Climate *pClim,
                                 unsigned int day, unsigned int hour,
                                 function<void(PipePair*, bool)> const& advectPipeSide);
    /**
     * @brief Advances the whole network over one synchronisation interval, coupling the nodes in
     * flow order: supply from the producers down to the consumers, then return back up.
     * @param advectPipeSide Advances one side of one pipe pair over the whole interval.
     * @param pairHeatTransportPrePass True when the heat transport equation couples the supply
     * and return parcels, so it has to be applied to the pair before either side is advected.
     */
    void advanceTransportInterval(unsigned int syncDt, float rho, float cp, Climate *pClim,
                                  unsigned int day, unsigned int hour,
                                  function<void(PipePair*, bool)> const& advectPipeSide,
                                  bool pairHeatTransportPrePass);

    // Interval-weighted averages of substation and thermal station powers over one hydraulic step.
    void resetTransportAverages();
    void accumulateTransportAverages(unsigned int dt);
    void applyTransportAverages();

    void initSubHourlyDiagnostic();
    void writeSubHourlyDiagnostic(unsigned int day, unsigned int hour, unsigned int elapsedSeconds);

    static float computeLaminarNusseltNumber(float const& length, float const& radius, float const& reynolds, float const& prandtl);
    static float computeTransitionNusseltNumber(float const& length, float const& radius, float const& reynolds, float const& prandtl);
    static float computeTurbulentNusseltNumber(float const& length, float const& radius, float const& reynolds, float const& prandtl);

public:
    ostream logStream;

    Network(TiXmlHandle hdl, DistrictEnergyCenter *pDEC, ostream *pLogStr);
    ~Network() { /*logStream << "Destructor of Network" << endl;*/ deleteNodesAndPipes(); }

    float getSoilkValue() { return soilkValue; }

    // Dubon start add
    float getAirHValue() { return airHValue; }
    float getRoughness() { return roughness; }

    DistrictEnergyCenter* getDEC() { return pDEC; }
    void addRegulatingElement(RegulatingElement* re) { regulatingElements.push_back(re); }

    void writeTHHeaderText(fstream& textFile);
    void writeTHResultsText(fstream& textFile, unsigned int i);
    int nbFloatsRecorded();

    void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp);
    void eraseRecords(unsigned int keepValue);
    void eraseRecords_back();


    /**
     * @return Sum of all thermal losses in the supply and return pipes of the network.
     */
    float computeTotalThermalLoss();

    NodePair* pointerOfNodeWithId(unsigned int nodeId);
    SubstationNodePair *pointerOfSubstationNodeWithId(unsigned int nodeId);
    ThermalStationNodePair* pointerOfThermalStationNodeWithId(unsigned int nodeId);
    ValveNodePair* pointerOfValveNodeId(unsigned int nodeId);

    vector<ThermalStationNodePair*> getThermalStationNodePairs(){return thermalStationNodePairs;} //Added by Max

    void checkThermalStationsFound();
    void substationsConstructionFinished();

    void initializeMassFlows();

    void propagateNetwork(float cp, Climate* pClim, bool isSupply, unsigned int day, unsigned hour); //Added by Max
    void advectNetworkWaterVolumes(Climate* pClim, unsigned int dtSeconds, unsigned int day, unsigned int hour); // Dubon

    ostream& print(ostream& stream);

    void computeEdgeIdx(map<PipePair*,size_t>& edgeIdxPipe, map<NodePairConnectingSupplyReturn*,size_t>& edgeIdxNode, size_t& nbEdges); // Changed by Max
    size_t nbLoops() { return loopMatrix.size(); }
    size_t nbEdges() { return loopMatrix[0].size(); }
    size_t nbRegPaths() { return regulatedPathMatrix.size(); }
    void findLoops(vector<vector<float>>& matB, map<PipePair*,size_t> const& edgeIdxPipe, map<NodePairConnectingSupplyReturn*,size_t> const& edgeIdxNode, size_t const& nbEdges); // Changed by Max
    void findRegulatedPaths(vector<vector<float>>& matR, map<NodePairConnectingSupplyReturn*,size_t> const& edgeIdxNode, size_t const& nbEdges); // Added by Max
    void computeJacobianRegulatedPathColumns(vector<vector<float>>& matJRPC, vector<vector<float>> const& matB, vector<vector<float>> const& matR, map<NodePairConnectingSupplyReturn*,size_t> const& edgeIdxNode, size_t const& nbEdges); // Changed by Max

    void matTransposedVecMult(vector<vector<float>> const& matTransposed, vector<float> const& vec, vector<float>& result);
    /**
     * @brief matVecMult Multiplies matrix*vector
     * @param mat Matrix to multiply.
     * @param vec Vector to multiply.
     * @param result Vector where product is stored (but can contain be bigger).
     * @param iResStart Start index in the result vector, to store the product.
     */
    void matVecMult(vector<vector<float>> const& mat, vector<float> const& vec, vector<double>& result, size_t const& iResStart);

    void matDiagMatTransposedMult(vector<vector<float>> const& mat, vector<float> const& diag, vector<vector<float>> const& matT, double *resultMat, int resultMat_n, size_t const& iResStart, size_t const& jResStart);

    /**
     * @brief computeDarcyFrictionFactor
     * @param massFlow Must be positive!
     * @param radius
     * @return
     */
    float computeDarcyFrictionFactor(float const& massFlow, float const& radius, float const& temp, float const& roughness);
    void computePressureDiff(vector<float> const& massFlows, float const& rho, vector<float>& deltaP, vector<float>& dDeltaP_dm);

    // Dubon add
    void dynamicForwardPassFixedDt(MCR* mcr, float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour);
    void dynamicForwardPassFixedDV(MCR* mcr, float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour);
    static float computeNusseltNumber(float const& massFlow, float const& innerRadius, float const& length, float const& waterTemp, float const& kWater, float const& muWater, float const& cpWater);

    void convergeToEquilibrium(MCR* mcr, float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour);
    void computeDerivative(float& deriv, vector<float>& prev, float const& curr, size_t const& idx);

    void updateDesiredMassFlow(float const& cp, Climate* pClim, unsigned int day, unsigned int hour);
    void updateOperationMode();
    float computeMassFlowSupplyToReturn();
    float avgPressureDiffSupplyReturn(Substation* apartFromMe);
};

struct DECSimulationConfig {
    enum BalanceDtMode {
        BALANCE_DT_MAX,
        BALANCE_DT_USER_CAP,
        BALANCE_DT_PERCENTILE
    };

    enum HeatTransportEquation {
        HEAT_TRANSPORT_THREE_NODE,
        HEAT_TRANSPORT_THREE_NODE_FULL
    };

    bool defined = false;
    string thermalModel = "steady-state";

    unsigned int steadyStepInSeconds = 3600u;
    unsigned int fixedDt = 0u;
    float fixedDV = 0.f;
    // massGrid: the cell mass is |m_dot| x massGridRefDt of the first hydraulic solution. It sets the
    // transport resolution independently of fixedDt, which only sets how often the heat is exchanged.
    unsigned int massGridRefDt = 5u;
    // pointerMerged: same immutable mass-coordinate plugs as massGrid, plus a transient CRS merge
    // view rebuilt before every advection.  The physical plugs are never coalesced.
    unsigned int pointerMergeRefDt = 5u;
    float pointerMergeTemperatureThreshold = 0.05f;
    bool hasSteadyStepInSeconds = false;
    bool hasFixedDt = false;
    bool hasFixedDV = false;
    bool hydraulicControlOuterLoop = true;

    BalanceDtMode balanceDtMode = BALANCE_DT_MAX;
    unsigned int balanceDtUserCapSeconds = 0u;
    float balanceDtPercentile = 75.f;
    HeatTransportEquation heatTransportEquation = HEAT_TRANSPORT_THREE_NODE;

    // Optional sub-hourly trace of node temperatures. The regular _TH.out is written once per
    // hour, which cannot resolve the transport phase the dynamic models are meant to capture.
    string subHourlyDiagnosticFile = "";        // Empty disables the trace.
    vector<unsigned int> subHourlyDiagnosticNodeIds; // Empty traces every substation and station.
    unsigned int subHourlyDiagnosticPeriod = 1u;     // Write one line every N synchronisation steps.
};

class DistrictEnergyCenter {
private:
    District* pDistrict; // pointer back to the district it belongs to

    unsigned int id;

    double cp, rho;
    static vector<string> muPossibilities; string mu = "0.0004"; // Dynamic viscocity of water [Pa*s]
    vector<ThermalStation*> thermalStations;
    Network* pipelineNetwork = nullptr;

    vector<float> totalThermalLossRecord; // Useful ? It can be computed just by summing up all the losses in the output file.

    MCR* mcr = nullptr;

public:
    ostream logStream;

    DistrictEnergyCenter(TiXmlHandle hdl, District* pDistrict, ostream *pLogStr);
    ~DistrictEnergyCenter() { /*logStream << "Destructor of DistrictEnergyCenter" << endl;*/ deleteDynAllocated(); }
    void deleteDynAllocated();

    District* getDistrict() {return pDistrict;}
    double getCp() { return cp; }
    double getRho() { return rho; }
    double getMu(double temp);
    double getkWater(double temp) { return computeWaterThermalConductivityIAPWS(temp); }
    unsigned int getId() { return id; }
    Network* getNetwork() { return pipelineNetwork; }
    float getInitTemp() { return 15.f; } // TODO Imporve this. Avg the thermal station init temperatures ?

    int nbFloatsRecorded() { return 1+pipelineNetwork->nbFloatsRecorded(); }

    float getTotalThermalLoss(unsigned int step) {return totalThermalLossRecord.at(step); }
    float getTotalThermalLoss() { return totalThermalLossRecord.back(); }
    float getYearlyTotalThermalLoss();
    size_t getnThermalStations() {return thermalStations.size();}
    ThermalStation* getThermalStation(unsigned int index) {return thermalStations[index];}

    void recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour);
    void eraseRecords(unsigned int keepValue);
    void eraseRecords_back();

    /**
     * Once the thermal power needs of all substation are computed (and set), this converges to the solution of temperatures and mass flows throughout the DEC. In particular it will compute the thermal power that the DEC will by able to provide to each substation.
     * @param pClim Pointer to climate
     * @param day Day of the year
     * @param hour Hour of the day
     */
    void convergeToEquilibrium(Climate* pClim, unsigned int day, unsigned int hour);

    /**
     * Searches if the DEC contains that node, if not returns false. If it does contain the node saves pointer to substation in a substations vector and in the node, and sets the substations pDEC and pNode, then returns true.
     * @param pSub Substation to add.
     * @param nodeId Id of the node that the substation is trying to connect to.
     * @return True if the DEC contains node of this id and the add was successful; false otherwise
     */
    bool addSubstationIfContainNode(Substation* pSub, unsigned int nodeId);

    void checkThermalStationsFound() { getNetwork()->checkThermalStationsFound(); }
    void substationsConstructionFinished() { getNetwork()->substationsConstructionFinished(); }

    float avgPressureDiffSupplyReturn(Substation* apartFromMe) { return pipelineNetwork->avgPressureDiffSupplyReturn(apartFromMe); }

    static float computeDynamicViscosityWaterIAPWS(float const& temp);
    static float computeWaterThermalConductivityIAPWS(float const& temp);

    // Fast table-lookup alternatives to the IAPWS formulas above.
    // Use 10°C-interval linear interpolation: error < 0.3% over the DHC range.
    // Intended for the per-pipe representative temperature in the Lagrangian inner loop,
    // where one call per pipe replaces one call per volume element.
    static float computeKFast(float temp);
    static float computeMuFast(float temp);
};






/*
class BoilerAndElectricElement : public EnergyConversionSystem {

    private:

        double boilerThermalPower, boilerThermalEfficiency;
        double electricElementPower;

    public:

        BoilerAndElectricElement(double boilerThermalPower, double boilerThermalEfficiency, double electricElementPower) {

            this->boilerThermalPower = boilerThermalPower;
            this->boilerThermalEfficiency = boilerThermalEfficiency;
            this->electricElementPower = electricElementPower;

        }

        friend ostream& operator<< (ostream& s, BoilerAndElectricElement& unit) {

            s << "\n\nBoiler:\nThermal Power: " << unit.boilerThermalPower << " W(th)\nEfficiency: " << unit.boilerThermalEfficiency << endl;
            s << "\nElectric element power: " << unit.electricElementPower << endl;

            return s;

        }

        void getMaxThermalPower(double thermalPower1Needed, double thermalPower2Needed, double &thermalPower1Available, double &thermalPower2Available, double sourceTemp) {

            if ( (thermalPower1Needed+thermalPower2Needed) <= (boilerThermalPower+electricElementPower) ) {

                thermalPower1Available = thermalPower1Needed;
                thermalPower2Available = thermalPower2Needed;

            }
            else {

                if ( boilerThermalPower + electricElementPower >= thermalPower1Needed ) {
                    // th1 satisfait
                    thermalPower1Available = thermalPower1Needed;
                    thermalPower2Available = boilerThermalPower + electricElementPower - thermalPower1Needed;
                }
                else {
                    // th1 meme pas satisfait, donne le max
                    thermalPower1Available = boilerThermalPower + electricElementPower;
                    thermalPower2Available = 0.0;
                }
            }
        }

        double getCO2Production(double time, double thermalPower1, double thermalPower2, double sourceTemp, double outputTemp1, double outputTemp2) {

            if (thermalPower1+thermalPower2 <= boilerThermalPower) return (time*(thermalPower1+thermalPower2)/boilerThermalEfficiency)*boilerCO2EmissionCoefficient;
            else return (time*boilerThermalPower/boilerThermalEfficiency)*boilerCO2EmissionCoefficient + (time*(thermalPower1+thermalPower2-boilerThermalPower)/1.0) * electricCO2EmissionCoefficient;

        }

        double getFuelConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) {

            if (thermalPower1+thermalPower2 <= boilerThermalPower) return (time*(thermalPower1+thermalPower2)/boilerThermalEfficiency);
            else return (time*boilerThermalPower/boilerThermalEfficiency);

        }

        double getElectricConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) {

            if (thermalPower1+thermalPower2 <= boilerThermalPower) return 0.0;
            else return (time*(thermalPower1+thermalPower2-boilerThermalPower)/1.0);

        }

        void getMaxThermalPower(vector<double> thermalPowerNeeded, vector<double> &thermalPowerAvailable, double sourceTemp) {

            if ( accumulate(thermalPowerNeeded.begin(), thermalPowerNeeded.end(), 0.0) <= (boilerThermalPower+electricElementPower) ) {

                thermalPowerAvailable = thermalPowerNeeded;

            }
            else { // accu
                // pas assez d'�nergie, regardons ce qu'il manque pour le boiler
                double accu = 0.0;

                for (unsigned int i=0; i<thermalPowerNeeded.size(); i++) {

                    accu += thermalPowerNeeded[i];
                    if ( electricElementPower + boilerThermalPower >= accu) {
                        // ce bout satisfait
                        thermalPowerAvailable.push_back(thermalPowerNeeded[i]);
                    }
                    else {
                        // ce bout pas satisfait
                        thermalPowerAvailable.push_back( max(electricElementPower + boilerThermalPower - (accu-thermalPowerNeeded[i]), 0.0) );
                    }

                }

            }
        }

        double getCO2Production(double time, vector<double> thermalPower, double sourceTemp, double outputTemp) {

            if (accumulate(thermalPower.begin(), thermalPower.end(), 0.0) <= boilerThermalPower) return (time*(accumulate(thermalPower.begin(), thermalPower.end(), 0.0))/boilerThermalEfficiency)*boilerCO2EmissionCoefficient;
            else return (time*boilerThermalPower/boilerThermalEfficiency)*boilerCO2EmissionCoefficient + (time*(accumulate(thermalPower.begin(), thermalPower.end(), 0.0)-boilerThermalPower)/1.0) * electricCO2EmissionCoefficient;
        }

};
*/


#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Woverloaded-virtual"
#endif
