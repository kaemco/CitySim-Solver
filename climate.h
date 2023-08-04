#ifndef CLIMATE
#define CLIMATE

#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <map>

#include "models.h"
#include "GENAngle.h"

using namespace std;

// *** Climate class, CitySim *** //
// *** jerome.kaempf@epfl.ch  *** //

class Climate {

private :

    // parent
    ostream logStream;

    // these values are from the climate file
    string filename;
    string location;
    float altitude, latitudeN, longitudeE;
    int meridianE; // in hours
    vector<float> Tout, Tground; // in celsius
    vector<float> windSpeed,windDirection;
    vector<float> relativeHumidity, Ibn, Idh, Prec;
    vector<float> cloudiness;

    // these values are calculated from the previous values
    vector<float> Td;                   // dew point temperature in celsius
    vector<float> meanDailyTemperature; // in celsius
    float meanAnnualTemperature;        // in celsius
    unsigned int coolDay, hotDay;       // index of the coldest and hottest days
    vector<float> Igh;                  // global horizontal irradiance
    vector<float> Igh_vis;              // global horizontal illuminance

    // values from the .cli2 if available
    map<unsigned int, vector<float> > KeCoeff1, KeCoeff2, KeCoeff3, Tla;

    // values from the .cli3 if available
    map<float,vector<float> > FF, DD, Ta;

public :

    Climate(ostream* logFileStream=nullptr);
    Climate(string filename, ostream* logFileStream=nullptr);
    void clear() {
        Igh.clear(); Igh_vis.clear();
    }
    void importPVGIS_TMY_CSV(string filename, int defaultCloudiness);
    void computeComplementaryValues();
    void exportCliFile(string filename);

    string getLocation() { return location; }

    float getToutCelsius(unsigned int day, unsigned int hour)   { return Tout.at((day-1)*24 + hour -1); }
    float getToutCelsius(unsigned int day, unsigned int hour, float height) {
        for (map<float,vector<float> >::iterator it=Ta.begin(); it!=Ta.end();++it) {
            if (height < it->first) return it->second.at((day-1)*24 + hour -1);
        }
        return NAN;
    }
    float getToutCelsius(int step)                              { return (step >= 0) ? Tout[step] : Tout[step + Tout.size()]; }

    float getTgroundCelsius(unsigned int day, unsigned int hour, float depth=0.f, float alpha=(0.25e-6f/(0.89f*1.6f))*24.f*3600.f, float depthMax=0.f); // default alpha value taken for clay (Arya, 2001) in m²/day

    float getWindSpeed(unsigned int day, unsigned int hour)     { return windSpeed[(day-1)*24 + hour -1]; }
    float getWindSpeed(unsigned int step)                       { return windSpeed.at(step); }
    float getWindSpeed(unsigned int day, unsigned int hour, float height) {
        for (map<float,vector<float> >::iterator it=FF.begin(); it!=FF.end();++it) {
            if (height < it->first) return it->second.at((day-1)*24 + hour -1);
        }
        return NAN;
    }
    float getWindDirection(unsigned int day, unsigned int hour) { return windDirection[(day-1)*24 + hour -1]; }
    float getWindDirection(unsigned int day, unsigned int hour, float height) {
        for (map<float,vector<float> >::iterator it=DD.begin(); it!=DD.end();++it) {
            if (height < it->first) return it->second.at((day-1)*24 + hour -1);
        }
        return NAN;
    }

    float getRelativeHumidity(unsigned int day, unsigned int hour) { return relativeHumidity[(day-1)*24 + hour -1]/100.f; } // relative humidity \in [0,1]
    float getRelativeHumidity(int step) { return (step >= 0) ? relativeHumidity[step]/100.f : relativeHumidity[step + relativeHumidity.size()]/100.f; } // relative humidity \in [0,1]

    double getPatm(unsigned int day,unsigned int hour) { return 101325.0*exp(-(28.97/1000)*9.81*altitude/(8.3145*(Tout[(day-1)*24 + hour -1]+273.15))); /*P=P0exp(Mgh/RT)*/ }
    double getPatm(unsigned int it) { return 101325.0*exp(-(28.97/1000)*9.81*altitude/(8.3145*(Tout[it]+273.15))); /*P=P0exp(Mgh/RT)*/ }
    float getLatitudeN()  { return latitudeN; }
    float getLongitudeE() { return longitudeE; }
    float getAltitude()  { return altitude; } // JK - added return altitude 16.02.09
    int getMeridian()    { return meridianE; } // JK - added meridian, the meridian refers to the local time difference due to the time zone (ex. GMT+1)

    float getAirDensity() { return 1.201385*exp(-1.219755e-4*altitude); /* equation from Appendix B - IEA BESTEST */ }
    static float getAirDensity(float myAltitude) { return 1.201385*exp(-1.219755e-4*myAltitude); /* equation from Appendix B - IEA BESTEST */ }

    float getTd(unsigned int step) { return Td[step]; }
    float getTd(unsigned int day, unsigned int hour) { return Td[(day-1)*24 + hour -1]; }

    float getClearness(unsigned int day, unsigned int hour) { return (8.f-cloudiness[(day-1)*24 + hour -1])/8.f; } // retourne la clearness p.r. a la n�bulosit� du fichier climatique en coeff et non en Oktas
    float getCloudCoverFraction(unsigned int day, unsigned int hour) { return cloudiness[(day-1)*24 + hour -1]/8.f; } // returns the cloud cover fraction between 0 and 1

    float getIdh(unsigned int day, unsigned int hour) { return Idh[(day-1)*24 + hour -1]; }
    float getIbn(unsigned int day, unsigned int hour) { return Ibn[(day-1)*24 + hour -1]; }
    float getIdh(unsigned int step){ return Idh[step]; }
    float getIbn(unsigned int step){ return Ibn[step]; }

    // defines the emissivity of the clear sky from Berdahl and Martin (1984), with Td dew point temperature (°C)
    float getEpsilon_clrsky(unsigned int day, unsigned int hour) { return 0.711 + 0.56*(getTd(day,hour)/100.f) + 0.73*pow((getTd(day,hour)/100.f),2) + 0.013*cos(2.*M_PI*(hour/24.)) + 0.00012*(getPatm(day,hour)*0.01-1013.25); } // hour being between 1 and 24
    // defines the emissivity of the cloudy sky from Eicker and Dalibard (2011), as a function of the cloud cover fraction and 0.9 is the default emissivity of the clouds
    float getEpsilon_sky(unsigned int day, unsigned int hour) { return getEpsilon_clrsky(day,hour) + (1.f-getEpsilon_clrsky(day,hour))*getCloudCoverFraction(day,hour)*0.9f; }
    // defined the TSky for output
    float getTskyCelsius(unsigned int day, unsigned int hour) { return pow(getEpsilon_sky(day,hour),0.25)*(getToutCelsius(day,hour)+273.15) - 273.15; }

    // variables for the ET model
    float getGamma(unsigned int day, unsigned int hour) { return 37.*getWindSpeed(day,hour)*694.5/((getToutCelsius(day,hour)+273.15)*(1.+0.34*getWindSpeed(day,hour))); } // Hourly determination of ET0
    static float getSaturatedVapourPressure(float T)           { return 0.611*exp(17.502*T/(T+240.97)); } // saturation vapour pressure kPa (Campbell and Norman, 1998)
    static float getSaturatedVapourPressureDerivative(float T) { return 240.97*17.502*0.611*exp(17.502*T/(T+240.97))/pow((T+240.97),2); } // slope of saturation vapour pressure (kPa/°C)
    float getVapourPressure(unsigned int day, unsigned int hour) { return getSaturatedVapourPressure(getToutCelsius(day,hour))*getRelativeHumidity(day,hour); } // Vapour pressure (kPa)
    float getVapourPressure(unsigned int step)                   { return getSaturatedVapourPressure(getToutCelsius(step))*getRelativeHumidity(step); } // Vapour pressure (kPa)
    float getVapourDeficit(unsigned int day, unsigned int hour)  { return getSaturatedVapourPressure(getToutCelsius(day,hour))-getVapourPressure(day,hour); } // Difference (kPa)

    float getAnnualMeanTemperature() { return meanAnnualTemperature; }
    vector<float> getMeanDailyTemperature() { return meanDailyTemperature; }
    float getSwingInMeanDailyTemperature() { return meanDailyTemperature[hotDay]-meanDailyTemperature[coolDay]; }
    unsigned int getDayWithMinMeanTemperature() { return coolDay; }

    // gets the average temperature of the last 24 hours
    float getDailyMeanTemperature(unsigned int day, unsigned int hour) {
        float Tdm = 0.f;
        for (unsigned int dayIndex = 0; dayIndex < 24; ++dayIndex) {
            Tdm += getToutCelsius(((static_cast<int>(day)-1)*24+static_cast<int>(hour)-1)+(dayIndex-24))/24.f;
        }
        return Tdm; // in celsius
    }

    bool getRainPresence(unsigned int day, unsigned int hour) { return (Prec[(day-1)*24 + hour -1] > 0.f) ? true : false; }

    // Maxwell model (Igh -> Ibn)
    static float dirint(float *g, float *z, float *td, int *doy, float *alt);
    // Perez model (Igh -> Ibn)
    static double dirint_(float *g, float *z__, float *td, float *alt, float *i0);

    // gets and erases the global horizontal irradiance and illuminance
    void setIgh(float value) { Igh.push_back(value); }
    void addIgh(float value) { if (Igh.empty()) Igh.push_back(value); else Igh.back()+=value; }
    float getIgh() { if (Igh.empty()) return 0.f; else return Igh.back(); }
    float getIgh(unsigned int step) { if (Igh.empty()) return 0.f; else return Igh.at(step); }
    float getIgh(unsigned int day, unsigned int hour) { if (Igh.empty()) return 0.f; else return Igh.at((day-1)*24 + hour -1); }
    void eraseIgh(unsigned int keepValue) { Igh.erase(Igh.begin(),Igh.end()-min(keepValue,(unsigned int)Igh.size())); }

    void setIgh_vis(float value) { Igh_vis.push_back(value); }
    bool isIgh_vis_empty() { return Igh_vis.empty(); }
    float getIgh_vis() { if (Igh_vis.empty()) return 0.f; else return Igh_vis.back(); }
    float getIgh_vis(unsigned int step) { if (Igh_vis.empty()) return 0.f; else return Igh_vis.at(step); }
    float getIgh_vis(unsigned int day, unsigned int hour) { if (Igh_vis.empty()) return 0.f; else return Igh_vis.at((day-1)*24 + hour -1); }
    void eraseIgh_vis(unsigned int keepValue) { Igh_vis.erase(Igh_vis.begin(),Igh_vis.end()-min(keepValue,(unsigned int)Igh_vis.size())); }

    // sets the values from the .cli2
    bool isEmptyCli2() { return KeCoeff1.empty(); }
    void loadCli2(string filename);
    float getKeCoeff1(unsigned int surfaceId, unsigned int day, unsigned int hour) { return KeCoeff1[surfaceId].at((day-1)*24 + hour -1); }
    float getKeCoeff2(unsigned int surfaceId, unsigned int day, unsigned int hour) { return KeCoeff2[surfaceId].at((day-1)*24 + hour -1); }
    float getKeCoeff3(unsigned int surfaceId, unsigned int day, unsigned int hour) { return KeCoeff3[surfaceId].at((day-1)*24 + hour -1); }
    float getTla(unsigned int surfaceId, unsigned int day, unsigned int hour) { return Tla[surfaceId].at((day-1)*24 + hour -1); }
    void loadCli3(string filename);

};
#endif
