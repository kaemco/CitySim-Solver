#ifndef SCENE_H
#define SCENE_H

#include <map>
#include <cmath>
#include <cassert>

#include "VFCViewFactorCalculation.h"
#include "DATARadiationScene.h"
#include "DATASurfaceDelegateABC.h"
#include "DATASurfaceIterator.h"
#include "DATAViewFactorSetSparse.h"
#include "DATASurface.h"
#include "DATAInsolationFactors.h"
#include "SKYSun.h"
#include "GEOMPolygonInfo.h"// to compute polygon area
#include "cPerezSkyModel.h" // includes sun.h, for the Perez sky model

#include "SKYTregenza.h" // to the the patches information

#include "climate.h"
#include "tinyxml.h"
#include "district.h"
#include <iostream>

#ifdef FMI
#include "FMILibrary/include/fmilib.h"
#endif

using namespace std;

// *** Scene class, CitySim  *** //
// *** jerome.kaempf@epfl.ch *** //

class Scene {

protected:

    // the input file (.rad or .xml)
    string inputFile;

    // the climate
    Climate *pClimate = nullptr;

    // the sun, the Tregenza vault and the Perez sky model
    SKYSun *pSun = nullptr;
    cPerezSkyModel sky;
    // creates the SKY in form of Tregenza Patches
    SKYTregenza tregenzaSky;
    // stores the actual radiance of the sky and ground (to be shared by different methods)
    vector<float> lv = vector<float>(tregenzaSky.PatchCount()/2, 0.f);
    float groundRadiance = 0.f;

    // for the VFC
	DATARadiationScene scene;
	VFCViewFactorCalculation v;

    // the reflection matrix is a Compressed Row Storage, which includes three vectors:
    vector<float> An;         // the values in a row wise fashion
    vector<unsigned int> Aj;  // the indices of the row of the elements
    vector<unsigned int> Ai;  // the indices of the first element of each row

    // define the number of inter-reflections in the irradiance calculation
    unsigned int mNbReflections = 2;

public:

    const int daysPerMonth[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    ostream logStream;
    // DP: This stream captures most of the informative log from the scene, including logs from the district, building, zone and surface objects.
    // It is by default directed on the console but can be used to direct this log to a file according to the need, independently of cout or cerr.
    // However it does not (yet?) capture logs from occupants, plants and other libraries objects such as DATARadiationScene or VFCDirect, which are still logged on "cerr"
    Scene():logStream(std::cout.rdbuf()) {}
    Scene(string inputFile, string climateFile);
    virtual ~Scene() {
        //logStream << "Destructor of Scene." << endl << flush;
        // NOTE: the destructor is virtual so that we can call delete directly on the Scene object and it will call the correct destructor for the derived class XmlScene
        if (pClimate) delete pClimate;
        if (pSun) delete pSun;
    }
    void clear() {
        pClimate->clear();
        lv.assign(tregenzaSky.PatchCount()/2, 0.f);
        groundRadiance = 0.f;
    }

    // get the pointer to the climate class
    Climate* getClimate() { return pClimate; }

    // get the reflection matrix components
    float getAn(unsigned int index) { return An[index]; }
    unsigned int getAj(unsigned int index) { return Aj[index]; }
    unsigned int getnAi() { return Ai.size(); }
    unsigned int getAi(unsigned int index) { if (index >= Ai.size()) return An.size(); else return Ai[index]; }
    VFCViewFactorCalculation* getV(){ return &v;}
    DATARadiationScene* getDATARadiationScene(){ return &scene;}

    // method to compute the view factors
    void computeViewFactors();
    void computeProjectedSolidAngles();

    // method the create the reflection matrices
    void buildSparseMatrix();
    void showViewFactors();
    // export files in a Radiance format to check
    void exportRadFile(string radFile="", bool triangulated=false);
    void exportInpFile(string radFile="", bool buildingsOnly=false);
    virtual void exportDXF(string fileName="");
    // computes the sky and ground (diffuse) radiance, in all generality
    void computeRadiance(const unsigned int& day, const float& Idh, const float& Ibn, const float& albedo=0.2f);
    // computes the cumulative sky and ground radiance
    void computeCumulativeRadiance(unsigned int beginDay=1, unsigned int endDay=365, float albedo=0.2f);
    // saves the current sky and ground radiance in a file for radiance
    void exportSkyAndGround(string radFile,float luminousEfficacy=1.f);
    // exports the cumulative sky file and the cumulative sun file
    virtual void exportCumulativeRadiance();
};

class Radscene : public Scene {

private:

    vector<float> irradiationDiffuseSky,irradiationDiffuseGround,irradiationBeam,irradiationReflection;
    // outData for the comparison with Radiance
    ostringstream outData;

public:

    Radscene(string inputFile, string climateFile);

    // clear the results
    void clearResults();
    // exports the results
    void exportSWFile(string filename);
    // computes the shortwave
    void computeShortWave(unsigned int day, unsigned int hour);
    void computeDaylight(unsigned int day, unsigned int hour);
    // export the Sun in a Radiance file
    void exportSunRadFile(string sunRadFile,unsigned int day,unsigned int hour);
    // export the diffuse sky in Tregenza zones in .cal and .rad files
    void exportSkyRadFile(string skyRadFile,unsigned int day,unsigned int hour);
    // export the diffuse ground
    void exportGroundRadFile(string groundRadFile,unsigned int day, unsigned int hour);
    // creates the Radiance files and start the simulation to finally save everything in a nice table
    void compareWithRadianceExternalIrradiance(unsigned int day, unsigned int hour);
    void withRadianceClimate(string climateFile);
    // reads results
    void readResults(string filename, vector<float> &results);
    // writes results
    void writeHeader(string filename);
    void writeResults(string filename);
    // subdivision of the surfaces in Triangles
    void gridTriangle(vector<GENPoint> triangle, double maxDetectorArea, vector<GENPoint> &grid);
    GENPoint centroid(vector<GENPoint> &shape);
    GENPoint middlePoint(GENPoint &point1, GENPoint &point2);
    // average irradiance on surfaces
    void averageIrradianceOnSurfaces(vector<unsigned int> count,vector<float> &irrad);

};

class XmlScene : public Scene {

private:

    string climateFile;

    // the begin and end of the simulation
    unsigned int beginDay=1, endDay=365;

    // simulation time steps
    unsigned int timeStepsSimulated = 0, timeSteps2Simulated = 0, simulationIndex = 0;
    // pre-conditioning period simulation
    unsigned int preTimeStepsSimulated = 0, preTimeSteps2Simulated = 0;

    // the district itself
    District* pDistrict;

    // stores the actual surfaces' irradiation from beam and from the sky and ground
    vector<float> irradiationBeam, irradiationDiffuseSky, irradiationDiffuseGround;

    // far field obstructions
    vector<double> farFieldOccludedPatchFraction;

    #ifdef FMI
    // the FMU
    fmi1_import_t* fmu;
    fmi_import_context_t* context;
    fmi1_status_t fmistatus;
    jm_status_enu_t jmstatus;
    #endif

#ifdef DEBUG
    stringstream ss_IAM, ss_IAM_irradiance;
#endif // DEBUG

public:

    XmlScene():Scene(),pDistrict(new District(this)) {}
    XmlScene(string inputFile, ostream* pLogFileStr=NULL, bool climateFileRequired=true);

    ~XmlScene();
    void clear() {
        pDistrict->clear();
        timeStepsSimulated = 0;
        timeSteps2Simulated = 0;
        simulationIndex = 0;
        preTimeStepsSimulated = 0;
        preTimeSteps2Simulated = 0;
        irradiationBeam.clear();
        irradiationDiffuseSky.clear();
        irradiationDiffuseGround.clear();
        Scene::clear();
    }

    void addAllSurfacesToScene();

    District* getDistrict() {return pDistrict;}
    void setInputFile(string inputFile) { this->inputFile = inputFile; }
    string getInputFile() { return inputFile; }
    string getInputFileNoExt() { return inputFile.substr(0,inputFile.size()-4); }
    string getInputFileNoExtNoPath() { return inputFile.substr(inputFile.find_last_of('/')+1,inputFile.size()-4); }

    unsigned int getBeginDay(){return beginDay;}
    unsigned int getEndDay(){return endDay;}
    unsigned int getPreTimeStepsSimulated(){ return preTimeStepsSimulated;}
    unsigned int getTimeStepsSimulated(){ return timeStepsSimulated;}
    unsigned int getSimulationIndex(){ return simulationIndex;}

    void exportXMLFile(string fileName="");
    void exportGML(string fileName="", const vector<double>& origin={0.,0.,0.});
    void exportDXF(string fileName="");
    void exportSTL(string fileName="") { exportSTL_binary(fileName); }
    void exportSTL_ascii(string fileName="");
    void exportSTL_binary(string fileName="");
    void exportRAD(string fileName="");

    string getClimateFile(){return climateFile;}
    void readClimate(string fileName);
    void importClimatePVGIS(string fileName, int defaultCloudiness);
    void importHorizonPVGIS(string fileName);
    void setHorizon();

    // compute the far field obstructions
    void initialiseFarField();
    void computeFarField();
    bool sunVisibleFarField(float sunAzimuth, float sunElevation);

    // computes the short wave exchanges
    void computeShortWave_Beam(unsigned int day, unsigned int hour);
    void computeShortWave_Diffuse();
    void computeShortWave_Interreflected();
    void computeShortWave(unsigned int day, unsigned int hour);
    // computes the internal illuminance levels (2 values)
    void computeDaylight(unsigned int day, unsigned int hour);
    // computes the longwave exchange
    void computeLongWave(unsigned int day, unsigned int hour);
    // computes the thermal in the buildings
    int computeWarmUp(); // number of days for the pre-conditioning period
    void computeThermal(unsigned int day, unsigned int hour);
    #ifdef FMI
    void initialiseThermal_EnergyPlus();
    fmi1_value_reference_t getValueReference(string variableName);
    void computeThermal_EnergyPlus(unsigned int day, unsigned int hour);
    void terminateThermal_EnergyPlus();
    #endif
    // computes the comfort indices outdoors
    void computeComfort(unsigned int day, unsigned int hour);
    // simulate
    void simulate();
    // simulate only the irradiation on the surfaces
    void simulateCumulativeIrradiance();
    void simulateRadiation(bool daylight=true);

    // decomposed simulation steps
    void initialiseSimulationParameters();
    void updateSimulationModelParameters();
    void simulateTimeStep(int day,int hour, bool preCond = false, bool radOnly = false, bool doDayLightSim=true);
    void clearPreConditioningResults();
    void checkMemoryUsage(double limit, bool radOnly=false);

    // erase the last result from memory in the thermal results
    void eraseResultsThermal_back();
    // erase the results from the memory (all but the last one)
    void eraseResults(unsigned int keepValue, bool eraseAllResults);
    void eraseResultsIrradiation(unsigned int keepValue);
    // writes the files at the end of the simulation
    unsigned int getColumnIndex(Surface* surface);
    void writeSWHeaderText(string fileOut, string unit="Irradiance(W/m2)");
    void writeSWResultsText(string fileOut); // in W/m2
    void writeSWResultsBinary(string fileOut);
    void writeSWvHeaderText(string fileOut);
    void writeSWvResultsText(string fileOut); // in lux
    void writeDLHeaderText(string fileOut);
    void writeDLResultsText(string fileOut); // in lux
    void writeLWHeaderText(string fileOut);
    void writeLWResultsText(string fileOut); // in W/m2
    void writeTHHeaderText(string fileOut);
    void writeTHResultsText(string fileOut); // exports temperature, heating and cooling
    void writeTHExplicitHeaderText(string fileOut);
    void writeTHExplicitResultsText(string fileOut);
    void writeTSHeaderText(string fileOut);
    void writeTSResultsText(string fileOut);
    void writeHCHeaderText(string fileOut);
    void writeHCResultsText(string fileOut);
    void writeETHeaderText(string fileOut);
    void writeETResultsText(string fileOut);
    void writeCMHeaderText(string fileOut);
    void writeCMResultsText(string fileOut);
    void writeClimaticDataText(string fileOut);
    void writeMonthlyResultsText(string fileOut);
    void writeYearlyResultsText(string fileOut);
    void writeYearlyResultsPerBuildingText(string fileOut);
//    void writeYearlyResultsPerDECText(string fileOut); //added by Dapeng // Cognet: Deleted it (never used?).
    map<unsigned int,vector<double> > getHeatingHourlyResultsPerBuilding();
    map<unsigned int,vector<double> > getHeatingMonthlyResultsPerBuilding();
    map<unsigned int,double> getHeatingYearlyResultsPerBuilding();
    map<unsigned int,vector<double> > getCoolingHourlyResultsPerBuilding();
    map<unsigned int,vector<double> > getCoolingMonthlyResultsPerBuilding();
    map<unsigned int,double> getCoolingYearlyResultsPerBuilding();
    map<unsigned int,vector<double> > getMachinePowerHourlyResultsPerBuilding();
    map<unsigned int,vector<double> > getMachinePowerMonthlyResultsPerBuilding();
    map<unsigned int,double> getMachinePowerYearlyResultsPerBuilding();
    void writeAreaText(string fileOut);
    void writeVFText(string fileOut);
    void writeInertiaText(string fileOut);
    void exportCumulativeRadiance();
    void exportHourlyRadiance();
    // comparison with Radiance for the internal illuminance
    void exportSkyRadFile(string skyRadFile, float *lv);
    // computes the memory usage for the XmlScene
    size_t memoryUsage();
    size_t memoryUsageIrradiation();

};
#endif
