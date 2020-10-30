//File: CitySimWebS.cpp
#include "soapH.h"
#include "CitySimWebS.nsmap"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <stdlib.h>
#include <ctime>
#include <map>
#include <vector>
#include <sys/stat.h>
#include <dirent.h>
#include "CitySimWebS.h"
#include "../../scene.h"

/*
Autor : Renaud Sauvain
Date : 20.03.2011
Content : Logique for the CitySimWebS gsoap web service
*/

/// NOTE: multi-threading code comes from www.cs.fsu.edu/~engelen/soapdoc2.html#sec:mt
/// NOTE2: as CitySim is not thread-safe yet, multi-threading code has been disabled

	/* Granularity definitions */
    enum ResultsGranularity { Monthly, Yearly, Hourly };
    static const double NA=-999.0;

	/* static datas such as the log file location & the output files */
	const std::string tempDirectory = "/var/tmp/dev/";
	const std::string logFile = tempDirectory + "citysimWS.log";
    ofstream WSLogStream;


	class ns__NormMaps
    {
    public:
        std::map<unsigned int, double> hotwater;
        std::map<unsigned int, double> electricity;
        ns__NormMaps(const ns__NormsW* nw);
        ns__NormMaps(){};
    };

    /* Function prototypes */
    int filesize(string filename);
    void cleanLogDir ();
	std::string buildNameFromTime(int sec);
	bool FileExists(string strFilename);
	std::string readFile(std::string fileURI);
	bool writeFile(std::string buildingXML,std::string fileURI);
	std::map<unsigned int,std::vector<double> > makeVectorMap(std::map<unsigned int,double> results);
	std::map<unsigned int, double> computeDemand(std::map<unsigned int, double> normValues, ns__DistrictW const& dw);
	void buildResults(std::map<unsigned int,std::vector<double> > heatingR, std::map<unsigned int,std::vector<double> > coolingR, std::map<unsigned int, double> hotwaterResults, std::map<unsigned int, double> electricityResults, ns__ResultsW &results);
    ns__NormMaps* buildNormMaps(const ns__NormsW* nw);
	int simulate(struct soap *soap, const ns__DistrictW* district, const ns__ClimateW* climate, const ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_, ResultsGranularity rGranularity, std::string infile);


    ns__NormMaps* buildNormMaps(const ns__NormsW* nw)
    {
        ns__NormMaps* n = new ns__NormMaps();
        for(unsigned int i = 0 ; i < nw->__sizeHotWater ; ++i)
            {
            //ns__NormPairW phw = nw->hotwater[i];
            //ns__NormPairW pe = nw->electricity[i];
            //unsigned int allocId = phw.allocId;
            n->hotwater.insert( std::pair<unsigned int, double> ( nw->hotwater[i].allocId, nw->hotwater[i].value) );
            n->electricity.insert( std::pair<unsigned int, double> ( nw->electricity[i].allocId, nw->electricity[i].value) );
        }
        return n;
    }

    ns__NormMaps::ns__NormMaps(const ns__NormsW* nw)
    {
        for(unsigned int i = 0 ; i < nw->__sizeHotWater ; ++i)
        {
            hotwater.insert( std::pair<unsigned int, double> ( nw->hotwater[i].allocId, nw->hotwater[i].value) );
            electricity.insert( std::pair<unsigned int, double> ( nw->electricity[i].allocId, nw->electricity[i].value) );
        }
    }

//	void* process_request(void* soap);

	/* Initialization of the gsoap web service */
	int main()
	{
		// create soap context and serve one CGI-based request: Uncomment if use of CGI mod instead of stand-alone
		// return soap_serve(soap_new());

        //redirect the error log to the log file
        WSLogStream.open(logFile.c_str(),ios_base::app);
        // Also capture the "cerr" logs in WSLogStream
        std::cerr.rdbuf(WSLogStream.rdbuf());
        //freopen(std::string(logFile).c_str(), "w", stderr);

		//comment all following if CGI mod
		struct soap soap;
		int m, s; // master and slave sockets
		soap_init(&soap);
		m = soap_bind(&soap, NULL , 9201, 100); // Null means local server URL and can be replaced by "lesopbpc27.epfl.ch"
		if (m < 0)
			soap_stream_fault(&soap, WSLogStream);
		else
		{
//			pthread_t tid;
			WSLogStream << "Socket connection successful: master socket = " << m << endl << flush;
			struct soap *tsoap;
			for (int i = 1; ; i++)
			{
				s = soap_accept(&soap);
				if (s < 0)
				{
					soap_stream_fault(&soap, WSLogStream);
					break;
				}
                // Limit logFile size
                if(filesize(logFile)>10000){
                    WSLogStream.close();
                    WSLogStream.open(logFile.c_str(),ios_base::app);
                }
                // Log time
                char buffer[BUFSIZ] = { '\0' };
                time_t now = time(NULL);
                struct tm * local_time=localtime( &now );
                strftime( buffer, BUFSIZ, "%x %X", local_time);
                WSLogStream << endl << "***** " << buffer << " *****" << endl << flush;

				WSLogStream << "Soap: accepted connection from IP= " << flush;
				WSLogStream << static_cast<int>((soap.ip >> 24)&0xFF) << "." << static_cast<int>((soap.ip >> 16)&0xFF) << "." << static_cast<int>((soap.ip >> 8)&0xFF) << "." << static_cast<int>(soap.ip&0xFF) << flush;
				WSLogStream << " socket=" << s << " (i=" << i << ")" << endl << flush;
//				tsoap = soap_copy(&soap); // make a safe copy
//				if(!tsoap) break;
//				pthread_create(&tid, NULL, (void*(*)(void*))process_request, (void*)tsoap);

				WSLogStream << "Soap: request is being processed" << endl << flush;
				if (soap_serve(&soap) != SOAP_OK) // process RPC request
					soap_stream_fault(&soap, WSLogStream); // print error
                else
                    WSLogStream << "Soap: request executed" << endl << flush;
				soap_destroy(&soap); // clean up class instance
                //WSLogStream << "Soap " << i << ": soap_destroy executed" << endl << flush;
				soap_end(&soap); // clean up everything and close socket
                WSLogStream << "Soap: soap_destroy and soap_end executed" << endl << flush;
			}
		}
		soap_done(&soap); // close master socket and detach context
	}



	int ns__MonthlySimulation(struct soap *soap, ns__DistrictW* district, ns__ClimateW* climate, ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_)
	{
        int result;
		cerr << "ns__MonthlySimulation: simulate" << endl;
        result = simulate(soap, district, climate, norms, userID, return_, Monthly,"");
		cerr << "ns__MonthlySimulation: return" << endl;
        return result;
	}
	int ns__YearlySimulation(struct soap *soap, ns__DistrictW* district, ns__ClimateW* climate, ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_)
	{
		return simulate(soap, district, climate, norms, userID, return_, Yearly,"");
	}
	int ns__HourlySimulation(struct soap *soap, ns__DistrictW* district, ns__ClimateW* climate, ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_)
	{
	    int res;
		cerr << "ns__HourlySimulation: simulate" << endl;
		res = simulate(soap, district, climate, norms, userID, return_, Hourly, "");
		cerr << "ns__HourlySimulation: return" << endl;
		return res;
	}

	/* -------------------------------- Simulation call methods -------------------------------- */

	int simulate(struct soap *soap, const ns__DistrictW* district, const ns__ClimateW* climate, const ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_, ResultsGranularity rGranularity, std::string infile)
	{
	    cleanLogDir ();

	    std::string simulationName;
	    std::string logFileName;
	    int sec=0;
	    do{
            //get a unique simulation name
            simulationName = buildNameFromTime(sec);
            //get a unique logFile name
            logFileName = tempDirectory+"CSsim_"+simulationName+".log";
            ++sec;
	    } while(FileExists(logFileName));

        //open a logStream
        ofstream simLogStream(logFileName.c_str());
        if(simLogStream.fail()){
            WSLogStream << "Error: Could not open simulation logStream" << endl << flush;
            exit;
        }

        // Write simulation log file name in WS log
	    WSLogStream << "simulate: CitySim simulation log file " << logFileName << endl << flush;

	    // log time
		char buffer[BUFSIZ] = { '\0' };
		time_t now = time(NULL);
		struct tm * local_time=localtime( &now );
		strftime( buffer, BUFSIZ, "%x %X", local_time);
		simLogStream << buffer << endl << flush;

		//transform user id to string
		stringstream st; st << userID; string sUserID = st.str();

	    simLogStream << "WS simulate start" << endl;
		 // if client supports gzip, enable the compression mod
		if (soap->zlib_out == SOAP_ZLIB_GZIP)
			soap_set_omode(soap, SOAP_ENC_ZLIB);
		soap_set_omode(soap, SOAP_IO_CHUNK); //enable HTTP chunking for the outputs messages

		// use normal console output

		//fprintf(stderr, "simulate \r\n");
		try{
		    // *************************************************************
		    // **** UNITS: energy results in Wh for heating and cooling ****
		    // ****                          kWh for hotwater and elec  ****
		    // *************************************************************

            // Estimating hotwater demand and electricity demand based on norms
			simLogStream << "WS Instancing Norms" << endl;
			ns__NormMaps * p_normMaps = new ns__NormMaps(norms);
			//ns__NormMaps * p_normMaps = buildNormMaps(norms);
	        WSLogStream << "simulate: Computing hotwater demand... " << flush;
			std::map<unsigned int, double> hotwaterR = computeDemand(p_normMaps->hotwater, *district);
			WSLogStream << "Done" << endl << flush;
            WSLogStream << "simulate: Computing electricity demand... " << flush;
			std::map<unsigned int, double> electricityR = computeDemand(p_normMaps->electricity, *district);
			WSLogStream << "Done" << endl << flush;
			simLogStream << "WS Norm results ok" << endl;

			// Create scene instance
			simLogStream << "WS Instancing Scene" << endl;
			XmlScene* p_xmlScene;
			if(infile != ""){
				std::string fileUri = tempDirectory+simulationName+ "_" + sUserID + "_.in";
				writeFile(infile, fileUri);
				p_xmlScene = new XmlScene(fileUri,&simLogStream);
			}
			else{
				//parameters initialisation
				//unsigned int beginDay=1, beginMonth=1, endDay=31, endMonth=12;
				p_xmlScene = new XmlScene(simulationName, climate, district, &electricityR, &simLogStream);
			}

			//Writing input for debuging
			p_xmlScene->exportXMLFile(tempDirectory+"simulationModel.xml");

			WSLogStream << "simulate: Computing heating and cooling demand... " << flush;
			simLogStream << "WS Launching Simulation" << endl << flush;
			p_xmlScene->simulate();
			WSLogStream << "Done" << endl << flush;

			simLogStream << "WS Simulate Done" << endl << "WS Getting Results" << endl << flush;

			//getting the outputs depending on the needed granularity
			std::map<unsigned int,std::vector<double> > heatingRV;
			std::map<unsigned int,std::vector<double> > coolingRV;

			switch (rGranularity)
			{
			case Monthly:
				heatingRV = p_xmlScene->getMachinePowerMonthlyResultsPerBuilding();
				coolingRV = p_xmlScene->getCoolingMonthlyResultsPerBuilding();
				break;
			case Yearly:
				heatingRV = makeVectorMap(p_xmlScene->getMachinePowerYearlyResultsPerBuilding());
				coolingRV = makeVectorMap(p_xmlScene->getCoolingYearlyResultsPerBuilding());
				break;
			case Hourly:
				heatingRV = p_xmlScene->getMachinePowerHourlyResultsPerBuilding();
				coolingRV = p_xmlScene->getCoolingHourlyResultsPerBuilding();
				break;
			}
			simLogStream << "WS CitySim results ok" << endl;

			delete p_normMaps;
			delete p_xmlScene;
			simLogStream << "WS Deleted scene" << endl;

			//construct output objects
            buildResults(heatingRV, coolingRV, hotwaterR, electricityR, return_);
			simLogStream << "WS return_ defined" << endl;
		}
		catch (std::exception &why) {
			std::ostringstream sst;//create a stringstream
			sst << "Caught exception : " << why.what() << "\r\n";
			WSLogStream << sst.str() << flush;
            simLogStream << sst.str() << flush;
			return_.comments = sst.str();
			return soap_receiver_fault(soap, sst.str().c_str(), NULL); // return fault to sender
		}
		catch(std::string str){
			WSLogStream << "WS Caught exception : " << str << endl << flush;
			simLogStream << "WS Caught exception : " << str << endl << flush;
			return_.comments = str;
			return soap_receiver_fault(soap, str.c_str(), NULL); // return fault to sender
		}
		simLogStream << "WS simulate return" << endl << flush;
		return SOAP_OK;
	}



/*
	void* process_request(void* soap)
	{
		pthread_detach(pthread_self());
		soap_serve((struct soap*)soap);
		soap_destroy((struct soap*)soap); // dealloc C++ data
		soap_end((struct soap*)soap); // dealloc data and clean up
		soap_done((struct soap*)soap); // detach soap struct
		free(soap);
		return NULL;
	}
*/

	/* -------------------------------- Tools -------------------------------- */
    /* gets the file size */
    int filesize(string filename)
    {
        std::ifstream in(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        in.seekg(0, std::ifstream::end);
        return (int) in.tellg();
    }

    /* clear old log files */
    void cleanLogDir ()
    {
        DIR *dp;
        struct dirent *dirp;
        string fileName;
        if((dp  = opendir(tempDirectory.c_str())) == NULL) {
            WSLogStream << "Error opening " << tempDirectory << endl << flush;
            return;
        }

		char buffer[BUFSIZ] = { '\0' };
        time_t yesterday = time(NULL)-60*60*24;
        struct tm * local_time = localtime ( &yesterday );
		strftime( buffer, BUFSIZ, "%Y%m%d%H%M%S", local_time);
		std::stringstream sst;//create a stringstream
		sst << buffer;
        string yesterdayLogName = sst.str();
        while ((dirp = readdir(dp)) != NULL) {
            fileName= string(dirp->d_name);
            if(fileName.substr(0,6)=="CSsim_" && fileName.substr(20,4)==".log" && fileName.substr(6,14) < yesterdayLogName){
                remove((tempDirectory+fileName).c_str());
            }
        }
    }

	/* make a string from current time */
	std::string buildNameFromTime(int sec=0)
	{
	    //cerr << "buildNameFromTime" << endl;
		char buffer[BUFSIZ] = { '\0' };
		time_t now = time(NULL)+sec;
		struct tm * local_time=localtime( &now );
		strftime( buffer, BUFSIZ, "%Y%m%d%H%M%S", local_time);
		std::stringstream sst;//create a stringstream
		sst << buffer;
		//delete local_time;
	    // << "buildNameFromTime" << endl;
		return sst.str();
	}

	/* check if a file is present or not */
	bool FileExists(string strFilename) {
	  struct stat stFileInfo;
	  bool blnReturn;
	  int intStat;

	  // Attempt to get the file attributes
	  intStat = stat(strFilename.c_str(),&stFileInfo);
	  if(intStat == 0) {
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	  } else {
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	  }
	  return(blnReturn);
	}

	/* read a file like the results */
	std::string readFile(std::string fileURI)
	{
		std::string result= "";
		std::string line;
		std::ifstream read (fileURI.data());//reading a file
		if (read.is_open())
		{
			while (! read.eof() )
			{
				getline (read,line);
				result += line;
			}
			read.close();
			return result;
		}
		else
		{
			fprintf(stderr, "Unable to open file");
			return "Unable to open file";
		}
	}

	/* file writing function */
	bool writeFile(std::string buildingXML,std::string fileURI)
	{
		fprintf(stderr, "Writing input file in %s \r\n",fileURI.data());
		std::ofstream write (fileURI.data());//writing to a file
		if (write.is_open())
		{
			write << buildingXML;
			write.close();
			return true;
		}
		else
		{
			return false;
		}
	}

	/* make a std::map<unsigned int,std::vector<double> > from a std::map<unsigned int,double> */
	std::map<unsigned int,std::vector<double> > makeVectorMap(std::map<unsigned int,double> results)
	{
		std::map<unsigned int,std::vector<double> > outputMap;
		int i = 0;
		std::map<unsigned int, double>::const_iterator iter;

		for (iter=results.begin(); iter != results.end(); ++iter) {
			std::vector<double> vect;
			vect.push_back(iter->second);
			outputMap.insert ( std::pair<unsigned int,std::vector<double> >(iter->first,vect) );
		}
		return outputMap;
	}

	/* make a <building_id, annual_demand> map given the norms value */
	/* can be used with both hotwater and electricity demand norms */
	/* Notice : this simulation part should be move inside the CitySim library itself */
	std::map<unsigned int, double> computeDemand(std::map<unsigned int, double> normValues, ns__DistrictW const& dw)
	{
		std::map<unsigned int, double> outputMap;
		for(unsigned int i = 0 ; i < dw.__sizeBuildings ; ++i)
		{
			if(!dw.buildings[i].ignoreResults){
                outputMap.insert( std::pair<unsigned int, double> (
                    dw.buildings[i].id,  (normValues[dw.buildings[i].mainAllocationId])*(dw.buildings[i].sre)
                ));
                //WSLogStream << "norm results: " << dw.buildings[i].id <<" "<< normValues[dw.buildings[i].mainAllocationId] << "*" << dw.buildings[i].sre << endl;
			}
		}
		return outputMap;
	}

	/* build the output class ns__ResultsW from a std::map<unsigned int,std::vector<double> > */
	void buildResults(std::map<unsigned int,std::vector<double> > heatingR, std::map<unsigned int,std::vector<double> > coolingR, std::map<unsigned int, double> hotwaterResults, std::map<unsigned int, double> electricityResults, ns__ResultsW &results)
	{
	    WSLogStream << "simulate: Building results object" << endl;
		//ns__ResultsW results;
		results.__size = hotwaterResults.size();
		results.comments = "";
		results.__ptr = new ns__BuildingResultW[results.__size];

        int i = 0;
		//matching heat & cold result
		if(results.__ptr == 0)
			fprintf(stderr, "Error: memory could not be allocated");
		else
		{
			for (std::map<unsigned int,double >::const_iterator iter=hotwaterResults.begin(); iter != hotwaterResults.end(); ++iter) {
			// hotwaterResults does not contain values for "ignoreResults" buildings -> ignored for heating and cooling too
				results.__ptr[i].building_ID = iter->first;
				results.__ptr[i].__size = heatingR[iter->first].size();
				ns__TCRW * ptr = new ns__TCRW[results.__ptr[i].__size];

				for(int k=0;k<results.__ptr[i].__size;++k)
				{
					ptr[k].c = coolingR[iter->first][k];
					ptr[k].h = heatingR[iter->first][k];
					//ptr[k].hw = hotwaterResults[iter->first];
					//ptr[k].e = electricityResults[iter->first];
				}
				results.__ptr[i].WhConsumptions = ptr;
				results.__ptr[i].yearlykWhConsumptionE = electricityResults[iter->first];
				results.__ptr[i].yearlykWhConsumptionHW = hotwaterResults[iter->first];
				//WSLogStream << "Building " << results.__ptr[i].building_ID  << " electricity demand: " << electricityResults[iter->first] << " " << results.__ptr[i].yearlykWhConsumptionE << endl;
				++i;
			}
		}
		WSLogStream << "simulate: " << i << " buildings' results in Results" << endl << flush;
		return; // results;
	}

/*  int ns__HourlyXMLSimulation(struct soap *soap, std::string buildingXML, unsigned int userID, ns__ResultsW &return_)
	{
		cerr << "ns__HourlyXMLSimulation call simulate" << endl;
		ns__DistrictW district;
		ns__ClimateW climate;
		ns__NormsW norms;
		return simulate(soap, &district, &climate, &norms, userID, return_, Hourly, buildingXML);
	}
	int ns__MonthlyXMLSimulation(struct soap *soap, std::string buildingXML, unsigned int userID, ns__ResultsW &return_)
	{
		cerr << "ns__MonthlyXMLSimulation call simulate" << endl;
		ns__DistrictW district;
		ns__ClimateW climate;
		ns__NormsW norms;
		cerr << "ns__MonthlyXMLSimulation return" << endl;
		return simulate(soap, &district, &climate, &norms, userID, return_, Monthly, buildingXML);
	}
	int ns__YearlyXMLSimulation(struct soap *soap, std::string buildingXML, unsigned int userID, ns__ResultsW &return_)
	{
		ns__DistrictW district;
		ns__ClimateW climate;
		ns__NormsW norms;
		return simulate(soap, &district, &climate, &norms, userID, return_, Yearly, buildingXML);
	}
*/
	// ns__constructeurs et ns__destructeurs pour gestion correcte de la mémoire et des pointeurs (defs dans CitySimWebS.h)

    ns__BuildingResultW::ns__BuildingResultW():WhConsumptions(NULL) {}

    void ns__BuildingResultW::copy(ns__BuildingResultW const& br){
        building_ID = br.building_ID;
		__size = br.__size;
		WhConsumptions = new ns__TCRW[__size];
		for(int i=0; i < __size ; ++i)
			WhConsumptions[i] = br.WhConsumptions[i];
    }
    ns__BuildingResultW::ns__BuildingResultW(const ns__BuildingResultW& s){
        cerr << "Copying ns__BuildingResultW: should not happen" << endl; // Basically this should not happen.
        copy(s);
    }
    ns__BuildingResultW& ns__BuildingResultW::operator=(ns__BuildingResultW const& s){
        cerr << "Equaling ns__BuildingResultW: should not happen" << endl;
        copy(s);
        return *this;
    }
    ns__BuildingResultW::~ns__BuildingResultW(){
        //cerr << "ns__BuildingResultW destructor" << endl;
        delete[] WhConsumptions;
    }


    ns__ResultsW::ns__ResultsW():__ptr(NULL) {}//cerr << "Default ns__ResultsW" << endl;}

    void ns__ResultsW::copy(ns__ResultsW const& r){
        comments = r.comments;
		__size = r.__size;
		__ptr = new ns__BuildingResultW[__size];
		for(int i=0; i < __size ; ++i)
			__ptr[i] = r.__ptr[i];
    }

    ns__ResultsW::ns__ResultsW(const ns__ResultsW& s){
        WSLogStream << "Copying ns__ResultsW: should not happen" << endl; // Basically this should not happen.
        copy(s);
    }

    ns__ResultsW& ns__ResultsW::operator=(ns__ResultsW const& s){
        WSLogStream << "Equaling ns__ResultsW: should not happen" << endl;
        copy(s);
        return *this;
    }
    ns__ResultsW::~ns__ResultsW(){
        //WSLogStream << "ns__ResultsW destructor" << endl;
        delete[] __ptr;
    }


   ns__NormsW::ns__NormsW(){
        //WSLogStream << "Default ns__NormsW" << endl;
    }
   ns__NormsW& ns__NormsW::operator=(ns__NormsW const& n){
        WSLogStream << "Equaling ns__NormsW: should not happen" << endl;
        return *this;
    }
   ns__NormsW::~ns__NormsW(){
        //WSLogStream << "ns__NormsW destructor" << endl;
    }

/* Apparently these CitySim object constructors need to be defined here, to avoid compilation incompatibilities */
/* ns__...W objects are not understood correctly in code compiled through CitySim's Makefile.webservice         */


    /* Building class builder, using the XML wrapper defined in CitySimWebS.h */

    Building::Building(const ns__BuildingW* buildingXMLWrapper, District* pDistrict, double yearlyElecCons, unsigned int heatingBeginDay, unsigned int heatingEndDay):pDistrict(pDistrict),logStream(pDistrict->logStream.rdbuf()) {
        // yearlyElecCons is in kWh

        // initialises the pointers to NULL and values
        nre = 0.f;
        gwp = 0.f;
        ubp = 0.f;
        blindsLambda = 0.017f; // shape of the sigmoid curve
        blindsIrradianceCutOff = buildingXMLWrapper->blindsIrrCutOff; // 1376 W/m2 = no cut-off; 100 W/m2 = Jérôme office cut-off

        // gets the buildings ID (must be unique)
        id = buildingXMLWrapper->id;
        key = "";
        logStream << "Building: " << id << endl << flush; //  << "\twith key: " << key

        // load the building's characteristics
        Tmax = buildingXMLWrapper->Tmax;
        Tmin = buildingXMLWrapper->Tmin;
        Ninf = buildingXMLWrapper->Ninf;

        float totBuildingOccupants= buildingXMLWrapper->occupantsNumber;
        float totBuildingSre= buildingXMLWrapper->sre;
        float totBuildingVolume=0;
        for (unsigned int zoneIndex = 0; zoneIndex < buildingXMLWrapper->__sizeZones; ++zoneIndex) {
            totBuildingVolume+=buildingXMLWrapper->zones[zoneIndex].Vi;
        }
        logStream << "   SRE: "<< totBuildingSre << ", Volume: " << totBuildingVolume << ", occupants: "<< totBuildingOccupants << endl << flush;

        double dailyWhBuildingMetabolicHeatGains = buildingXMLWrapper->occupantsSensibleHeat*buildingXMLWrapper->meanDailyPresenceHours*totBuildingOccupants; // Wh / day
        double dailyWhBuildingElecCons = yearlyElecCons*1000/365.; // Wh / day
        double dailyWhBuildingInternalGains = dailyWhBuildingElecCons + dailyWhBuildingMetabolicHeatGains;

        logStream << "  Heat gains: dailyWhBuildingMetabolicHeatGains " << dailyWhBuildingMetabolicHeatGains << ", dailyWhBuildingElecCons " << dailyWhBuildingElecCons << endl << flush;

        vector<float> onekWhHeatGainProfile;
        logStream << "              onekWhHeatGainProfile [" ;
        for(int h=0; h<24; ++h)
        {
            onekWhHeatGainProfile.push_back(
                    (dailyWhBuildingElecCons*buildingXMLWrapper->OneKWhElectricityUseProfile[h] +
                     dailyWhBuildingMetabolicHeatGains*buildingXMLWrapper->OneHourPresenceProfile[h] ) /
                    dailyWhBuildingInternalGains);
            logStream << onekWhHeatGainProfile[h];
            if (h<23)
                logStream << ", ";
        }
        logStream << "]" << endl << flush;
        pDistrict->getOccupancyProfiles()->addDayProfile(id, onekWhHeatGainProfile); // dayProfile id = yearProfile id = building id
        pDistrict->getOccupancyProfiles()->addYearProfile(id, vector<unsigned int> (365,id)); // same profile every day

        // no HVAC system
        HVACpresence = false;
        // no tank means only the demands are needed
        // but a heating unit is necessary to consider an off period during summer
        heatStock = new Tank(/*V*/0.01, /*phi*/20, /*rho*/1000, /*Cp*/4180, /*Tmin*/20, /*Tmax*/60);
        coldStock = new Tank(/*V*/0.01, /*phi*/20, /*rho*/1000, /*Cp*/4180, /*Tmin*/5, /*Tmax*/20);

        // loading the heating properties
        //heatingUnit = NULL;
        heatingUnit = new Boiler(/*Pmax*/ 1e12,/*eff*/1,heatingBeginDay,heatingEndDay);

        // loading the cooling properties
        coolingUnit = NULL;

        // THERMAL MODEL
        // creates thermal zone(s)
        for (unsigned int zoneIndex = 0; zoneIndex < buildingXMLWrapper->__sizeZones; ++zoneIndex) {
            // containers of the elements
            vector<Wall*>    zoneWalls;
            vector<Roof*>    zoneRoofs;
            vector<Surface*> zoneSurfaces;
            vector<Floor*>   zoneFloors;

            // creates the OCCUPANTS for the Zone, representing both metabolic heat and electricity use gains
            // n_occupants * occupantsSensibleHeat = zoneDailyWhHeatGains ; profile for oneWh
            // (see setOccupantsSensibleHeat(1), setOccupantsLatentHeat)
            double zoneDailyWhHeatGains = dailyWhBuildingInternalGains * buildingXMLWrapper->zones[zoneIndex].zoneSre /  buildingXMLWrapper->sre;
            //Occupants* occupants =
              //new DeterministicOccupantsPresence(zoneDailyWhHeatGains,onekWhHeatGainProfile,onekWhHeatGainProfile,onekWhHeatGainProfile);

            // read the walls
            unsigned int wallIndex = 0;
            for (wallIndex = 0; wallIndex < buildingXMLWrapper->zones[zoneIndex].__sizeWalls; ++wallIndex) {

                // creation of the pointer to the WallType defined for this wall
                WallType* pWallType;
                if (to<float>(toString(buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].type))!=0) {
                    // get the type according to the value in "type"
                    pWallType = pDistrict->getWallType(toString(buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].type));
                }
                else {
                    // get or create a new WallType in the building for this special wall
                    pWallType = pDistrict->getUvalueWallType(to<float>(toString(buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].wallUValue)));

                }

                // add a new wall to the zoneWalls vector
                zoneWalls.push_back(new Wall(buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].id,
                                             pWallType,
                                             buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].ShortWaveReflectance,
                                             buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].GlazingRatio,
                                             buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].GlazingGValue,
                                             buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].GlazingUValue,
                                             buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].GlazingOpenableRatio,&logStream));

                // loop on the vertices to get what is needed
                for (unsigned int vertexIndex = 0; vertexIndex < buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].__sizeVertexs; ++vertexIndex) {
                    zoneWalls.back()->pushVertex(buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].vertexs[vertexIndex].x,
                                                 buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].vertexs[vertexIndex].y,
                                                 buildingXMLWrapper->zones[zoneIndex].walls[wallIndex].vertexs[vertexIndex].z);
                }
                zoneWalls.back()->computeNormalAndArea();

                if (zoneWalls.back()->getArea() > 0.f) {
                    // computes for this wall the eco-indicators and add them to the whole building's values
                    nre += zoneWalls.back()->getNRE();
                    gwp += zoneWalls.back()->getGWP();
                    ubp += zoneWalls.back()->getUBP();
                }
                else {
                    logStream << "(Warning) Wall id=" << zoneWalls.back()->getId() << " has a too small surface area, removing it" << endl << flush;
                    delete zoneWalls.back();
                    zoneWalls.pop_back();
                }
            }
            // end of the reading of the walls
            // read the roofs
            for (unsigned int roofIndex = 0; roofIndex < buildingXMLWrapper->zones[zoneIndex].__sizeRoofs; ++roofIndex) {
                // add the parameters to the new roof
                zoneRoofs.push_back(new Roof(buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].id,
                                             buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].Uvalue,
                                             buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].ShortWaveReflectance,
                                             buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].GlazingRatio,
                                             buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].GlazingGValue,
                                             buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].GlazingUValue,
                                             buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].GlazingOpenableRatio,&logStream));

                // loop on the vertices to get what is needed
                for (unsigned int vertexIndex = 0; vertexIndex < buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].__sizeVertexs; ++vertexIndex) {
                        zoneRoofs.back()->pushVertex(buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].vertexs[vertexIndex].x,
                                                     buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].vertexs[vertexIndex].y,
                                                     buildingXMLWrapper->zones[zoneIndex].roofs[roofIndex].vertexs[vertexIndex].z);
                }
                zoneRoofs.back()->computeNormalAndArea();

            }
            // end of the loop on roofs
            // reads the floor area and conductance
            for (unsigned int floorIndex = 0; floorIndex < buildingXMLWrapper->zones[zoneIndex].__sizeFloors; ++floorIndex) {
                // reading of the zoneFloors' id
                zoneFloors.push_back(new Floor(buildingXMLWrapper->zones[zoneIndex].floors[floorIndex].id,0.f,&logStream));
                // reads the vertices
                for (unsigned int vertexIndex = 0; vertexIndex < buildingXMLWrapper->zones[zoneIndex].floors[floorIndex].__sizeVertexs; ++vertexIndex) {
                    zoneFloors.back()->pushVertex(buildingXMLWrapper->zones[zoneIndex].floors[floorIndex].vertexs[vertexIndex].x,
                                                  buildingXMLWrapper->zones[zoneIndex].floors[floorIndex].vertexs[vertexIndex].y,
                                                  buildingXMLWrapper->zones[zoneIndex].floors[floorIndex].vertexs[vertexIndex].z);
                }
                zoneFloors.back()->computeNormalAndArea();

                // get the conductance to the ground
                zoneFloors.back()->setKground(buildingXMLWrapper->zones[zoneIndex].floors[floorIndex].Kground);
            }
            // end of floor loop

            // updates the link matrix in a sparse format
            // a variable to know which one is the first one of the row
            linksAi.push_back( linksAn.size() );
            for (unsigned int zoneSurfaceIndex = 0; zoneSurfaceIndex < buildingXMLWrapper->zones[zoneIndex].__sizeZoneSurfaces; ++zoneSurfaceIndex) {

                float Swall = buildingXMLWrapper->zones[zoneIndex].zoneSurfaces[zoneSurfaceIndex].Area;
                string zoneSurfaceType = toString(buildingXMLWrapper->zones[zoneIndex].zoneSurfaces[zoneSurfaceIndex].type);

                float Kwall;
                if ( buildingXMLWrapper->zones[zoneIndex].zoneSurfaces[zoneSurfaceIndex].Vertical ) {
                    Kwall = Swall/(1.f/pDistrict->getWallType(zoneSurfaceType)->getConductance()
                                   + 1.f/3.f + 1.f/3.f );
                }
                else { // by default a horizontal surface (buoyant 4.3 and stratified 1.5)
                    Kwall = Swall/(1.f/pDistrict->getWallType(zoneSurfaceType)->getConductance()
                                   + 1.f/1.5f + 1.f/4.3f );
                }

                float Cwall = Swall*pDistrict->getWallType(zoneSurfaceType)->getCapacitance();

                // reading the area and the connected zone
                linksAn.push_back( pair<float,float>(Kwall,Cwall) );
                linksAj.push_back( buildingXMLWrapper->zones[zoneIndex].zoneSurfaces[zoneSurfaceIndex].linkZone );

                // putting back the eco-indicators for the intermediary walls, half for each zone
                nre += Swall*pDistrict->getWallType(zoneSurfaceType)->getNRE()/2.f;
                gwp += Swall*pDistrict->getWallType(zoneSurfaceType)->getGWP()/2.f;
                ubp += Swall*pDistrict->getWallType(zoneSurfaceType)->getUBP()/2.f;

            }

            /// TODO: these things will have to change if zoning per floor appears
            // creates the associated zone (NOT only 1 Zone per Building !)
            //unsigned int id  = buildingXMLWrapper->id;
            float zoneVolume = buildingXMLWrapper->zones[zoneIndex].Vi;
            float zoneSre    = buildingXMLWrapper->zones[zoneIndex].zoneSre;
            float Kpsi       = buildingXMLWrapper->zones[zoneIndex].Psi;
            logStream << "zoneVolume: " << zoneVolume << ", zoneSre: " << zoneSre << endl << flush;
            // hypothesis, the building has only one zone and is @ ground floor
            bool groundFloor = true;

            // need to check the consistency of this zone, to define its type (1N, 2N or even 3N)
            if (zoneWalls.empty()) {
                // no external walls, create a zone with only one node
                zones.push_back(new Zone1N(zoneIndex,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,NULL,this));
            }
            else {
                // one or more external walls
                WallType* pWallType = zoneWalls[0]->getWallType();
                if (!pWallType->getLayers()->empty()) {
                    // walls not given by Uvalues, then they should have the same wallType
                    for (wallIndex=1;wallIndex<zoneWalls.size();++wallIndex) {
                        if (zoneWalls[wallIndex]->getWallType() != pWallType)
                            throw(string("In building: ") + toString(id)
                                  + ", all the walls must have the same wall type.");
                    }
                    // check complete, create the 2N Zone
                    zones.push_back(new Zone2N(zoneIndex,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,NULL,this));
                }
                else {
                    // walls are given by Uvalues, create a 1N Zone
                    zones.push_back(new Zone1N(zoneIndex,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,NULL,this));
                }
            }
            // sets the zone's thermal bridges
            zones.back()->setKpsi(Kpsi);
            zones.back()->setOccupantsSensibleHeat(1);
            zones.back()->setOccupantsLatentHeat(1);
            zones.back()->setOccupantsNumber(zoneDailyWhHeatGains);
            zones.back()->setOccupantsYearProfileID(id);

        }

        // adds the capacitance value using addZoneC to the corresponding nodes
        for (unsigned int i=0;i<getnLinksAi();++i) {
            // loop on the elements of the rows
            for (unsigned int index=getLinksAi(i);index<getLinksAi(i+1);++index) {

                // prints on the screen some outputs
                logStream << "Zone i: " << i << "\tlinked with zone id: " << getLinksAj(index) << "\tindex: " << getZoneIndexFromId(getLinksAj(index)) << "\tCwall: " << getLinksAn(index).second << endl << flush;

                // adds the capacitance to the correct zone and node
                if (zones[i]->getnNodes() == 2 ) {
                  addZoneC(i, 0, getLinksAn(index).second/4.);
                  addZoneC(i, 1, getLinksAn(index).second/4.);
                }
                else if (zones[i]->getnNodes() == 3) {
                  addZoneC(i, 2, getLinksAn(index).second/2.);
                }
                else throw ("Links not defined for more than 3 nodes.");

            }
        }
        // end of the creation of the zone and the links between the zones
        logStream << "Zones loaded: " << zones.size() << endl << flush;

        // output of the eco-indicators
        logStream << "NRE: " << nre << " MJ\tGWP: " << gwp << " kgCO2\tUBP: " << ubp << " pts" << endl << flush;

        // creates the three matrices matrices that describe the problem thermally
        int NP = getnNodes();
        // matrix of the capacitances C and constant conductances G1
        C = new Mat_DP(NP,NP);
        G1 = new Mat_DP(NP,NP);
        initialiseMatrix(C);
        initialiseMatrix(G1);
        // the initialisation of the matrices C and G is done with the fixed values of the conductances/capacitances
        for (unsigned int i=0; i<zones.size(); ++i) {
            for (unsigned int j=0; j<zones[i]->getnNodes(); ++j) {
                // put the values in the C matrix
                (*C)[getMatrixPosition(i)+j][getMatrixPosition(i)+j] = zones[i]->getC(j);
                for (unsigned int k=0; k<zones[i]->getnNodes(); ++k) {
                    // put the values in the G1 matrix
                    (*G1)[getMatrixPosition(i)+j][getMatrixPosition(i)+k] = zones[i]->getFixedMatrixElement(j,k);
                }
            }
        }
        // loads the matrix of the conductances of the links between the zones (addition to the present fixed conductances matrix)
        loadLinkSparse(G1);

        // end of the building constructor
        return;
    }

    /* District class builder, using the XML wrapper defined in CitySimWebS.h */

    District::District(const ns__DistrictW* districtXMLWrapper, XmlScene* pScene, const std::map<unsigned int, double>* electricityR, unsigned int heatingBeginDay, unsigned int heatingEndDay):pScene(pScene),logStream(pScene->logStream.rdbuf()) {

        // Far Field Obstructions loading in the District
        logStream << "Loading Far Field Obstruction profile: ";
        if (districtXMLWrapper->__sizeFarFieldObstructionPoints > 0) {
            for (unsigned int index=0; index < districtXMLWrapper->__sizeFarFieldObstructionPoints; ++index) {
                farFieldObstructions.push_back(pair<float,float>( districtXMLWrapper->farFieldObstructionPoints[index].phi,
                                                                  districtXMLWrapper->farFieldObstructionPoints[index].theta ));
            }
            logStream << farFieldObstructions.size() << " loaded." << endl << flush;
        }
        else logStream << "no profile given." << endl << flush;

        // Wall Types loading in the District
        logStream << "Loading wall type(s): ";
        if (districtXMLWrapper->__sizeWallTypes > 0) {
            for (unsigned int index=0; index < districtXMLWrapper->__sizeWallTypes; ++index) {
                wallTypes.insert( pair<string,WallType*>(toString(districtXMLWrapper->wallTypes[index].id), new WallType(&districtXMLWrapper->wallTypes[index],&logStream)));
                logStream << "Added walltype with id " << districtXMLWrapper->wallTypes[index].id << endl << flush;
            }
            logStream << wallTypes.size() << " loaded." << endl << flush;
        }
        else logStream << "Warning: no WallType defined." << endl << flush;

        logStream << "Loading geometry of building(s)." << endl << flush;
        // loop on the buildings, if at least one building exists
        if (districtXMLWrapper->__sizeBuildings > 0) {
            unsigned int b_id;
            for (unsigned int index=0; index < districtXMLWrapper->__sizeBuildings; ++index) {
                b_id = districtXMLWrapper->buildings[index].id;
                // check if the building should be simulated or not (in the building tag itself)
                if (!districtXMLWrapper->buildings[index].Simulate) {
                    // if simulate=false, then easy, we just add surfaces to the ground for the SW calculation only
                    logStream << "Building: " << b_id;
                    logStream << "\tnot thermally simulated." << endl << flush;
                    // loop on thermal zones withing the building itself
                    for (unsigned int zoneIndex=0; zoneIndex < districtXMLWrapper->buildings[index].__sizeZones; ++zoneIndex) {
                        // loading the wall vertices and put them in a vector
                        for (unsigned int wallIndex = 0; wallIndex < districtXMLWrapper->buildings[index].zones[zoneIndex].__sizeWalls; ++wallIndex) {
                            // puts the points in the groundSurfaces for no thermal simulation
                            groundSurfaces.push_back(new Surface(districtXMLWrapper->buildings[index].zones[zoneIndex].walls[wallIndex].id,
                                                                 districtXMLWrapper->buildings[index].zones[zoneIndex].walls[wallIndex].ShortWaveReflectance,
                                                                 0.f,0.f,0.f,0.f,&logStream));

                            for (unsigned int pointIndex=0; pointIndex < districtXMLWrapper->buildings[index].zones[zoneIndex].walls[wallIndex].__sizeVertexs; ++pointIndex) {
                                // put back the vertex of each point
                                groundSurfaces.back()->pushVertex(districtXMLWrapper->buildings[index].zones[zoneIndex].walls[wallIndex].vertexs[pointIndex].x,
                                                                  districtXMLWrapper->buildings[index].zones[zoneIndex].walls[wallIndex].vertexs[pointIndex].y,
                                                                  districtXMLWrapper->buildings[index].zones[zoneIndex].walls[wallIndex].vertexs[pointIndex].z);
                            }
                            groundSurfaces.back()->computeNormalAndArea();
                        }
                        // loading the roof vertices
                        for (unsigned int roofIndex = 0; roofIndex < districtXMLWrapper->buildings[index].zones[zoneIndex].__sizeRoofs; ++roofIndex) {
                            // puts the points in the groundSurfaces for no thermal simulation
                            groundSurfaces.push_back(new Surface(districtXMLWrapper->buildings[index].zones[zoneIndex].roofs[roofIndex].id,
                                                                 districtXMLWrapper->buildings[index].zones[zoneIndex].roofs[roofIndex].ShortWaveReflectance,
                                                                 0.f,0.f,0.f,0.f,&logStream));

                            for (unsigned int pointIndex=0; pointIndex < districtXMLWrapper->buildings[index].zones[zoneIndex].roofs[roofIndex].__sizeVertexs; ++pointIndex) {
                                // put back the vertex of each point
                                groundSurfaces.back()->pushVertex(districtXMLWrapper->buildings[index].zones[zoneIndex].roofs[roofIndex].vertexs[pointIndex].x,
                                                                  districtXMLWrapper->buildings[index].zones[zoneIndex].roofs[roofIndex].vertexs[pointIndex].y,
                                                                  districtXMLWrapper->buildings[index].zones[zoneIndex].roofs[roofIndex].vertexs[pointIndex].z);
                            }
                            groundSurfaces.back()->computeNormalAndArea();
                        }
                    }
                    // end of the loop on all zones belonging to a building
                }
                else {
                    // the building should be simulated, we then add it to the vector of buildings
                    if(electricityR!=NULL) {// electricity heat gains are provided
                        logStream << "Electricity demand of building: " << electricityR->find(b_id)->second << endl << flush;
                        buildings.push_back(new Building(&districtXMLWrapper->buildings[index],this,electricityR->find(b_id)->second,heatingBeginDay,heatingEndDay));
                    }
                    else {
                        logStream << "No electricity demand values" << endl << flush;
                        buildings.push_back(new Building(&districtXMLWrapper->buildings[index],this,0));
                    }
                }
            }
            // enf of the loop on the buildings
            logStream << buildings.size() << " loaded." << endl << flush;
        }
        else throw(string("Error in XML file: no Buildings."));

        logStream << "Loading ground surface(s): ";
        for (unsigned int index=0; index < districtXMLWrapper->__sizeGroundSurfaces; ++index) {

            // reads the vertices and creates a ground, with an id
            groundSurfaces.push_back(new Ground(districtXMLWrapper->groundSurfaces[index].id,
                                                districtXMLWrapper->groundSurfaces[index].ShortWaveReflectance,&logStream));

            for (unsigned int pointIndex=0; pointIndex < districtXMLWrapper->groundSurfaces[index].__sizeVertexs; ++pointIndex) {
                groundSurfaces.back()->pushVertex(districtXMLWrapper->groundSurfaces[index].vertexs[pointIndex].x,
                                                  districtXMLWrapper->groundSurfaces[index].vertexs[pointIndex].y,
                                                  districtXMLWrapper->groundSurfaces[index].vertexs[pointIndex].z);
            }
            groundSurfaces.back()->computeNormalAndArea();
        }
        logStream << groundSurfaces.size() << " loaded." << endl << flush;

        logStream << "District created." << endl << flush;

    }

    /* WallType class builder, using the XML wrapper defined in CitySimWebS.h */

    WallType::WallType(ns__WallTypeW* wallTypeXMLWrapper, ostream* pLogStr):wallTypeId(wallTypeXMLWrapper->id),logStream(std::cout.rdbuf()) {

        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        for (unsigned int index=0; index < wallTypeXMLWrapper->__sizeLayers; ++index) {
            float nre = 0.f;
            float gwp = 0.f;
            float ubp = 0.f;
            vLayer.push_back(Layer(wallTypeXMLWrapper->layers[index].Thickness,
                                   wallTypeXMLWrapper->layers[index].Conductivity,
                                   wallTypeXMLWrapper->layers[index].Cp,
                                   wallTypeXMLWrapper->layers[index].Density,
                                   nre, gwp, ubp));
        }
        // computes the standard Uvalue according to inside and outside layers of conductance 8 and 25
        Uvalue = 1.f/(0.125f + getResistance() + 0.04f);
        // screen output
        logStream << "Adding a layered walltype with id " << wallTypeXMLWrapper->id << endl << flush;
    }

    /* Climate class builder, using the XML wrapper defined in CitySimWebS.h */

    Climate::Climate(const ns__ClimateW* climateXMLWrapper, ostream* pLogFileStream):logStream(std::cout.rdbuf()) {

        // logStream is directed by default to the "cout" streambuf
        if(pLogFileStream!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogFileStream->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        location = climateXMLWrapper->location;
        latitudeN = climateXMLWrapper->latitudeN;
        longitudeE = climateXMLWrapper->longitudeE;
        altitude = climateXMLWrapper->altitude;
        meridianE = climateXMLWrapper->meridianE;

        for (unsigned int index=0; index<climateXMLWrapper->__sizeIbn; ++index)
            Ibn.push_back(climateXMLWrapper->ibn[index]);

        for (unsigned int index=0; index<climateXMLWrapper->__sizeIdh; ++index)
            Idh.push_back(climateXMLWrapper->idh[index]);

        // store the Igh if necessary
        for (unsigned int index=0; index<climateXMLWrapper->__sizeIgh; ++index)
            Igh.push_back(climateXMLWrapper->igh[index]);

        for (unsigned int index=0; index<climateXMLWrapper->__sizeTout; ++index)
            Tout.push_back(climateXMLWrapper->tout[index]);

        // Define t_ground temperature at 8°C (SIA norm) to avoid bad model
        for (unsigned int index=0; index<climateXMLWrapper->__sizeTout; ++index)
            Tground.push_back(8);

        for (unsigned int index=0; index<climateXMLWrapper->__sizeWindSpeed; ++index)
            windSpeed.push_back(climateXMLWrapper->windSpeed[index]);

        for (unsigned int index=0; index<climateXMLWrapper->__sizeWindDirection; ++index)
            windDirection.push_back(climateXMLWrapper->windDirection[index]);

        for (unsigned int index=0; index<climateXMLWrapper->__sizeRelativeHumidity; ++index)
            relativeHumidity.push_back(climateXMLWrapper->relativeHumidity[index]);

        for (unsigned int index=0; index<climateXMLWrapper->__sizePrec; ++index)
            Prec.push_back(climateXMLWrapper->prec[index]);

        for (unsigned int index=0; index<climateXMLWrapper->__sizeCloudiness; ++index)
            cloudiness.push_back(climateXMLWrapper->cloudiness[index]);

        logStream << "Climate constructor: read all data" << endl << flush;

        // calculation of several variables that we use in the radiation model
        meanAnnualTemperature = accumulate(Tout.begin(),Tout.end(),0.f)/float(Tout.size());
        //calculate the temperature for a day and set the coolest day in the year and a vector of daily mean temperature
        meanDailyTemperature.assign(365,0.f);
        for(unsigned int i=0;i<365;++i) {
            meanDailyTemperature[i] = accumulate(Tout.begin()+24*i,Tout.begin()+24*(i+1),0.f)/24.f;
            // note: accumulate does not take into account the last element pointed by end: [first,last)
            //logStream << *(Tout.begin()+24*i) << "\t" << *(Tout.begin()+24*(i+1)) << endl << flush;
        }
        hotDay = distance(meanDailyTemperature.begin(), max_element(meanDailyTemperature.begin(), meanDailyTemperature.end()));
        coolDay = distance(meanDailyTemperature.begin(), min_element(meanDailyTemperature.begin(), meanDailyTemperature.end()));

        // determine dew point temperature
        Td = new float[Tout.size()];
        for(unsigned int i=0;i<Tout.size();++i){
    //        double a = 17.27;
    //        double b = 237.7;
    //        double gamma = a * Tout[i] /(b+Tout[i])+ log(Hr[i]/100.0);
    //        Td[i] = b*gamma /(a-gamma);
    //        if(Td[i] < 0 || Tout[i] < 0)
    //            Td[i] = Tout[i];
            Td[i] = 1.f / ( 1.f/(Tout[i] + 273.15f) - 1.85e-4f * log(relativeHumidity[i]/100.f) ) - 273.15f;
        }

        // If only global irradiance available (nothing for Ibn and Igh), use model in dirint.c to compute Idh and Ibn
        if (Ibn.empty() && Idh.empty()) {

            #ifdef DEBUG
            ofstream file("maxwellIbnIdh.txt");
            file << "#Ibn\tIdh" << endl << flush;
            #endif

            // creates a sun for the current location
            SKYSiteLocation location(GENAngle::Degrees(latitudeN), GENAngle::Degrees(longitudeE), GENAngle::Degrees(meridianE*360.f/24.f), GENAngle::Degrees(0.f));
            SKYSun sunClimate(location);
            int dayOfYear;

            float Ib;
            float z[3], Ig[3]; // previous, current and next data
            z[0]=NA;
            Ig[0]=NA;
            z[2] = (90.f - sunClimate.GetPosition().Altitude().degrees())/(180./M_PI); // angle entre le zénith et le soleil, initialiser avec la première valeur

            Ig[2] = Igh[0];

            for (size_t j=0; j<Igh.size()-1; ++j) {

                Ig[1] = Ig[2];
                Ig[2] = Igh[j+1];
                z[1] = z[2];

                dayOfYear = (j+1)/24+1;
                sunClimate.SetDay(dayOfYear);
                sunClimate.SetClockTime(j%24);
                z[2] = (90.f - sunClimate.GetPosition().Altitude().degrees())/(180./M_PI);	// angle entre le zénith et le soleil (RADIANS), valeur suivante

                /*
                 ** float dirint(float *g, float *z, float *td, int *doy, float *alt)
                 ** 'g' is the address of someplace that has 3 floats representing
                 ** 	global irradiance in watts. The 3 refer to the previous, current
                 **	and next observed data
                 ** 'z' is the address of someplace with 3 floats representing
                 **	the zenith angle in RADIANS of the three observations.
                 ** 'td' points to the dewpoint temperature in degrees C.
                 ** 'doy' points to the day of year (from 1).
                 ** 'alt' points to the altitude of the site, in meters.
                 **
                 ** The returned value is the modelled direct normal irradiance in watts.
                */

                Ib = dirint(&Ig[0], &z[0], &Td[j], &dayOfYear, &altitude);
                // float io = 1367.f; // beam normal extraterrestrial irradiance
                // Ib = dirint_(&Ig[0], &z[0], &Td[j], &altitude, &i0);
                Ibn.push_back(Ib);
                Idh.push_back((Ig[1] - Ib*cos(z[1])));

                #ifdef DEBUG
                file << Ibn.back() << "\t" << Idh.back() << endl;
                #endif

                /*logStream << "G_h : " << setprecision(3) << setw(5) << Ig[1]
                    << ", G_bn : " << setprecision(3) << setw(5) << Ib
                    << ", G_dh : " << setprecision(3) << setw(5) << Ig[1] - Ib*cos(z[1])
                    << ", dy : " << setprecision(3) << setw(3) << dayOfYear
                    << ", h : " << setprecision(3) << setw(3) << (j+1)%24
                    << ", dT : " << setprecision(3) << setw(6) <<  Td[j]
                    << ", z : " << setprecision(3) << setw(5) << z[0] << " " << setprecision(3) << setw(5) << z[1] << " " << setprecision(3) << setw(5) << z[2]
                    << ", Ig : " << setprecision(3) << setw(5) << Ig[0] << " " << setprecision(3) << setw(5) << Ig[1] << " " << setprecision(3) << setw(5) << Ig[2] << endl << flush;
    */
                Ig[0] = Ig[1];
                z[0] = z[1];
            }

            Ig[1] = Ig[2];
            Ig[2]=NA;
            z[1] = z[2];
            z[2]=NA;

            Ib = dirint(&Ig[0], &z[0], &Td[8759], &dayOfYear, &altitude);
            // float io = 1367.f; // beam normal extraterrestrial irradiance
            // Ib = dirint_(&Ig[0], &z[0], &Td[8759], &altitude, &i0);
            Ibn.push_back(Ib);
            Idh.push_back(Ig[1] - Ib*cos(z[1]));

            // clear the useless Igh vector
            Igh.clear();
        }
    }

    /* XmlScene class builder, using the XML wrappers defined in CitySimWebS.h */

    XmlScene::XmlScene(string simulationName, const ns__ClimateW* climateXMLWrapper, const ns__DistrictW* districtXMLWrapper, const std::map<unsigned int, double>* electricityR, ostream* pLogFileStr){
		// , unsigned int beginDay, unsigned int beginMonth, unsigned int endDay, unsigned int endMonth
        // logStream is directed by default to the "cout" streambuf
        // If a log file stream is provided, redirect logStream to the file stream.
        if(pLogFileStr!=NULL)
            logStream.rdbuf(pLogFileStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        // initialise the simulated number of time steps
        preTimeStepsSimulated = 0;
        preTimeSteps2Simulated = 0;
        timeStepsSimulated = 0;
        timeSteps2Simulated = 0;
        simulationIndex = 0;

        // copy the inputFile
        this->inputFile = simulationName;

        // reads the Radiance file description and put the facades in the surfaceVector
        logStream << "Reading XML wrapper..." << endl << flush;

        logStream << "Checking ns__ClimateW object in XmlScene constructor: "<< climateXMLWrapper->location <<" "<< climateXMLWrapper->latitudeN <<" "<< climateXMLWrapper->longitudeE<<" "<<climateXMLWrapper->altitude<<endl<<flush;
        // load the climate
        pClimate = new Climate(climateXMLWrapper,&logStream);

        // shows the information about the location
        logStream << "Location: " << pClimate->getLocation();
        logStream << "\t(Latitude: " << pClimate->getLatitudeN() << " N, Longitude: " << pClimate->getLongitudeE();
        logStream << " E, Altitude: " << pClimate->getAltitude() << " m, Meridian: " << pClimate->getMeridian() << " h)" << endl << flush;

        // gets the information from the climate file
        float latN = pClimate->getLatitudeN();
        float longE = pClimate->getLongitudeE();
        float merE = pClimate->getMeridian()*360./24.;
        float northW = 0.f;

        // initilisation of the location for the sun and the scene
        SKYSiteLocation location(GENAngle::Degrees(latN), GENAngle::Degrees(longE), GENAngle::Degrees(merE), GENAngle::Degrees(northW));
        pSun = new SKYSun(location);
        scene.SetLocation(location);

		// Year simulation
        this->beginDay = 1;
        this->endDay   = 365;

		// Set a heating period excluding a period in summer
		// The excluded period is centered around the 1st of August and lasts the number of day with a mean T > 13°C
		// (beginDay > endDay -> off period during the summer)
		vector<float> meanDailyTemp = pClimate->getMeanDailyTemperature();
		int offPeriodLength = std::count_if(meanDailyTemp.begin(), meanDailyTemp.end(), bind2nd(greater<float>(), 13));

		unsigned int heatingBeginDay = 212 + offPeriodLength/2;
        unsigned int heatingEndDay   = 212 - offPeriodLength/2 - offPeriodLength%2;

        // sets the seed for the random number generator (must be a uint32_t)
        zigset( static_cast<uint32_t>(26041978) );
        logStream << "sprng seed: " << static_cast<uint32_t>(26041978) << endl << flush;

        // shows the simulation time span
        logStream << "begin day: " << beginDay << "\tend day: " << endDay << endl << flush;

        // creates the district
        pDistrict = new District(districtXMLWrapper, this, electricityR, heatingBeginDay, heatingEndDay);

        // browse the district to create buildings surfaces
        for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) { // loop on all buildings
            for (unsigned int j=0; j<pDistrict->getBuilding(i)->getnZones(); ++j) { // loop on all zones in the building
                // loop for the walls on this zone
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnWalls(); ++k) {
                    logStream << "Building " << i << "\tZone: " << j << "\tWall " << k << endl << flush;
                    // add the surface to the Buildings surfaces (daylight calculation)
                    scene.AddBuildingSurface(GENHandle<Wall>(pDistrict->getBuilding(i)->getZone(j)->getWall(k)));
                }
                // loop for the roofs on this zone
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnRoofs(); ++k) {
                    logStream << "Building " << i << "\tZone: " << j << "\tRoof " << k << endl << flush;
                    // add the surface to the Ground surfaces (meaning NO daylight calculation)
                    scene.AddGroundSurface(GENHandle<Roof>(pDistrict->getBuilding(i)->getZone(j)->getRoof(k)));
                }
                // loop on the obstructing surfaces on this zone
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnSurfaces(); ++k) {
                    logStream << "Building " << i << "\tZone: " << j << "\tSurface " << k << endl << flush;
                    // add the surface to the Ground surfaces (meaning NO daylight calculation)
                    scene.AddGroundSurface(GENHandle<Surface>(pDistrict->getBuilding(i)->getZone(j)->getSurface(k)));
                }
            }
        }
        logStream << "Buildings' surfaces added to the scene." << endl << flush;

        // browse the district to create ground surfaces
        for (unsigned int i=0; i<pDistrict->getnGroundSurfaces(); i++) {
            // adds the ground surfaces to the scene
            scene.AddGroundSurface(GENHandle<Surface>(pDistrict->getGroundSurface(i)));
        }
        logStream << "Ground surfaces added to the scene." << endl << flush;

        // calculates the view factors of the scene
        // N.B.: the view factor calculation starts the direct, diffuse and daylight calculations in sequence
        // the direct calculation needs the Site Location in order to compute all sun positions
        v.CalculateViewFactors(scene);
        logStream << "View factors calculated." << endl << flush;

        // initialise some of the important parameters
        //
        // define the number of surfaces to play with
        m_NbSurf = scene.SurfaceCount();
        mNbReflections = 2;
        logStream << "Building inter-reflection matrix." << endl << flush;
        logStream << "Number of surfaces: " << m_NbSurf << endl << flush;
        logStream << "Number of reflections: " << mNbReflections << endl << flush;
        buildSparseMatrix();
        // initialises the Far Field obstructions vector
        initialiseFarField();

    }
