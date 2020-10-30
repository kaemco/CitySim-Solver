#ifndef CITYSIMWEBS_H
#define CITYSIMWEBS_H
//File: CitySimWebS.h
//gsoap ns service name: CitySimWebS
//gsoap ns service namespace: urn:CitySimWebS
//gsoap ns service location: http://lesopbpc27.epfl.ch:9201

/*
Autor : Renaud Sauvain
Date : 20.03.2011
Content : header & presentation for CitySimWebS gsoap web service
*/

	class ns__FarFieldObstructionPointW {
	public:
		float phi, theta;
	};

	class ns__LayerW {
	public:
		std::string name;
		float Thickness, Conductivity, Cp, Density;
	};

	class ns__WallTypeW {
	public:
		unsigned int id;
		std::string name;

		unsigned int __sizeLayers; // number of elements pointed to
		ns__LayerW * layers;// points to array elements
        float wallUValue;

        // Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        //ns__WallTypeW();
        //ns__WallTypeW(ns__WallTypeW const&);
        //ns__WallTypeW& operator=(ns__WallTypeW const&);
		//~ns__WallTypeW();
    //private:
        //void copy(ns__WallTypeW const&);
	};

	class ns__VertexW {
	public:
		double x,y,z;
	};

	class ns__SurfaceW {
    public:
        unsigned int __sizeVertexs; // number of elements pointed to
        ns__VertexW * vertexs;// points to array elements

        unsigned int id;

        // Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        //ns__SurfaceW();
        //ns__SurfaceW(ns__SurfaceW const&);
        //ns__SurfaceW& operator=(ns__SurfaceW const&);
		//~ns__SurfaceW();
    //private:
        //void copy(ns__SurfaceW const&);
	};

	class ns__WallW : public ns__SurfaceW{
	public:
		 unsigned int type;
		 // JK - each wall has a GlazingRatio & a GlazingOpenableRatio
		 float ShortWaveReflectance, GlazingRatio, GlazingGValue, GlazingUValue, GlazingOpenableRatio;
		 float wallUValue;
	};

	class ns__RoofW : public ns__SurfaceW{
	public:
		 float Uvalue, ShortWaveReflectance, GlazingRatio, GlazingGValue, GlazingUValue, GlazingOpenableRatio;
	};

	class ns__FloorW : public ns__SurfaceW{
	public:
		 float Kground;
	};

	class ns__ZoneSurfaceW{
	public:
		unsigned int id, type, linkZone;
		float Area;
		bool Vertical;
	};

	class ns__GroundSurfaceW : public ns__SurfaceW{
	public:
		float ShortWaveReflectance;
	};

	class ns__ZoneW{
	public:
        // the volume of the zone inside the building and its parameters
		float Psi, Vi, zoneSre;

		unsigned int __sizeWalls; // number of elements pointed to
		ns__WallW * walls;// points to array elements

		// JK - modifié pour avoir notation vectorielle de roofs et floors
		unsigned int __sizeRoofs; // number of elements pointed to
		ns__RoofW * roofs;// points to array elements

		unsigned int __sizeFloors; // number of elements pointed to
		ns__FloorW * floors;// points to array elements

		unsigned int __sizeZoneSurfaces; // number of elements pointed to
		ns__ZoneSurfaceW * zoneSurfaces;// points to array elements

		// Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        //ns__ZoneW();
        //ns__ZoneW(ns__ZoneW const&);
        //ns__ZoneW& operator=(ns__ZoneW const&);
		//~ns__ZoneW();
    //private:
        //void copy(ns__ZoneW const&);
	};

	class ns__BuildingW{
	public:
		unsigned int id;
		float Ninf, Tmin, Tmax;
		float blindsIrrCutOff;
		bool Simulate;
		bool ignoreResults;

		// the vector containing the building Zones, it should be created by the constructor of building
		unsigned int __sizeZones; // number of elements pointed to
		ns__ZoneW * zones;// points to array elements
		unsigned short int mainAllocationId; // id of allocation type
		float sre; // surface of the building

        float occupantsNumber;
        float meanDailyPresenceHours; // daily presence in h
        float occupantsSensibleHeat; // power in W
        //float annualElectricityGainsKWh=0; // total yearly kWh

		unsigned int __sizeOneKWhElectricityUseProfile;
		float * OneKWhElectricityUseProfile;
		unsigned int __sizeOneHourPresenceProfile;
		float * OneHourPresenceProfile;

		// Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        //ns__BuildingW();
        //ns__BuildingW(ns__BuildingW const&);
        //ns__BuildingW& operator=(ns__BuildingW const&);
		//~ns__BuildingW();
    //private:
        //void copy(ns__BuildingW const&);
	};

	class ns__DistrictW{
	public:
		// elements in the District
		unsigned int __sizeBuildings; // number of elements pointed to
		ns__BuildingW * buildings;// points to array elements

		unsigned int __sizeGroundSurfaces; // number of elements pointed to
		ns__GroundSurfaceW * groundSurfaces;// points to array elements

		unsigned int __sizeFarFieldObstructionPoints; // number of elements pointed to
		ns__FarFieldObstructionPointW * farFieldObstructionPoints;// points to array elements

		unsigned int __sizeWallTypes; // number of elements pointed to
		ns__WallTypeW * wallTypes;// points to array elements

        // Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        //ns__DistrictW();
        //ns__DistrictW(ns__DistrictW const&);
        //ns__DistrictW& operator=(ns__DistrictW const&);
		//~ns__DistrictW();
    //private:
        //void copy(ns__DistrictW const&);
	};

	class ns__NormPairW{
	public:
		unsigned int allocId;
		double value;
	};

	class ns__NormsW{
	public:
		// elements in a norm map a main allocation id to a norm value
		unsigned int __sizeHotWater;
		ns__NormPairW * hotwater;

		unsigned int __sizeElectricity;
		ns__NormPairW * electricity;

        // Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        ns__NormsW();
        ns__NormsW(ns__NormsW const&);
        ns__NormsW& operator=(ns__NormsW const&);
		~ns__NormsW();

    private:
        void copy(ns__NormsW const&);
	};

	class ns__ClimateW
	{
    public:
        std::string location;
        float altitude, latitudeN, longitudeE;
        int meridianE;

        unsigned int __sizeTout;
        float * tout;

        unsigned int __sizeWindSpeed;
        float * windSpeed;

        unsigned int __sizeWindDirection;
        float * windDirection;

        unsigned int __sizeRelativeHumidity;
        float * relativeHumidity;

        unsigned int __sizeIbn;
        float * ibn;

        unsigned int __sizeIdh;
        float * idh;

        unsigned int __sizeIgh;
        float * igh;

        unsigned int __sizePrec;
        float * prec;

        unsigned int __sizeCloudiness;
        float * cloudiness;

        // Soap does not like its received objects being deleted, so no destructor will be defined and
        // to avoid memory leaks, copies are always superficial
        // (In the best case scenario, the ns__ object should simply not be copied, which should be the case now)
        //ns__ClimateW();
        //ns__ClimateW(ns__ClimateW const&);
        //ns__ClimateW& operator=(ns__ClimateW const&);
		//~ns__ClimateW();
    //private:
        //void copy(ns__ClimateW const&);
	};

	/* TimedConsumptionResults */ //limité à des abréviations pour diminuer la charge lors de demande horaire
	class ns__TCRW
	{
	public:
		/* heating */
		float h;
		/* cooling */
		float c;
		/* hotwater */
		//double hw;
		/* electricity */
		//double e;
	};

	class ns__BuildingResultW
	{
	public:
		unsigned int building_ID;
		unsigned int __size; // number of elements pointed to
		ns__TCRW * WhConsumptions;// points to array elements

		float yearlykWhConsumptionE, yearlykWhConsumptionHW;

		ns__BuildingResultW();
        ns__BuildingResultW(ns__BuildingResultW const&);
        ns__BuildingResultW& operator=(ns__BuildingResultW const&);
		~ns__BuildingResultW();

    private:
        void copy(ns__BuildingResultW const&);
	};

	class ns__ResultsW
	{
	public:
		std::string comments;
		unsigned int __size; // number of elements pointed to
		ns__BuildingResultW * __ptr;// points to array elements

		ns__ResultsW();
		ns__ResultsW(ns__ResultsW const&);
        ns__ResultsW& operator=(ns__ResultsW const&);
		~ns__ResultsW();

    private:
        void copy(ns__ResultsW const&);
	};

/*	int ns__HourlyXMLSimulation(std::string buildingXML, unsigned int userID, ns__ResultsW &return_);
	int ns__MonthlyXMLSimulation(std::string buildingXML, unsigned int userID, ns__ResultsW &return_);
	int ns__YearlyXMLSimulation(std::string buildingXML, unsigned int userID, ns__ResultsW &return_);
*/
	int ns__HourlySimulation(ns__DistrictW* district, ns__ClimateW* climate, ns__NormsW* norms,  unsigned int userID, ns__ResultsW &return_);
	int ns__MonthlySimulation(ns__DistrictW* district, ns__ClimateW* climate, ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_);
	int ns__YearlySimulation(ns__DistrictW* district, ns__ClimateW* climate, ns__NormsW* norms, unsigned int userID, ns__ResultsW &return_);

#endif
