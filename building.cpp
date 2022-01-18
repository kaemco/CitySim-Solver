#include "district.h"
#include "building.h"
#include "zone.h"
#include "scene.h"
#include "util.h"

// *** Building class, CitySim *** //
// *** jerome.kaempf@epfl.ch   *** //

XmlScene* Building::getScene(){ return getDistrict()->getScene();}

Building::Building(TiXmlHandle hdl, District* pDistrict):pDistrict(pDistrict),logStream(pDistrict->logStream.rdbuf()) {

    // default XML pointing element
    TiXmlElement* xmlElem = NULL;

    // gets the buildings ID (used?)
    id = to<unsigned int>(hdl.ToElement()->Attribute("id"));
    if (hdl.ToElement()->Attribute("key")) key = hdl.ToElement()->Attribute("key");
    logStream << "Building: " << id << "\twith key: " << key << endl << flush;

    // gets the fmu file to simulate this building
    if (hdl.ToElement()->Attribute("fmu")) fmuFile = hdl.ToElement()->Attribute("fmu");
    if (hdl.ToElement()->Attribute("tmp")) tmpPath = hdl.ToElement()->Attribute("tmp");

    // load the building's characteristics
    //float Vi = 0.f; // temporary building's volume for one zone creation
    //if (hdl.ToElement()->Attribute("Vi")) Vi=to<float>(hdl.ToElement()->Attribute("Vi"));
    if (hdl.ToElement()->Attribute("BlindsLambda")) blindsLambda = to<float>(hdl.ToElement()->Attribute("BlindsLambda"));
    if (hdl.ToElement()->Attribute("BlindsIrradianceCutOff")) blindsIrradianceCutOff = to<float>(hdl.ToElement()->Attribute("BlindsIrradianceCutOff"));
    if (hdl.ToElement()->Attribute("Simulate")) simulateEP = (hdl.ToElement()->Attribute("Simulate") == string("ep"));

    // load the MRT specifications
    if (hdl.ToElement()->Attribute("mrt")) mrt = (hdl.ToElement()->Attribute("mrt") == string("true"));
    if (hdl.ToElement()->Attribute("mrtEpsilon")) mrtEpsilon = to<float>(hdl.ToElement()->Attribute("mrtEpsilon"));

    // load the parameters for the HVAC system
    if (hdl.ChildElement("HVAC", 0).ToElement()) {
        HVACpresence = true;
        coileff    = to<float>(hdl.ChildElement("HVAC", 0).ToElement()->Attribute("coilEfficiency"));
        TminSupply = to<float>(hdl.ChildElement("HVAC", 0).ToElement()->Attribute("TminSupply"));
        TmaxSupply = to<float>(hdl.ChildElement("HVAC", 0).ToElement()->Attribute("TmaxSupply"));
        deltaT     = to<float>(hdl.ChildElement("HVAC", 0).ToElement()->Attribute("deltaT"));
        evaporativeCooling =  (hdl.ChildElement("HVAC", 0).ToElement()->Attribute("evaporativeCooling") == string("true"));
        coilHTWeff = to<float>(hdl.ChildElement("HVAC", 0).ToElement()->Attribute("hydroThermalWheel"));
    }
    else HVACpresence = false;

    // loading and creating the tanks
    if(hdl.ChildElement("HeatTank", 0).ToElement())
        heatStock = new Tank(hdl.ChildElement("HeatTank", 0));
    else throw(string("Error in the XML file: no HeatTank tag."));
    if(hdl.ChildElement("DHWTank", 0).ToElement())
        dhwStock = new Tank(hdl.ChildElement("DHWTank", 0));
    if(hdl.ChildElement("CoolTank", 0).ToElement())
        coldStock = new Tank(hdl.ChildElement("CoolTank", 0));
    else throw(string("Error in the XML file: no CoolTank tag."));

    // Cognet: Start of added code.
    heatingNeeds = 1.0; coolingNeeds = 1.0;
    HS_needs = 1.0; DHW_needs = 1.0; CS_needs = 1.0;
    HS_SolPp = 1.0; DHW_SolPp = 1.0; CS_SolPp = 1.0;
    SolTherFracLeft = 1.0;
    HS_Pp = 1.0; DHW_Pp = 1.0; CS_Pp = 1.0;
    VdotUsed = 0.01;
    Tamb = 10.0;
    // Cognet: End of added code.


    // loading the heating properties
    if(hdl.FirstChildElement("HeatSource").ToElement()) {
        TiXmlHandle heatSource = hdl.FirstChildElement("HeatSource");
        // gets the working period of the heatSource
        unsigned int beginDay, endDay;
        if ( heatSource.ToElement()->Attribute("beginDay") == NULL ) beginDay = 1;
        else beginDay = to<unsigned int>(heatSource.ToElement()->Attribute("beginDay"));
        if ( heatSource.ToElement()->Attribute("endDay") == NULL )   endDay   = 365;
        else endDay   = to<unsigned int>(heatSource.ToElement()->Attribute("endDay"));

        bool heatSourceFound = false; // Cognet: Added this, to detect error if multiple heatsources present in xml.
        string errorMsg = "Error in the XML file: a Building has multiple units in HeatSource, there must be only one."; // Cognet: Added this

        // Load source of heat
        if (heatSource.FirstChildElement("Boiler").ToElement()) {
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true; // Cognet: Added this.
            heatingUnit = new Boiler(heatSource.FirstChildElement("Boiler"),beginDay,endDay,&(this->logStream));
        }
        if (heatSource.FirstChildElement("HeatPump").ToElement()) { // Cognet: Added this, to detect possible errors in xml.
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true; // Cognet: Added this.
            heatingUnit = new HeatPump(heatSource.FirstChildElement("HeatPump"),beginDay,endDay,&(this->logStream));
        }
        if (heatSource.FirstChildElement("CHP").ToElement()) { // Cognet: Added this, to detect possible errors in xml.
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true; // Cognet: Added this.
            heatingUnit = new CoGeneration(heatSource.FirstChildElement("CHP"),beginDay,endDay,&(this->logStream));
        }
        if (heatSource.FirstChildElement("CHP-HP").ToElement()) { // Cognet: Added this, to detect possible errors in xml.
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true; // Cognet: Added this.
            heatingUnit = new CoGenerationHeatPump(heatSource.FirstChildElement("CHP-HP"),beginDay,endDay,&(this->logStream));
        }
        //beginning of contents added by Dapeng // Cognet: Made modifications.
        if(heatSource.FirstChildElement("Substation").ToElement() ){
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true;
            heatingUnit = Substation::createNewSubstation(heatSource.FirstChildElement("Substation"),  this, beginDay, endDay,  &(this->logStream));
        }
        if(heatSource.FirstChildElement("SubstationHP").ToElement() ){
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true;
            heatingUnit = new SubstationHeatPump(heatSource.FirstChildElement("SubstationHP"),  this, beginDay, endDay,  &(this->logStream));
        }
        if(heatSource.FirstChildElement("SubstationHP2stages").ToElement() ){
            if ( heatSourceFound ) { throw errorMsg; } heatSourceFound = true;
            heatingUnit = new SubstationHeatPump2stages(heatSource.FirstChildElement("SubstationHP2stages"),  this, beginDay, endDay,  &(this->logStream));
        }
        //end of contents added by Dapeng, later changed by Cognet
        if ( not heatSourceFound ) { throw string("Error in the XML file: a Building has a HeatSource containing no units, there must be exactly one."); } // Cognet: Added this.
    }

    // loading the cooling properties
    //beginning of contents added by Dapeng
    if(hdl.FirstChildElement("CoolSource").ToElement()) {
        TiXmlHandle coolSource = hdl.FirstChildElement("CoolSource");
        // gets the working period of the coolSource
        unsigned int beginDay, endDay;
        if ( hdl.FirstChildElement("CoolSource").ToElement()->Attribute("beginDay") == NULL ) beginDay = 1;
        else beginDay = to<unsigned int>(hdl.FirstChildElement("CoolSource").ToElement()->Attribute("beginDay"));
        if ( hdl.FirstChildElement("CoolSource").ToElement()->Attribute("endDay") == NULL )   endDay = 365;
        else endDay = to<unsigned int>(hdl.FirstChildElement("CoolSource").ToElement()->Attribute("endDay"));

        bool coolSourceFound = false; // Cognet: Added this, to detect error if multiple coolsources present in xml.
        string errorMsg = "Error in the XML file: a Building has multiple units in CoolSource, there must be only one."; // Cognet: Added this

        if (coolSource.FirstChildElement("HeatPump").ToElement()) {
            // the air conditioning system
            if ( coolSourceFound ) { throw errorMsg; } coolSourceFound = true; // Cognet: Added this.
            coolingUnit = new HeatPump(coolSource.FirstChildElement("HeatPump"),beginDay,endDay,&(this->logStream));
        }
        //beginning of contents added by Dapeng // Cognet: Made modifications.
        if(coolSource.FirstChildElement("Substation").ToElement() ){ // Cognet: Added this, to detect possible errors in xml.
            if ( coolSourceFound ) { throw errorMsg; } coolSourceFound = true; // Cognet: Added this.
            coolingUnit = Substation::createNewSubstation(coolSource.FirstChildElement("Substation"),  this, beginDay, endDay,  &(this->logStream));
            // DP: simplified compared to previous version; one could check the "supplemental source" is a heat pump (i.e. not a possible cold producing system)
        }
        //end of contents added by Dapeng
        if ( not coolSourceFound ) { throw string("Error in the XML file: a Building has a CoolSource containing no units, there must be exactly one."); } // Cognet: Added this.
    }

    // Cognet: Start added content to impose the heat demanded to the heating unit.
    if(hdl.FirstChildElement("ImposedHeatDemand").ToElement()) {

        TiXmlAttribute* attrib = hdl.FirstChildElement("ImposedHeatDemand").ToElement()->FirstAttribute();
        while ( attrib!=NULL ) {
            string dayHour = attrib->Name();
            double demand;
            if ( attrib->QueryDoubleValue(&demand) ) { throw string("Error in the XML file: a Building has ImposedHeatDemand with attribute"+dayHour+" that isn't a double."); }
            if (imposedHeatDemand.count(dayHour) != 0) {
                throw string("Error in the XML file: a Building has on ImposedHeatDemand multiple times the attribute "+dayHour+", there must be at most one.");
            } else {
                imposedHeatDemand[dayHour] = demand;
            }
            attrib = attrib->Next();
        }
    }
    // Cognet: End added content to impose the heat demanded to the heating unit.


    /* DP : Previous version of EnergyConversionSystem reading
    // loading the heating properties
    if(hdl.FirstChildElement("HeatSource").ToElement()) {
        TiXmlHandle heatSource = hdl.FirstChildElement("HeatSource");
        // gets the working period of the heatSource
        unsigned int beginDay, endDay;
        if ( heatSource.ToElement()->Attribute("beginDay") == NULL ) beginDay = 1;
        else beginDay = to<unsigned int>(heatSource.ToElement()->Attribute("beginDay"));
        if ( heatSource.ToElement()->Attribute("endDay") == NULL )   endDay   = 365;
        else endDay   = to<unsigned int>(heatSource.ToElement()->Attribute("endDay"));
        // different sources of heat
        if (heatSource.FirstChildElement("Boiler").ToElement()) {
            heatingUnit = new Boiler(to<double>(heatSource.FirstChildElement("Boiler").ToElement()->Attribute("Pmax")),
                                     to<double>(heatSource.FirstChildElement("Boiler").ToElement()->Attribute("eta_th")),
                                     beginDay,endDay);
        }
        if (heatSource.FirstChildElement("HeatPump").ToElement()) {

            logStream << "Heat Pump." << endl << flush;

            heatingUnit = new HeatPump(to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("Pmax")),
                                       to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("eta_tech")),
                                       to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("Ttarget")),
                                       beginDay,endDay);

            if (string(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("Tsource"))==string("ground")) {

                logStream << "Ground source." << endl << flush;

                if (heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("z1")) {
                    // vertical pipes in the ground
                    heatingUnit->setGround(to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("depth")),
                                           to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("z1")),
                                           to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("alpha")));
                }
                else {
                    // horizontal pipes in the ground
                    heatingUnit->setGround(to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("depth")),
                                           to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("depth")),
                                           to<double>(heatSource.FirstChildElement("HeatPump").ToElement()->Attribute("alpha")));
                }
            }
            else logStream << "Air source." << endl << flush;

        }
        if (heatSource.FirstChildElement("CHP").ToElement()) {
            logStream << "Cogen." << endl << flush;
            heatingUnit = new CoGeneration(to<double>(heatSource.FirstChildElement("CHP").ToElement()->Attribute("Pmax")),
                                           to<double>(heatSource.FirstChildElement("CHP").ToElement()->Attribute("eta_el")),
                                           to<double>(heatSource.FirstChildElement("CHP").ToElement()->Attribute("eta_th")),
                                           to<double>(heatSource.FirstChildElement("CHP").ToElement()->Attribute("minPartLoadCoeff")),
                                           beginDay,endDay);
        }
        if (heatSource.FirstChildElement("CHP-HP").ToElement()) {
            logStream << "Cogen + HP." << endl << flush;
            heatingUnit = new CoGenerationHeatPump(to<double>(heatSource.FirstChildElement("CHP-HP").ToElement()->Attribute("Pmax")),
                                                   to<double>(heatSource.FirstChildElement("CHP-HP").ToElement()->Attribute("eta_el")),
                                                   to<double>(heatSource.FirstChildElement("CHP-HP").ToElement()->Attribute("eta_th")),
                                                   to<double>(heatSource.FirstChildElement("CHP-HP").ToElement()->Attribute("minPartLoadCoeff")),
                                                   to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("eta_tech")),
                                                   to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("Ttarget")),
                                                   beginDay,endDay);

            if (string(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("Tsource"))==string("ground")) {

                logStream << "Ground source." << endl << flush;
                if (heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("z1")) {
                    // vertical pipes in the ground
                    heatingUnit->setGround(to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("depth")),
                                           to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("z1")),
                                           to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("alpha")));
                }
                else {
                    // horizontal pipes in the ground
                    heatingUnit->setGround(to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("depth")),
                                           to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("depth")),
                                           to<double>(heatSource.FirstChildElement("CHP-HP").FirstChildElement("HeatPump").ToElement()->Attribute("alpha")));
                }
            }
            else logStream << "Air source." << endl << flush;

        }

        //beginning of contents added by Dapeng
        if(heatSource.FirstChildElement("substation").ToElement() ){
            logStream <<"heating substation" <<endl << flush;
            linkedDEC_ForHeat = to<int>(heatSource.FirstChildElement("substation").ToElement()->Attribute("linked"));

            config = "HX";
            TiXmlHandle supplementalUnit = heatSource.FirstChildElement("substation");
            if(supplementalUnit.FirstChildElement("HeatPump").ToElement()) {
                heatPumpElectricPower = to<double>(supplementalUnit.FirstChildElement("HeatPump").ToElement()->Attribute("Pmax") );
                heatPumpEtaTech = to<double>(supplementalUnit.FirstChildElement("HeatPump").ToElement()->Attribute("eta_tech") );
                targetTemp = to<double>(supplementalUnit.FirstChildElement("HeatPump").ToElement()->Attribute("Ttarget") );
                config = "HX+HeatPump";
            }
            else {
                heatPumpElectricPower = 0;
                heatPumpEtaTech = 0;
                targetTemp = 0;
            }
//            logStream <<"test output: " <<heatPumpElectricPower<<" "<<heatPumpEtaTech<<" "<<targetTemp<<endl << flush;
            if(supplementalUnit.FirstChildElement("Boiler").ToElement() ) {
                boilerThermalPower = to<double>(supplementalUnit.FirstChildElement("Boiler").ToElement()->Attribute("Pmax") );
                boilerThermalEfficiency = to<double>(supplementalUnit.FirstChildElement("Boiler").ToElement()->Attribute("eta_th") );
                config = "HX+Boiler";
            }
            else {
                boilerThermalPower = 0;
                boilerThermalEfficiency = 0;
            }
            if(heatSource.FirstChildElement("substation").ToElement()->Attribute("epsilon")==NULL) {
                epsilon=0.0;
                config = "HeatPump";
            }
            else epsilon = to<double>(heatSource.FirstChildElement("substation").ToElement()->Attribute("epsilon"));

            heatDesignThermalPower = to<double>(heatSource.FirstChildElement("substation").ToElement()->Attribute("designThermalPower"));
            heatDesignTempDifference = to<double>(heatSource.FirstChildElement("substation").ToElement()->Attribute("designTempDifference"));

            heatingUnit = new Substation(heatDesignThermalPower, heatDesignTempDifference,
                                         config, epsilon,
                                         boilerThermalPower, boilerThermalEfficiency,
                                         heatPumpElectricPower, heatPumpEtaTech, targetTemp,
                                         beginDay, endDay);
        }
        //end of contents added by Dapeng

    }

    // loading the cooling properties
    //beginning of contents added by Dapeng
    if(hdl.FirstChildElement("CoolSource").ToElement()) {
        TiXmlHandle coolSource = hdl.FirstChildElement("CoolSource");
        // gets the working period of the coolSource
        unsigned int beginDay, endDay;
        if ( hdl.FirstChildElement("CoolSource").ToElement()->Attribute("beginDay") == NULL ) beginDay = 1;
        else beginDay = to<unsigned int>(hdl.FirstChildElement("CoolSource").ToElement()->Attribute("beginDay"));
        if ( hdl.FirstChildElement("CoolSource").ToElement()->Attribute("endDay") == NULL )   endDay = 365;
        else endDay = to<unsigned int>(hdl.FirstChildElement("CoolSource").ToElement()->Attribute("endDay"));
        if(coolSource.FirstChildElement("HeatPump").ToElement() ) {
        // the air conditioning system
            logStream << "AC." << endl << flush;
            coolingUnit = new HeatPump(to<double>(coolSource.FirstChildElement("HeatPump").ToElement()->Attribute("Pmax")),
                                       to<double>(coolSource.FirstChildElement("HeatPump").ToElement()->Attribute("eta_tech")),
                                       to<double>(coolSource.FirstChildElement("HeatPump").ToElement()->Attribute("Ttarget")),
                                       beginDay,endDay);
        }
        if(coolSource.FirstChildElement("substation").ToElement() ) {
            logStream <<"cooling substation" <<endl << flush;
            linkedDEC_ForCool = to<int>(coolSource.FirstChildElement("substation").ToElement()->Attribute("linked"));
            config = "HX";

            TiXmlHandle supplementalUnit = coolSource.FirstChildElement("substation");
            if(supplementalUnit.FirstChildElement("HeatPump").ToElement() ) {
                heatPumpElectricPower = to<double>(supplementalUnit.FirstChildElement("HeatPump").ToElement()->Attribute("Pmax") );
                heatPumpEtaTech = to<double>(supplementalUnit.FirstChildElement("HeatPump").ToElement()->Attribute("eta_tech") );
                targetTemp = to<double>(supplementalUnit.FirstChildElement("HeatPump").ToElement()->Attribute("Ttarget") );
                config = "HX+HeatPump";
            }
            else {
                heatPumpElectricPower = 0;
                heatPumpEtaTech = 0;
                targetTemp = 0;
            }
            if(coolSource.FirstChildElement("substation").ToElement()->Attribute("epsilon")==NULL) {
                epsilon=0.0;
                config = "HeatPump";
            }
            else epsilon = to<double>(coolSource.FirstChildElement("substation").ToElement()->Attribute("epsilon"));

            coolDesignThermalPower = to<double>(coolSource.FirstChildElement("substation").ToElement()->Attribute("designThermalPower")),
            coolDesignTempDifference = to<double>(coolSource.FirstChildElement("substation").ToElement()->Attribute("designTempDifference")),
            coolingUnit = new Substation(coolDesignThermalPower, coolDesignTempDifference,
                                         config, epsilon,
                                         heatPumpElectricPower, heatPumpEtaTech, targetTemp,
                                         beginDay, endDay);
        }
    }
*/

    // THERMAL MODEL
    // creates thermal zone(s)
    if ( hdl.FirstChildElement("Zone").ToElement() ) {
        // zoning exists
        unsigned int zoneIndex = 0;
        while ( hdl.ChildElement("Zone",zoneIndex).ToElement() ) {
            // check if only one zone per building is defined for the E+ simulation
            if (simulateEP && (zoneIndex > 0)) throw(string("Only one Thermal Zone per building allowed when co-simulating with EnergyPlus."));

            // get some details about zone
            if (!hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("id"))
                throw (string("No id attribute in the Zone."));
            unsigned int zoneId = to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("id"));
            if (!hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("volume"))
                throw (string("Zone id ")+toString(zoneId)+string(": no volume attribute."));
            float zoneVolume = to<float>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("volume"));
            // sets if the ground floor is there
            bool groundFloor = false;
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("groundFloor"))
                groundFloor = ( hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("groundFloor") == string("true") );

            // containers of the elements
            vector<Wall*>    zoneWalls;
            vector<Roof*>    zoneRoofs;
            vector<Surface*> zoneSurfaces;
            vector<Floor*>   zoneFloors;

            // read the walls
            unsigned int wallIndex = 0;
            while ( hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).ToElement() ) {

                // add a new wall to the zoneWalls vector
                zoneWalls.push_back(new Wall(hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex), this, &logStream));

                /*TiXmlElement *wallElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).ToElement();
                // read the openableRatio or if not defined, leave it to 0
                float glazingOpenableRatio = 0.f;
                if (wallElem->Attribute("OpenableRatio"))
                    glazingOpenableRatio = to<float>(wallElem->Attribute("OpenableRatio"));
                // creation of the pointer to the WallType defined for this wall
                Composite* pWallType;
                if (wallElem->Attribute("type")) {
                    // gets the type according to the value in "type"
                    pWallType = pDistrict->getComposite(wallElem->Attribute("type"));
                }
                else if (wallElem->Attribute("Uvalue")) {
                    // gets or creates a new WallType in the building for this special wall
                    pWallType = pDistrict->getUvalueComposite(to<float>(wallElem->Attribute("Uvalue")));
                }
                else throw(string("In Wall id: ")+wallElem->Attribute("id")+", no type or Uvalue given.");

                // add a new wall to the zoneWalls vector
                zoneWalls.push_back(new Wall(to<unsigned int>(wallElem->Attribute("id")),pWallType, to<float>(wallElem->Attribute("ShortWaveReflectance")) ,to<float>(wallElem->Attribute("GlazingRatio")),to<float>(wallElem->Attribute("GlazingGValue")),to<float>(wallElem->Attribute("GlazingUValue")),glazingOpenableRatio,&logStream));
                TiXmlElement* pvElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).FirstChildElement("PV").ToElement();
                if (pvElem) {
                    PhotoVoltaic* pv = new PhotoVoltaic(hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).FirstChildElement("PV"));
                    zoneWalls.back() -> setPVPanel(to<float>(pvElem->Attribute("pvRatio")),pv);
                }
                TiXmlElement* pvThermalElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).FirstChildElement("SolarHeater").ToElement();
                if(pvThermalElem){// Add Solar thermal
                    SolarHeater* sh = new SolarHeater(hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).FirstChildElement("SolarHeater"));
                    zoneWalls.back() -> setSolarThermalPanel(to<float>(pvThermalElem->Attribute("solarHeaterRatio")), sh);
                }

                // add the emissivity of present
                if (wallElem->Attribute("LongWaveEmissivity"))
                    zoneWalls.back()->setLongWaveEmissivity(to<float>(wallElem->Attribute("LongWaveEmissivity")));

                // add the ep_id
                if (wallElem->Attribute("ep_id"))
                    zoneWalls.back()->setEp_id(wallElem->Attribute("ep_id"));

                // add the vertices
                wallElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).FirstChildElement("V0").ToElement();
                unsigned int vertexIndex = 0;
                while(wallElem) {
                    // check the presence of coordinates in the vertex tag
                    if (wallElem->Attribute("x")==NULL || wallElem->Attribute("y")==NULL || wallElem->Attribute("z")==NULL) break;
                    zoneWalls.back()->pushVertex(to<float>(wallElem->Attribute("x")),to<float>(wallElem->Attribute("y")),to<float>(wallElem->Attribute("z")));
                    string vertexj = "V" + toString(++vertexIndex);
                    wallElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Wall", wallIndex).FirstChildElement(vertexj.c_str()).ToElement();
                }
                if (vertexIndex<3) throw(string("Wall id=" + toString(zoneWalls.back()->getId()) + " has only " + toString(vertexIndex) + " vertices"));
                zoneWalls.back()->computeNormalAndArea();
*/
                if ((zoneWalls.back()->getArea() > 0.f) && (zoneWalls.back()->getRadius() > 0.f)) {
                    //logStream << "Wall surface loaded";
                    // computes for this wall the eco-indicators and add them to the whole building's values
                    nre += zoneWalls.back()->getNRE();
                    gwp += zoneWalls.back()->getGWP();
                    ubp += zoneWalls.back()->getUBP();
                    // outputs of the calculation
                    //logStream << "Wall id: " << zoneWalls.back()->getId() << "\tNRE: " << zoneWalls.back()->getNRE() << endl << flush;
                }
                else {
                    logStream << "(Warning) Wall id=" << zoneWalls.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                    delete zoneWalls.back();
                    zoneWalls.pop_back();
                }
                // increment the wall index
                ++wallIndex;
            }
            // read the roofs
            unsigned int roofIndex = 0; bool roofUvalueOnly = false;
            while ( hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement() ) {

                // add a new roof to the zoneRoofs vector
                zoneRoofs.push_back(new Roof(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex), this, &logStream));

                // test if the composite is fully described with layers
                if (pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement()->Attribute("type"))->getnLayers()==0) roofUvalueOnly = true;

                /*TiXmlElement* roofElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement();
                // takes care of the glazing openable ratio
                float glazingOpenableRatio = 0.f;
                if (roofElem->Attribute("OpenableRatio"))
                    glazingOpenableRatio = to<float>(roofElem->Attribute("OpenableRatio"));
                // takes care of the type of Roof
                if (!roofElem->Attribute("type"))
                    throw(string("No type given for Roof id=")+roofElem->Attribute("id"));
                // add a new roof to the zoneRoofs vector
                zoneRoofs.push_back(new Roof(to<unsigned int>(roofElem->Attribute("id")), pDistrict->getComposite(roofElem->Attribute("type")), to<float>(roofElem->Attribute("ShortWaveReflectance")), to<float>(roofElem->Attribute("GlazingRatio")),to<float>(roofElem->Attribute("GlazingGValue")),to<float>(roofElem->Attribute("GlazingUValue")),glazingOpenableRatio,&logStream));
                // look for PV element
                TiXmlElement* pvElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).FirstChildElement("PV").ToElement();
                if (pvElem) {
                    PhotoVoltaic* pv = new PhotoVoltaic(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).FirstChildElement("PV"));
                    zoneRoofs.back() -> setPVPanel(to<float>(pvElem->Attribute("pvRatio")),pv);
                }
                TiXmlElement* pvThermalElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).FirstChildElement("SolarHeater").ToElement();
                if(pvThermalElem){// Add Solar thermal
                    SolarHeater* sh = new SolarHeater(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).FirstChildElement("SolarHeater"));
                    logStream << "New SolarHeater built !" << endl << flush;
                    zoneRoofs.back() -> setSolarThermalPanel(to<float>(pvThermalElem->Attribute("solarHeaterRatio")), sh);
                }


                // add the emissivity
                if (roofElem->Attribute("LongWaveEmissivity"))
                    zoneRoofs.back()->setLongWaveEmissivity(to<float>(roofElem->Attribute("LongWaveEmissivity")));

                // add the ep_id
                if (roofElem->Attribute("ep_id"))
                    zoneRoofs.back()->setEp_id(roofElem->Attribute("ep_id"));

                // add the vertices
                roofElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).FirstChildElement("V0").ToElement();
                unsigned int vertexIndex = 0;
                while (roofElem) {
                    // check the presence of coordinates in the vertex tag
                    if (roofElem->Attribute("x")==NULL || roofElem->Attribute("y")==NULL || roofElem->Attribute("z")==NULL) break;
                    zoneRoofs.back()->pushVertex(to<float>(roofElem->Attribute("x")),to<float>(roofElem->Attribute("y")),to<float>(roofElem->Attribute("z")));
                    string vertexj = "V" + toString(++vertexIndex);
                    roofElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).FirstChildElement(vertexj.c_str()).ToElement();
                }
                if (vertexIndex<3) throw(string("Roof id=" + toString(zoneRoofs.back()->getId()) + " has only " + toString(vertexIndex) + " vertices"));
                zoneRoofs.back()->computeNormalAndArea();
                */
                if ((zoneRoofs.back()->getArea() > 0.f) && (zoneRoofs.back()->getRadius() > 0.f)) {
                    // computes for this wall the eco-indicators and add them to the whole building's values
                    if (hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement()->Attribute("type") != NULL) {
                        nre += zoneRoofs.back()->getArea()*(1.f-zoneRoofs.back()->getGlazingRatio())*pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement()->Attribute("type"))->getNRE();
                        gwp += zoneRoofs.back()->getArea()*(1.f-zoneRoofs.back()->getGlazingRatio())*pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement()->Attribute("type"))->getGWP();
                        ubp += zoneRoofs.back()->getArea()*(1.f-zoneRoofs.back()->getGlazingRatio())*pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement()->Attribute("type"))->getUBP();
                        //logStream << "Roof id: " << zoneRoofs.back()->getId() << "\tNRE: " << zoneRoofs.back()->getArea()*(1.f-zoneRoofs.back()->getGlazingRatio())*pDistrict->getType(hdl.ChildElement("Zone",zoneIndex).ChildElement("Roof", roofIndex).ToElement()->Attribute("type"))->getNRE() << endl << flush;
                    }
                }
                else {
                    logStream << "(Warning) Roof id=" << zoneRoofs.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                    delete zoneRoofs.back();
                    zoneRoofs.pop_back();
                }
                ++roofIndex;
            }

            // read the surfaces
            unsigned int surfaceIndex = 0;
            while (hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).ToElement()) {

                // add a new surface to the zoneSurfaces vector
                zoneSurfaces.push_back(new Surface(hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex), this, &logStream));


/*                TiXmlElement* surfaceElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).ToElement();
                zoneSurfaces.push_back(new Surface(to<unsigned int>(surfaceElem->Attribute("id")),to<float>(surfaceElem->Attribute("ShortWaveReflectance")),0.f,0.f,0.f,0.f,&logStream));

                TiXmlElement* pvElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).FirstChildElement("PV").ToElement();
                if (pvElem) {
                    PhotoVoltaic* pv = new PhotoVoltaic(hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).FirstChildElement("PV"));
                    zoneSurfaces.back() -> setPVPanel(to<float>(pvElem->Attribute("pvRatio")),pv);
                }
                TiXmlElement* pvThermalElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).FirstChildElement("SolarHeater").ToElement();
                if(pvThermalElem){// Add Solar thermal
                    SolarHeater* sh = new SolarHeater(hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).FirstChildElement("SolarHeater"));
                    zoneSurfaces.back() -> setSolarThermalPanel(to<float>(pvThermalElem->Attribute("solarHeaterRatio")), sh);
                }

                // add the emissivity
                if (surfaceElem->Attribute("LongWaveEmissivity"))
                    zoneSurfaces.back()->setLongWaveEmissivity(to<float>(surfaceElem->Attribute("LongWaveEmissivity")));

                // add the vertices
                surfaceElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).FirstChildElement("V0").ToElement();
                unsigned int vertexIndex = 0;
                while(surfaceElem) {
                    // check the presence of coordinates in the vertex tag
                    if (surfaceElem->Attribute("x")==NULL || surfaceElem->Attribute("y")==NULL || surfaceElem->Attribute("z")==NULL) break;
                    zoneSurfaces.back()->pushVertex(to<float>(surfaceElem->Attribute("x")),to<float>(surfaceElem->Attribute("y")),to<float>(surfaceElem->Attribute("z")));
                    string vertexj = "V" + toString(++vertexIndex);
                    surfaceElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Surface", surfaceIndex).FirstChildElement(vertexj.c_str()).ToElement();
                }
                if (vertexIndex<3) throw(string("Surface id=" + toString(zoneSurfaces.back()->getId()) + " has only " + toString(vertexIndex) + " vertices"));
                zoneSurfaces.back()->computeNormalAndArea();

*/
                if ((zoneSurfaces.back()->getArea() > 0.f) && (zoneSurfaces.back()->getRadius() > 0.f)) {
                    // that's OK
                }
                else {
                    logStream << "(Warning) Surface id=" << zoneSurfaces.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                    delete zoneSurfaces.back();
                    zoneSurfaces.pop_back();
                }
                ++surfaceIndex;
            }

            // reads the floor area and conductance
            unsigned int floorIndex = 0; bool floorUvalueOnly = false;
            while (hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()) {

                // add a new floor to the zoneFloors vector
                zoneFloors.push_back(new Floor(hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex), this, &logStream));

                // test is the floor is fully described
                if (pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()->Attribute("type"))->getnLayers()==0) floorUvalueOnly = true;

                /*
                TiXmlElement *floorElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement();
                // add a new floor element
                zoneFloors.push_back(new Floor(to<unsigned int>(floorElem->Attribute("id")),pDistrict->getComposite(floorElem->Attribute("type")),&logStream));
                floorElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).FirstChildElement("V0").ToElement();
                unsigned int vertexIndex = 0;
                while(floorElem){
                    // check the presence of coordinates in the vertex tag
                    if (floorElem->Attribute("x")==NULL || floorElem->Attribute("y")==NULL || floorElem->Attribute("z")==NULL) break;
                    zoneFloors.back()->pushVertex(to<float>(floorElem->Attribute("x")),to<float>(floorElem->Attribute("y")),to<float>(floorElem->Attribute("z")));
                    string vertexj = "V" + toString(++vertexIndex);
                    floorElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).FirstChildElement(vertexj.c_str()).ToElement();
                }
                if (vertexIndex<3) throw(string("Floor id=" + toString(zoneFloors.back()->getId()) + " has only " + toString(vertexIndex) + " vertices"));
                zoneFloors.back()->computeNormalAndArea();

                */
                if ((zoneFloors.back()->getArea() > 0.f) && (zoneFloors.back()->getRadius() > 0.f)) {
                    // computes for this wall the eco-indicators and add them to the whole building's values
                    if (hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()->Attribute("type") != NULL) {
                        nre += zoneFloors.back()->getArea()*pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()->Attribute("type"))->getNRE();
                        gwp += zoneFloors.back()->getArea()*pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()->Attribute("type"))->getGWP();
                        ubp += zoneFloors.back()->getArea()*pDistrict->getComposite(hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()->Attribute("type"))->getUBP();
                        //logStream << "Floor id: " << zoneFloors.back()->getId() << "\tNRE: " << zoneFloors.back()->getArea()*pDistrict->getType(hdl.ChildElement("Zone",zoneIndex).ChildElement("Floor",floorIndex).ToElement()->Attribute("type"))->getNRE() << endl << flush;
                    }
                }
                else {
                    logStream << "(Warning) Floor id=" << zoneFloors.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                    delete zoneFloors.back();
                    zoneFloors.pop_back();
                }
                // increments the floor index
                ++floorIndex;
            }

            // some outputs for verification
            logStream << "Zone: " << zoneIndex << "\tWalls: " << wallIndex << "\tRoofs: " << roofIndex << "\tSurfaces: " << surfaceIndex << "\tFloors: " << floorIndex << endl << flush;

            // updates the link matrix in a sparse format
            // a variable to know which one is the first one of the row
            linksAi.push_back( linksAn.size() );
            unsigned int zoneSurfaceIndex = 0;
            while ( hdl.ChildElement("Zone",zoneIndex).ChildElement("ZoneSurface", zoneSurfaceIndex).ToElement() ) {

                TiXmlElement* zoneSurfaceElem = hdl.ChildElement("Zone",zoneIndex).ChildElement("ZoneSurface", zoneSurfaceIndex).ToElement();

                float Swall = to<float>( zoneSurfaceElem->Attribute("Area") );

                float Kwall;
                if ( zoneSurfaceElem->Attribute("Vertical") == string("true") ) {
                    Kwall = Swall/(1.f/pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getConductance()
                                   + 1.f/3.f + 1.f/3.f );
                }
                else { // by default a horizontal surface (buoyant 4.3 and stratified 1.5)
                    Kwall = Swall/(1.f/pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getConductance()
                                   + 1.f/1.5f + 1.f/4.3f );
                }

                float Cwall = Swall*pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getCapacitance();

                // reading the area and the connected zone
                linksAn.push_back( pair<float,float>(Kwall,Cwall) );
                linksAj.push_back( to<unsigned int>( zoneSurfaceElem->Attribute("LinkZone") ) );

                // putting back the eco-indicators for the intermediary walls, half for each zone
                nre += Swall*pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getNRE()/2.f;
                gwp += Swall*pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getGWP()/2.f;
                ubp += Swall*pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getUBP()/2.f;
                logStream << "Slab: " << zoneSurfaceIndex << "\tNRE: " << Swall*pDistrict->getComposite(zoneSurfaceElem->Attribute("type"))->getNRE()/2.f << endl;

                // increment the zoneSurfaceIndex
                ++zoneSurfaceIndex;

            }

            // creates the OCCUPANTS for the Zone
            Occupants *occupants;
            if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()) {

                float occupantsNumber = to<float>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("n"));
                logStream << "Occupants defined: " << occupantsNumber << endl << flush;

                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Stochastic").ToElement()) {
                    // stochastic profile
                    Model::setThermalExplicit(true); // chooses to use the thermal explicit
                    // loads the different arguments for the stochastic model
                    string windowModel = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Models").ToElement()->Attribute("Windows");
                    // start off with the stochastic presence profile
                    vector<float> pMon, pTue, pWed, pThu, pFri, pSat, pSun;
                    TiXmlElement* buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Mon").ToElement();
                    TiXmlAttribute* attrib = buildingElem->FirstAttribute();
                    do {
                        pMon.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Tue").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        pTue.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Wed").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        pWed.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Thu").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        pThu.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Fri").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        pFri.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Sat").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        pSat.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Sun").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        pSun.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );

                    // start off with P01arr
                    vector<float> P01arr, P10arr, P01int, P10int, P01dep, P10dep, duropen;
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("P01arr").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        P01arr.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("P10arr").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        P10arr.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("P01int").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        P01int.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("P10int").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        P10int.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("P01dep").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        P01dep.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    // P10dep
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("P10dep").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        P10dep.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    // duropen
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("duropen").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        duropen.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );

                    // loads the different arguments for the stochastic blinds model
                    vector<float> Plowerarr_UB, Praisearr_UB, Plowerint_UB, Praiseint_UB, Pfulllower_UB, Pfullraise_UB;
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Plowerarr_UB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Plowerarr_UB .push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Praisearr_UB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Praisearr_UB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Plowerint_UB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Plowerint_UB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Praiseint_UB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Praiseint_UB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Pfulllower_UB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Pfulllower_UB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Pfullraise_UB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Pfullraise_UB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );

                    vector<float> Plowerarr_LB, Praisearr_LB, Plowerint_LB, Praiseint_LB, Pfulllower_LB, Pfullraise_LB, DistFrac_LB;
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Plowerarr_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Plowerarr_LB .push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Praisearr_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Praisearr_LB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Plowerint_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Plowerint_LB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Praiseint_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Praiseint_LB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Pfulllower_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Pfulllower_LB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Pfullraise_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Pfullraise_LB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("DistFrac_LB").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        DistFrac_LB.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );

                    // loads the different arguments for the stochastic lights model
                    vector<float> Ponarr, Ponint;
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Ponarr").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Ponarr.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );
                    buildingElem = hdl.FirstChildElement("Occupants").FirstChildElement("Stochastic").FirstChildElement("Ponint").ToElement();
                    attrib = buildingElem->FirstAttribute();
                    do {
                        Ponint.push_back(attrib->DoubleValue());
                    } while ( (attrib = attrib->Next()) );

                    // creates the classes that contains the parameters
                    StochasticPresenceParameters *presParam = new StochasticPresenceParameters(pMon,pTue,pWed,pThu,pFri,pSat,pSun);
                    StochasticWindowParameters *winParam = new StochasticWindowParameters(windowModel, P01arr, P10arr, P01int, P10int, P01dep, P10dep, duropen);
                    StochasticBlindsParameters *blindsParam = new StochasticBlindsParameters(Plowerarr_UB, Praisearr_UB, Plowerint_UB, Praiseint_UB, Pfulllower_UB, Pfullraise_UB, Plowerarr_LB, Praisearr_LB, Plowerint_LB, Praiseint_LB, Pfulllower_LB, Pfullraise_LB, DistFrac_LB);
                    StochasticLightsParameters *lightsParam = new StochasticLightsParameters(Ponarr, Ponint);

                    // creates the occupants
                    occupants = new StochasticOccupantsPresence(occupantsNumber,presParam,winParam,blindsParam,lightsParam);

                }
                else {

                    vector<float> occupantsWeekday, occupantsSaturday, occupantsSunday;

                    TiXmlElement* buildingElem = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Weekday").ToElement();
                    if (buildingElem) {
                        logStream << "Day type=Weekday ";
                        TiXmlAttribute* attrib = buildingElem->FirstAttribute();
                        do {
                            logStream << attrib->DoubleValue() << " ";
                            occupantsWeekday.push_back(attrib->DoubleValue());
                        } while ( (attrib = attrib->Next()) );
                        logStream << endl << flush;
                    }
                    buildingElem = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Saturday").ToElement();
                    if (buildingElem) {
                        logStream << "Day type=Saturday ";
                        TiXmlAttribute* attrib = buildingElem->FirstAttribute();
                        do {
                            logStream << attrib->DoubleValue() << " ";
                            occupantsSaturday.push_back(attrib->DoubleValue());
                        } while ( (attrib = attrib->Next()) );
                        logStream << endl << flush;
                    }
                    buildingElem = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Sunday").ToElement();
                    if (buildingElem) {
                        logStream << "Day type=Sunday ";
                        TiXmlAttribute* attrib = buildingElem->FirstAttribute();
                        do {
                            logStream << attrib->DoubleValue() << " ";
                            occupantsSunday.push_back(attrib->DoubleValue());
                        } while ( (attrib = attrib->Next()) );
                        logStream << endl << flush;
                    }
                    // create the occupants deterministically
                    occupants = new DeterministicOccupantsPresence(occupantsNumber, occupantsWeekday, occupantsSaturday, occupantsSunday);
                }
            }
            else { // no occupants defined, empty
                occupants = new Occupants(0.f);
            }

            // need to check the consistency of this zone, to define its type (1N, 2N, 3N or even 4N)
            if ((roofUvalueOnly && floorUvalueOnly)||(zoneRoofs.empty() && zoneFloors.empty()))
                zones.push_back(new Zone2N(zoneId,this,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,occupants));
            else if ((roofUvalueOnly)||zoneRoofs.empty())
                zones.push_back(new Zone3N_floor(zoneId,this,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,occupants));
            else if ((floorUvalueOnly)||zoneFloors.empty())
                zones.push_back(new Zone3N(zoneId,this,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,occupants));
            else if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("detailedSimulation")
                     && hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("detailedSimulation")==string("true"))
                zones.push_back(new ZoneN(zoneId,this,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,occupants));
            else
                zones.push_back(new Zone4N(zoneId,this,groundFloor,zoneVolume,zoneWalls,zoneRoofs,zoneSurfaces,zoneFloors,occupants));

            // adds the thermal bridges to the zone if they exist
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("psi")) {
                zones.back()->setKpsi(to<float>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("psi")));
            }

            // adds the Tmin and Tmax to the zone if they exist in the Tag Zone or take it from the building itself
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("Tmin")) {
                zones.back()->setTmin(to<float>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("Tmin")));
            }
            else zones.back()->setTmin(to<float>(hdl.ToElement()->Attribute("Tmin")));
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("Tmax")) {
                zones.back()->setTmax(to<float>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("Tmax")));
            }
            else zones.back()->setTmax(to<float>(hdl.ToElement()->Attribute("Tmax")));

            // adds the night ventilation
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("nightVentilationBegin")) {
                zones.back()->setNightVentilationBegin(to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("nightVentilationBegin")));
            }
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("nightVentilationEnd")) {
                zones.back()->setNightVentilationEnd(to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("nightVentilationEnd")));
            }

            // adds the ep_id if it exists
            if (hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("ep_id")) {
                    zones.back()->setEp_id(hdl.ChildElement("Zone",zoneIndex).ToElement()->Attribute("ep_id"));
            }
            else zones.back()->setEp_id("SINGLE_ZONE");

            // setting of the occupants
            if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()) {
                logStream << "Zone occupants" << endl;
                // *** stochastic ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("stochastic"))
                    zones.back()->setOccupantsStochastic(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("stochastic") == string("true"));
                // *** number ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("n"))
                    zones.back()->setOccupantsNumber(to<float>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("n")));
                // *** sensible heat ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("sensibleHeat"))
                    zones.back()->setOccupantsSensibleHeat(to<float>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("sensibleHeat")));
                // *** sensible heat radiant fraction ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("sensibleHeatRadiantFraction"))
                    zones.back()->setOccupantsSensibleHeatRadiantFraction(to<float>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("sensibleHeatRadiantFraction")));
                // *** latent heat ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("latentHeat"))
                    zones.back()->setOccupantsLatentHeat(to<float>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("latentHeat")));
                // *** activity type ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("activityType"))
                    zones.back()->setActivityType(to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("activityType")));
                // *** year profile ID ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("type")) {
                    unsigned int occupantsYearId = to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("type"));
                    if (pDistrict->getOccupancyProfiles()->getYearProfileFromXmlId(occupantsYearId) == nullptr)
                        throw("Occupants' Year Profile type=" + toString(occupantsYearId) + " not found in Zone id=" + toString(zones.back()->getId()));
                    else
                        zones.back()->setOccupantsYearProfile(pDistrict->getOccupancyProfiles()->getYearProfileFromXmlId(occupantsYearId));
                }
                else {
                    logStream << "Old occupants definition\t";
                    // add the day profiles
                    vector<float> dayProfile;
                    DayProfile* weekday=NULL;
                    DayProfile* saturday=NULL;
                    DayProfile* sunday=NULL;
                    // weekday
                    if ( (xmlElem = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Weekday").ToElement()) ) {
                        logStream << "Day zoneId= " << zoneId*10+1 << ": ";
                        unsigned int hourIndex=1;
                        do {
                            string attrib = *(xmlElem->Attribute("p"+toString(hourIndex)));
                            logStream << attrib << " ";
                            dayProfile.push_back(to<float>(attrib));
                        } while ( xmlElem->Attribute("p"+toString(++hourIndex)) );
                        logStream << endl << flush;
                        weekday = new DayProfile(zoneId*10+1,"Day profile weekday "+toString(zoneId*10+1),dayProfile);
                        pDistrict->getOccupancyProfiles()->addDayProfile(weekday);
                        //pDistrict->getOccupancyProfiles()->addDayProfileOld(weekday->getId(),*weekday);
                    }
                    dayProfile.clear();
                    // saturday
                    if ( (xmlElem = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Saturday").ToElement()) ) {
                        logStream << "Day zoneId= " << zoneId*10+2 << ": ";
                        unsigned int hourIndex=1;
                        do {
                            string attrib = *(xmlElem->Attribute("p"+toString(hourIndex)));
                            logStream << attrib << " ";
                            dayProfile.push_back(to<float>(attrib));
                        } while ( xmlElem->Attribute("p"+toString(++hourIndex)) );
                        logStream << endl << flush;
                        saturday = new DayProfile(zoneId*10+2,"Day profile saturday "+toString(zoneId*10+2),dayProfile);
                        pDistrict->getOccupancyProfiles()->addDayProfile(saturday);
                        //pDistrict->getOccupancyProfiles()->addDayProfileOld(saturday->getId(),*saturday);
                    }
                    dayProfile.clear();
                    // sunday
                    if ( (xmlElem = hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").FirstChildElement("Sunday").ToElement()) ) {
                        logStream << "Day zoneId= " << zoneId*10+3 << ": ";
                        unsigned int hourIndex=1;
                        do {
                            string attrib = *(xmlElem->Attribute("p"+toString(hourIndex)));
                            logStream << attrib << " ";
                            dayProfile.push_back(to<float>(attrib));
                        } while ( xmlElem->Attribute("p"+toString(++hourIndex)) );
                        logStream << endl << flush;
                        sunday = new DayProfile(zoneId*10+3,"Day profile sunday "+toString(zoneId*10+3),dayProfile);
                        //pDistrict->getOccupancyProfiles()->addDayProfileOld(sunday->getId(),*sunday);
                    }

                    vector<DayProfile*> yearProfile;
                    logStream << "Occupants Year zoneId= " << zoneId*10 << ": ";
                    for (unsigned int day=1;day<=365;++day) {
                        unsigned int dayOfTheWeek = (day-1) % 7;
                        //logStream << "dayOfTheWeek: " << dayOfTheWeek;
                        if (dayOfTheWeek < 5) { /*logStream << "\toccupancy: " << occupantsWeekday[hour-1] << endl << flush;*/
                            yearProfile.push_back(weekday);
                        }
                        else if (dayOfTheWeek < 6) { /*logStream << "\toccupancy: " << occupantsSaturday[hour-1] << endl << flush;*/
                            yearProfile.push_back(saturday);
                        }
                        else { /*logStream << "\toccupancy: " << occupantsSunday[hour-1] << endl << flush;*/
                            yearProfile.push_back(sunday);
                        }
                        logStream << yearProfile.back()->getId() << " ";
                    }
                    logStream << endl << flush;
                    YearProfile* yp = new YearProfile(zoneId*10,"Year profile "+toString(zoneId*10),yearProfile);
                    pDistrict->getOccupancyProfiles()->addYearProfile(yp);

                    zones.back()->setOccupantsYearProfile(yp);
                }
                // *** DHW year profile ID ***
                if (hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("DHWType")) {
                    logStream << "DHW Year id=" << to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("DHWType")) << endl << flush;
                    zones.back()->setDHWYearProfile(pDistrict->getDHWProfiles()->getYearProfileFromXmlId(to<unsigned int>(hdl.ChildElement("Zone",zoneIndex).FirstChildElement("Occupants").ToElement()->Attribute("DHWType"))));
                }

            }
            // *** output to the screen ***
            logStream << "Occupants number: " << zones.back()->getOccupantsNumber()
                      << "\tsensible heat convective: " << zones.back()->getOccupantsSensibleHeatConvective() << " W " << "/radiative: " << zones.back()->getOccupantsSensibleHeatRadiative() << " W\n"
                      << "\tlatent heat: " << zones.back()->getOccupantsLatentHeat() << " W" << endl;

            // increment the zone index
            ++zoneIndex;
        }
/* // DP: Moved to updateSimulationModelParameters()
        // adds the capacitance value using addZoneC to the corresponding nodes
        for (unsigned int i=0;i<getnLinksAi();++i) {
            // loop on the elements of the rows
            for (unsigned int index=getLinksAi(i);index<getLinksAi(i+1);++index) {

                // prints on the screen some outputs
                logStream << "Zone i: " << i << "\tlinked with zone id: " << getLinksAj(index) << "\tindex: " << getZoneIndexFromId(getLinksAj(index)) << "\tCwall: " << getLinksAn(index).second << endl << flush;

                // adds the capacitance to the correct zone and node
                addZoneC(i, getLinksAn(index).second/2.);
            }
        }
        // end of the creation of the zone and the links between the zones
        logStream << "Zones loaded: " << zones.size() << endl << flush;*/
    }
    else throw(string("No thermal zone given for Building id=")+hdl.ToElement()->Attribute("id"));

    // sets the infiltration rate for all zones
    setNinf(to<float>(hdl.ToElement()->Attribute("Ninf")));

    // output of the eco-indicators
    logStream << "NRE: " << nre << " MJ\tGWP: " << gwp << " kgCO2\tUBP: " << ubp << " pts" << endl;

    update();

    computeHasSolarThermal();

    // end of the building constructor
    return;
}

Building::Building(vector<Wall*> walls, vector<Roof*> roofs, vector<Floor*> floors, vector<Surface*> surfaces, District *pDistrict):pDistrict(pDistrict),id(pDistrict->getnBuildings()),logStream(std::cout.rdbuf()){

    /* Default energy systems for buildings */
    heatStock = new Tank(0.01, 20.0, 1000.0, 4180.0, 20.0, 35.0);
    coldStock = new Tank(0.01, 20.0, 1000.0, 4180.0, 5.0, 20.0);
    heatingUnit = new Boiler(1.0e7, 0.95, 258, 135);
    coolingUnit = new HeatPump(1.0e7, 0.3, 5.0, 136, 257);

    // creates the UNIQUE thermal zone for that building
    float Vi = 0.f;
    // computes the volume under the surfaces
    for (vector<Wall*>::iterator it=walls.begin();it!=walls.end();++it) Vi += (*it)->getVolume();
    for (vector<Roof*>::iterator it=roofs.begin();it!=roofs.end();++it) Vi += (*it)->getVolume();
    for (vector<Floor*>::iterator it=floors.begin();it!=floors.end();++it) Vi += (*it)->getVolume();

    // need to check the consistency of this zone, to define its type (1N, 2N, 3N or even 4N)
    zones.push_back(new Zone4N(id,this,true,Vi,walls,roofs,surfaces,floors,nullptr));

    update();

    computeHasSolarThermal();

    return;
}

Building::~Building() {
    //cout << "delete building " << id << endl;
    // deleting all pointers
    for (vector<Zone*>::iterator it=zones.begin();it!=zones.end();++it) delete *it;
    //cout << "zones deleted" << endl;
    if(heatStock != NULL) delete heatStock;
    if(coldStock != NULL) delete coldStock;
    if (heatingUnit != NULL) delete heatingUnit;
    if (coolingUnit != NULL) delete coolingUnit;
    //cout << "ecs deleted" << endl;
}

void Building::update() {

    // update the zones
    for (vector<Zone*>::iterator it=zones.begin();it!=zones.end();++it) (*it)->update();

    // DP: Add capacitance to zones based on connections with other zones ?
    // adds the capacitance value using addZoneC to the corresponding nodes
    for (unsigned int i=0;i<getnLinksAi();++i) {
        // loop on the elements of the rows
        for (unsigned int index=getLinksAi(i);index<getLinksAi(i+1);++index) {

            // prints on the screen some outputs
            logStream << "Zone i: " << i << "\tlinked with zone id: " << getLinksAj(index) << "\tindex: " << getZoneIndexFromId(getLinksAj(index)) << "\tCwall: " << getLinksAn(index).second << endl << flush;

            // adds the capacitance to the correct zone and node
            addZoneC(i, getLinksAn(index).second/2.);
        }
    }

    // end of the creation of the zone and the links between the zones
    logStream << "Zones loaded: " << zones.size() << endl << flush;

    // creates the three matrices matrices that describe the problem thermally
    const size_t NP = getnNodes();
    // matrix of the capacitances C and constant conductances G1
    C = vector<vector<float>>(NP, vector<float>(NP,0.f));
    G1 = vector<vector<float>>(NP, vector<float>(NP,0.f));
    // the initialisation of the matrices C and G is done with the fixed values of the conductances/capacitances
    for (unsigned int i=0; i<zones.size(); ++i) {
        for (unsigned int j=0; j<zones[i]->getnNodes(); ++j) {
            // put the values in the C matrix
            C[getMatrixPosition(i)+j][getMatrixPosition(i)+j] = zones[i]->getC(j);
            for (unsigned int k=0; k<zones[i]->getnNodes(); ++k) {
                // put the values in the G1 matrix
                G1[getMatrixPosition(i)+j][getMatrixPosition(i)+k] = zones[i]->getFixedMatrixElement(j,k);
            }
        }
    }
    // loads the matrix of the conductances of the links between the zones (addition to the present fixed conductances matrix)
    loadLinkSparse(G1);
    //logStream << "G1: " << *G1 << endl;
}

void Building::writeXML(ofstream& file, string tab){
    file << tab << "<Building id=\"" << id << "\" key=\"" << key << "\" Vi=\"" << getVolume();

    streamsize ss = file.precision();
    file.precision(6); // use max 6 significant numbers for floats...
    file << "\" Ninf=\""<< getNinf() << "\" BlindsLambda=\"" << blindsLambda << "\" BlindsIrradianceCutOff=\"" << blindsIrradianceCutOff;
    file.precision(ss); // restore default
    file << "\" Simulate=\"true\">" << endl;

    // Energy system not considered
    string subtab=tab+"\t";
    heatStock->writeXML(file,"HeatTank",subtab);
    coldStock->writeXML(file,"CoolTank",subtab);
    if(heatingUnit!=NULL){
        file << subtab << "<HeatSource beginDay=\"" << heatingUnit->getBeginDay() << "\" endDay=\"" << heatingUnit->getEndDay() << "\">" << endl;
        heatingUnit->writeXML(file,subtab+"\t");
        file << subtab << "</HeatSource>" << endl;
    }
    if(coolingUnit!=NULL){
        file << subtab << "<CoolSource beginDay=\"" << coolingUnit->getBeginDay() << "\" endDay=\"" << coolingUnit->getEndDay() << "\">" << endl;
        coolingUnit->writeXML(file,subtab+"\t");
        file << subtab << "</CoolSource>" << endl;
    }
    for (unsigned int i=0; i < zones.size(); ++i)
        zones[i]->writeXML(file,subtab);
    file << tab << "</Building>" << endl;
}

void Building::writeGML(ofstream& file, string tab) {

    // output of the PV panel
    for (size_t i=0; i < zones.size(); ++i) {
        // writes the different surface elements
        for (size_t j=0; j < zones[i]->getnWalls(); ++j) {
            if (zones[i]->getWall(j)->getPVRatio() > 0.f) {
                file << tab << "<core:cityObjectMember>\n"
                     << tab << "\t<energy:PhotovoltaicSystem gml:id=\"PV_1\">\n"
                     << tab << "\t\t<energy:nominalEfficiency uom=\"ratio\">" << zones[i]->getWall(j)->getPVPanel()->getMaxPowerEfficiency(800.,20.) << "</energy:nominalEfficiency>\n"
                     << tab << "\t\t<energy:collectorSurface uom=\"m2\">" << zones[i]->getWall(j)->getPVRatio()*zones[i]->getWall(j)->getArea() << "</energy:collectorSurface>\n"
                     << tab << "\t\t<gen:measureAttribute name=\"panelAzimuth\">\n"
                     << tab << "\t\t\t<gen:value uom=\"deg\">" << zones[i]->getWall(j)->getAzimuth() << "</gen:value>\n"
                     << tab << "\t\t</gen:measureAttribute>\n"
                     << tab << "\t\t<gen:measureAttribute name=\"panelInclination\">\n"
                     << tab << "\t\t\t<gen:value uom=\"deg\">\n" << zones[i]->getWall(j)->getAltitude() << "</gen:value>\n"
                     << tab << "\t\t</gen:measureAttribute>\n"
                     << tab << "\t\t<energy:installedOnBoundarySurface xlink:href=\"#";
                    if (zones[i]->getWall(j)->getKey().empty())
                        file << "Wall_" << zones[i]->getWall(j)->getId() << "\">" << endl;
                    else
                        file << zones[i]->getWall(j)->getKey() << "\">" << endl;
                file << tab << "\t</energy:PhotovoltaicSystem>\n"
                     << tab << "</core:cityObjectMember>\n" << flush;
            }
        }
        for (size_t j=0; j < zones[i]->getnRoofs(); ++j) {
            if (zones[i]->getRoof(j)->getPVRatio() > 0.f) {
                file << tab << "<core:cityObjectMember>\n"
                     << tab << "\t<energy:PhotovoltaicSystem gml:id=\"PV_1\">\n"
                     << tab << "\t\t<energy:nominalEfficiency uom=\"ratio\">" << zones[i]->getRoof(j)->getPVPanel()->getMaxPowerEfficiency(800.,20.) << "</energy:nominalEfficiency>\n"
                     << tab << "\t\t<energy:collectorSurface uom=\"m2\">" << zones[i]->getRoof(j)->getPVRatio()*zones[i]->getRoof(j)->getArea() << "</energy:collectorSurface>\n"
                     << tab << "\t\t<gen:measureAttribute name=\"panelAzimuth\">\n"
                     << tab << "\t\t\t<gen:value uom=\"deg\">" << zones[i]->getRoof(j)->getAzimuth() << "</gen:value>\n"
                     << tab << "\t\t</gen:measureAttribute>\n"
                     << tab << "\t\t<gen:measureAttribute name=\"panelInclination\">\n"
                     << tab << "\t\t\t<gen:value uom=\"deg\">" << zones[i]->getRoof(j)->getAltitude() << "</gen:value>\n"
                     << tab << "\t\t</gen:measureAttribute>\n"
                     << tab << "\t\t<energy:installedOnBoundarySurface xlink:href=\"#";
                    if (zones[i]->getRoof(j)->getKey().empty())
                        file << "Roof_" << zones[i]->getRoof(j)->getId() << "\">" << endl;
                    else
                        file << zones[i]->getRoof(j)->getKey() << "\">" << endl;
                file << tab << "\t</energy:PhotovoltaicSystem>\n"
                     << tab << "</core:cityObjectMember>\n" << flush;
            }
        }
    }

    file << tab << "<core:cityObjectMember>" << endl;
    file << tab << "\t<bldg:Building gml:id=\"";
    if (key.empty())
        file << "Bldg-" << id << "\">" << endl;
    else
        file << key << "\">" << endl;

    string subtab=tab+"\t\t";

    // add the energy demand computed by CitySim
    if (pDistrict->getScene()->getTimeStepsSimulated() > 0) {
        // add the energy demand
        file << subtab << "<energy:demands>" << endl;
        // heating
        file << subtab << tabs(1) << "<energy:EnergyDemand>" << endl;
        file << subtab << tabs(2) << "<energy:energyAmount>" << endl;
        file << subtab << tabs(3) << "<energy:RegularTimeSeries>" << endl;
        // properties to define the source of data
        file << subtab << tabs(4) << "<energy:variableProperties>" << endl;
        file << subtab << tabs(5) << "<energy:TimeValuesProperties>" << endl;
        file << subtab << tabs(6) << "<energy:acquisitionMethod>simulation</energy:acquisitionMethod>" << endl;
        file << subtab << tabs(6) << "<energy:interpolationType>averageInSucceedingInterval</energy:interpolationType>" << endl;
        file << subtab << tabs(6) << "<energy:source>CitySim</energy:source>" << endl;
        file << subtab << tabs(6) << "<energy:thematicDescription>Heating energy</energy:thematicDescription>" << endl;
        file << subtab << tabs(5) << "</energy:TimeValuesProperties>" << endl;
        file << subtab << tabs(4) << "</energy:variableProperties>" << endl;
        // energy amount for heating
        file << subtab << tabs(4) << "<energy:temporalExtent></energy:temporalExtent>" << endl;
        file << subtab << tabs(4) << "<energy:timeInterval unit=\"hour\">1</energy:timeInterval>" << endl;
        // loop on the number of time steps
        file << subtab << tabs(4) << "<energy:values uom=\"Wh\">";
        for (unsigned int i=0; i<pDistrict->getScene()->getTimeStepsSimulated(); ++i)
             file << getHeating(i) << " ";
        file << "\n" << subtab << tabs(4) << "</energy:values>" << endl;
        // close the tab EnergyDemand
        file << subtab << tabs(3) << "</energy:RegularTimeSeries>" << endl;
        file << subtab << tabs(2) << "</energy:energyAmount>" << endl;
        file << subtab << tabs(2) << "<energy:endUse>spaceHeating</energy:endUse>" << endl;
//        // adds the ECS
//        file << subtab << tabs(2) << "<energy:isProvidedBy>" << endl;
//        if (heatingUnit) heatingUnit->writeGML(file,subtab+tabs(3));
//        file << subtab << tabs(2) << "</energy:isProvidedBy>" << endl;
        file << subtab << tabs(1) << "</energy:EnergyDemand>" << endl;
        file << subtab << "</energy:demands>" << endl;
        // add the cooling
        file << subtab << "<energy:demands>" << endl;
        file << subtab << tabs(1) << "<energy:EnergyDemand>" << endl;
        file << subtab << tabs(2) << "<energy:energyAmount>" << endl;
        file << subtab << tabs(3) << "<energy:RegularTimeSeries>" << endl;
        // properties to define the source of data
        file << subtab << tabs(4) << "<energy:variableProperties>" << endl;
        file << subtab << tabs(5) << "<energy:TimeValuesProperties>" << endl;
        file << subtab << tabs(6) << "<energy:acquisitionMethod>simulation</energy:acquisitionMethod>" << endl;
        file << subtab << tabs(6) << "<energy:interpolationType>averageInSucceedingInterval</energy:interpolationType>" << endl;
        file << subtab << tabs(6) << "<energy:source>CitySim</energy:source>" << endl;
        file << subtab << tabs(6) << "<energy:thematicDescription>Cooling energy</energy:thematicDescription>" << endl;
        file << subtab << tabs(5) << "</energy:TimeValuesProperties>" << endl;
        file << subtab << tabs(4) << "</energy:variableProperties>" << endl;
        // energy amount for cooling
        file << subtab << tabs(4) << "<energy:temporalExtent></energy:temporalExtent>" << endl;
        file << subtab << tabs(4) << "<energy:timeInterval unit=\"hour\">1</energy:timeInterval>" << endl;
        // loop on the number of time steps
        file << subtab << tabs(4) << "<energy:values uom=\"Wh\">";
        for (unsigned int i=0; i<pDistrict->getScene()->getTimeStepsSimulated(); ++i)
            file << max(-getCooling(i),0.) << " ";
        file << "\n" << subtab << tabs(4) << "</energy:values>" << endl;
        file << subtab << tabs(3) << "</energy:RegularTimeSeries>" << endl;
        file << subtab << tabs(2) << "</energy:energyAmount>" << endl;
        file << subtab << tabs(2) << "<energy:endUse>spaceCooling</energy:endUse>" << endl;
//        // adds the ECS
//        file << subtab << tabs(2) << "<energy:isProvidedBy>" << endl;
//        if (coolingUnit) coolingUnit->writeGML(file,subtab+tabs(3));
//        file << subtab << tabs(2) << "</energy:isProvidedBy>" << endl;
        file << subtab << tabs(1) << "</energy:EnergyDemand>" << endl;
        file << subtab << "</energy:demands>" << endl;
    }

    // lod2Solid constructions from the GML
    file << tab << "\t<bldg:lod2Solid>\n"
         << tab << "\t\t<gml:Solid>\n"
         << tab << "\t\t\t<gml:exterior>\n"
         << tab << "\t\t\t\t<gml:CompositeSurface>\n";

    subtab=tab+"\t";

    for (size_t i=0; i < zones.size(); ++i) {
        // writes the different surface elements
        for (size_t j=0; j < zones[i]->getnWalls(); ++j) {
            file << subtab << "\t\t\t\t<gml:surfaceMember xlink:href=\"#";
            if (zones[i]->getWall(j)->getKey().empty())
                file << "b" << id << "_p_w_" << zones[i]->getWall(j)->getId() << "\">" << endl;
            else
                file << zones[i]->getWall(j)->getKey() << "\">" << endl;
        }
        for (size_t j=0; j < zones[i]->getnRoofs(); ++j) {
            file << subtab << "\t\t\t\t<gml:surfaceMember xlink:href=\"#";
            if (zones[i]->getRoof(j)->getKey().empty())
                file << "b" << id << "_p_r_" << zones[i]->getRoof(j)->getId() << "\">" << endl;
            else
                file << zones[i]->getRoof(j)->getKey() << "\">" << endl;
        }
        for (size_t j=0; j < zones[i]->getnFloors(); ++j) {
            file << subtab << "\t\t\t\t<gml:surfaceMember xlink:href=\"#";
            if (zones[i]->getFloor(j)->getKey().empty())
                file << "b" << id << "_p_g_" << zones[i]->getFloor(j)->getId() << "\">" << endl;
            else
                file << zones[i]->getFloor(j)->getKey() << "\">" << endl;
        }
    }

    file << tab << "\t\t\t\t</gml:CompositeSurface>\n"
         << tab << "\t\t\t</gml:exterior>\n"
         << tab << "\t\t</gml:Solid>\n"
         << tab << "\t</bldg:lod2Solid>" << endl;

    // add the bounded by tags, in order to recognize Wall, Roof and Floor
    subtab=tab+"\t";
    for (size_t i=0; i < zones.size(); ++i) {
        // writes the different surface elements
        for (size_t j=0; j < zones[i]->getnWalls(); ++j) {
            file << subtab << "<bldg:boundedBy>\n"
                 << subtab << "\t<bldg:WallSurface gml:id=\"Wall_" << zones[i]->getWall(j)->getId() << "\">\n"
                 << subtab << "\t\t<bldg:lod2MultiSurface>\n"
                 << subtab << "\t\t\t<gml:MultiSurface>\n"
                 << subtab << "\t\t\t\t<gml:surfaceMember>\n"
                 << subtab << "\t\t\t\t\t<gml:Polygon ";
                if (zones[i]->getWall(j)->getKey().empty())
                    file << "gml:id=\"b" << id << "_p_w_" << zones[i]->getWall(j)->getId() << "\">" << endl;
                else
                    file << "gml:id=\"" << zones[i]->getWall(j)->getKey() << "\">" << endl;
                zones[i]->getWall(j)->writeGML(file,subtab+"\t");
                file << subtab << "\t\t\t\t\t</gml:Polygon>\n"
                 << subtab << "\t\t\t\t</gml:surfaceMember>\n"
                 << subtab << "\t\t\t</gml:MultiSurface>\n"
                 << subtab << "\t\t</bldg:lod2MultiSurface>\n"
// removed in version 1.0 (to be checked)
//                 << subtab << "\t\t<energy:globalSolarIrradiance>\n"
//                 << subtab << "\t\t\t<energy:RegularTimeSeriesFile>\n"
//                 << subtab << "\t\t\t\t<energy:uom uom=\"W/m2\"/>\n"
//                 << subtab << "\t\t\t\t<energy:file>" << pDistrict->getScene()->getInputFileNoExtNoPath() << "_SW.tsv" << "</energy:file>\n"
//                 << subtab << "\t\t\t\t<energy:temporalExtent></energy:temporalExtent>\n"
//                 << subtab << "\t\t\t\t<energy:timeInterval unit=\"hour\">1</energy:timeInterval>\n"
//                 << subtab << "\t\t\t\t<energy:numberOfHeaderLines>1</energy:numberOfHeaderLines>\n"
//                 << subtab << "\t\t\t\t<energy:valueColumnNumber>" << pDistrict->getScene()->getColumnIndex(zones[i]->getWall(j)) << "</energy:valueColumnNumber>\n"
//                 << subtab << "\t\t\t\t<energy:fieldSeparator>\\t</energy:fieldSeparator>\n"
//                 << subtab << "\t\t\t</energy:RegularTimeSeriesFile>\n"
//                 << subtab << "\t\t</energy:globalSolarIrradiance>\n"
                 << flush;
            file << subtab << "\t</bldg:WallSurface>\n"
                 << subtab << "</bldg:boundedBy>" << endl;
        }
        for (size_t j=0; j < zones[i]->getnRoofs(); ++j) {
            file << subtab << "<bldg:boundedBy>\n"
                 << subtab << "\t<bldg:RoofSurface gml:id=\"Roof_" << zones[i]->getRoof(j)->getId() << "\">\n"
                 << subtab << "\t\t<bldg:lod2MultiSurface>\n"
                 << subtab << "\t\t\t<gml:MultiSurface>\n"
                 << subtab << "\t\t\t\t<gml:surfaceMember>\n"
                 << subtab << "\t\t\t\t\t<gml:Polygon ";
                if (zones[i]->getRoof(j)->getKey().empty())
                    file << "gml:id=\"b" << id << "_p_r_" << zones[i]->getRoof(j)->getId() << "\">" << endl;
                else
                    file << "gml:id=\"" << zones[i]->getRoof(j)->getKey() << "\">" << endl;
                zones[i]->getRoof(j)->writeGML(file,subtab+"\t");
                file << subtab << "\t\t\t\t\t</gml:Polygon>\n"
                 << subtab << "\t\t\t\t</gml:surfaceMember>\n"
                 << subtab << "\t\t\t</gml:MultiSurface>\n"
                 << subtab << "\t\t</bldg:lod2MultiSurface>\n"
//                 << subtab << "\t\t<energy:globalSolarIrradiance>\n"
//                 << subtab << "\t\t\t<energy:RegularTimeSeriesFile>\n"
//                 << subtab << "\t\t\t\t<energy:uom uom=\"W/m2\"/>\n"
//                 << subtab << "\t\t\t\t<energy:file>" << pDistrict->getScene()->getInputFileNoExtNoPath() << "_SW.tsv" << "</energy:file>\n"
//                 << subtab << "\t\t\t\t<energy:temporalExtent></energy:temporalExtent>\n"
//                 << subtab << "\t\t\t\t<energy:timeInterval unit=\"hour\">1</energy:timeInterval>\n"
//                 << subtab << "\t\t\t\t<energy:numberOfHeaderLines>1</energy:numberOfHeaderLines>\n"
//                 << subtab << "\t\t\t\t<energy:valueColumnNumber>" << pDistrict->getScene()->getColumnIndex(zones[i]->getRoof(j)) << "</energy:valueColumnNumber>\n"
//                 << subtab << "\t\t\t\t<energy:fieldSeparator>\\t</energy:fieldSeparator>\n"
//                 << subtab << "\t\t\t</energy:RegularTimeSeriesFile>\n"
//                 << subtab << "\t\t</energy:globalSolarIrradiance>\n"
                 << flush;
            file << subtab << "\t</bldg:RoofSurface>\n"
                 << subtab << "</bldg:boundedBy>" << endl;
        }
        for (size_t j=0; j < zones[i]->getnFloors(); ++j) {
            file << subtab << "<bldg:boundedBy>\n"
                 << subtab << "\t<bldg:GroundSurface gml:id=\"Floor_" << zones[i]->getFloor(j)->getId() << "\">\n"
                 << subtab << "\t\t<bldg:lod2MultiSurface>\n"
                 << subtab << "\t\t\t<gml:MultiSurface>\n"
                 << subtab << "\t\t\t\t<gml:surfaceMember>\n"
                 << subtab << "\t\t\t\t\t<gml:Polygon ";
                if (zones[i]->getFloor(j)->getKey().empty())
                    file << "gml:id=\"b" << id << "_p_g_" << zones[i]->getFloor(j)->getId() << "\">" << endl;
                else
                    file << "gml:id=\"" << zones[i]->getFloor(j)->getKey() << "\">" << endl;
                zones[i]->getFloor(j)->writeGML(file,subtab+"\t");
                file << subtab << "\t\t\t\t\t</gml:Polygon>\n"
                 << subtab << "\t\t\t\t</gml:surfaceMember>\n"
                 << subtab << "\t\t\t</gml:MultiSurface>\n"
                 << subtab << "\t\t</bldg:lod2MultiSurface>\n"
                 << subtab << "\t</bldg:GroundSurface>\n"
                 << subtab << "</bldg:boundedBy>" << endl;
        }
    }

    // energy adds-on - V1.0 volume specification
    file << subtab << "<energy:volume>\n"
         << subtab << "\t<energy:VolumeType>\n"
         << subtab << "\t\t<energy:type>grossVolume</energy:type>\n"
         << subtab << "\t\t<energy:value uom=\"m3\">" << getVolume() << "</energy:value>\n"
         << subtab << "\t</energy:VolumeType>\n"
         << subtab << "</energy:volume>\n" << endl;

    // adds the ThermalZones
    file << subtab << "<energy:thermalZone>" << endl;
    for (size_t i=0; i < zones.size(); ++i) {
        file << subtab << "\t<energy:ThermalZone gml:id=\"TZ_" << zones.at(i)->getId() << "\">" << endl;
        file << subtab << "\t\t<energy:additionalThermalBridgeUValue uom=\"W/(m2K)\">" << zones.at(i)->getKpsi()/(zones.at(i)->getSwi()+zones.at(i)->getWallArea()+zones.at(i)->getRoofArea()+zones.at(i)->getFloorArea()) << "</energy:additionalThermalBridgeUValue>" << endl;
        file << subtab << "\t\t<energy:infiltrationRate uom=\"1/h\">" << zones.at(i)->getNinf() << "</energy:infiltrationRate>" << endl;
        file << subtab << "\t\t<energy:isCooled>" << (this->getCoolingUnit()?"true":"false") << "</energy:isCooled>" << endl;
        file << subtab << "\t\t<energy:isHeated>" << (this->getHeatingUnit()?"true":"false") << "</energy:isHeated>" << endl;
        // writes the different surface elements
        for (size_t j=0; j < zones[i]->getnWalls(); ++j) {

            // the zone is bounded by different elements
            file << subtab << "\t\t<energy:boundedBy>" << endl;
            file << subtab << "\t\t\t<energy:ThermalBoundary>" << endl;
            file << subtab << "\t\t\t\t<energy:thermalBoundaryType>outerWall</energy:thermalBoundaryType>" << endl;

            zones[i]->getWall(j)->writeGML_composedOf(file,subtab+"\t\t\t\t");

            // the partOf is made if you have one surface that belongs to two different zones (e.g. ZoneSurface)
            file << subtab << "\t\t\t\t<energy:delimits xlink:href=\"#TZ_" << zones.at(i)->getId() << "\"/>" << endl;

            // this is where the link is given to the boundedBy Surface
            file << subtab << tabs(4) << "<energy:relatesTo xlink:href=\"#Wall_" << zones[i]->getWall(j)->getId() << "\"/>" << endl;

            file << subtab << tabs(3) << "</energy:ThermalBoundary>" << endl;
            file << subtab << tabs(2) << "</energy:boundedBy>" << endl;
        }
        for (size_t j=0; j < zones[i]->getnRoofs(); ++j) {

            // the zone is bounded by different elements
            file << subtab << tabs(2) << "<energy:boundedBy>" << endl;
            file << subtab << tabs(3) << "<energy:ThermalBoundary>" << endl;
            file << subtab << tabs(4) << "<energy:thermalBoundaryType>roof</energy:thermalBoundaryType>" << endl;

            zones.at(i)->getRoof(j)->writeGML_composedOf(file,subtab+"\t\t\t\t");

            // the partOf is made if you have one surface that belongs to two different zones (e.g. ZoneSurface)
            file << subtab << "\t\t\t\t<energy:delimits xlink:href=\"TZ_" << zones.at(i)->getId() << "\"/>" << endl;

            // this is where the link is given to the boundedBy Surface
            file << subtab << tabs(4) << "<energy:relatesTo xlink:href=\"#Roof_" << zones[i]->getRoof(j)->getId() << "\"/>" << endl;
            file << subtab << tabs(3) << "</energy:ThermalBoundary>" << endl;
            file << subtab << tabs(2) << "</energy:boundedBy>" << endl;
        }
        for (size_t j=0; j < zones[i]->getnFloors(); ++j) {

            // the zone is bounded by different elements
            file << subtab << tabs(2) << "<energy:boundedBy>" << endl;
            file << subtab << tabs(3) << "<energy:ThermalBoundary>" << endl;
            file << subtab << tabs(4) << "<energy:thermalBoundaryType>groundSlab</energy:thermalBoundaryType>" << endl;

            zones[i]->getFloor(j)->writeGML_composedOf(file,subtab+"\t\t\t\t");

            // the partOf is made if you have one surface that belongs to two different zones (e.g. ZoneSurface)
            file << subtab << "\t\t\t\t<energy:delimits xlink:href=\"TZ_" << zones.at(i)->getId() << "\"/>" << endl;

            // this is where the link is given to the boundedBy Surface
            file << subtab << tabs(4) << "<energy:relatesTo xlink:href=\"#Floor_" << zones[i]->getFloor(j)->getId() << "\"/>" << endl;
            file << subtab << tabs(3) << "</energy:ThermalBoundary>" << endl;
            file << subtab << tabs(2) << "</energy:boundedBy>" << endl;
        }
        //file << subtab << tabs(2) << "<energy:relates xlink:href=\"#ID_" << zones.at(i)->getId() << "\"/>" << endl;
        file << subtab << tabs(1) << "</energy:ThermalZone>" << endl;
    }
    file << subtab << "</energy:thermalZone>" << endl;

    // add the occupancy through the UsageZone
    file << subtab << "<energy:usageZone>" << endl;
    for (size_t i=0; i < zones.size(); ++i) {
        file << subtab << tabs(1) << "<energy:UsageZone gml:id=\"UZ_" << zones.at(i)->getId() << "\">" << endl;
        // cooling schedule
        file << subtab << tabs(2) << "<energy:coolingSchedule>" << endl;
        file << subtab << tabs(3) << "<energy:ConstantValueSchedule>" << endl;
        file << subtab << tabs(4) << "<energy:averageValue uom=\"celsius\">" << zones.at(i)->getTmax() << "</energy:averageValue>" << endl;
        file << subtab << tabs(3) << "</energy:ConstantValueSchedule>" << endl;
        file << subtab << tabs(2) << "</energy:coolingSchedule>" << endl;
        // heating schedule
        file << subtab << tabs(2) << "<energy:heatingSchedule>" << endl;
        file << subtab << tabs(3) << "<energy:ConstantValueSchedule>" << endl;
        file << subtab << tabs(4) << "<energy:averageValue uom=\"celsius\">" << zones.at(i)->getTmin() << "</energy:averageValue>" << endl;
        file << subtab << tabs(3) << "</energy:ConstantValueSchedule>" << endl;
        file << subtab << tabs(2) << "</energy:heatingSchedule>" << endl;
        // specific attributes
        file << subtab << tabs(2) << "<energy:usageZoneType>residential</energy:usageZoneType>" << endl;
        // ventilation schedule
        file << subtab << tabs(2) << "<energy:ventilationSchedule>" << endl;
        file << subtab << tabs(3) << "<energy:ConstantValueSchedule>" << endl;
        file << subtab << tabs(4) << "<energy:averageValue uom=\"1/h\">" << zones.at(i)->getNinf() << "</energy:averageValue>" << endl;
        file << subtab << tabs(3) << "</energy:ConstantValueSchedule>" << endl;
        file << subtab << tabs(2) << "</energy:ventilationSchedule>" << endl;
        // floor area
        file << subtab << tabs(2) << "<energy:floorArea>" << endl;
        file << subtab << tabs(3) << "<energy:FloorArea>" << endl;
        file << subtab << tabs(4) << "<energy:type>grossFloorArea</energy:type>" << endl;
        file << subtab << tabs(4) << "<energy:value uom=\"m2\">" << this->getFloorArea() << "</energy:value>" << endl;
        file << subtab << tabs(3) << "</energy:FloorArea>" << endl;
        file << subtab << tabs(2) << "</energy:floorArea>" << endl;
        // occupancy in terms of people
        file << subtab << tabs(2) << "<energy:occupiedBy>" << endl;
        file << subtab << tabs(3) << "<energy:Occupants>" << endl;
        file << subtab << tabs(4) << "<energy:heatDissipation>" << endl;
        file << subtab << tabs(5) << "<energy:HeatExchangeType>" << endl;
        file << subtab << tabs(6) << "<energy:convectiveFraction uom=\"ratio\">" << zones.at(i)->getOccupantsSensibleHeat()*(1.f-zones.at(i)->getOccupantsSensibleHeatRadiantFraction())/(zones.at(i)->getOccupantsSensibleHeat()+zones.at(i)->getOccupantsLatentHeat()) << "</energy:convectiveFraction>" << endl;
        file << subtab << tabs(6) << "<energy:latentFraction uom=\"ratio\">" << zones.at(i)->getOccupantsLatentHeat()/(zones.at(i)->getOccupantsSensibleHeat()+zones.at(i)->getOccupantsLatentHeat()) << "</energy:latentFraction>" << endl;
        file << subtab << tabs(6) << "<energy:radiantFraction uom=\"ratio\">" << zones.at(i)->getOccupantsSensibleHeat()*zones.at(i)->getOccupantsSensibleHeatRadiantFraction()/(zones.at(i)->getOccupantsSensibleHeat()+zones.at(i)->getOccupantsLatentHeat()) << "</energy:radiantFraction>" << endl;
        file << subtab << tabs(6) << "<energy:totalValue uom=\"W\">" << zones.at(i)->getOccupantsNumber()*(zones.at(i)->getOccupantsSensibleHeat()+zones.at(i)->getOccupantsLatentHeat()) << "</energy:totalValue>" << endl;
        file << subtab << tabs(5) << "</energy:HeatExchangeType>" << endl;
        file << subtab << tabs(4) << "</energy:heatDissipation>" << endl;
        file << subtab << tabs(4) << "<energy:numberOfOccupants>" << round(zones.at(i)->getOccupantsNumber()) << "</energy:numberOfOccupants>" << endl;
        /* -> compute occupancyRate according to stochastic and deterministic profiles ? */
        file << subtab << tabs(4) << "<energy:occupancyRate>" << endl;
        file << subtab << tabs(5) << "<energy:TimeSeriesSchedule>" << endl;
        file << subtab << tabs(6) << "<energy:timeDependingValues>" << endl;
        file << subtab << tabs(7) << "<energy:RegularTimeSeries>" << endl;
        file << subtab << tabs(8) << "<energy:variableProperties>" << endl;
        file << subtab << tabs(9) << "<energy:TimeValuesProperties>" << endl;
        file << subtab << tabs(10) << "<energy:acquisitionMethod>estimation</energy:acquisitionMethod>" << endl;
        file << subtab << tabs(10) << "<energy:interpolationType>averageInSucceedingInterval</energy:interpolationType>" << endl;
        file << subtab << tabs(10) << "<energy:source>" << zones.at(i)->getOccupantsYearProfile()->getName() << "</energy:source>" << endl;
        file << subtab << tabs(9) << "</energy:TimeValuesProperties>" << endl;
        file << subtab << tabs(8) << "</energy:variableProperties>" << endl;
        file << subtab << tabs(8) << "<energy:temporalExtent></energy:temporalExtent>" << endl;
        file << subtab << tabs(8) << "<energy:timeInterval unit=\"hour\">1</energy:timeInterval>" << endl;
        // loop on the number of time steps
        file << subtab << tabs(8) << "<energy:values uom=\"none\">";
        for (unsigned int day=1;day<=365;++day) {
            for (unsigned int hour=1;hour<=24;++hour) {
                file << zones.at(i)->getOccupantsYearProfile()->getDayProfile(day)->getHourValue(hour) << " ";
            }
        }
        file << "</energy:values>" << endl;
        file << subtab << tabs(7) << "</energy:RegularTimeSeries>" << endl;
        file << subtab << tabs(6) << "</energy:timeDependingValues>" << endl;
        file << subtab << tabs(5) << "</energy:TimeSeriesSchedule>" << endl;
        file << subtab << tabs(4) << "</energy:occupancyRate>" << endl;
        file << subtab << tabs(3) << "</energy:Occupants>" << endl;
        file << subtab << tabs(2) << "</energy:occupiedBy>" << endl;
        // appliances model
        unsigned int deviceType;
        if (zones.at(i)->getActivityType() != numeric_limits<unsigned int>::signaling_NaN()) { // activity is defined
            for(size_t j=0; j < pDistrict->getActivityType(zones.at(i)->getActivityType())->getnActivities(); ++j) {
                // loop on all devices listed in this activity
                deviceType = pDistrict->getActivityType(zones.at(i)->getActivityType())->getActivityDeviceType(j);
                for (size_t k = 0; k < pDistrict->getDeviceType(deviceType)->getnDevices(); ++k) {
                    // write the electrical appliance information
        file << subtab << tabs(2) << "<energy:has>" << endl;
        file << subtab << tabs(3) << "<energy:ElectricalAppliances>" << endl;
        file << subtab << tabs(4) << "<energy:heatDissipation>" << endl;
        file << subtab << tabs(5) << "<energy:HeatExchangeType>" << endl;
        file << subtab << tabs(6) << "<energy:convectiveFraction uom=\"none\">" << pDistrict->getDeviceType(deviceType)->getDeviceConvectiveFraction(k) << "</energy:convectiveFraction>" << endl;
        file << subtab << tabs(6) << "<energy:latentFraction uom=\"\">0</energy:latentFraction>" << endl;
        file << subtab << tabs(6) << "<energy:radiantFraction uom=\"\">" << pDistrict->getDeviceType(deviceType)->getDeviceRadiativeFraction(k) << "</energy:radiantFraction>" << endl;
        file << subtab << tabs(6) << "<energy:totalValue uom=\"\">" << pDistrict->getDeviceType(deviceType)->getDeviceAvgPower(k) << "</energy:totalValue>" << endl;
        file << subtab << tabs(5) << "</energy:HeatExchangeType>" << endl;
        file << subtab << tabs(4) << "</energy:heatDissipation>" << endl;
        file << subtab << tabs(4) << "<energy:operationSchedule>" << endl;
        file << subtab << tabs(5) << "<energy:DailyPatternSchedule>" << endl;
        file << subtab << tabs(6) << "<energy:dailySchedule>" << endl;
        file << subtab << tabs(7) << "<energy:DailySchedule>" << endl;
        file << subtab << tabs(8) << "<energy:dayType>WeekDay</energy:dayType>" << endl;
        file << subtab << tabs(8) << "<energy:schedule>" << endl;
        file << subtab << tabs(9) << "<energy:RegularTimeSeries>" << endl;
        file << subtab << tabs(10) << "<energy:temporalExtent></energy:temporalExtent>" << endl;
        file << subtab << tabs(10) << "<energy:timeInterval unit=\"hour\">1</energy:timeInterval>" << endl;
        file << subtab << tabs(10) << "<energy:values uom=\"none\">";
        for (size_t hour=1;hour<=24;++hour) {
            file << pDistrict->getActivityType(zones.at(i)->getActivityType())->getActivityProbability(j,hour)*pDistrict->getDeviceType(deviceType)->getDeviceProbability(k,hour) << " ";
        }
        file << "</energy:values>" << endl;
        file << subtab << tabs(9) << "</energy:RegularTimeSeries>" << endl;
        file << subtab << tabs(8) << "</energy:schedule>" << endl;
        file << subtab << tabs(7) << "</energy:DailySchedule>" << endl;
        file << subtab << tabs(6) << "</energy:dailySchedule>" << endl;
        file << subtab << tabs(5) << "</energy:DailyPatternSchedule>" << endl;
        file << subtab << tabs(4) << "</energy:operationSchedule>" << endl;
        file << subtab << tabs(4) << "<energy:electricalPower uom=\"W\">" << pDistrict->getDeviceType(deviceType)->getDeviceAvgPower(k) << "</energy:electricalPower>" << endl;
        file << subtab << tabs(3) << "</energy:ElectricalAppliances>" << endl;
        file << subtab << tabs(2) << "</energy:has>" << endl;
                }
            }
        }
        file << subtab << tabs(1) << "</energy:UsageZone>" << endl;
    }
    file << subtab << "</energy:usageZone>" << endl;

    // close the tag building
    file << subtab << "\t</bldg:Building>" << endl;
    file << tab << "</core:cityObjectMember>" << endl;
}

void Building::computeVolume() {

    // computes the volume in each thermal zone
    for (size_t i=0; i<zones.size();++i) zones[i]->computeVolume();

}

unsigned int Building::getnNodes()
{
    unsigned int n=0;
    for (unsigned int i=0; i<zones.size(); i++) n += zones[i]->getnNodes();
    return n;
}

unsigned int Building::getZonenNodes(unsigned int i) { return zones[i]->getnNodes(); }

unsigned int Building::getZoneIndexFromId(unsigned int id)
{

    for (unsigned int index=0; index<zones.size(); ++index) if ( zones[index]->getId() == id ) return index;
    throw string("Zone id=" + toString(id) + " not found in the zones of building id: " + toString(this->getId()));

}

double Building::getZoneC(unsigned int zoneIndex, unsigned int nodeIndex) { return zones[zoneIndex]->getC(nodeIndex); }

void Building::addZoneC(unsigned int zoneIndex, float value) { zones[zoneIndex]->addC(value); }

double Building::getZoneT(unsigned int i, unsigned int j) { return zones[i]->getT(j); }

void Building::setZoneT(unsigned int i, unsigned int j, double value) { zones[i]->setT(j,value); }

bool Building::getHasSolarThermal(){
    return hasSolarThermal;
}

void Building::computeHasSolarThermal(){
    for (unsigned int i=0; i<getnZones(); i++) {
        for (unsigned int j=0; j<getZone(i)->getnRoofs(); ++j) {
            hasSolarThermal = hasSolarThermal or getZone(i)->getRoof(j)->getHasSolarThermal();
        }
        for (unsigned int j=0; j<getZone(i)->getnWalls(); ++j) {
            hasSolarThermal = hasSolarThermal or getZone(i)->getWall(j)->getHasSolarThermal();
        }
        for (unsigned int j=0; j<getZone(i)->getnSurfaces(); ++j) {
            hasSolarThermal = hasSolarThermal or getZone(i)->getSurface(j)->getHasSolarThermal();
        }
    }
}

void Building::deterministicShadingAction(/*unsigned int day*/) {

    // the shading state depends solely on the irradiance on the facade and roof, but also on the day of year
    float irradiance = 0.f;
    // loop on all the zones
    for (size_t i=0;i<zones.size();++i) {
        // loop on all walls
        for (size_t j=0;j<zones[i]->getnWalls();++j) {
            irradiance = zones[i]->getWall(j)->getShortWaveIrradiance(); // in W/m≤
            // sets the lower blinds shading state for each wall
            /// TODO: change here to allow a strategy based on the temperature (e.g. BSol)
            //if (heatingUnit!=NULL && heatingUnit->isWorking(day))
            //    zones[i]->getWall(j)->setLowerShadingState(1.f);
            //else
                zones[i]->getWall(j)->setLowerShadingState(Model::deterministicShadingAction(irradiance,blindsIrradianceCutOff,blindsLambda));

        }
        // loop on all roofs
        for (size_t j=0;j<zones[i]->getnRoofs();++j) {
            irradiance = zones[i]->getRoof(j)->getShortWaveIrradiance(); // in W/m≤
            // sets the lower blinds shading state for each wall
            /// TODO: change here to allow a strategy based on the temperature (e.g. BSol)
            //if (heatingUnit!=NULL && heatingUnit->isWorking(day))
            //    zones[i]->getRoof(j)->setShadingState(1.f);
            //else
                zones[i]->getRoof(j)->setShadingState(Model::deterministicShadingAction(irradiance,blindsIrradianceCutOff,blindsLambda));

        }
    }

}

double Building::getMachinePower(unsigned int day, unsigned int hour) {
    return machinePower.at((day-1)*24 + hour -1 +getScene()->getPreTimeStepsSimulated());
}


bool Building::hasImposedHeatDemand(unsigned int day, unsigned int hour, float& retValue) {
    string dayHour = "d"+toString(day)+"h"+toString(hour);
    if (imposedHeatDemand.count(dayHour.c_str())==0) {
        return false;
    } else {
        retValue = imposedHeatDemand[dayHour];
        return true;
    }
}

bool Building::hasImposedHeatDemand(unsigned int day, unsigned int hour) {
    string dayHour = "d"+toString(day)+"h"+toString(hour);
    return (imposedHeatDemand.count(dayHour.c_str())!=0);
}



Tree::Tree(TiXmlHandle hdl, ostream* pLogStream):logStream(std::cout.rdbuf()) {

    // the read buffers are associated
    associate(pLogStream,logStream);

    // loading the diverse attributes
    if (hdl.ToElement()->Attribute("id")) id = to<unsigned int>(hdl.ToElement()->Attribute("id"));
    if (hdl.ToElement()->Attribute("name")) key = hdl.ToElement()->Attribute("name");
    if (hdl.ToElement()->Attribute("key")) key = hdl.ToElement()->Attribute("key");
    if (hdl.ToElement()->Attribute("leafAreaIndex")) layers = to<unsigned int>(hdl.ToElement()->Attribute("leafAreaIndex"));
    if (hdl.ToElement()->Attribute("leafWidth")) leafWidth = to<float>(hdl.ToElement()->Attribute("leafWidth"));
    if (hdl.ToElement()->Attribute("leafDistance")) layersDistance = to<float>(hdl.ToElement()->Attribute("leafDistance"));
    if (hdl.ToElement()->Attribute("deciduous")) deciduous = (hdl.ToElement()->Attribute("leafDistance") == string("true"));
    if (hdl.ToElement()->Attribute("class")) leafClass = hdl.ToElement()->Attribute("class");

    // read the leaves and create subleaves
    unsigned int surfaceIndex = 0;
    while (hdl.ChildElement("Leaf", surfaceIndex).ToElement()) {
        leaves.push_back(new Surface(hdl.ChildElement("Leaf", surfaceIndex++), NULL, &logStream));
        if (leaves.back()->getArea() <= 0.f) {
            logStream << "Leaf id: " << leaves.back()->getId() << " has a too small surface, removing it." << endl;
            delete leaves.back();
            leaves.pop_back();
        }
        else {
            // create the subleaves
            subLeaves.push_back(new Surface(*leaves.back()));
            subLeaves.back()->reverseOrientation();
            for (unsigned int index = 1; index < layers; ++index) {
                subLeaves.push_back(new Surface(*leaves.back()));
                subLeaves.back()->translate(GENPoint::Cartesian(0.f,0.f,-static_cast<float>(index)*layersDistance));
                subLeaves.push_back(new Surface(*subLeaves.back()));
                subLeaves.back()->reverseOrientation();
            }
        }
    }

    // read the trunc
    surfaceIndex = 0;
    while (hdl.ChildElement("Trunc", surfaceIndex).ToElement()) {
        trunc.push_back(new Surface(hdl.ChildElement("Trunc", surfaceIndex++), NULL, &logStream));
        if (trunc.back()->getArea() <= 0.f) {
            logStream << "Trunc id: " << trunc.back()->getId() << " has a too small surface, removing it." << endl;
            delete trunc.back();
            trunc.pop_back();
        }
    }

    logStream << "Tree with key: " << key << endl;

}

void Tree::writeXML(ofstream& file, string tab){
    file << tab << "<Tree id=\"" << id << "\" name=\"" << name << "\" key=\"" << key << "\" leafAreaIndex=\"" << layers << "\" leafWidth=\"" << leafWidth << "\" leafDistance=\"" << layersDistance << "\" deciduous=\"" << (deciduous?"true":"false") << "\" class=\"" << leafClass << "\">" << endl;
    string subtab=tab+"\t";
    // leaves
    for (size_t i=0; i < leaves.size(); ++i) {
        file << subtab << "<Leaf id=\"" << leaves.at(i)->getId() << "\" ShortWaveReflectance=\"" << leaves.at(i)->getShortWaveReflectance() << "\" LongWaveEmissivity=\"" << leaves.at(i)->getLongWaveEmissivity() << "\">" << endl;
        leaves.at(i)->writeXML(file,subtab+"\t");
        file << subtab << "</Leaf>" << endl;
    }
    // trunc
    for (size_t i=0; i < trunc.size(); ++i) {
        file << subtab << "<Trunc id=\"" << trunc.at(i)->getId() << "\" ShortWaveReflectance=\"" << trunc.at(i)->getShortWaveReflectance() << "\" LongWaveEmissivity=\"" << trunc.at(i)->getLongWaveEmissivity() << "\">" << endl;
        trunc.at(i)->writeXML(file,subtab+"\t");
        file << subtab << "</Trunc>" << endl;
    }
    file << tab << "</Tree>" << endl;
}
