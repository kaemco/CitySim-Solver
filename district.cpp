#include "district.h"
#include "scene.h"

// *** District class, CitySim *** //
// *** jerome.kaempf@epfl.ch   *** //

District::District(TiXmlHandle XMLHandler, XmlScene* pScene):pScene(pScene),occupancyProfiles(&(pScene->logStream)),dhwProfiles(&(pScene->logStream)),logStream(pScene->logStream.rdbuf()) {

    try { // Cognet: Added this to delete dynamically allocated objects if exception is thrown in the constructor.

        if (!XMLHandler.FirstChild("District").ToElement()) throw(string("Error in XML file: no District tag."));

        unsigned int i=0;//,j=0;
        TiXmlElement *elem;

        // Far Field Obstructions loading in the District
        logStream << "Loading Far Field Obstruction profile: ";
        if (XMLHandler.FirstChild("District").FirstChildElement("FarFieldObstructions").ToElement()) {
            // the element FarFieldObstructions exists
            while (XMLHandler.FirstChild("District").FirstChild("FarFieldObstructions").ChildElement("Point",i).ToElement()) {
                farFieldObstructions.push_back(pair<float,float>(
                to<float>(XMLHandler.FirstChild("District").FirstChild("FarFieldObstructions").ChildElement("Point",i).ToElement()->Attribute("phi")),
                to<float>(XMLHandler.FirstChild("District").FirstChild("FarFieldObstructions").ChildElement("Point",i).ToElement()->Attribute("theta"))));
                ++i;
            }
            logStream << farFieldObstructions.size() << " loaded." << endl << flush;
        }
        else logStream << "no profile given." << endl << flush;

        // Composites loading in the District, they can be named: Composite, WallType and GroundType
        logStream << "Loading composite(s):\n";
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("Composite",i).ToElement()) {
            composites.insert( pair<string,Composite*>(XMLHandler.FirstChild("District").ChildElement("Composite",i).ToElement()->Attribute("id"), new Composite(XMLHandler.FirstChild("District").ChildElement("Composite",i),&logStream)));
            ++i;
        }
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("WallType",i).ToElement()) {
            composites.insert( pair<string,Composite*>(XMLHandler.FirstChild("District").ChildElement("WallType",i).ToElement()->Attribute("id"), new Composite(XMLHandler.FirstChild("District").ChildElement("WallType",i),&logStream)));
            composites[XMLHandler.FirstChild("District").ChildElement("WallType",i).ToElement()->Attribute("id")]->setCategory(Surface::STypeNames[Surface::WALL]);
            ++i;
        }
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("GroundType",i).ToElement()) {
            composites.insert( pair<string,Composite*>(XMLHandler.FirstChild("District").ChildElement("GroundType",i).ToElement()->Attribute("id"), new Composite(XMLHandler.FirstChild("District").ChildElement("GroundType",i),&logStream)));
            composites[XMLHandler.FirstChild("District").ChildElement("GroundType",i).ToElement()->Attribute("id")]->setCategory(Surface::STypeNames[Surface::GROUND]);
            ++i;
        }
        logStream << composites.size() << " loaded." << endl << flush;

        // Occupancy profiles loading in the District
        occupancyProfiles.clearIdMap();
        occupancyProfiles.readFromXml(XMLHandler.FirstChild("District"));
        // Domestic Hot Water (DHW) profiles loading in the District
        dhwProfiles.clearIdMap();
        dhwProfiles.readFromXml(XMLHandler.FirstChild("District"));
        // Temperature profiles loading in the District
        temperatureProfiles = new TemperatureProfiles(XMLHandler.FirstChild("District"));

        // Device type profiles loading in the District
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("DeviceType",i).ToElement()) {
            deviceTypes.push_back(new DeviceType(XMLHandler.FirstChild("District").ChildElement("DeviceType",i++), &logStream));
        }
        // Activity type profiles loading in the District
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("ActivityType",i).ToElement()) {
            activityTypes.push_back(new ActivityType(XMLHandler.FirstChild("District").ChildElement("ActivityType",i++), &logStream));
        }

        // Cognet: Start of code moved to here, to initialize the DistrictEnergyCenters before the Substations.
        logStream << "Loading district energy center(s) : " << endl << flush;
        i=0;
        if (XMLHandler.FirstChild("District").ChildElement("DistrictEnergyCenter",i).ToElement()) {
            logStream << "DistrictEnergyCenter tag: " << endl << flush;
        }
        else {
            logStream << "There are no DistrictEnergyCenters." << endl << flush;
        }

        while (XMLHandler.FirstChild("District").ChildElement("DistrictEnergyCenter",i).ToElement() ){
            districtEnergyCenters.push_back(new DistrictEnergyCenter(XMLHandler.FirstChild("District").ChildElement("DistrictEnergyCenter", i), this, &(this->logStream)) ); // Cognet: Added logStream.
            ++i;
        }
        // Cognet: End of moved code.


        // loading the geometry of the buildings
        logStream << "Loading geometry of building(s)." << endl << flush;
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("Building",i).ToElement()) {

            // check if the building should be simulated or not (in the building tag itself)
            if (string(XMLHandler.FirstChild("District").ChildElement("Building",i).ToElement()->Attribute("Simulate"))==string("false")) {

                // if simulate=false, then easy, we just add surfaces to the ground for the SW calculation only
                logStream << "Building: " << string(XMLHandler.FirstChild("District").ChildElement("Building",i).ToElement()->Attribute("id")) << "\twith key: " << string(XMLHandler.FirstChild("District").ChildElement("Building",i).ToElement()->Attribute("key"));
                logStream << "\tnot thermally simulated." << endl << flush;

                // loop on the zones
                TiXmlElement *zoneElem = XMLHandler.FirstChild("District").ChildElement("Building",i).FirstChildElement("Zone").ToElement();
                while (zoneElem) { // zoning exists

                    // loading the wall vertices and put them in a vector
                    TiXmlHandle zoneHdl = zoneElem;
                    elem = zoneHdl.FirstChildElement("Wall").ToElement();
                    while(elem){

                        // reads the vertices and creates a shading, with an id
                        surfaces.push_back(new Surface(elem,nullptr,&logStream));

                        // checks the size of the surface
                        if ((surfaces.back()->getArea() <= 0.f) || (surfaces.back()->getRadius() <= 0.f)) {
                            logStream << "(Warning) Wall id=" << surfaces.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                            delete surfaces.back();
                            surfaces.pop_back();
                        }
                        elem = elem->NextSiblingElement("Wall");
                    }

                    // loading the roof vertices
                    elem = zoneHdl.FirstChildElement("Roof").ToElement();
                    while(elem){

                        // reads the vertices and creates a shading, with an id
                        surfaces.push_back(new Surface(elem,nullptr,&logStream));

                        // checks the size of the surface
                        if ((surfaces.back()->getArea() <= 0.f) || (surfaces.back()->getRadius() <= 0.f)) {
                            logStream << "(Warning) Roof id=" << surfaces.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                            delete surfaces.back();
                            surfaces.pop_back();
                        }
                        elem = elem->NextSiblingElement("Roof");
                    }

                    // loading the floor vertices
                    elem = zoneHdl.FirstChildElement("Floor").ToElement();
                    while(elem){

                        // reads the vertices and creates a shading, with an id
                        surfaces.push_back(new Surface(elem,nullptr,&logStream));

                        // checks the size of the surface
                        if ((surfaces.back()->getArea() <= 0.f) || (surfaces.back()->getRadius() <= 0.f)) {
                            logStream << "(Warning) Roof id=" << surfaces.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                            delete surfaces.back();
                            surfaces.pop_back();
                        }
                        elem = elem->NextSiblingElement("Floor");
                    }

                    // loading the surface vertices and put them in a vector (obstructing surfaces)
                    elem = zoneHdl.FirstChildElement("Surface").ToElement();
                    while(elem){

                        surfaces.push_back(new Surface(elem,nullptr,&logStream));

                        // checks the size of the surface
                        if ((surfaces.back()->getArea() <= 0.f) || (surfaces.back()->getRadius() <= 0.f)) {
                            logStream << "(Warning) Surface id=" << surfaces.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                            delete surfaces.back();
                            surfaces.pop_back();
                        }
                        elem = elem->NextSiblingElement("Surface");
                    }

                    // go to the next zone
                    zoneElem = zoneElem->NextSiblingElement("Zone");
                }
            }
            else {
                // the building should be simulated, we then add it to the vector of buildings
                buildings.push_back(new Building(XMLHandler.FirstChild("District").ChildElement("Building",i),this));
            }
            ++i;
        }
        logStream << buildings.size() << " loaded." << endl << flush;

        for(auto dec : districtEnergyCenters) { dec->substationsConstructionFinished(); } // To check that substations have all been found, find loops and regulated paths...


        logStream << "Loading wind turbine(s): ";
        i=0;
        while (XMLHandler.FirstChild("District").ChildElement("WindTurbine",i).ToElement()) {
            /// TODO: load the data of the WindTurbines
            /// save the output of the model
            //windTurbines.push_back(new WindTurbine());
        }
        logStream << windTurbines.size() << " loaded." << endl;

        // loads the Shading surfaces
        logStream << "Loading shading surface(s): " << flush;

        elem = XMLHandler.FirstChild("District").FirstChild("ShadingSurface").FirstChild("Surface").ToElement();
        while (elem) {

            // reads the vertices and creates a shading, with an id
            surfaces.push_back(new Surface(elem,NULL,&logStream));

            // checks the size of the surface
            if ((surfaces.back()->getArea() <= 0.f) || (surfaces.back()->getRadius() <= 0.f)) {
                logStream << "(Warning) Shading id=" << surfaces.back()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                delete surfaces.back();
                surfaces.pop_back();
            }
            elem = elem->NextSiblingElement();
        }
        logStream << surfaces.size() << " loaded." << endl;

        // loads the Trees
        logStream << "Loading tree(s): ";
        elem = XMLHandler.FirstChild("District").FirstChild("Trees").FirstChild("Tree").ToElement();
        while (elem) {
            // pushes back the whole tree
            trees.push_back(new Tree(elem));
            elem = elem->NextSiblingElement();
        }
        logStream << trees.size() << " loaded." << endl;

        // loads the Ground surfaces
        logStream << "Loading ground surface(s): " << flush;
        if (XMLHandler.FirstChild("District").FirstChild("GroundSurface").ToElement())
            if (XMLHandler.FirstChild("District").FirstChild("GroundSurface").ToElement()->Attribute("ShortWaveReflectance"))
                groundAlbedo = to<float>(XMLHandler.FirstChild("District").FirstChild("GroundSurface").ToElement()->Attribute("ShortWaveReflectance"));

        elem = XMLHandler.FirstChild("District").FirstChild("GroundSurface").FirstChild("Ground").ToElement();
        while (elem) {
            Composite *groundType = NULL;
            if (elem->Attribute("type")) {
                groundType = getComposite(elem->Attribute("type"));
            }

            // reads the vertices and creates a ground, with an id
            grounds.push_front(new Ground(elem,groundType,&logStream));

            // checks the size of the Ground surface
            if ((grounds.front()->getArea() <= 0.f) || (grounds.front()->getRadius() <= 0.f)) {
                logStream << "(Warning) Ground id=" << grounds.front()->getId() << " has a too small surface area or radius, removing it" << endl << flush;
                delete grounds.front();
                grounds.pop_front();
            }
            elem = elem->NextSiblingElement();
        }
        logStream << getnGrounds() << " loaded." << endl;

        logStream << "District created." << endl;
    } catch(...) { // Cognet: Added this to delete dynamically allocated objects if exception is thrown in the constructor.
        deleteDynamicallyAllocated();
        throw;
    }

}

District::District(XmlScene* pScene):pScene(pScene),occupancyProfiles(&(pScene->logStream)),logStream(pScene->logStream.rdbuf()){}

void District::addBuilding(Building *b) {

    // adds the building in the building vector
    buildings.push_back(b);

    // adds the building in the DATARadiationScene
    for (unsigned int j=0; j<buildings.back()->getnZones(); ++j) { // loop on all zones in the building
        // loop for the walls on this zone
        for (unsigned int k=0; k<buildings.back()->getZone(j)->getnWalls(); ++k) {
            // add the surface to the Buildings surfaces (daylight calculation)
            getScene()->getDATARadiationScene()->AddBuildingSurface(GENHandle<Wall>(buildings.back()->getZone(j)->getWall(k)));
        }
        // loop for the roofs on this zone
        for (unsigned int k=0; k<buildings.back()->getZone(j)->getnRoofs(); ++k) {
            // add the surface to the Ground surfaces (meaning NO daylight calculation)
            getScene()->getDATARadiationScene()->AddGroundSurface(GENHandle<Roof>(buildings.back()->getZone(j)->getRoof(k)));
        }
        // loop on the obstructing surfaces on this zone
        for (unsigned int k=0; k<buildings.back()->getZone(j)->getnSurfaces(); ++k) {
            // add the surface to the Ground surfaces (meaning NO daylight calculation)
            getScene()->getDATARadiationScene()->AddGroundSurface(GENHandle<Surface>(buildings.back()->getZone(j)->getSurface(k)));
        }
    }

    return;
}

void District::addTree(Tree *t) {

    // adds the Tree to the vector
    trees.push_back(t);

    // adds the tree surfaces to the DATARadiationScene
    for (size_t j=0; j<trees.back()->getnSurfaces(); ++j)
        getScene()->getDATARadiationScene()->AddGroundSurface(GENHandle<Surface>(trees.back()->getSurface(j)));

    return;
}

void District::addSurface(Surface *s) {

    // adds the surface to the vector
    surfaces.push_back(s);

    // adds the ground surfaces to the DATARadiationScene
    getScene()->getDATARadiationScene()->AddGroundSurface(GENHandle<Surface>(surfaces.back()));

    return;
}

void District::addGround(Ground *g) {

    // adds the ground to the vector
    grounds.push_front(g);

    // adds the ground surfaces to the DATARadiationScene
    getScene()->getDATARadiationScene()->AddGroundSurface(GENHandle<Ground>(grounds.front()));

    return;
}

void District::writeXML(ofstream& file, string tab){
    file << tab << "<District>" << endl;
    string subtab = tab+"\t";
    file << subtab << "<FarFieldObstructions>" << endl;
    for (unsigned int i=0; i<farFieldObstructions.size();++i){
        file << subtab << "\t<Point phi=\"" << farFieldObstructions[i].first << "\"  theta=\"" << farFieldObstructions[i].second << "\"/>" << endl;
    }
    file << subtab << "</FarFieldObstructions>" << endl;

    // Get the list of all used YearProfile and Composite
    set<YearProfile*> usedProfiles;
    set<YearProfile*> usedDHWProfiles;
    set<Composite*> usedComposites;

    // shading surfaces
    for (vector<Surface*>::iterator itSurfaces=getSurfaces()->begin(); itSurfaces != getSurfaces()->end(); ++itSurfaces) {
        if((*itSurfaces)->getComposite()!=nullptr) usedComposites.insert((*itSurfaces)->getComposite());
    }
    // ground surfaces
    for (forward_list<Ground*>::iterator itSurfaces=getGrounds()->begin(); itSurfaces != getGrounds()->end(); ++itSurfaces) {
        if((*itSurfaces)->getComposite()!=nullptr) usedComposites.insert((*itSurfaces)->getComposite());
    }
    // Buildings
    for (vector<Building*>::iterator itBuildings = getBuildings()->begin(); itBuildings != getBuildings()->end(); ++itBuildings) {
        for (vector<Zone*>::iterator itZones = (*itBuildings)->getZones()->begin(); itZones != (*itBuildings)->getZones()->end(); ++itZones) {
            if((*itZones)->getOccupantsYearProfile()!=nullptr) usedProfiles.insert((*itZones)->getOccupantsYearProfile());
            if((*itZones)->getDHWYearProfile()!=nullptr) usedDHWProfiles.insert((*itZones)->getDHWYearProfile());
            // walls
            for (vector<Wall*>::iterator itSurfaces = (*itZones)->getWalls()->begin(); itSurfaces != (*itZones)->getWalls()->end(); ++itSurfaces) {
                if((*itSurfaces)->getComposite()!=nullptr) usedComposites.insert((*itSurfaces)->getComposite());
            }
            // roofs
            for (vector<Roof*>::iterator itSurfaces = (*itZones)->getRoofs()->begin(); itSurfaces != (*itZones)->getRoofs()->end(); ++itSurfaces) {
                if((*itSurfaces)->getComposite()!=nullptr) usedComposites.insert((*itSurfaces)->getComposite());
            }
            // floors
            for (vector<Floor*>::iterator itSurfaces = (*itZones)->getFloors()->begin(); itSurfaces != (*itZones)->getFloors()->end(); ++itSurfaces) {
                if((*itSurfaces)->getComposite()!=nullptr) usedComposites.insert((*itSurfaces)->getComposite());
            }
            // surfaces
            for (vector<Surface*>::iterator itSurfaces = (*itZones)->getSurfaces()->begin(); itSurfaces != (*itZones)->getSurfaces()->end(); ++itSurfaces) {
                if((*itSurfaces)->getComposite()!=nullptr) usedComposites.insert((*itSurfaces)->getComposite());
            }
        }
    }

    // Write Composites
    // Make sure all Composites used in this scene have different ids before writing xml
    set<unsigned int> usedCompositeIds;
    unsigned int newId=0;
    for (set<Composite*>::const_iterator iter=usedComposites.begin(); iter != usedComposites.end(); ++ iter){
        while(usedCompositeIds.count(newId)!=0){
            ++newId;
        }
        if(usedCompositeIds.count((*iter)->getId())!=0){
            (*iter)->setId(newId);
        }
        usedCompositeIds.insert((*iter)->getId());
        (*iter)->writeXML(file, subtab);
    }
    file << endl;

    // Write OccupancyProfiles
    occupancyProfiles.writeXML(file, usedProfiles, subtab);
    // write DHW profiles
    dhwProfiles.writeXML(file, usedDHWProfiles, subtab);
    // write TemperatureProfiles
    if (temperatureProfiles) temperatureProfiles->writeXML(file, subtab);

    // Device type profiles
    for (vector<DeviceType*>::iterator it=deviceTypes.begin(); it!=deviceTypes.end();++it)
        (*it)->writeXML(file, subtab);
    // Activity type profiles
    for (vector<ActivityType*>::iterator it=activityTypes.begin(); it!=activityTypes.end();++it)
        (*it)->writeXML(file, subtab);

    // Write Buildings
    for (unsigned int i=0; i<buildings.size(); ++ i){
        buildings[i]->writeXML(file, subtab);
    }

    // Write Shading
    file << subtab << "<ShadingSurface>" << endl;
    for (size_t i=0; i<surfaces.size(); ++i){
        file << subtab << "\t<Surface id=\"" << surfaces[i]->getId() << "\" ";
        // write the key if it exists
        if (!surfaces[i]->getKey().empty()) file << "key=\"" << surfaces[i]->getKey() << "\" ";
        file << "ShortWaveReflectance=\"" << surfaces[i]->getShortWaveReflectance();
        file << "\">" << endl;
        surfaces[i]->writeXML(file,subtab+"\t\t");
        file << subtab << "\t</Surface>" << endl;
    }
    file << subtab << "</ShadingSurface>" << endl;

    // write Trees
    file << subtab << "<Trees>" << endl;
    for (size_t i=0; i<trees.size(); ++i){
        trees.at(i)->writeXML(file,subtab+"\t");
    }
    file << subtab << "</Trees>" << endl;

    // Write Ground
    file << subtab << "<GroundSurface ShortWaveReflectance=\"" << groundAlbedo << "\">" << endl;
    for (forward_list<Ground*>::iterator it=grounds.begin(); it!=grounds.end(); ++it){
        (*it)->writeXML(file,subtab+"\t");
    }
    file << subtab << "</GroundSurface>" << endl;

    file << tab << "</District>" << endl;
}

void District::writeGML(ofstream& file, string tab, const vector<double>& origin) {

    // write the composites
    for (map<string,Composite*>::iterator it=composites.begin();it!=composites.end();++it) {
        file << tab << "<gml:featureMember>" << endl;
        it->second->writeGML(file, tab+"\t");
        file << tab << "</gml:featureMember>" << endl;
    }

    // write the buildings
    for (size_t i=0; i<buildings.size(); ++i) {

        buildings[i]->writeGML(file, tab, origin);

    }
}

District::~District() {

    //logStream << "Destructor of District." << endl << flush;
    for (vector<Building*>::iterator it=buildings.begin();it!=buildings.end();++it) delete *it;
    //NB: the deletion of Ground pointers is done by GENHandle when destroying the scene
    // DP : not sure what this means... but deleting here create seg fault.
    //for (vector<Ground*>::iterator it=grounds.begin();it!=grounds.end();it++) delete *it;
    for (map<string,Composite*>::iterator it=composites.begin();it!=composites.end();++it) delete it->second;
    for (vector<EnergyConversionSystem*>::iterator it=plants.begin();it!=plants.end();++it) delete *it;
    //beginning of contents added by Dapeng
    for (vector<DistrictEnergyCenter*>::iterator it = districtEnergyCenters.begin(); it!=districtEnergyCenters.end(); ++it) delete *it;
    //ending of contents added by Dapeng
    if (temperatureProfiles) delete temperatureProfiles;

}

void District::deleteDynamicallyAllocated() { // Cognet: Added this, to delete objects that have been created if an exception is thrown from the constructor.
    for (auto pair : composites) { delete pair.second; } composites.clear();
    while(!deviceTypes.empty()) { delete deviceTypes.back(); deviceTypes.pop_back(); }
    while(!activityTypes.empty()) { delete activityTypes.back(); activityTypes.pop_back(); }
    while(!districtEnergyCenters.empty()) { delete districtEnergyCenters.back(); districtEnergyCenters.pop_back(); }
    while(!buildings.empty()) { delete buildings.back(); buildings.pop_back(); }
    while(!surfaces.empty()) { delete surfaces.back(); surfaces.pop_back(); }
    while(!trees.empty()) { delete trees.back(); trees.pop_back(); }
    while(!grounds.empty()) { delete grounds.front(); grounds.pop_front(); }
    if (temperatureProfiles) delete temperatureProfiles;
}

void District::readFarField(string fileName){

    // If no fileName to load -> go back directly
    if (fileName.empty()) {
        logStream << "WARNING: Horizon file name empty." << endl << flush;
        return;
    }

    farFieldObstructions.clear();

    // variables used to store the elements
    //char cbuffer[200];
    string buffer;

    // Climate filename opening
    fstream input(fileName.c_str(), ios::in | ios::binary);
    if (!input.is_open()) throw(string("Error opening climate file: " + fileName));
    //logStream << "Loading: " << fileName << endl << flush;

    float azimut;
    float elevation;

    input >> buffer; // azimut

    while ( !input.eof() ) {

        azimut = atof(buffer.c_str());
        if (azimut < 0) azimut += 360;

        input >> buffer; // elevation
        elevation = atof(buffer.c_str());

        //cout << "azimut=" << azimut << ", elevation="<< elevation << endl;
        farFieldObstructions.push_back(pair<float,float>(azimut, elevation));

        input >> buffer; // next azimut, if not end of file
    }
}
