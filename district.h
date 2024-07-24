#ifndef DISTRICT_H
#define DISTRICT_H

#include <map>
#include <forward_list>

#include "building.h"
#include "surface.h"

//beginning of contents added by Dapeng
#include "plant.h"
//end of contents added by Dapeng

using namespace std;

// *** District class, CitySim *** //
// *** jerome.kaempf@epfl.ch   *** //

class District {

private:

    // parent
    XmlScene* pScene;

    // elements in the District
    vector<Building*> buildings;
    vector<Surface*> surfaces; //!< obstructing surfaces
    vector<Tree*> trees; //!< trees in the district

    float groundAlbedo = 0.2f; //!< albedo of the lower hemisphere
    forward_list<Ground*> grounds; //!< ground surfaces
    vector<pair<float,float> > farFieldObstructions;
    map<string,Composite*> composites;
    vector<EnergyConversionSystem*> plants;
    vector<WindTurbine> windTurbines;
    OccupancyProfiles occupancyProfiles;
    DHWProfiles dhwProfiles;
    TemperatureProfiles* temperatureProfiles = nullptr;
    vector<DeviceType*> deviceTypes;
    vector<ActivityType*> activityTypes;

    vector<DistrictEnergyCenter*> districtEnergyCenters;

public:

    ostream logStream;

    District(TiXmlHandle XMLHandler, XmlScene* pScene);
    District(XmlScene* pScene); // Incomplete constructor for DXF reading, do not use without completing:
    ~District();
    void deleteDynamicallyAllocated(); // Cognet: Added this.
    void clear() {
        for (vector<Building*>::iterator it=buildings.begin();it!=buildings.end();++it) (*it)->clear();
        for (vector<Surface*>::iterator it=surfaces.begin();it!=surfaces.end();++it) (*it)->clear();
        for (vector<Tree*>::iterator it=trees.begin();it!=trees.end();++it) (*it)->clear();
        for (forward_list<Ground*>::iterator it=grounds.begin();it!=grounds.end();++it) (*it)->clear();
    }

    void writeXML(ofstream& file, string tab="");
    void writeGML(ofstream& file, string tab="", const vector<double>& origin={0.,0.,0.});

    void readFarField(string fileName);

    // gets the parent
    XmlScene* getScene() { return pScene; }

    size_t getnBuildings() { return buildings.size(); }
    Building* getBuilding(unsigned int i) { return buildings[i]; }
    vector<Building*>* getBuildings() { return &buildings; }
    void addBuilding(Building* b);
    vector<Tree*>* getTrees() { return &trees; }
    size_t getnTrees() { return trees.size(); }
    Tree* getTree(size_t i) { return trees.at(i); }
    void addTree(Tree* t);

    OccupancyProfiles* getOccupancyProfiles() { return &occupancyProfiles; }
    DHWProfiles* getDHWProfiles() { return &dhwProfiles; }
    size_t getnActivityTypes() { return activityTypes.size(); }
    ActivityType* getActivityType(unsigned int id) { for (size_t i=0; i<activityTypes.size();++i) { if (activityTypes.at(i)->getId()==id) return activityTypes.at(i); } throw(string("ActivityType id not found.")); }
    size_t getnDeviceTypes() { return deviceTypes.size(); }
    DeviceType* getDeviceType(unsigned int id) { for (size_t i=0; i<deviceTypes.size();++i) { if (deviceTypes.at(i)->getId()==id) return deviceTypes.at(i); } throw(string("DeviceType id not found.")); }
    TemperatureProfiles* getTemperatureProfiles() { return temperatureProfiles; }

    unsigned int getnSurfaces() { return surfaces.size(); }
    Surface* getSurface(unsigned int i) { return surfaces.at(i); }
    vector<Surface*>* getSurfaces() { return &surfaces;}
    Surface* getSurface() { if (surfaces.empty()) return NULL; else return surfaces.back(); }
    void addSurface(Surface* s);

    float getGroundAlbedo() { return groundAlbedo; }

    unsigned int getnGrounds() { return distance(grounds.begin(),grounds.end()); }
    forward_list<Ground*>* getGrounds() { return &grounds;}
    Ground* getGround() { if (grounds.empty()) return NULL; else return grounds.front(); }
    void addGround(Ground* g);
    unsigned int getLargestGroundId() { if (grounds.empty()) return 0; else return (*max_element(grounds.begin(),grounds.end()))->getId(); }

    unsigned int getnFarFieldObstructions() { return farFieldObstructions.size(); }
    vector<pair<float,float> > getFarFieldObstructions() { return farFieldObstructions; }
    pair<float,float> getFarFieldObstructions(unsigned int i) { return farFieldObstructions[i]; }
    void addFarFieldObstructions(float phi, float theta) { farFieldObstructions.push_back(pair<float,float>(phi,theta)); }
    void clearFarFieldObstructions() { farFieldObstructions.clear(); }

    Composite* getComposite(string type) { map<string,Composite*>::iterator it = composites.find(type); if (it!=composites.end()) return it->second; else throw(string("Composite: ") + type + string(" not found in XML file.")); }
    int getMaxWallTypeId() {
        int maxId=0;
        for (map<string,Composite*>::iterator it = composites.begin(); it!=composites.end(); ++it) {
            if (it->second->getId() > maxId)
                maxId=it->second->getId();
        }
        return maxId;
    }
    Composite* getUvalueComposite(float Uvalue) {
        int maxId=0;
        string base = "Simple U-value="+toString(Uvalue)+" composite";
        string name = base;
        for (map<string,Composite*>::iterator it = composites.begin(); it!=composites.end(); ++it) {
            if (it->second->getnLayers()==0 && it->second->getUvalue()==Uvalue)
                return it->second;
            if (it->second->getId() > maxId)
                maxId=it->second->getId();
        }
        ++maxId;
        if(composites.count(name)==1){
            name = base + " "+ toString(maxId);
        }
        while(composites.count(name)==1){
            ++maxId;
            name = base + " "+ toString(maxId);
        }
        composites.insert(pair<string,Composite*>(name, new Composite(Uvalue,maxId,name,&logStream)));
        return getComposite(name);
    }

    map<string,Composite*>* getComposites(){
        return &composites;
    }

    //beginning of contents added by Dapeng
    unsigned int getnDECs() { return districtEnergyCenters.size(); }
    DistrictEnergyCenter* getDEC(unsigned int i) {return districtEnergyCenters[i]; }
    //end of contents added by Dapeng

};

#endif
