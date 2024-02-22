#ifndef OCCUPANTS
#define OCCUPANTS

#include <map>
#include <set>
#include "models.h"
#include "tinyxml.h"
#include "util.h"

// *** Occupants class, CitySim  *** //
// *** jerome.kaempf@epfl.ch     *** //

class StochasticPresenceParameters {

private:

    vector<float> pMon, pTue, pWed, pThu, pFri, pSat, pSun;

public:
    StochasticPresenceParameters(vector<float> pMon, vector<float> pTue, vector<float> pWed, vector<float> pThu, vector<float> pFri, vector<float> pSat, vector<float> pSun)
        : pMon(pMon), pTue(pTue), pWed(pWed), pThu(pThu), pFri(pFri), pSat(pSat), pSun(pSun) {}

    float getpMon(unsigned int index) { return pMon.at(index); }
    float getpTue(unsigned int index) { return pTue.at(index); }
    float getpWed(unsigned int index) { return pWed.at(index); }
    float getpThu(unsigned int index) { return pThu.at(index); }
    float getpFri(unsigned int index) { return pFri.at(index); }
    float getpSat(unsigned int index) { return pSat.at(index); }
    float getpSun(unsigned int index) { return pSun.at(index); }

};

class StochasticWindowParameters {

private:

    string windowModel;
    vector<float> P01arr, P10arr, P01int, P10int, P01dep, P10dep;
    vector<float> duropen;

public:
    StochasticWindowParameters(string windowModel, vector<float> P01arr, vector<float> P10arr, vector<float> P01int, vector<float> P10int, vector<float> P01dep, vector<float> P10dep, vector<float> duropen)
        :windowModel(windowModel), P01arr(P01arr), P10arr(P10arr), P01int(P01int), P10int(P10int), P01dep(P01dep), P10dep(P10dep), duropen(duropen) {}

    string getWindowModel() { return windowModel; }

    float getP01arr(unsigned int index) { return P01arr.at(index); }
    float getP10arr(unsigned int index) { return P10arr.at(index); }
    float getP01int(unsigned int index) { return P01int.at(index); }
    float getP10int(unsigned int index) { return P10int.at(index); }
    float getP01dep(unsigned int index) { return P01dep.at(index); }
    float getP10dep(unsigned int index) { return P10dep.at(index); }

    float getDuropen(unsigned int index) { return duropen.at(index); }

};

class StochasticBlindsParameters {

private:
    vector<float> Plowerarr_UB, Praisearr_UB, Plowerint_UB, Praiseint_UB, Pfulllower_UB, Pfullraise_UB;
    vector<float> Plowerarr_LB, Praisearr_LB, Plowerint_LB, Praiseint_LB, Pfulllower_LB, Pfullraise_LB, DistFrac_LB;

public:
    StochasticBlindsParameters(vector<float> Plowerarr_UB, vector<float> Praisearr_UB, vector<float> Plowerint_UB, vector<float> Praiseint_UB, vector<float> Pfulllower_UB, vector<float> Pfullraise_UB,
        vector<float> Plowerarr_LB, vector<float> Praisearr_LB,  vector<float> Plowerint_LB, vector<float> Praiseint_LB, vector<float> Pfulllower_LB, vector<float> Pfullraise_LB, vector<float> DistFrac_LB)
        : Plowerarr_UB(Plowerarr_UB), Praisearr_UB(Praisearr_UB), Plowerint_UB(Plowerint_UB), Praiseint_UB(Praiseint_UB), Pfulllower_UB(Pfulllower_UB), Pfullraise_UB(Pfullraise_UB),
          Plowerarr_LB(Plowerarr_LB), Praisearr_LB(Praisearr_LB), Plowerint_LB(Plowerint_LB), Praiseint_LB(Praiseint_LB), Pfulllower_LB(Pfulllower_LB), Pfullraise_LB(Pfullraise_LB), DistFrac_LB(DistFrac_LB) {}

    // upper blinds model
    float getPlowerarr_UB(unsigned int index) { return Plowerarr_UB.at(index); }
    float getPraisearr_UB(unsigned int index) { return Praisearr_UB.at(index); }
    float getPlowerint_UB(unsigned int index) { return Plowerint_UB.at(index); }
    float getPraiseint_UB(unsigned int index) { return Praiseint_UB.at(index); }
    float getPfulllower_UB(unsigned int index) { return Pfulllower_UB.at(index); }
    float getPfullraise_UB(unsigned int index) { return Pfullraise_UB.at(index); }

    // lower blinds model
    float getPlowerarr_LB(unsigned int index) { return Plowerarr_LB.at(index); }
    float getPraisearr_LB(unsigned int index) { return Praisearr_LB.at(index); }
    float getPlowerint_LB(unsigned int index) { return Plowerint_LB.at(index); }
    float getPraiseint_LB(unsigned int index) { return Praiseint_LB.at(index); }
    float getPfulllower_LB(unsigned int index) { return Pfulllower_LB.at(index); }
    float getPfullraise_LB(unsigned int index) { return Pfullraise_LB.at(index); }
    float getDistFrac_LB(unsigned int index) { return DistFrac_LB.at(index); }

};

class StochasticLightsParameters {

private:

    vector<float> Ponarr, Ponint;

public:
    StochasticLightsParameters(vector<float> Ponarr, vector<float> Ponint)
        :Ponarr(Ponarr), Ponint(Ponint) {}

    float getPonarr(unsigned int index) { return Ponarr.at(index); }
    float getPonint(unsigned int index) { return Ponint.at(index); }

};

// Note: currently used only through derived class StochasticOccupantsPresence (virtual class)
class Occupants {

private:

    float occupantsNumber;

protected:

    // presence parameters
    StochasticPresenceParameters *presParam;
    // window parameters
    StochasticWindowParameters *winParam;
    // blinds parameters
    StochasticBlindsParameters *blindsParam;
    // lights parameters
    StochasticLightsParameters *lightsParam;

public:

    Occupants(float occupantsNumber) { this->occupantsNumber = occupantsNumber; presParam = NULL; winParam = NULL; blindsParam = NULL; lightsParam = NULL; }
    virtual ~Occupants() {}

#pragma GCC diagnostic ignored "-Wunused-parameter"
    float getOccupantsNumber() { return occupantsNumber; }
    virtual float getOccupantsFraction(unsigned int day, unsigned int hour, int fracHour=0) { return 1.f; }
    virtual float getCurrentDuration(unsigned int day, unsigned int hour, unsigned int fracHour) { return 0.f; }
    virtual float getFutureDuration(unsigned int day, unsigned int hour, unsigned int fracHour) { return 0.f; }
#pragma GCC diagnostic warning "-Wunused-parameter"

    // gets the pointer to the window parameters
    StochasticWindowParameters* getStochasticWindowParameters() { return winParam; }
    // gets the pointer to the blinds parameters
    StochasticBlindsParameters* getStochasticBlindsParameters() { return blindsParam; }
    // gets the pointer to the lights parameters
    StochasticLightsParameters* getStochasticLightsParameters() { return lightsParam; }

};

class DayProfile{
private:
    unsigned int id;
    string name;
    float waterConsumption = 0.f;
    vector<float> profile;
    ostream logStream;

public:
    DayProfile(int id, string name, vector<float> profile):id(id),name(name),profile(profile),logStream(std::cout.rdbuf()){}
    DayProfile():DayProfile(0,"empty day profile",vector<float>(24,0)) {}
    DayProfile(TiXmlHandle hdl, ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        TiXmlElement *elem = hdl.ToElement();

        id = to<unsigned int>(elem->Attribute("id"));
        if (elem->Attribute("name")){
            elem->QueryStringAttribute("name",&name);
        }
        else{
            elem->QueryStringAttribute("id",&name);
            name = "Day profile "+name;
        }
        // gets the water consumption if it exists
        if (elem->Attribute("waterConsumption")) waterConsumption = to<float>(elem->Attribute("waterConsumption"));

        logStream << "Day id= " << id << ": ";

        int j=1;
        do {
            string attrib;
            elem->QueryStringAttribute(string("p"+toString(j)).c_str(), &attrib);
            logStream << attrib << " ";
            profile.push_back(to<float>(attrib));

        } while ( elem->Attribute("p"+toString(++j)) );

        logStream << endl;
    }
    DayProfile(const DayProfile& dp):logStream(std::cout.rdbuf()){
        logStream.rdbuf(dp.logStream.rdbuf());
        id = dp.id;
        name = dp.name;
        profile = dp.profile;
    }

    unsigned int getId(){ return id;}
    void setId(unsigned int i) {id = i; }

    float getHourValue(unsigned int h){ return profile.at(h-1); }
    string getName(){ return name; }
    void setName(string n) { name = n; }
    float getWaterConsumption() { return waterConsumption; }

    void writeXML(ofstream& file, string tab="", string type="Occupancy"){
        file << tab << "<" << type << "DayProfile id=\""<< id << "\" name=\"" << name << "\" ";
        for (unsigned int i=0; i<profile.size(); ++i){
            file << "p" << i+1 << "=\"" << profile[i] << "\" ";
        }
        file << "/>" << endl;
    }
};

class YearProfile{
private:
    unsigned int id;
    string name;
    vector<DayProfile*> profile;
    ostream logStream;
    DayProfile* emptyDayProfile;

public:
    YearProfile(int id, string name, vector<DayProfile*> profile):id(id),name(name),profile(profile),logStream(std::cout.rdbuf()){}
    YearProfile(TiXmlHandle hdl, map<unsigned int, DayProfile*> mapId2DayProfile, ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

        // create an empty day profile
        emptyDayProfile = new DayProfile();

        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        TiXmlElement *elem = hdl.ToElement();

        id = to<unsigned int>(elem->Attribute("id"));
        if (elem->Attribute("name")){
            elem->QueryStringAttribute("name",&name);
        }
        else{
            elem->QueryStringAttribute("id",&name);
            name = "Year profile "+name;
        }

        logStream << "Year id= " << id << ": ";

        int j=1;
        do {
            unsigned int attrib;
            elem->QueryUnsignedAttribute(("d"+toString(j)).c_str(),&attrib);
            logStream << attrib << " ";
            if(mapId2DayProfile.count(attrib)==1){
                profile.push_back(mapId2DayProfile.at(attrib));
            }
            else if (attrib == 0) {
                // default empty day profile
                profile.push_back(emptyDayProfile);
            }
            else{
                throw("Year Profile id=" + toString(id) + " at Day d=" + toString(j) + ": Day Profile id=" + toString(attrib) + " does not exist.");
            }

        } while ( elem->Attribute("d"+toString(++j)) );

        logStream << endl;
    }
    YearProfile(const YearProfile& yp):logStream(std::cout.rdbuf()){
        logStream.rdbuf(yp.logStream.rdbuf());
        id = yp.id;
        name = yp.name;
        profile = yp.profile;
    }

    ~YearProfile() {
        if (emptyDayProfile) delete emptyDayProfile;
    }

    void fillNull(DayProfile* d){
        for (unsigned int i=0; i<profile.size(); ++i){
            if (profile[i]==NULL){
                profile[i] = d;
            }
        }
    }

    DayProfile* getDayProfile(unsigned int d){
        return profile.at((d-1)%365); // remain within 1-365 for the days, JK - 21.06.2015
    }

    unsigned int getNDayProfiles(){
        return profile.size();
    }

    unsigned int getId(){ return id;}
    void setId(unsigned int i) {id = i; }
    string getName(){ return name; }
    void setName(string n) { name = n; }

    void writeXML(ofstream& file, string tab="", string type="Occupancy"){
        file << tab << "<" << type << "YearProfile id=\""<< id << "\" name=\"" << name << "\" ";
        for (unsigned int i=0; i<profile.size(); ++i){
            file << "d" << i+1 << "=\"" << profile[i]->getId() << "\" ";
        }
        file << "/>" << endl;
    }
};

class OccupancyProfiles {

protected:

    map<string,DayProfile*> dayProfiles;
    map<string,YearProfile*> yearProfiles;
    map<unsigned int, DayProfile*> xmlId2DayProfile;
    map<unsigned int, YearProfile*> xmlId2YearProfile;
    ostream logStream;

public:

    static DayProfile emptyDay;
    static YearProfile emptyYear;
    OccupancyProfiles(ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        addDayProfile(&emptyDay);
        addYearProfile(&emptyYear);
    }

    virtual ~OccupancyProfiles(){
        for (map<string,DayProfile*>::iterator it=dayProfiles.begin();it!=dayProfiles.end();++it) {
            if (it->second == NULL){
                cout << "ERROR: null day profile with name " << it->first;
            }
            else if (it->second != &emptyDay){
                delete (*it).second;
            }
        }
        for (map<string,YearProfile*>::iterator it=yearProfiles.begin();it!=yearProfiles.end();++it) {
            if (it->second == NULL){
                cout << "ERROR: null year profile with name " << it->first;
            }
            else if (it->second != &emptyYear){
                delete (*it).second;
            }
        }
    }

    virtual void readFromXml(TiXmlHandle hdl) {
        int i=0;
        while (hdl.ChildElement("OccupancyDayProfile",i).ToElement()) {

            DayProfile* d = new DayProfile(hdl.ChildElement("OccupancyDayProfile",i),&logStream);
            xmlId2DayProfile.insert(pair<unsigned int, DayProfile*>(d->getId(),d));
            addDayProfile(d);
            ++i;
        }
        logStream << "Occupancy DAY profiles: " << dayProfiles.size()-1 << " loaded." << endl;

        logStream << "Loading occupancy YEAR profiles" << endl << flush;

        i=0;

        while (hdl.ChildElement("OccupancyYearProfile",i).ToElement()) {

            YearProfile* y = new YearProfile(hdl.ChildElement("OccupancyYearProfile",i), xmlId2DayProfile, &logStream);
            addYearProfile(y);
            xmlId2YearProfile.insert(pair<unsigned int, YearProfile*>(y->getId(),y));

            ++i;
        }
        logStream << "Occupancy YEAR profiles: " << yearProfiles.size()-1 << " loaded." << endl;
    }

    void clearIdMap(){
        xmlId2DayProfile.clear();
        xmlId2YearProfile.clear();
    }

    YearProfile* getYearProfileFromXmlId(unsigned int id) {
        if (xmlId2YearProfile.find(id)!=xmlId2YearProfile.end()) return xmlId2YearProfile.at(id);
        else return nullptr;
    }

    void addDayProfile(DayProfile* profile) {
        string name = profile->getName();
        string originalName = name;
        int i = 2;
        while(dayProfiles.count(name)==1){
            name = originalName+" ("+toString(i)+")";
            ++i;
        }
        profile->setName(name);
        dayProfiles.insert(pair<string,DayProfile*>(profile->getName(),profile));
        logStream << "Added day profile with name " << profile->getName() << endl;
    }

    const map<string, DayProfile*>* getDayProfiles(){ return &dayProfiles;}

    void addYearProfile(YearProfile* profile) {
        string name = profile->getName();
        string originalName = name;
        int i = 2;
        while(yearProfiles.count(name)==1){
            name = originalName+" ("+toString(i)+")";
            ++i;
        }
        profile->setName(name);
        profile->fillNull(dayProfiles["empty day profile"]);
        yearProfiles.insert(pair<string,YearProfile*>(profile->getName(),profile));
        logStream << "Added year profile with name " << profile->getName() << endl;
    }

    const map<string, YearProfile*>* getYearProfiles(){ return &yearProfiles;}

    virtual void print(){
        logStream << endl << "Occupancy profiles" << endl << "Day profiles: " << endl << flush;
        for (map<string, DayProfile*>::iterator it = dayProfiles.begin(); it != dayProfiles.end(); ++it){
            logStream << it->second->getName() << "(" << it->first <<"), id=" << it->second->getId() << endl;
        }
        logStream << endl << "Year profiles: " << endl;
        for (map<string, YearProfile*>::iterator it = yearProfiles.begin(); it != yearProfiles.end(); ++it){
            logStream << it->second->getName() << "(" << it->first <<"), id=" << it->second->getId() << endl;
        }
    }

    virtual void writeXML(ofstream& file, set<YearProfile*> usedYearProfiles, string tab=""){
        set<unsigned int> usedIds;
        unsigned int newId=0;
        set<DayProfile*> usedDayProfiles;

        // Check used YearProfiles: set unique ids and get used DayProfiles
        for (set<YearProfile*>::const_iterator iter=usedYearProfiles.begin(); iter != usedYearProfiles.end(); ++ iter){
            // find the next unused id
            while(usedIds.count(newId)!=0){
                ++newId;
            }
            if(usedIds.count((*iter)->getId())!=0){
                // set a free id for this YearProfile
                (*iter)->setId(newId);
                usedIds.insert(newId);
            }
            // Get all used DayProfiles
            for (unsigned int i=0; i<(*iter)->getNDayProfiles(); ++i){
                usedDayProfiles.insert((*iter)->getDayProfile(i+1));
            }
        }

        // Set unique ids and write used DayProfiles
        usedIds.clear();
        newId = 0;
        for (set<DayProfile*>::const_iterator iter=usedDayProfiles.begin(); iter != usedDayProfiles.end(); ++ iter){
            // find the next unused id
            while(usedIds.count(newId)!=0){
                ++newId;
            }
            if(usedIds.count((*iter)->getId())!=0){
                // set a free id for this DayProfile
                (*iter)->setId(newId);
                usedIds.insert(newId);
            }
            (*iter)->writeXML(file, tab);
        }
        file << endl;

        // Write used YearProfiles
        for (set<YearProfile*>::const_iterator iter=usedYearProfiles.begin(); iter != usedYearProfiles.end(); ++ iter){
            (*iter)->writeXML(file, tab);
        }
        file << endl;
    }
};

class DHWProfiles: public OccupancyProfiles {

public:

    DHWProfiles(ostream* pLogStr=NULL):OccupancyProfiles(pLogStr) {}

    void readFromXml(TiXmlHandle hdl) {

        logStream << "Loading Domestic Hot Water (DHW) DAY profiles" << endl;

        int i=0;
        while (hdl.ChildElement("DHWDayProfile",i).ToElement()) {

            DayProfile* d = new DayProfile(hdl.ChildElement("DHWDayProfile",i),&logStream);
            xmlId2DayProfile.insert(pair<unsigned int, DayProfile*>(d->getId(),d));
            addDayProfile(d);
            ++i;
        }
        logStream << "DHW DAY profiles: " << dayProfiles.size()-1 << " loaded." << endl;

        logStream << "Loading Domestic Hot Water (DHW) YEAR profiles" << endl << flush;

        i=0;
        while (hdl.ChildElement("DHWYearProfile",i).ToElement()) {

            YearProfile* y = new YearProfile(hdl.ChildElement("DHWYearProfile",i), xmlId2DayProfile, &logStream);
            addYearProfile(y);
            xmlId2YearProfile.insert(pair<unsigned int, YearProfile*>(y->getId(),y));
            ++i;
        }
        logStream << "DHW YEAR profiles: " << yearProfiles.size()-1 << " loaded." << endl;
    }

    void print(){
        logStream << endl << "Domestic Hot Water (DHW) profiles" << endl << "Day profiles: " << endl << flush;
        for (map<string, DayProfile*>::iterator it = dayProfiles.begin(); it != dayProfiles.end(); ++it){
            logStream << it->second->getName() << "(" << it->first <<"), id=" << it->second->getId() << endl;
        }
        logStream << endl << "Year profiles: " << endl;
        for (map<string, YearProfile*>::iterator it = yearProfiles.begin(); it != yearProfiles.end(); ++it){
            logStream << it->second->getName() << "(" << it->first <<"), id=" << it->second->getId() << endl;
        }
    }

    void writeXML(ofstream& file, set<YearProfile*> usedYearProfiles, string tab=""){
        set<unsigned int> usedIds;
        unsigned int newId=0;
        set<DayProfile*> usedDayProfiles;

        // Check used YearProfiles: set unique ids and get used DayProfiles
        for (set<YearProfile*>::const_iterator iter=usedYearProfiles.begin(); iter != usedYearProfiles.end(); ++ iter){
            // find the next unused id
            while(usedIds.count(newId)!=0){
                ++newId;
            }
            if(usedIds.count((*iter)->getId())!=0){
                // set a free id for this YearProfile
                (*iter)->setId(newId);
                usedIds.insert(newId);
            }
            // Get all used DayProfiles
            for (unsigned int i=0; i<(*iter)->getNDayProfiles(); ++i){
                usedDayProfiles.insert((*iter)->getDayProfile(i+1));
            }
        }

        // Set unique ids and write used DayProfiles
        usedIds.clear();
        newId = 0;
        for (set<DayProfile*>::const_iterator iter=usedDayProfiles.begin(); iter != usedDayProfiles.end(); ++ iter){
            // find the next unused id
            while(usedIds.count(newId)!=0){
                ++newId;
            }
            if(usedIds.count((*iter)->getId())!=0){
                // set a free id for this DayProfile
                (*iter)->setId(newId);
                usedIds.insert(newId);
            }
            (*iter)->writeXML(file, tab, "DHW");
        }
        file << endl;

        // Write used YearProfiles
        for (set<YearProfile*>::const_iterator iter=usedYearProfiles.begin(); iter != usedYearProfiles.end(); ++ iter){
            (*iter)->writeXML(file, tab, "DHW");
        }
        file << endl;
    }

};

class TemperatureProfiles {

private:

    class DayProfile {
    private:
        unsigned int id;
        string name;
        vector<float> profile_Tmin;
        vector<float> profile_Tmax;
        ostream logStream;

    public:
        DayProfile(TiXmlHandle hdl, ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

            // logStream is directed by default to the "cout" streambuf
            if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
                logStream.rdbuf(pLogStr->rdbuf());
            if (!logStream.good())
                throw(string("Unable to define correctly the logStream."));

            TiXmlElement* elem = hdl.ToElement();

            id = to<unsigned int>(elem->Attribute("id"));
            if (elem->Attribute("name")){
                elem->QueryStringAttribute("name",&name);
            }
            else{
                elem->QueryStringAttribute("id",&name);
                name = "Day profile "+name;
            }

            logStream << "Day id= " << id << " - Tmin: ";
            // reads the Tmin profile
            elem = hdl.FirstChildElement("Tmin").ToElement();
            if (elem) {
                int j=1;
                do {
                    string attrib;
                    elem->QueryStringAttribute(string("h"+toString(j)).c_str(), &attrib);
                    logStream << attrib << " ";
                    profile_Tmin.push_back(to<float>(attrib));

                } while ( elem->Attribute("h"+toString(++j)) );
                logStream << endl;
            }
            else logStream << "none." << endl;

            logStream << "Day id= " << id << " - Tmax: ";
            // reads the Tmax profile
            elem = hdl.FirstChildElement("Tmax").ToElement();
            if (elem) {
                int j=1;
                do {
                    string attrib;
                    elem->QueryStringAttribute(string("h"+toString(j)).c_str(), &attrib);
                    logStream << attrib << " ";
                    profile_Tmax.push_back(to<float>(attrib));

                } while ( elem->Attribute("h"+toString(++j)) );
                logStream << endl;
            }
            else logStream << "none." << endl;

        }

        unsigned int getId() { return id; }
        string getName() { return name; }

        float getHourValue_Tmin(unsigned int h){ return profile_Tmin.at(h-1); }
        float getHourValue_Tmax(unsigned int h){ return profile_Tmax.at(h-1); }

        void writeXML(ofstream& file, string tab="") {
            file << tab << "<TemperatureDayProfile id=\""<< id << "\" name=\"" << name << "\" ";
            if (!profile_Tmin.empty()) {
                file << tab << "\t" << "<Tmin ";
                for (unsigned int i=0; i<profile_Tmin.size(); ++i) {
                    file << "h" << i+1 << "=\"" << profile_Tmin[i] << "\" ";
                }
                file << "/>" << endl;
            }
            if (!profile_Tmax.empty()) {
                file << tab << "\t" << "<Tmax ";
                for (unsigned int i=0; i<profile_Tmax.size(); ++i) {
                    file << "h" << i+1 << "=\"" << profile_Tmax[i] << "\" ";
                }
                file << "/>" << endl;
            }
        }
    };
    class YearProfile{
    private:
        unsigned int id;
        string name;
        vector<DayProfile*> profile;
        ostream logStream;

    public:
        YearProfile(TiXmlHandle hdl, map<unsigned int, DayProfile*> mapId2DayProfile, ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

            // logStream is directed by default to the "cout" streambuf
            if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
                logStream.rdbuf(pLogStr->rdbuf());
            if (!logStream.good())
                throw(string("Unable to define correctly the logStream."));

            TiXmlElement *elem = hdl.ToElement();

            id = to<unsigned int>(elem->Attribute("id"));
            if (elem->Attribute("name")){
                elem->QueryStringAttribute("name",&name);
            }
            else{
                elem->QueryStringAttribute("id",&name);
                name = "Year profile "+name;
            }

            logStream << "Year id= " << id << ": ";

            int j=1;
            do {
                unsigned int attrib;
                elem->QueryUnsignedAttribute(("d"+toString(j)).c_str(),&attrib);
                logStream << attrib << " ";
                if(mapId2DayProfile.count(attrib)==1){
                    profile.push_back(mapId2DayProfile.at(attrib));
                }
                else{
                    throw("Year Profile id=" + toString(id) + " at Day d=" + toString(j) + ": Day Profile id=" + toString(attrib) + " does not exist.");
                }

            } while ( elem->Attribute("d"+toString(++j)) );

            logStream << endl;
        }

        DayProfile* getDayProfile(unsigned int d) {
            return profile.at((d-1)%365); // remain within 1-365 for the days, JK - 21.06.2015
        }

        unsigned int getId() { return id;}
        string getName(){ return name; }

        void writeXML(ofstream& file, string tab="", string type="Occupancy"){
            file << tab << "<" << type << "TemperatureYearProfile id=\""<< id << "\" name=\"" << name << "\" ";
            for (unsigned int i=0; i<profile.size(); ++i){
                file << "d" << i+1 << "=\"" << profile[i]->getId() << "\" ";
            }
            file << "/>" << endl;
        }
    };
    map<unsigned int, DayProfile*> dayProfiles;
    map<unsigned int, YearProfile*> yearProfiles;
    ostream logStream;

public:

    TemperatureProfiles(TiXmlHandle hdl):logStream(std::cout.rdbuf()) {
        int i=0;
        while (hdl.ChildElement("TemperatureDayProfile",i).ToElement()) {
            DayProfile* d = new DayProfile(hdl.ChildElement("TemperatureDayProfile",i),&logStream);
            dayProfiles.insert(pair<unsigned int, DayProfile*>(d->getId(),d));
            ++i;
        }
        logStream << "Temperature DAY profiles: " << dayProfiles.size() << " loaded." << endl;

        logStream << "Loading temperature YEAR profiles" << endl << flush;

        i=0;

        while (hdl.ChildElement("TemperatureYearProfile",i).ToElement()) {
            YearProfile* y = new YearProfile(hdl.ChildElement("TemperatureYearProfile",i), dayProfiles, &logStream);
            yearProfiles.insert(pair<unsigned int, YearProfile*>(y->getId(),y));
            ++i;
        }
        logStream << "Temperature YEAR profiles: " << yearProfiles.size() << " loaded." << endl;
    }

    ~TemperatureProfiles() {
        for (map<unsigned int,DayProfile*>::iterator it=dayProfiles.begin();it!=dayProfiles.end();++it) {
            if (it->second) delete (*it).second;
        }
        for (map<unsigned int,YearProfile*>::iterator it=yearProfiles.begin();it!=yearProfiles.end();++it) {
            if (it->second) delete (*it).second;
        }
    }

    YearProfile* getYearProfile(unsigned int id) {
        if (yearProfiles.find(id)!=yearProfiles.end()) return yearProfiles.at(id);
        else return nullptr;
    }

    void print() {
        logStream << endl << "Temperature profiles" << endl << "Day profiles: " << endl << flush;
        for (map<unsigned int, DayProfile*>::iterator it = dayProfiles.begin(); it != dayProfiles.end(); ++it) {
            logStream << it->second->getName() << "(" << it->first <<"), id=" << it->second->getId() << endl;
        }
        logStream << endl << "Year profiles: " << endl;
        for (map<unsigned int, YearProfile*>::iterator it = yearProfiles.begin(); it != yearProfiles.end(); ++it) {
            logStream << it->second->getName() << "(" << it->first <<"), id=" << it->second->getId() << endl;
        }
    }

    void writeXML(ofstream& file, string tab="") {
        // write the dayProfiles first
        for (map<unsigned int, DayProfile*>::iterator it = dayProfiles.begin(); it != dayProfiles.end(); ++it) {
            (it->second)->writeXML(file, tab+"\t");
        }
        // write the yearProfiles afterwards
        for (map<unsigned int, YearProfile*>::iterator it = yearProfiles.begin(); it != yearProfiles.end(); ++it) {
            (it->second)->writeXML(file, tab+"\t");
        }
    }
};

class DeviceType {

private:
    unsigned int id;
    string name;
    class Device {
    public:
        string name;
        float avgPower;
        float convectiveFraction = 0.8f, radiativeFraction = 0.2f;
        vector<float> profile;

        Device(TiXmlHandle hdl) {

            if (hdl.ToElement()->Attribute("name")){
                hdl.ToElement()->QueryStringAttribute("name",&name);
            }
            if (to<float>(hdl.ToElement()->Attribute("avgPower")))
                avgPower = to<float>(hdl.ToElement()->Attribute("avgPower"));
            else throw(string("Missing avgPower in Device: " + name));
            if (to<float>(hdl.ToElement()->Attribute("convectiveFraction")))
                convectiveFraction = to<float>(hdl.ToElement()->Attribute("convectiveFraction"));
            if (to<float>(hdl.ToElement()->Attribute("radiativeFraction")))
                radiativeFraction  = to<float>(hdl.ToElement()->Attribute("radiativeFraction"));

            int j=1;
            do {
                string attrib;
                hdl.ToElement()->QueryStringAttribute(string("p"+toString(j)).c_str(), &attrib);
                profile.push_back(to<float>(attrib));

            } while ( hdl.ToElement()->Attribute("p"+toString(++j)) );
            if (j<24) throw(string("Missing profile, only " + toString(j) + "/24 fractions in Device: " + name));

        }
        void writeXML(ofstream& file, string tab="") {
            file << tab << "<Device name=\"" << name << "\" avgPower=\"" << avgPower << "\" convectiveFraction=\"" << convectiveFraction << "\" radiativeFraction=\"" << radiativeFraction << "\"";
            for (size_t i=0;i<profile.size();++i) file << " p" << i+1 << "=\"" << profile.at(i) << "\"";
            file << "/>" << endl;
        }
    };
    vector<Device> devices;
    ostream logStream;

public:
    DeviceType(TiXmlHandle hdl, ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        id = to<unsigned int>(hdl.ToElement()->Attribute("id"));
        if (hdl.ToElement()->Attribute("name")){
            hdl.ToElement()->QueryStringAttribute("name",&name);
        }

        // load the devices
        int i=0;
        while (hdl.ChildElement("Device",i).ToElement()) {
            devices.emplace_back(hdl.ChildElement("Device",i++));
        }
        logStream << "Device type: " << id << "(" << name << "): " << devices.size() << " loaded." << endl;

    }
    unsigned int getId() { return id; }
    size_t getnDevices() { return devices.size(); }
    float getDeviceAvgPower(size_t device_index) { return devices.at(device_index).avgPower; }
    float getDeviceConvectiveFraction(size_t device_index) { return devices.at(device_index).convectiveFraction; }
    float getDeviceRadiativeFraction(size_t device_index)  { return devices.at(device_index).radiativeFraction; }
    float getDeviceProbability(size_t device_index, size_t hour) { return devices.at(device_index).profile.at(hour-1); } // hour \in [1,24]
    void writeXML(ofstream& file, string tab="") {
        file << tab << "<DeviceType id=\"" << id << "\" name=\"" << name << "\">" << endl;
        for (vector<Device>::iterator it=devices.begin(); it!=devices.end(); ++it) (*it).writeXML(file, tab+"\t");
        file << tab << "</DeviceType>" << endl;
    }

};

class ActivityType {

private:
    unsigned int id;
    string name;
    class Activity {
    public:
        string name;
        unsigned int deviceType;
        vector<float> profile;

        Activity(TiXmlHandle hdl) {

            if (hdl.ToElement()->Attribute("name")){
                hdl.ToElement()->QueryStringAttribute("name",&name);
            }
            deviceType = to<unsigned int>(hdl.ToElement()->Attribute("deviceType"));

            int j=1;
            do {
                string attrib;
                hdl.ToElement()->QueryStringAttribute(string("p"+toString(j)).c_str(), &attrib);
                profile.push_back(to<float>(attrib));

            } while ( hdl.ToElement()->Attribute("p"+toString(++j)) );

        }
        void writeXML(ofstream& file, string tab="") {
            file << tab << "<Activity name=\"" << name << "\" deviceType=\"" << deviceType << "\"";
            for (size_t i=0;i<profile.size();++i) file << " p" << i+1 << "=\"" << profile.at(i) << "\"";
            file << "/>" << endl;
        }

    };
    vector<Activity> activites;
    ostream logStream;

public:
    ActivityType(TiXmlHandle hdl, ostream* pLogStr=NULL):logStream(std::cout.rdbuf()) {

        // logStream is directed by default to the "cout" streambuf
        if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
            logStream.rdbuf(pLogStr->rdbuf());
        if (!logStream.good())
            throw(string("Unable to define correctly the logStream."));

        id = to<unsigned int>(hdl.ToElement()->Attribute("id"));
        if (hdl.ToElement()->Attribute("name")){
            hdl.ToElement()->QueryStringAttribute("name",&name);
        }

        // load the devices
        int i=0;
        while (hdl.ChildElement("Activity",i).ToElement()) {
            activites.emplace_back(hdl.ChildElement("Activity",i++));
        }
        logStream << "Activity type: " << id << "(" << name << "): " << activites.size() << " loaded." << endl;

    }
    unsigned int getId() { return id; }
    size_t getnActivities() { return activites.size(); }
    unsigned int getActivityDeviceType(size_t activity_index) { return activites.at(activity_index).deviceType; }
    float getActivityProbability(size_t activity_index, size_t hour) { return activites.at(activity_index).profile.at(hour-1); } // hour \in [1,24]
    void writeXML(ofstream& file, string tab="") {
        file << tab << "<ActivityType id=\"" << id << "\" name=\"" << name << "\">" << endl;
        for (vector<Activity>::iterator it=activites.begin(); it!=activites.end(); ++it) (*it).writeXML(file, tab+"\t");
        file << tab << "</ActivityType>" << endl;
    }

};

#pragma GCC diagnostic ignored "-Wunused-parameter"
// Note: not used anymore !!!
class DeterministicOccupantsPresence : public Occupants {

    private:

    vector<float> occupantsWeekday,occupantsSaturday,occupantsSunday;

    public:

    DeterministicOccupantsPresence(float occupantsNumber, vector<float> occupantsWeekday, vector<float> occupantsSaturday, vector<float> occupantsSunday)
     : Occupants(occupantsNumber), occupantsWeekday(occupantsWeekday), occupantsSaturday(occupantsSaturday), occupantsSunday(occupantsSunday) {}

    float getOccupantsFraction(unsigned int day, unsigned int hour, int fracHour=0) {
        unsigned int dayOfTheWeek = (day-1) % 7;
        //cout << "dayOfTheWeek: " << dayOfTheWeek;
        if (dayOfTheWeek < 5) { /*cout << "\toccupancy: " << occupantsWeekday[hour-1] << endl;*/ return occupantsWeekday[hour-1]; }
        else if (dayOfTheWeek < 6) { /*cout << "\toccupancy: " << occupantsSaturday[hour-1] << endl;*/ return occupantsSaturday[hour-1]; }
        else { /*cout << "\toccupancy: " << occupantsSunday[hour-1] << endl;*/ return occupantsSunday[hour-1]; }
    }

};

class StochasticOccupantsPresence : public Occupants {

private:

    vector<float> profile;
    vector<float> current_dur;
    vector<float> future_dur;

public:

    StochasticOccupantsPresence(float occupantsNumber, StochasticPresenceParameters *presParam, StochasticWindowParameters *winParam, StochasticBlindsParameters *blindsParam, StochasticLightsParameters *lightsParam);
    ~StochasticOccupantsPresence() { delete presParam; delete winParam; delete blindsParam; delete lightsParam; }

    float getOccupantsFraction(unsigned int day, unsigned int hour, int fracHour);
    float getCurrentDuration(unsigned int day, unsigned int hour, unsigned int fracHour);
    float getFutureDuration(unsigned int day, unsigned int hour, unsigned int fracHour);

    static double getT01(double pcurr, double pnext, double shuff);
    static double getT11(double pcurr, double pnext, double shuff);
    vector<bool> presence();

};

class Behaviour {

  public:

    virtual double Nvent(double Ta, double Tout) { return 0.0; }

};

class DeterministicBehaviour: public Behaviour {

  private:

    double NventMax;

  public:

    DeterministicBehaviour(double ach) { NventMax = ach; }

    double Nvent(double Ta, double Tout) {

        double a = -6.22;
        double b = 0.23;
        if (Ta < 26.0) return NventMax*exp(a+b*Ta)/(1.0+exp(a+b*Ta)); else return 0.0;

    }

};

#pragma GCC diagnostic warning "-Wunused-parameter"
#endif
