#include "surface.h"
#include "building.h"
#include "district.h"
#include "scene.h"
#include "GEOMBoundingSphereCalc.h"

Material::Material(TiXmlHandle hdl, ostream* pLogStr):logStream(std::cout.rdbuf()) {

    // logStream is directed by default to the "cout" streambuf
    if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    id = to<int>(hdl.ToElement()->Attribute("id"));
    if((hdl.ToElement())->Attribute("name")!=NULL){
        name = hdl.ToElement()->Attribute("name");
    }
    if (hdl.ToElement()){
        TiXmlElement *elem = hdl.ToElement();
        name = (elem->Attribute("name"));
        conductivity = to<float>(elem->Attribute("Conductivity"));
        cp = to<float>(elem->Attribute("Cp"));
        density = to<float>(elem->Attribute("Density"));
        nre = (elem->Attribute("NRE") == NULL) ? 0.f : to<float>(elem->Attribute("NRE"));
        gwp = (elem->Attribute("GWP") == NULL) ? 0.f : to<float>(elem->Attribute("GWP"));
        ubp = (elem->Attribute("UBP") == NULL) ? 0.f : to<float>(elem->Attribute("UBP"));
    }
    else logStream << "Could not read a Material correctly";
}

Composite::Composite(TiXmlHandle hdl, ostream* pLogStr):insulationLayer(NULL),logStream(std::cout.rdbuf()) {

    // logStream is directed by default to the "cout" streambuf
    if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    id = to<int>((hdl.ToElement())->Attribute("id"));
    if((hdl.ToElement())->Attribute("name")!=NULL){
        hdl.ToElement()->QueryStringAttribute("name",&name);
    }
    if((hdl.ToElement())->Attribute("category")!=NULL){
        hdl.ToElement()->QueryStringAttribute("category",&category);
    }
    if (hdl.FirstChildElement().ToElement()){
        TiXmlElement *elem = hdl.FirstChildElement().ToElement();
        while(elem){
            float nre = (elem->Attribute("NRE") == NULL) ? 0.f : to<float>(elem->Attribute("NRE"));
            float gwp = (elem->Attribute("GWP") == NULL) ? 0.f : to<float>(elem->Attribute("GWP"));
            float ubp = (elem->Attribute("UBP") == NULL) ? 0.f : to<float>(elem->Attribute("UBP"));
            vLayer.push_back(Layer(to<float>(elem->Attribute("Thickness")),
                                   to<float>(elem->Attribute("Conductivity")),
                                   to<float>(elem->Attribute("Cp")),
                                   to<float>(elem->Attribute("Density")),
                                   nre, gwp, ubp));
            elem = elem->NextSiblingElement();
        }
        // computes the standard Uvalue according to inside and outside layers of conductance 8 and 25
        Uvalue = 1.f/(0.125f + getResistance() + 0.04f);
    }
    else Uvalue = to<float>((hdl.ToElement())->Attribute("Uvalue"));
    logStream << "Walltype " << id << ", Uvalue = " << Uvalue << endl << flush;
}

/**
 * @brief Composite::Composite: Create a Uvalue-only composite with given id
 * @param Uvalue
 * @param wtId
 * @param name
 * @param pLogStr
 */
Composite::Composite(float Uvalue, int wtId, string name, ostream* pLogStr):name(name),Uvalue(Uvalue),id(wtId),logStream(std::cout.rdbuf()) {
    // logStream is directed by default to the "cout" streambuf
    if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    logStream << "Adding a Uvalue only composite with id " << wtId << endl << flush;
}


Composite::Composite(TiXmlHandle hdl, map<string,Material*> materials, ostream* pLogStr):insulationLayer(NULL),logStream(std::cout.rdbuf()) {

    // logStream is directed by default to the "cout" streambuf
    if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    id = to<int>((hdl.ToElement())->Attribute("id"));
    if((hdl.ToElement())->Attribute("name")!=NULL){
        hdl.ToElement()->QueryStringAttribute("name",&name);
        //name = to<string>((hdl.ToElement())->Attribute("name"));
    }
    if((hdl.ToElement())->Attribute("category")!=NULL){
        hdl.ToElement()->QueryStringAttribute("category",&category);
    }
    if (hdl.FirstChildElement().ToElement()){
        TiXmlElement *elem = hdl.FirstChildElement().ToElement();
        while(elem){
            try{
                Material* m = materials.at(elem->Attribute("materialID"));

                addLayer(to<float>(elem->Attribute("thickness")), m,
                        /*isInsulationLayer:*/(elem->Attribute("insulation")!=NULL && to<string>(elem->Attribute("insulation"))=="true"));
            }
            catch(exception e){
                logStream << "ERROR: could not find material with id=" << elem->Attribute("materialID") << ". " << e.what() << endl;
            }
            elem = elem->NextSiblingElement();
        }
        // computes the standard Uvalue according to inside and outside layers of conductance 8 and 25
        Uvalue = 1.f/(0.125f + getResistance() + 0.04f);
    }
    else Uvalue = to<float>((hdl.ToElement())->Attribute("Uvalue"));
    logStream << "Composite " << id << ", Uvalue = " << Uvalue << endl << flush;
}

<<<<<<< HEAD
Composite::Composite(Composite const& c, float insulationTh):logStream(c.logStream.rdbuf()) {

    id = c.id;
    name = c.name;
    category = c.category;
    Uvalue = c.Uvalue;
    for(unsigned int i=0; i<c.vLayer.size();++i){
        vLayer.push_back(c.vLayer[i]);
        if(&(c.vLayer[i])==c.insulationLayer){
            insulationLayer = &(vLayer.back());
            if (insulationTh != -1){
                insulationLayer->setxx(insulationTh);
            }
        }
    }
}

=======
>>>>>>> cf042395e1dcb11627fd9bf5653a3121bccaadb0
void Composite::getSimplifiedNode(float& Cw, float& Kw1, float& Kw2) {

    // creation of C and R (of size vLayer.size())
    vector<double> C(vLayer.size(), 0.), R(vLayer.size(), 0.);

    // Capacitance calculation
    Cw=0.f;
    for (unsigned int i=0;i<vLayer.size();++i) {
        C[i]=vLayer.at(i).getrho()*vLayer.at(i).getCp()*vLayer.at(i).getxx();
        Cw+=C[i];
    }

    // Resistor calculation, between outside and node i
    for (unsigned int i=0;i<vLayer.size();++i) {
      R[i]=0.5*(vLayer.at(i).getxx()/vLayer.at(i).getkw());
      for (unsigned int j=0;j<i;++j) R[i]+=(vLayer.at(j).getxx()/vLayer.at(j).getkw());
    }

    // temporary variables for Resistor calculation
    double Rw1=0., Rtot=0.;

    for (unsigned int i=0;i<vLayer.size();i++) {
      Rw1+=R[i]*C[i]/Cw;
      Rtot+=(vLayer.at(i).getxx()/vLayer.at(i).getkw());
    }

    // Conductance calculation
    Kw1=1./Rw1;
    Kw2=1./(Rtot-Rw1);

    return;
}

void Composite::addLayer(float thickness, Material* m, bool isInsulationLayer){
    vLayer.push_back(Layer(thickness,
                       m->getConductivity(),
                       m->getCp(),
                       m->getDensity(),
                       m->getNRE(),
                       m->getGWP(),
                       m->getUBP()));
    if (isInsulationLayer){
        insulationLayer = &(vLayer[vLayer.size()-1]);
    }
}

void Composite::writeXML(ofstream& file, string tab=""){
    if(!vLayer.empty()){
        file << tab << "<Composite id=\"" << id << "\" name=\"" << name << "\" category=\"" << category << "\">" << endl;
        for (unsigned int i=0; i<vLayer.size(); ++i)
            vLayer[i].writeXML(file, &vLayer[i]==insulationLayer, tab+"\t");
        file << tab << "</Composite>" << endl;
    }
    else{
        file << tab << "<Composite id=\"" << id << "\" Uvalue =\"" << Uvalue << "\"/>" << endl;
    }
}

void Composite::writeGML(ofstream& file, string tab=""){
    // writes the opaque construction op_#id
     file << tab << "<energy:Construction gml:id=\"op_" << id << "\">" << endl;
     file << tab << tabs(1) << "<energy:uValue uom=\"W/(m2K)\">" << this->getUvalue() << "</energy:uValue>" << endl;
     if(!vLayer.empty()){
         for (size_t i=0; i<vLayer.size(); ++i)
             vLayer.at(i).writeGML(file, tab+tabs(1));
     }
     file << tab << "</energy:Construction>" << endl;
}

bool Composite::equals(const Composite& c) const {
    if (c.name != name) return false;
    //if (c->id != id) return false;
    if (c.Uvalue != Uvalue) return false;
    if (vLayer.size() != c.vLayer.size()) return false;
    for (unsigned int i = 0; i<vLayer.size(); ++i){
        if (c.vLayer[i] != vLayer[i]) return false;
        if ((&(c.vLayer[i])==c.insulationLayer) != (&(vLayer[i])==insulationLayer)){
            cout << "Composites " << name << " and " << c.name << " have a different insulation layer !!!" << endl;
            return false;
        }
    }
    return true;
}

const string Surface::STypeNames[]={"Wall","Roof","Floor","Ground","Surface"};

Surface::Surface(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr):b(pBuilding),logStream(std::cout.rdbuf()) {

    // logStream is directed by default to the "cout" streambuf
    if(pLogStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    id = to<uint64_t>(hdl.ToElement()->Attribute("id"));

    if (hdl.ToElement()->Attribute("ShortWaveReflectance"))
        shortWaveReflectance = to<float>(hdl.ToElement()->Attribute("ShortWaveReflectance"));

    if (hdl.ToElement()->Attribute("GlazingRatio"))
        glazingRatio = to<float>(hdl.ToElement()->Attribute("GlazingRatio"));

    if (hdl.ToElement()->Attribute("GlazingGValue"))
        glazingGvalue = to<float>(hdl.ToElement()->Attribute("GlazingGValue"));

    if (hdl.ToElement()->Attribute("GlazingUValue"))
        glazingUvalue = to<float>(hdl.ToElement()->Attribute("GlazingUValue"));

    // read the openableRatio or if not defined, leave it to 0
    if (hdl.ToElement()->Attribute("OpenableRatio"))
        glazingOpenableRatio = to<float>(hdl.ToElement()->Attribute("OpenableRatio"));

    // add the emissivity of present
    if (hdl.ToElement()->Attribute("LongWaveEmissivity"))
        longWaveEmissivity = to<float>(hdl.ToElement()->Attribute("LongWaveEmissivity"));

    // add the vertices
    TiXmlElement* vertexElem;
    vertexElem = hdl.FirstChildElement("V0").ToElement();
    unsigned int vertexIndex = 0;
    while(vertexElem) {
        // check the presence of coordinates in the vertex tag
        if (vertexElem->Attribute("x")==NULL || vertexElem->Attribute("y")==NULL || vertexElem->Attribute("z")==NULL) break;
        // check if not repeated vertex
        if (!vertices.empty()
            && to<float>(vertexElem->Attribute("x"))==vertices.back()[0]
            && to<float>(vertexElem->Attribute("y"))==vertices.back()[1]
            && to<float>(vertexElem->Attribute("z"))==vertices.back()[2])
        {
            logStream << "Warning: In Surface id=" << toString(id) << " repeated vertex V" << vertexIndex << ", removing it." << endl;
        }
        else pushVertex(to<float>(vertexElem->Attribute("x")),to<float>(vertexElem->Attribute("y")),to<float>(vertexElem->Attribute("z")));
        string vertexj = "V" + toString(++vertexIndex);
        vertexElem = hdl.FirstChildElement(vertexj.c_str()).ToElement();
    }

    if (vertices.size()<3) throw(string("Surface id=" + toString(id) + " has only " + toString(vertices.size()) + " vertices"));

    computeNormalAndArea();

    // Get pv and solar thermal information
    TiXmlElement* pvElem = hdl.FirstChildElement("PV").ToElement();
    if (pvElem) {
        PhotoVoltaic* pv = new PhotoVoltaic(hdl.FirstChildElement("PV"));
        setPVPanel(to<float>(pvElem->Attribute("pvRatio")),pv);
    }
    TiXmlElement* pvThermalElem = hdl.FirstChildElement("ST").ToElement();
    if(pvThermalElem){ // Add Solar thermal
        stPanel = new SolarThermal(hdl.FirstChildElement("ST"));
        stRatio = to<float>(pvThermalElem->Attribute("stRatio"));
    }
    // Get pv and solar thermal information
    TiXmlElement* pvtElem = hdl.FirstChildElement("PVT").ToElement();
    if (pvtElem) {
        pvtPanel = new SolarHybrid(hdl.FirstChildElement("PVT"));
        pvtRatio = to<float>(pvtElem->Attribute("pvtRatio"));
    }

    glazingGvalue_hemispherical = computeGlazingGvalueHemispherical();
}

Surface::Surface(unsigned int id, TiXmlElement* pPosList):id(id),logStream(std::cout.rdbuf()) {
    // creates the surface from a posList in text format read in an XML file
    string value;
    stringstream ss;
    float px,py,pz;
    ss << pPosList->GetText();
    while (ss >> value) {
        px = to<float>(value);
        ss >> value;
        py = to<float>(value);
        ss >> value;
        pz = to<float>(value);
        //cout << "Vertex: " << px << " " << py << " " << pz << endl;
        pushVertex(px,py,pz);
    }
    removeRepeatedVertices();
    // checks the quality of the surface
    if (vertexCount() > 2) {
        computeNormalAndArea();
        // check if positive area
        if (getArea() > 0.f && getRadius() > 0.f) {
            // everything is OK
        }
        else {
            throw("Surface id=" + toString(getId()) + " has a too small surface area, not creating it.");
        }
    }
    else {
        throw("Surface id=" + toString(getId()) + " has less than 3 vertices, not creating it.");
    }
}

float Surface::getRadius() const {

    GEOMBoundingSphereCalc boundingSphereCalculator;
    GEOMSphere s=boundingSphereCalculator.Calculate(vertices.begin(),vertices.end());
    return s.Radius();

}

float Surface::getGlazingGvalue(float theta, float p, float q) {

    if (theta >= 90) return 0.f;
    else return glazingGvalue*(1-8*pow((double)(theta/90.f),(double)(5.2+0.7*q))-0.25/q*pow((double)(theta/90.f),(double)2.f)-(1-8-0.25/q)*pow((double)(theta/90.f),(double)((5.26+0.06*p)+(0.73+0.04*p)*q)));

}

float Surface::computeGlazingGvalueHemispherical() {

    // numerical quadrature
    double value = 0.; unsigned int div = 9;
    for (unsigned int i=1;i<div;++i) {
        value += getGlazingGvalue(i*90./div)*cos(i*M_PI/(2.*div))*sin(i*M_PI/(2*div))*(M_PI/(2*div));
    }
    return 2.*value; // integral of g(theta)*cos(theta)*sin(theta) dtheta between 0 and 90∞ with numerical quadrature using triangular approximation with 9 divisions (first one and last one are zero)

}

Wall::Wall(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr):Surface(hdl, pBuilding, pLogStr) {

    // Get the composite
    if (hdl.ToElement()->Attribute("type")) {
        // gets the type according to the value in "type"
        composite = b->getDistrict()->getComposite(hdl.ToElement()->Attribute("type"));
<<<<<<< HEAD

        if (hdl.ToElement()->Attribute("insulationThickness"))
            insulationThickness = to<float>(hdl.ToElement()->Attribute("insulationThickness"));
=======
>>>>>>> cf042395e1dcb11627fd9bf5653a3121bccaadb0
    }
    else if (hdl.ToElement()->Attribute("Uvalue")) {
        // gets or creates a new WallType in the building for this special wall
        composite = b->getDistrict()->getUvalueComposite(to<float>(hdl.ToElement()->Attribute("Uvalue")));
    }
    else throw(string("In Wall id: ")+toString(id)+", no composite type nor U-value given.");

    // specific for the Wall
    if (hdl.ToElement()->Attribute("ep_id")) ep_id = hdl.ToElement()->Attribute("ep_id"); // add the ep_id
    else ep_id = "Wall" + toString(id);
}

Floor::Floor(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr):Surface(hdl, pBuilding, pLogStr) {

    // Get the composite
    if (hdl.ToElement()->Attribute("type")) {
        // gets the type according to the value in "type"
        composite = b->getDistrict()->getComposite(hdl.ToElement()->Attribute("type"));
<<<<<<< HEAD

        if (hdl.ToElement()->Attribute("insulationThickness"))
            insulationThickness = to<float>(hdl.ToElement()->Attribute("insulationThickness"));
=======
>>>>>>> cf042395e1dcb11627fd9bf5653a3121bccaadb0
    }
    else if (hdl.ToElement()->Attribute("Uvalue")) {
        // gets or creates a new WallType in the building for this special wall
        composite = b->getDistrict()->getUvalueComposite(to<float>(hdl.ToElement()->Attribute("Uvalue")));
    }
    else throw(string("In Floor id: ")+toString(id)+", no composite type nor U-value given.");

    // Ensure some parameters are 0 ?
    shortWaveReflectance = 0;
    glazingGvalue = 0;
    glazingUvalue = 0;
    glazingOpenableRatio = 0;
    glazingRatio = 0;

}

Roof::Roof(TiXmlHandle hdl, Building* pBuilding, ostream* pLogStr):Surface(hdl, pBuilding, pLogStr) {

    // Get the composite
    if (hdl.ToElement()->Attribute("type")) {
        // gets the type according to the value in "type"
        composite = b->getDistrict()->getComposite(hdl.ToElement()->Attribute("type"));
<<<<<<< HEAD

        if (hdl.ToElement()->Attribute("insulationThickness"))
            insulationThickness = to<float>(hdl.ToElement()->Attribute("insulationThickness"));
=======
>>>>>>> cf042395e1dcb11627fd9bf5653a3121bccaadb0
    }
    else if (hdl.ToElement()->Attribute("Uvalue")) {
        // gets or creates a new WallType in the building for this special wall
        composite = b->getDistrict()->getUvalueComposite(to<float>(hdl.ToElement()->Attribute("Uvalue")));
    }
    else throw(string("In Roof id: ")+toString(id)+", no composite type nor U-value given.");

    // adds the micro wind turbine if it exists
    TiXmlElement* microWindElem = hdl.FirstChildElement("MicroWindTurbine").ToElement();
    if (microWindElem) microWindTurbine = new MicroWindTurbine(to<float>(microWindElem->Attribute("cutInSpeed")),to<float>(microWindElem->Attribute("ratedSpeed")),
                                                               to<float>(microWindElem->Attribute("cutOutSpeed")),to<float>(microWindElem->Attribute("c1")),
                                                               to<float>(microWindElem->Attribute("c2")), to<float>(microWindElem->Attribute("c3")),
                                                               to<float>(microWindElem->Attribute("airDensity")),to<float>(microWindElem->Attribute("Height")),
                                                               to<float>(microWindElem->Attribute("Alpha")),to<float>(microWindElem->Attribute("Gamma")));
    // adds the centroid information
    centroid = computeCentroid();

    // add the ep_id
    if (hdl.ToElement()->Attribute("ep_id")) ep_id = hdl.ToElement()->Attribute("ep_id");
    else ep_id = "Roof" + toString(id);

    // adds the kFactor for Evapotranspiration
    if (hdl.ToElement()->Attribute("kFactor")) kFactor = to<float>(hdl.ToElement()->Attribute("kFactor"));

}

Roof::~Roof() {
    //cout << "deleting Roof ";
    if (microWindTurbine) delete microWindTurbine;
    //cout << "done."<< endl;
}

Ground::Ground(TiXmlHandle hdl, Composite* c, ostream* pLogStream):Surface(hdl, NULL, pLogStream) {

    // sets the composite
    composite = c;

    // add the kfactor for Evapotranspiration by GU
    if (hdl.ToElement()->Attribute("kFactor")) {
       kFactor = to<float>(hdl.ToElement()->Attribute("kFactor"));
    }

    // adds the detailed simulation attribute
    if (hdl.ToElement()->Attribute("detailedSimulation")) {
       detailedSimulation = (hdl.ToElement()->Attribute("detailedSimulation") == string("true"));
    }

    initialiseModel();

}
