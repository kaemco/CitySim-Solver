#include "zone.h"

#include "district.h"
#include "building.h"
#include "models.h"
#include "occupants.h"
#include "scene.h"

// *** Zone class, CitySim   *** //
// *** jerome.kaempf@epfl.ch *** //

Zone::Zone(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)
    :id(id),pBuilding(pBuilding),walls(walls),roofs(roofs),surfaces(surfaces),floors(floors),Vi(Vi),groundFloor(groundFloor),pOccupants(pOccupants),logStream(pBuilding->logStream.rdbuf()) {

    logStream << "Zone creator." << endl << flush;

    // the zone cannot have a volume of zero
    if (Vi<=0.f) throw(string("In Zone id ") + toString(id) + string(": Cannot create a zone without a positive volume value."));


    for (unsigned int i=0; i<walls.size(); ++i) {
        walls[i]->setBuildingRef(pBuilding);
    }
    for (unsigned int i=0; i<roofs.size(); ++i) {
        roofs[i]->setBuildingRef(pBuilding);
    }
    for (unsigned int i=0; i<floors.size(); ++i) {
        floors[i]->setBuildingRef(pBuilding);
    }
    for (unsigned int i=0; i<surfaces.size(); ++i) {
        surfaces[i]->setBuildingRef(pBuilding);
    }

    update(true);

}

/**
 * @brief Zone::updateSimulationModelParameters: update all "secondary" parameters, i.e.
 *  parameters which are based on primary parameters that can be modified in the interface CitySimPro
 */
void Zone::update(bool constructor){
    // note: the computation of Qsun1 and Qsun2 is made when a call to the method getQsun1 and getQsun2 is made
    if(!constructor){
        Swa=0.;
        Swi=0.;
        SwiO=0.;
        Sro=0.f;
        Kroof = 0.f;
        Kwindow=0.f;
    }
    // loop on the surfaces to get the surface of walls (Swa) and the surface of windows (Swi)
    for (unsigned int i=0; i<walls.size(); ++i) {
        logStream << "Wall #" << i << " Glazing Ratio: " << walls[i]->getGlazingRatio() << "\tArea: " << walls[i]->getArea() << endl << flush;
        Swa  += walls[i]->getWallArea();
        Swi  += walls[i]->getGlazingArea();
        SwiO += walls[i]->getGlazingArea()*walls[i]->getGlazingOpenableRatio();
        Kwindow += walls[i]->getGlazingUvalue()*walls[i]->getGlazingArea();
        walls[i]->updateGlazingGvalueHemispherical();
    }
    for (unsigned int i=0; i<roofs.size(); ++i) {
        logStream << "Roof #" << i << " Glazing Ratio: " << roofs[i]->getGlazingRatio() << "\tArea: " << roofs[i]->getArea() << endl << flush;
        Sro += roofs[i]->getRoofArea();
        Swi += roofs[i]->getGlazingArea();
        SwiO += roofs[i]->getGlazingArea()*roofs[i]->getGlazingOpenableRatio();
        Kroof += roofs[i]->getComposite()->getUvalue()*roofs[i]->getRoofArea();
        Kwindow += roofs[i]->getGlazingUvalue()*roofs[i]->getGlazingArea();
        roofs[i]->updateGlazingGvalueHemispherical();
    }
    logStream << "Total wall surface (m^2): " << Swa << "\tTotal window surface (m^2): " << Swi << "\tTotal roof surface (m^2): " << Sro << endl << flush;
    logStream << "Total window conductance (W/K): " << Kwindow << "\tTotal roof conductance (W/K): " << Kroof << endl << flush;

    // computation of Ww and Wa, as a first approximation, the ratio of the surfaces to the total surface, JK - 18.04.2009
    //if ( (Swa+Swi) > 0. ) { Ww=Swa/(Swa+Swi); Wa=Swi/(Swa+Swi); }
    //logStream << "Ww: " << Ww << "\tWa: " << Wa << endl << flush;

    if (pBuilding->getDistrict()->getScene()->getClimate()!=NULL)
        setAirDensity(pBuilding->getDistrict()->getScene()->getClimate()->getAirDensity());
    else
        setAirDensity(Climate::getAirDensity(0.f)); // if no Climate make an assumption that the altitude of the situation is 0 m (sea level)

    // output of the values
    logStream << "Vi: " << Vi << endl << flush;
}

vector<Surface*> Zone::getAllSurfaces(){
    vector<Surface*> v = surfaces;
    v.insert(v.begin(),walls.begin(),walls.end());
    v.insert(v.begin(),floors.begin(),floors.end());
    v.insert(v.begin(),roofs.begin(),roofs.end());
    return v;
}

void Zone::computeVolume() {

    // initialize the air volume
    Vi = 0.f;
    // computes the volume under the surfaces
    for (vector<Wall*>::iterator it=walls.begin();it!=walls.end();++it) Vi += (*it)->getVolume();
    for (vector<Roof*>::iterator it=roofs.begin();it!=roofs.end();++it) Vi += (*it)->getVolume();
    for (vector<Floor*>::iterator it=floors.begin();it!=floors.end();++it) Vi += (*it)->getVolume();

}

double Zone::getKappa3() {

    double kappa3 = 0.;
    for (size_t i=0; i<roofs.size(); ++i){
        kappa3 += roofs[i]->getKappa()*roofs[i]->get_hr()*roofs[i]->getRoofArea();
        //cout << "Zone::getKappa3: Kappa="<< roofs[i]->getKappa()<<" hr="<<roofs[i]->get_hr() << " area="<<roofs[i]->getRoofArea()<<endl;
    }
    return kappa3;

}


float Zone::getQsun1() {
    // computes Qsun1 (what goes on the enveloppe of the building) from the values on the facades
    float value1 = 0.f;
    for (size_t i=0; i<walls.size(); ++i) { // loop on the external walls
        value1 += ( walls[i]->getWallArea()
                    *(walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                      +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature()) );
    }
    return value1;
}

float Zone::getQsun2() {
    // computes Qsun2 (what goes inside the building with shading devices in use) from the values on the facades
    float value2 = 0.f;
    for (size_t i=0; i<walls.size(); ++i) { // loop on the external walls
        //logStream << "Irradiation: " << walls[i]->getShortWaveIrradiation()/walls[i]->getArea()
        value2 += ( walls[i]->getGlazingRatio()*walls[i]->getArea()*
                    (walls[i]->getGlazingGvalue(walls[i]->getBeamAngle())*walls[i]->getBeamIrradiance()
                     +walls[i]->getGlazingGvalueHemispherical()*(walls[i]->getShortWaveIrradiance()-walls[i]->getBeamIrradiance()))*
                    (0.15f+0.85f*walls[i]->getLowerShadingState()) );

        //cout << "Wall: " << walls[i]->getId() << " Irradiance: " << walls[i]->getShortWaveIrradiation()/walls[i]->getArea() << " shading state: " << walls[i]->getLowerShadingState() << "\tshading: " << (0.15f+0.85f*walls[i]->getLowerShadingState()) << endl;
    }
    for (size_t i=0; i<roofs.size(); ++i) { // loop on the external roofs
        value2 += ( roofs[i]->getGlazingRatio()*roofs[i]->getArea()*
                    (roofs[i]->getGlazingGvalue(roofs[i]->getBeamAngle())*roofs[i]->getBeamIrradiance()
                     +roofs[i]->getGlazingGvalueHemispherical()*(roofs[i]->getShortWaveIrradiance()-roofs[i]->getBeamIrradiance()))*
                    (0.15f+0.85f*roofs[i]->getShadingState()) );

        //cout << "Roof: " << roofs[i]->getId() << " Irradiance: " << roofs[i]->getShortWaveIrradiation()/roofs[i]->getArea() << " shading state: " << roofs[i]->getShadingState() << "\tshading: " << (0.15f+0.85f*roofs[i]->getShadingState()) << endl;
    }
    return value2;
}

float Zone::getQsun3() {
    // computes Qsun3 (what goes on the envelope of the building) from the values on the roofs
    float Qsun3 = 0.f;
    for (size_t i=0; i<roofs.size(); ++i) {
        Qsun3 += roofs[i]->getKappa()*roofs[i]->getRoofArea()
                 *( roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())
                    +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() );
    }
    return Qsun3;
}

double Zone::getTa(unsigned int day, unsigned int hour) {
    // pre-simulation results are not erased, but part of the vector might have been saved already
    return Ta.at((day-1)*24 + hour -1 + pBuilding->getDistrict()->getScene()->getPreTimeStepsSimulated() - pBuilding->getDistrict()->getScene()->getSimulationIndex());
}

void Zone::setOccupantsCountAndActivity(unsigned int day, unsigned int hour)
{
    // initialises the sensible and radiative gains Lc and Lr
    Lc = 0.f;
    Lr = 0.f;

    // two different loops, if stochastic or not
    if (occupantsStochastic) {
        // loop to determine the occupants' count
        occupantsCount = 0.f; // has to be a float due to deterministic procedures
        for (unsigned int i=0; i<static_cast<unsigned int>(round(occupantsNumber)); ++i) {
            if (randomUniform(0,1) <= occupantsYearProfile->getDayProfile(day)->getHourValue(hour)) occupantsCount+=1.f;
        }
        // loop to determine the occupants' behavior
        float activity, randomNumber;
        unsigned int deviceType;
        for (size_t n=0; n < occupantsCount; ++n) {
            activity = 0.f;
            randomNumber = randomUniform(0,1); // draw one random number per occupant to determine the activity performed
            if (activityType != numeric_limits<unsigned int>::signaling_NaN()) { // activity is defined
                for(size_t j=0; j < pBuilding->getDistrict()->getActivityType(activityType)->getnActivities(); ++j) {
                    activity += pBuilding->getDistrict()->getActivityType(activityType)->getActivityProbability(j,hour);
                    if (randomNumber <= activity) {
                        // execute the activity
                        deviceType = pBuilding->getDistrict()->getActivityType(activityType)->getActivityDeviceType(j);
                        // loop on all devices listed in this activity
                        for (size_t k = 0; k < pBuilding->getDistrict()->getDeviceType(deviceType)->getnDevices(); ++k) {
                            if (randomUniform(0,1) <= pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceProbability(k,hour)) {
                                // this device is used, add its electricity consumption to the building
                                pBuilding->addElectricConsumption(float(Model::dt)*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k));
                                // adds the sensible and convective gains due to the device use
                                Lc += pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k)*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceConvectiveFraction(k);
                                Lr += pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k)*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceRadiativeFraction(k);
                            }
                        }
                        // once the activity is accomplished, break the loop on activities
                        break;
                    }
                }
            }
        }

    }
    else {
        // define the occupants count
        occupantsCount = occupantsNumber * occupantsYearProfile->getDayProfile(day)->getHourValue(hour);
        // loop to determine occupants behaviour
        unsigned int deviceType;
        float electricity;
        if (activityType != numeric_limits<unsigned int>::signaling_NaN()) { // activity is defined
            for(size_t j=0; j < pBuilding->getDistrict()->getActivityType(activityType)->getnActivities(); ++j) {
                // loop on all devices listed in this activity
                deviceType = pBuilding->getDistrict()->getActivityType(activityType)->getActivityDeviceType(j);
                for (size_t k = 0; k < pBuilding->getDistrict()->getDeviceType(deviceType)->getnDevices(); ++k) {
                    // this device is used, add its electricity consumption to the building
                    electricity = occupantsCount
                                  *pBuilding->getDistrict()->getActivityType(activityType)->getActivityProbability(j,hour)
                                  *pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceProbability(k,hour)
                                  *pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceAvgPower(k);
                    pBuilding->addElectricConsumption(float(Model::dt)*electricity);
                    // adds the sensible and convective gains due to the device use
                    Lc += electricity*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceConvectiveFraction(k);
                    Lr += electricity*pBuilding->getDistrict()->getDeviceType(deviceType)->getDeviceRadiativeFraction(k);
                }
            }
        }
    }

    // sensible and radiative gains due to the occupants
    Lc += occupantsCount*occupantsSensibleHeat*(1.f-occupantsSensibleHeatRadiantFraction);
    Lr += occupantsCount*occupantsSensibleHeat*occupantsSensibleHeatRadiantFraction;

}

float Zone::getDHWConsumption(unsigned int day, unsigned int hour)
{
    if (occupantsStochastic) {
        float DHWconsumption = 0.f;
        for (unsigned int i=0; i<occupantsCount; ++i) {
            if (randomUniform(0,1) <= dhwYearProfile->getDayProfile(day)->getHourValue(hour)) DHWconsumption+=dhwYearProfile->getDayProfile(day)->getWaterConsumption()*dhwYearProfile->getDayProfile(day)->getHourValue(hour); // Cognet: Deleted the division by 24 and times "getHourValue", so that the expected consumed volume per person after one day is "getWaterConsumption". Basically each person probabilistically consumes the whole "getWaterConsumption" liters in one go. TODO check this.
        }
        return DHWconsumption;
    }
    else return occupantsCount*dhwYearProfile->getDayProfile(day)->getWaterConsumption()*dhwYearProfile->getDayProfile(day)->getHourValue(hour); // Cognet: Deleted the division by 24, so that the probabilities are normalized. TODO: check that this is correct.
}

void Zone::writeXML(ofstream& file, string tab=""){
    streamsize ss = file.precision();
    file.precision(6); // use max 6 significant numbers for floats...
    file << tab << "<Zone id=\"" << id << "\" volume=\"" << Vi << "\" psi=\"" << Kpsi
         << "\" Tmin=\"" << Tmin << "\" Tmax=\"" << Tmax << "\" groundFloor=\"";

    if (groundFloor)
        file << "true";
    else
        file << "false";
    file << "\" nightVentilationBegin=\"" << nightVentilationBegin << "\" nightVentilationEnd=\"" << nightVentilationEnd;
    file << "\">" << endl;
    string subtab =tab+"\t";
    file << subtab << "<Occupants n=\"" << occupantsNumber << "\" sensibleHeat=\"" << getOccupantsSensibleHeat()
                                                           << "\" sensibleHeatRadiantFraction=\"" << getOccupantsSensibleHeatRadiantFraction()
                                                           << "\" latentHeat=\"" << getOccupantsLatentHeat();
    file << "\" type=\"" << occupantsYearProfile->getId() << "\"";
    if (activityType != numeric_limits<unsigned int>::signaling_NaN()) file << " activityType=\"" << activityType << "\"";
    file << " DHWType=\"" << getDHWYearProfile()->getId() << "\"";
    file << "/>" << endl;
    file.precision(ss); // restore default

    for (unsigned int i=0; i<walls.size(); ++i){
        walls[i]->writeXML(file,subtab);
    }
    for (unsigned int i=0; i<roofs.size(); ++i){
        roofs[i]->writeXML(file,subtab);
    }
    for (unsigned int i=0; i<floors.size(); ++i){
        floors[i]->writeXML(file,subtab);
    }
    // Write Shading Surfaces inside the building
    for (size_t i=0; i<surfaces.size(); ++i) {
        file << subtab << "<Surface id=\"" << surfaces[i]->getId();
        streamsize ss = file.precision();
        file.precision(6); // for fixed format, two decimal p
        file << "\" ShortWaveReflectance=\"" << surfaces[i]->getShortWaveReflectance();
        file.precision(ss); // restore default
        file << "\">" << endl;
        surfaces[i]->writeXML(file,subtab+"\t");
        file << subtab << "</Surface>" << endl;
    }

    file << tab << "</Zone>" << endl;
}

// Zone2N
Zone2N::Zone2N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)
        :Zone(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {

    // intialisation
    nNodes=2;
    Tw=15.f;

    update(true);

}

/**
 * @brief Zone2N::updateSimulationModelParameters: update all "secondary" parameters, i.e.
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro
 */
void Zone2N::update(bool constructor){
    if(!constructor){
        Zone::update();
    }
    // If constructor, this was already called when the parent Zone was initialised

    if (walls.empty()) { // no external wall in this zone, temperature decoupled from the model
        Cw = 1.f;
        Kw1 = 1.f;
        Kw2 = 1.f;
    }
    else {
        walls[0]->getComposite()->getSimplifiedNode(Cw, Kw1, Kw2);
        // multiply with the total walls area
        Cw *= Swa;
    }

    // outputs the values
    //logStream << "Cw: " << Cw << "\tKw1: " << Kw1*Swa << "\tKw2: " << Kw2*Swa << endl;
}

double Zone2N::getKappa1() {

    // if no walls, returns a dummy value to couple with the outdoor environment
    if (walls.empty()) return 1.;

    double kappa1 = 0.;
    for (size_t i=0; i<walls.size(); ++i)
        kappa1 += Kw1*walls[i]->getWallArea()*(walls[i]->get_hc()+walls[i]->get_hr())/(Kw1+walls[i]->get_hc()+walls[i]->get_hr());
    return kappa1;

}

void Zone2N::setTos(float Tout) {
    // define the surface temperature per surface
    for (size_t i=0; i<walls.size(); ++i) {
        float wallTemperature = ( Kw1*Tw + walls[i]->get_hc()*Tout
                                  +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                                  +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature() )
                                / (Kw1 + walls[i]->get_hc() + walls[i]->get_hr());
        walls[i]->setTemperature(wallTemperature);
    }
    for (size_t i=0; i<roofs.size(); ++i) {
        float roofTemperature = ( roofs[i]->getKr()*Ta.back() + roofs[i]->get_hc()*Tout
                                  +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())
                                  +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() )
                                / (roofs[i]->getKr() + roofs[i]->get_hc() + roofs[i]->get_hr());
        roofs[i]->setTemperature(roofTemperature);
    }
}

// Zone3N
Zone3N::Zone3N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)
        :Zone2N(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {

    // intialisation
    nNodes=3;
    Tr=15.f;

    update(true);
}

/**
 * @brief Zone3N::updateSimulationModelParameters: update all "secondary" parameters, i.e.
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro
 */
void Zone3N::update(bool constructor){
    if(!constructor){
        Zone2N::update();
    }
    // If constructor, this was already called when the parent Zone2N was initialised

    if (roofs.empty()) { // no roofs in this zone, temperature decoupled from the model
        Cr = 1.f;
        Kr1 = 1.f;
        Kr2 = 1.f;
    }
    else {
        roofs[0]->getComposite()->getSimplifiedNode(Cr, Kr1, Kr2);
        // multiply with the total roofs area
        Cr *= Sro;
    }

    // outputs the values
    //logStream << "Cr: " << Cr << "\tKr1: " << Kr1*Sro << "\tKr2: " << Kr2*Sro << endl;
}

double Zone3N::getKappa3() {

    // if no roofs, returns a dummy value to couple with the outdoor environment (avoid singular matrix)
    if (roofs.empty()) return 1.;

    double kappa3 = 0.;
    for (size_t i=0; i<roofs.size(); ++i){
        kappa3 += Kr1*roofs[i]->getRoofArea()*(roofs[i]->get_hc()+roofs[i]->get_hr()+roofs[i]->get_X())/(Kr1+roofs[i]->get_hc()+roofs[i]->get_hr()+roofs[i]->get_X());
        //cout << "Zone3N::getKappa3: hc="<< roofs[i]->get_hc()<<" hr="<<roofs[i]->get_hr() << " area="<<roofs[i]->getRoofArea()<<endl<<flush;
    }
    return kappa3;

}

void Zone3N::setTos(float Tout) {
    // define the surface temperature per surface
    for (size_t i=0; i<walls.size(); ++i) {
        float wallTemperature = ( Kw1*Tw + walls[i]->get_hc()*Tout
                                  +walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                                  +walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature() )
                                / (Kw1 + walls[i]->get_hc() + walls[i]->get_hr());
        if (wallTemperature>100) {
            logStream << "Wall: " + toString(walls[i]->getId()) + "(" + toString(walls[i]->getKey()) + "), Temperature: " + toString(wallTemperature)
                      << "\tTw: " << Tw << "\tTout: " << Tout << "\tIrradiance: " << walls[i]->getShortWaveIrradiance() << "\tIrradiance absorbed: " << walls[i]->getShortWaveIrradiance()*(1.-walls[i]->getShortWaveReflectance_opaque())
                      << "\tEnvironmental Temperature: " << walls[i]->get_hr()*walls[i]->getEnvironmentalTemperature() << endl;
        }
        walls[i]->setTemperature(wallTemperature);
    }
    for (size_t i=0; i<roofs.size(); ++i) {
        float roofTemperature = ( Kr1*Tr + roofs[i]->get_hc()*Tout
                                  +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())
                                  +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature()
                                  -roofs[i]->get_Y() )
                                / (Kr1 + roofs[i]->get_hc() + roofs[i]->get_hr() + roofs[i]->get_X());
        if (roofTemperature>100) {
            logStream << "Roof: " + toString(roofs[i]->getId()) + "(" + toString(roofs[i]->getKey()) + "), Temperature: " + toString(roofTemperature)
                      << "Tr: " << Tr << "\tTout: " << Tout << "\tIrradiance: " << roofs[i]->getShortWaveIrradiance() << "\tIrradiance absorbed: " << roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())
                      << "\tEnvironmental Temperature: " << roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() << endl;
        }
        roofs[i]->setTemperature(roofTemperature);
    }
}

// Zone3N
Zone3N_floor::Zone3N_floor(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)
        :Zone2N(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {

    // intialisation
    nNodes=3;
    Tf=15.f;

    update(true);

}

/**
 * @brief Zone3N_floor::updateSimulationModelParameters: update all "secondary" parameters, i.e.
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro
 */
void Zone3N_floor::update(bool constructor){
    if(!constructor){
        Zone2N::update();
    }
    // If constructor, this was already called when the parent Zone2N was initialised

    if (floors.empty()) { // no floor in this zone
        throw(string("Creation of a thermal zone without floor using the four nodes model."));
        //Cw = 1.;
        //Kw1 = 1.;
        //Kw2 = 1.;
        //Ki = 0.; // the computation of the wall temperature is completely disconnected from the air node temperature
    }
    else floors[0]->getComposite()->getSimplifiedNode(Cf, Kf1, Kf2);

    // multiply with the total floor area
    Sf = 0.f;
    for (size_t i=0; i<floors.size(); ++i) Sf += floors[i]->getArea();
    Cf *= Sf;

    // outputs the values
    //logStream << "Cf: " << Cf << "\tKf1: " << Kf1*Sf << "\tKf2: " << Kf2*Sf << endl;
}

// Zone4N
Zone4N::Zone4N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)
        :Zone3N(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {

    logStream << "Zone4N with 4 thermal nodes." << endl;

    // intialisation
    nNodes=4;
    Tf=15.f;

    update(true);
}

Zone4N::Zone4N(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, Surface* floor, float elevation)
    :Zone3N(id,pBuilding,groundFloor,Vi,vector<Wall*>(),vector<Roof*>(),vector<Surface*>(),vector<Floor*>(),NULL){
    // intialisation
    nNodes=4;
    Tf=15.f;

    floors.push_back(new Floor(*floor));
    roofs.push_back(new Roof(*floor, elevation)); // will get floor id + 1
    vector<GENPoint>* vertices = floor->getVertices();
    for(unsigned int i=0; i<vertices->size()-1; ++i){
        walls.push_back(new Wall(floor->getId()+1+i,(*vertices)[i],(*vertices)[i+1],elevation));
    }
    walls.push_back(new Wall(floor->getId()+vertices->size(),(*vertices)[vertices->size()-1],(*vertices)[0],elevation));
    delete floor;
}

Zone4N::Zone4N(Building* pBuilding, bool groundFloor):Zone3N(pBuilding->getId(),pBuilding,groundFloor){
    // Incomplete constructor for DXF reading, do not use without completing the geometry...
    nNodes=4;
    Tf=15.f;
}

void Zone4N::addSurface(Surface* s){
    //cout << "Zone::addSurface, normal altitude: " << s->normal().Altitude().degrees() << endl;
    if(s->normal().Altitude().degrees() > 1.){ // roof
        //cout << "Zone::addSurface -> roof" << endl;
        roofs.push_back(new Roof(*s));
    }
    else if(s->normal().Altitude().degrees() < -1.){ // floor
        //cout << "Zone::addSurface -> floor" << endl;
        floors.push_back(new Floor(*s));
    }
    else{ // all other surfaces added through this function are considered walls
        //cout << "Zone::addSurface -> wall" << endl;
        walls.push_back(new Wall(*s));
    }
    delete s;
}

/**
 * @brief Zone4N::updateSimulationModelParameters: update all "secondary" parameters, i.e.
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro
 */
void Zone4N::update(bool constructor){
    if(!constructor){
        Zone3N::update();
    }
    // If constructor, this was already called when the parent Zone3N was initialised

    if (floors.empty()) { // no floor in this zone, temperature decoupled from the model
        Cf = 1.f;
        Kf1 = 1.f;
        Kf2 = 1.f;
    }
    else  {
        floors[0]->getComposite()->getSimplifiedNode(Cf, Kf1, Kf2);
        // multiply with the total floor area
        Sf = 0.f;
        for (size_t i=0; i<floors.size(); ++i) Sf += floors[i]->getArea();
        Cf *= Sf;
    }

    // outputs the values
    //logStream << "Cf: " << Cf << "\tKf1: " << Kf1*Sf << "\tKf2: " << Kf2*Sf << endl;
}

double Zone4N::getKappa5() {

    // if no floors, returns a dummy value to couple with the outdoor environment (avoid singular matrix)
    if (floors.empty()) return 1.;
    else return Kf1*Sf;

}

// ZoneN
ZoneN::ZoneN(unsigned int id, Building* pBuilding, bool groundFloor, float Vi, vector<Wall*> walls, vector<Roof*> roofs, vector<Surface*> surfaces, vector<Floor*> floors, Occupants *pOccupants)
        :Zone(id,pBuilding,groundFloor,Vi,walls,roofs,surfaces,floors,pOccupants) {

    update(true);
}

/**
 * @brief ZoneN::updateSimulationModelParameters: update all "secondary" parameters, i.e.
 *  parameters which are based on other parameters that can be modified in the interface CitySimPro
 */
void ZoneN::update(bool constructor){
    if(!constructor){
        Zone::update();
    }
    // If constructor, this was already called when the parent ZoneN was initialised

    // check the potential errors
    if (walls.empty()) { // no external wall in this zone
        throw(string("Creation of a thermal zone without external walls using the n-nodes model."));
        //Cw = 1.;
        //Kw1 = 1.;
        //Kw2 = 1.;
        //Ki = 0.; // the computation of the wall temperature is completely disconnected from the air node temperature
    }
    else {
        // intialisation
        nNodes=walls[0]->getComposite()->getnLayers()+1; // plus one node for the internal air
    }

    Tw.assign(nNodes-1, 15.f);
    if (Model::thermalExplicit) TwExpl.assign(nNodes-1, 15.f);

    // initialisation of the Ki (using a fixed value of hc=3.0)
    Ki=3.0*Swa;
}

void ZoneN::setTos(float Tout) {
    // clamp the wall temperature, at minimum the outside air temperature, at maximum 80∞C (upper limit value taken from CIBSE guide)
    //float wallTemperature = min(max((Kw1*(Tw.back()) + getKe()*(Tout) + getQsun1())/(Kw1 + getKe()), Tout), 80.);
    float wallTemperature = (getG0()*Tw[0] + getKe()*Tout + getQsun1())/(getG0() + getKe() + getHr());
    for (size_t i=0; i<walls.size(); ++i) {
        // sets the same wall temperature for all walls as only one temperature for the walls
        walls[i]->setTemperature(wallTemperature);
    }
    for (size_t i=0; i<roofs.size(); ++i) {
        // clamp the roof temperature, at minimum the outside air temperature, at maximum 80∞C (upper limit value taken from CIBSE guide)
        //float roofTemperature = min(max(( Kr*Ta.back() + (getKe()/getSwa())*Tout + getQsun3(i)/roofs[i]->getRoofArea())/(Kr + 23.f), Tout), 80.);
        float roofTemperature = ( roofs[i]->getKr()*Ta.back() + roofs[i]->get_hc()*Tout
                                  +roofs[i]->getShortWaveIrradiance()*(1.-roofs[i]->getShortWaveReflectance_opaque())
                                  +roofs[i]->get_hr()*roofs[i]->getEnvironmentalTemperature() )
                                / (roofs[i]->getKr() + roofs[i]->get_hc() + roofs[i]->get_hr());
        roofs[i]->setTemperature(roofTemperature);
    }
}
