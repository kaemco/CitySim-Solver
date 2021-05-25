#include "occupants.h"

#include <limits>
#include <cassert>

#include "util.h"

// *** Occupants class, CitySim  *** //
// *** jerome.kaempf@epfl.ch     *** //

DayProfile OccupancyProfiles::emptyDay(0,"empty day profile",vector<float>(365,0));
YearProfile OccupancyProfiles::emptyYear = YearProfile(0,"empty year profile",vector<DayProfile*>(365,&OccupancyProfiles::emptyDay));

StochasticOccupantsPresence::StochasticOccupantsPresence(float occupantsNumber, StochasticPresenceParameters *presParam, StochasticWindowParameters *winParam, StochasticBlindsParameters *blindsParam, StochasticLightsParameters *lightsParam) : Occupants(occupantsNumber)
{
    // copy the pointers
    this->presParam=presParam;
    this->winParam=winParam;
    this->blindsParam=blindsParam;
    this->lightsParam=lightsParam;

    cerr << "Creation of the stochastic profile." << endl;

    // creation of the average presence profile
    profile.assign((Model::dt/Model::dt2)*24*365,0.f);
    for (unsigned int personIndex=0; personIndex < occupantsNumber; ++personIndex) {
        vector<bool> pres = presence();
        for (size_t timeIndex=0;timeIndex<pres.size();++timeIndex) {
            cerr << "profile.size(): " << profile.size() << "\tpres.size(): " << pres.size() << endl;
            assert(profile.size()==pres.size());
            if (pres[timeIndex]) profile[timeIndex] += 1.f/occupantsNumber;
        }
    }
    cerr << "Stochastic profile created: " << profile.size() << " elements." << endl;
    // saves the average profile
    save("averagePresenceProfile.txt",profile);

    // creates from the presence profile the current duration and future duration vectors
    current_dur.push_back(0.f);
    for (unsigned int i=1; i<profile.size(); ++i) {
        if ( profile[i] == profile[i-1] ) current_dur.push_back(current_dur[i-1] + Model::dt2);
        else current_dur.push_back(0.f);
    }

    future_dur.assign(profile.size(), std::numeric_limits<float>::quiet_NaN());
    for (unsigned int i=1; i<profile.size(); ++i) {
        if (profile[i] == 1 && profile[i+1]==0) {
            unsigned int j = i+1;
            while (profile[j]==0 && j < profile.size()) ++j;
            future_dur[i] = current_dur[j-1];
        }
    }

}

float StochasticOccupantsPresence::getOccupantsFraction(unsigned int day, unsigned int hour, int fracHour)
{
    return profile[ (fracHour < 0 && day == 1 && hour == 1) ?
                    static_cast<int>(Model::dt/Model::dt2)*((static_cast<int>(day)-1)*24 + static_cast<int>(hour) -1) + fracHour + static_cast<int>(profile.size()) :
                    Model::dt/Model::dt2*((day-1)*24 + hour -1) + fracHour ];
}

float StochasticOccupantsPresence::getCurrentDuration(unsigned int day, unsigned int hour, unsigned int fracHour)
{
    return current_dur[ Model::dt/Model::dt2*((day-1)*24 + hour -1) + fracHour ];
}

float StochasticOccupantsPresence::getFutureDuration(unsigned int day, unsigned int hour, unsigned int fracHour)
{
    return future_dur[ Model::dt/Model::dt2*((day-1)*24 + hour -1) + fracHour ];
}

double StochasticOccupantsPresence::getT01(double pcurr, double pnext, double shuff) {

    // This function returns the transition probabilities T01
    // pcurr: current step occupancy probability
    // pnext: next step occupancy probability
    // shuff: shuffling parameter
    // beta:  adjusted value of shuff

    double T01; // Probability to leave the space
    double beta=shuff; // default: no adjustment needed

    if (pnext == 0.) T01=0.;
    else if (pnext == 1.) T01=1.;
    else {
        if (pcurr == 1.) T01=0.;
        else if (pcurr == 0.) T01=pnext;
        else if (pcurr == pnext) {
            if (pcurr+pnext>1.) {
                if (shuff>1./(2.*pcurr-1.)) beta=1./(2.*pcurr-1.);
                else beta=shuff;
            }
            else if (pcurr+pnext<1.) {
                if (shuff>1./(1.-2.*pcurr)) beta=1./(1.-2.*pcurr);
                else beta=shuff;
            }
            else beta=shuff;

            T01=2.*beta*pcurr/(beta+1.);
        }
        else if (pcurr<pnext) {
            if (shuff < (pnext-pcurr)/(2.-(pnext+pcurr))) beta=(pnext-pcurr)/(2.-(pnext+pcurr));
            else {
                if ((pcurr+pnext>1.) && (shuff>(pcurr-pnext+1.)/(pnext+pcurr-1.)))
                    { beta=(pcurr-pnext+1.)/(pnext+pcurr-1.); }
                else if ((pcurr+pnext<1.) && (shuff>(1.-pcurr+pnext)/(1.-pcurr-pnext)))
                    { beta=(1.-pcurr+pnext)/(1.-pcurr-pnext); }
                else beta=shuff;
            }
            T01=pnext+pcurr*(beta-1.)/(beta+1.);
        }
        else { // Case of (pcurr>pnext)
            if (shuff<(pcurr-pnext)/(pnext+pcurr)) beta=(pcurr-pnext)/(pnext+pcurr);
            else {
                if ((pcurr+pnext>1.) && (shuff>(pcurr-pnext+1.)/(pnext+pcurr-1.)))
                    { beta=(pcurr-pnext+1.)/(pnext+pcurr-1.); }
                else if ((pcurr+pnext<1.) && (shuff>(1.-pcurr+pnext)/(1.-pcurr-pnext)))
                    { beta=(1.-pcurr+pnext)/(1.-pcurr-pnext); }
                else beta=shuff;
            }
            T01=pnext+pcurr*(beta-1.)/(beta+1.);
        }
    }

    return T01;

}

double StochasticOccupantsPresence::getT11(double pcurr, double pnext, double shuff) {

    // This function returns the transition probabilities T01 and T11
    // pcurr: current step occupancy probability
    // pnext: next step occupancy probability
    // shuff: shuffling parameter
    // beta:  adjusted value of shuff

    double T11; // Probability to stay in the space
    double beta = shuff; // default: no adjustment needed

    if (pnext == 0.) { T11=0.; }
    else if (pnext == 1.) { T11=1.; }
    else {
        if (pcurr == 1.) { T11=pnext; }
        else if (pcurr == 0.) { T11=0.; }
        else if (pcurr == pnext) {
            if (pcurr+pnext>1.) {
                if (shuff>1./(2.*pcurr-1.)) { beta=1./(2.*pcurr-1.); }
                else beta=shuff;
            }
            else if (pcurr+pnext<1.) {
                if (shuff>1./(1.-2.*pcurr)) { beta=1./(1.-2.*pcurr); }
                else { beta=shuff; }
            }
            else { beta=shuff; }

            T11=1.-(1.-pcurr)*getT01(pcurr,pnext,beta)/pcurr;
        }
        else if (pcurr<pnext) {
            if (shuff < (pnext-pcurr)/(2.-(pnext+pcurr))) { beta=(pnext-pcurr)/(2.-(pnext+pcurr)); }
            else {
                if ((pcurr+pnext>1.) && (shuff>(pcurr-pnext+1.)/(pnext+pcurr-1.)))
                    { beta=(pcurr-pnext+1.)/(pnext+pcurr-1.); }
                else if ((pcurr+pnext<1.) && (shuff>(1.-pcurr+pnext)/(1.-pcurr-pnext)))
                    { beta=(1.-pcurr+pnext)/(1.-pcurr-pnext); }
                else { beta=shuff; }
            }
            T11=1./pcurr*(pnext-(1.-pcurr)*getT01(pcurr,pnext,beta));
        }
        else {// Case of (pcurr>pnext)
            if (shuff<(pcurr-pnext)/(pnext+pcurr))
                { beta=(pcurr-pnext)/(pnext+pcurr); }
            else {
                if ((pcurr+pnext>1.) && (shuff>(pcurr-pnext+1.)/(pnext+pcurr-1.)))
                    { beta=(pcurr-pnext+1.)/(pnext+pcurr-1.); }
                else if ((pcurr+pnext<1.) && (shuff>(1.-pcurr+pnext)/(1.-pcurr-pnext)))
                    { beta=(1.-pcurr+pnext)/(1.-pcurr-pnext); }
                else { beta=shuff; }
            }
            T11=1./pcurr*(pnext-(1.-pcurr)*getT01(pcurr,pnext,beta));
        }
    }

    return T11;

}

vector<bool> StochasticOccupantsPresence::presence() {

    // Model for the prediction of presence derived by J. Page
    // Reference: J. Page, D. Robinson, N. Morel, J.-L. Scartezzini, A generalised stochastic
    // model for the simulation of occupant presence, Energy and Buildings 40(2), 83-98 (2008).

    vector<bool> occ;
    occ.push_back(false);
    double shuff=0.11; // Mean observed value for mobility parameter

    double LongAbsCurrentDuration=0.;

    for (unsigned int day=1;day<=365;++day) {
        // determination of the day of the week
        unsigned int dayOfTheWeek = (day-1) % 7;
        for (unsigned int hour=1;hour<=24;++hour) {

            double pHour = 0.;
            double pNextHour= 0.;

            if (dayOfTheWeek == 0)      { pHour=presParam->getpMon(hour-1); pNextHour=((hour==24) ? presParam->getpTue(0) : presParam->getpMon(hour)); }
            else if (dayOfTheWeek == 1) { pHour=presParam->getpTue(hour-1); pNextHour=((hour==24) ? presParam->getpWed(0) : presParam->getpTue(hour)); }
            else if (dayOfTheWeek == 2) { pHour=presParam->getpWed(hour-1); pNextHour=((hour==24) ? presParam->getpThu(0) : presParam->getpWed(hour)); }
            else if (dayOfTheWeek == 3) { pHour=presParam->getpThu(hour-1); pNextHour=((hour==24) ? presParam->getpFri(0) : presParam->getpThu(hour)); }
            else if (dayOfTheWeek == 4) { pHour=presParam->getpFri(hour-1); pNextHour=((hour==24) ? presParam->getpSat(0) : presParam->getpFri(hour)); }
            else if (dayOfTheWeek == 5) { pHour=presParam->getpSat(hour-1); pNextHour=((hour==24) ? presParam->getpSun(0) : presParam->getpSat(hour)); }
            else                        { pHour=presParam->getpSun(hour-1); pNextHour=((hour==24) ? presParam->getpMon(0) : presParam->getpSun(hour)); }

            for (unsigned int fracHour=0;fracHour<Model::dt/Model::dt2;++fracHour) {

                // probabilities of the current fracHour and the next fracHour
                double pcurr = ( (static_cast<double> (Model::dt/Model::dt2-fracHour)) * pHour
                                + (static_cast<double> (fracHour)) * pNextHour ) / (static_cast<double> (Model::dt/Model::dt2));
                double pnext = ( (static_cast<double> (Model::dt/Model::dt2-(fracHour+1))) * pHour
                                + (static_cast<double> (fracHour+1)) * pNextHour ) / (static_cast<double> (Model::dt/Model::dt2));

                double ProbLongAbsence=0.001; // To be adjusted later, lack of calibration data.

                // --- 1. If a long absence is ongoing -----------------------
                if (LongAbsCurrentDuration > 0.) {
                    LongAbsCurrentDuration = max(LongAbsCurrentDuration-1., 0.);
                    occ.push_back(false);
                }

                // --- 2. If there is no long absence ------------------------
                else {
                    // --- 2a. If a long absence starts ------------
                    if (randomUniform(0.f,1.f) < ProbLongAbsence)
                        { occ.push_back(false); LongAbsCurrentDuration = 1000.; } // 1000 constante arbitraire, longueur de longue absence dans XML

                    // --- 2b. If a long absence does not start ----
                    else {

                        // Room currently not occupied
                        if (!occ.back()) {
                            if ( randomUniform(0.f,1.f) < getT01(pcurr,pnext,shuff) ) { occ.push_back(true); }
                            else { occ.push_back(false); }
                        }
                        // Room currently occupied
                        else {
                            if ( randomUniform(0.f,1.f) < (1.-getT11(pcurr,pnext,shuff)) ) { occ.push_back(false); }
                            else { occ.push_back(true); }
                        }
                    }
                }
            }
        }
    }

    // remove the first element of the vector (which was just for starting the process)
    occ.erase(occ.begin());

    return occ;

}
