#include "plant.h"
#include "district.h" //added by Dapeng

// *** Plant class, CitySim  *** //
// *** jerome.kaempf@epfl.ch *** //

#pragma GCC diagnostic ignored "-Wunused-parameter"


PhotoVoltaic::PhotoVoltaic(TiXmlHandle hdl, ostream* pLogStr):logStream(std::cout.rdbuf()) {

    // logStream is directed by default to the "cout" streambuf
    if(pLogStr)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    //logStream << "PhotoVoltaic constructor from xml " << endl << flush;

    if (hdl.ToElement()->Attribute("Pmp") && hdl.ToElement()->Attribute("Ac")) etampref = to<float>(hdl.ToElement()->Attribute("Pmp"))/(to<float>(hdl.ToElement()->Attribute("Ac"))*1000.0);
    else if (hdl.ToElement()->Attribute("Etampref")) etampref = to<float>(hdl.ToElement()->Attribute("Etampref"));
    else cout << "ERROR PhotoVoltaic: Etampref parameter not defined (nor Pmp & Ac)";

    if (hdl.ToElement()->Attribute("Tref")) tref = to<float>(hdl.ToElement()->Attribute("Tref"));
    else cout << "ERROR PhotoVoltaic: Tref parameter not defined";

    if (hdl.ToElement()->Attribute("Tcnoct")) tcnoct = to<float>(hdl.ToElement()->Attribute("Tcnoct"));
    else cout << "ERROR PhotoVoltaic: Tcnoct parameter not defined";

    if (hdl.ToElement()->Attribute("muVoc")) muvoc = to<float>(hdl.ToElement()->Attribute("muVoc"));
    else cout << "ERROR PhotoVoltaic: muVoc parameter not defined";

    if (hdl.ToElement()->Attribute("Vmp")) vmp = to<float>(hdl.ToElement()->Attribute("Vmp"));
    else cout << "ERROR PhotoVoltaic: vmp parameter not defined";

    if (hdl.ToElement()->Attribute("name")){ hdl.ToElement()->QueryStringAttribute("name",&name);}

    if (hdl.ToElement()->Attribute("IAM")) iam = to<float>(hdl.ToElement()->Attribute("IAM"));

}

void PhotoVoltaic::writeXML(ofstream& file, float ratio, string tab){
    file << tab << "<PV pvRatio=\"" << ratio << "\" name=\"" << name <<
            "\" Etampref=\"" << etampref << "\" Vmp=\"" << vmp << "\" muVoc=\"" << muvoc <<
            "\" Tcnoct=\"" << tcnoct << "\" Tref=\"" << tref << "\"/>" << endl;
}

double PhotoVoltaic::getMaxPowerEfficiency(double gt, double tout) {
    //cout << "getMaxPowerEfficiency: etampref=" << etampref <<" muvoc="<<muvoc <<" /vmp="<<vmp<<" tout="<<tout << " gt="<< gt <<" tcnoct="<<tcnoct <<" toutsoc="<<toutsoc <<" /gtsoc=" <<gtsoc << " tref=" << tref << endl;
    return etampref*(1.0 + (muvoc/vmp)*(tout + gt*(tcnoct-toutsoc)/gtsoc - tref))/(1.0 + etampref*gt*(tcnoct-toutsoc)/gtsoc*(1.0/0.9)*(muvoc/vmp));
}

SolarThermal::SolarThermal(TiXmlHandle hdl, ostream* pLogStr):logStream(std::cout.rdbuf()) {

    // logStream is directed by default to the "cout" streambuf
    if(pLogStr)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    //logStream << "SolarHeater constructor from xml " << endl << flush;

    if (hdl.ToElement()->Attribute("eta0")) eta0 = to<float>(hdl.ToElement()->Attribute("eta0"));
    else if (hdl.ToElement()->Attribute("etaOptical")) eta0 = to<float>(hdl.ToElement()->Attribute("etaOptical"));
    else throw("Solar Thermal: eta0 parameter not defined");

    if (hdl.ToElement()->Attribute("a1")) a1 = to<float>(hdl.ToElement()->Attribute("a1"));
    else if (hdl.ToElement()->Attribute("heatLoss1")) a1 = to<float>(hdl.ToElement()->Attribute("heatLoss1"));
    else throw("Solar Thermal: a1 parameter not defined");

    if (hdl.ToElement()->Attribute("a2")) a2 = to<float>(hdl.ToElement()->Attribute("a2"));
    else if (hdl.ToElement()->Attribute("heatLoss2")) a2 = to<float>(hdl.ToElement()->Attribute("heatLoss2"));
    else throw("Solar Thermal: a2 parameter not defined");

    if (hdl.ToElement()->Attribute("name")){ hdl.ToElement()->QueryStringAttribute("name",&name);}
}

void SolarThermal::writeXML(ofstream& file, float ratio, string tab){
    file << tab << "<ST stRatio=\"" << ratio << "\" name=\"" << name <<
            "\" eta0=\"" << eta0 << "\" a1=\"" << a1 << "\" a2=\"" << a2 << "\"/>" << endl;
}

SolarHybrid::SolarHybrid(TiXmlHandle hdl, ostream* pLogStr) : PhotoVoltaic(hdl,pLogStr) {

    logStream << "SolarHybrid constructor from xml " << endl << flush;

    if (hdl.ToElement()->Attribute("Pth")) Pth = to<float>(hdl.ToElement()->Attribute("Pth"));
    else cout << "ERROR PhotoVoltaic: Pth parameter not defined";

    // if given by the user, shadow the default value
    if (hdl.ToElement()->Attribute("massFlowRate")) massFlowRate = to<float>(hdl.ToElement()->Attribute("massFlowRate"));

}

float SolarHybrid::getThermalSurfacePowerDensity(float gt, float tout, float windspeed, float Tsky, float Tin) {

    // Convection naturelle hc_n :
    float hc_n = 5.67+3.86*windspeed;                 // [W/m2.K]

    // define the ground temperature (back of the panel)
    float Tground = tout + 2.f;

    #ifdef DEBUG
        ofstream file("debug_PVT.out", ios::app);
        file << gt << "\t" << tout << "\t" << hc_n << "\t" << Tsky << "\t" << Tground << "\t" << Tin << "\t";
    #endif


    // Rayonnement entre ciel et la vitre solAc hrad_sky :
    float c_boltz = 5.67e-8;
    float hrad_sky = 4.f*c_boltz*epsilon_g*pow((Tg + (Tsky+273.15))/2.,3);

    // Verre solAc / cellules pv :
    // Conduction :
    float hc_g = lambda_g/x_g;
    float hc_air = lambda_air/x_air;
    float hc_eva = lambda_eva/x_eva;

    // Rayonnement entre verre solAc / cellule pv
    float E = 1./(1./x_pv + 1./x_g - 1.);
    float Tmoy = (Tpv + Tg)/2.;
    float hrad_g_pv = 4.*c_boltz*E*pow(Tmoy,3);                     // [W/m2.K]

    // On a 3 résistances en parallèles :
    float h1 = hc_g + hc_air + hrad_g_pv;                     // [W/m2.K]

    // Puis 2 résistances en série :
    float hpv_g = 1./(1./h1 + 1./hc_eva);                               // [W/m2.K]

    // Cellules pv/absorbeur :
    // Conduction :
    float hc_bs = lambda_bs/x_bs;                                       // [W/m2.K]
    float hc_ab = lambda_ab/x_ab;                                       // [W/m2.K]
    float hc_pv = lambda_pv/x_pv;                                       // [W/m2.K]

    float  hab = 1./(1./hc_eva + 1./hc_bs + 1./hc_ab + 1./hc_pv); // [W/m2.K]

    // Absorbeur/air ambiant
    // Pertes derriere l'échangeur radiation
    float hrad_loss = 4.*c_boltz*epsilon_ab*pow((Tmw + (Tground+273.15))/2.,3); // [W/m2.K]

    // Perte totale arrière
    float hloss = hc_n + hrad_loss; // [W/m2.K]

    // Absorbeur/fluide
    // hab_f :
    float hab_f = (eff_ab*hloss)/(1.-eff_ab);         // [W/m2.K]

    float Ac=1.58f;

    // initialise the B vector
    double b[4];
    b[0] = gt*alpha_g + (tout+273.15)*hc_n + (Tsky+273.15f)*hrad_sky;
    b[1] = gt*FF*tau_g*(alpha_pv-etampref+etampref*(muvoc/vmp)*(tref+273.15));
    b[2] = 0.;
    b[3] = ((Tin+273.15f)*(massFlowRate*Cp_w*eta_exch))/Ac + (tout+273.15)*hloss;

    // initialise the A matrix
    double A[4*4];
    A[0+4*0] = hc_n +hrad_sky+hpv_g;
    A[0+4*1] = -hpv_g;
    A[0+4*2] = 0.;
    A[0+4*3] = 0.;

    A[1+4*0] = -hpv_g;
    A[1+4*1] = hpv_g + hab + gt*FF*etampref*tau_g*(muvoc/vmp);
    A[1+4*2] = -hab;
    A[1+4*3] = 0.;

    A[2+4*0] = 0.;
    A[2+4*1] = hab;
    A[2+4*2] = hab_f - hab;
    A[2+4*3] = -hab_f;

    A[3+4*0] = 0.;
    A[3+4*1] = 0.;
    A[3+4*2] = (massFlowRate*Cp_w*eta_exch)/Ac - hab_f;
    A[3+4*3] = hab_f + hloss;

//    cout << "matrice A: " << endl;
//
//    for (int i=0;i<A.nrows();++i) {
//        for (int j=0;j<A.ncols();++j) {
//            cout << A[i][j] << "\t";
//        }
//        cout << endl;
//    }
//
//    cout << "vecteur b: " << endl;
//
//    for (int i=0;i<b.size();++i) {
//            cout << b[i] << "\t";
//    }
//    cout << endl;

    solve_Ax_equal_b(A,b,4); // solution A * x = b put in b

    Tg = b[0];
    Tpv = b[1];
    Tab = b[2];
    Tmw = b[3];

    // Calcul puissance instantannée électrique
    //float P_el = Pmp*(gt/1000.)*(1-(Tpv-tref)*0.0044);  // [W] - 1000 W/m2 STC
    ///float rendement = etampref*(1+(muvoc/vmp)*(Tpv-(tref+273.15f)));  // [W]

    // Avec l'onduleur : P_el_reel = P_el*rend_ond;

    // Calcul puissance instantannée thermique
    float Tout_ab = (2.*Tmw)-(Tin+273.15f);

    float DT = Tout_ab-(Tin+273.15f);
    if (DT<0.) DT=0.;

    #ifdef DEBUG
    file << hloss << "\t" << hab_f << "\t" << Tg << "\t" << Tpv << "\t" << Tab << "\t" << Tmw << "\t" << Tout_ab << "\t" << DT << "\t" << massFlowRate*Cp_w*(DT) << endl;
    file.close();
    #endif

    return massFlowRate*Cp_w*(DT)/Ac;

}

double Equations::solarHeaterEfficiency(SolarThermal *panel, double gt, double xsi) {

    return panel->getEta0() - panel->getA1()*xsi - panel->getA2()*gt*xsi*xsi;

}

double Equations::windTurbinePower(WindTurbine *turbine, double v) {

    if ( v <= turbine->getvi() ) return 0.0;
    else if ( v > turbine->getvm() ) return 0.0;
    else if ( v > turbine->getvr() ) return turbine->getPr();
    else {

        double a = turbine->getPr()*pow(turbine->getvi(), turbine->getc())/(pow(turbine->getvi(), turbine->getc()) - pow(turbine->getvr(), turbine->getc()));
        double b = turbine->getPr()/(pow(turbine->getvr(), turbine->getc()) - pow(turbine->getvi(), turbine->getc()));

        return a + b * pow(v, turbine->getc());

    }

}

double Equations::windSpeedRatio(int type, double height, int typeRef, double heightRef) {

    double alpha, alphaprime, gamma, gammaprime;

    switch ( type ) {

        case 1:

            gamma = 0.10;
            alpha = 1.30;
            break;

        case 2:

            gamma = 0.15;
            alpha = 1.00;
            break;

        case 3:

            gamma = 0.2;
            alpha = 0.85;
            break;

        case 4:

            gamma = 0.25;
            alpha = 0.67;
            break;

        case 5:

            gamma = 0.35;
            alpha = 0.47;
            break;

        default:

            gamma = 0.35;
            alpha = 0.47;
            break;

    }

    switch ( typeRef ) {

        case 1:

            gammaprime = 0.10;
            alphaprime = 1.30;
            break;

        case 2:

            gammaprime = 0.15;
            alphaprime = 1.00;
            break;

        case 3:

            gammaprime = 0.2;
            alphaprime = 0.85;
            break;

        case 4:

            gammaprime = 0.25;
            alphaprime = 0.67;
            break;

        case 5:

            gammaprime = 0.35;
            alphaprime = 0.47;
            break;

        default:

            gammaprime = 0.35;
            alphaprime = 0.47;
            break;

    }

    return (alpha*pow( height/10.0, gamma)) / (alphaprime*pow( heightRef/10.0, gammaprime));


}

double Tank::temperature(double t, double VdotUsed, double Pp2, double Pup2, double T0, double Tinlet, double Tamb) {

    return (Pp2 + Pup2 + getCp()*Tinlet*VdotUsed*getRho() + Tamb*getPhi() -
    (Pp2 + Pup2 - getCp()*T0*VdotUsed*getRho() + getCp()*Tinlet*VdotUsed*getRho() - T0*getPhi() + Tamb*getPhi())/
      std::exp((t*(getCp()*VdotUsed*getRho() + getPhi()))/(getCp()*getVolume()*getRho())))/(getCp()*VdotUsed*getRho() + getPhi());

}

double Tank::power(double t, double Tf, double VdotUsed, double Pup2, double T0, double Tinlet, double Tamb) {

    return (Pup2 - std::exp(((getCp()*VdotUsed*getRho() + getPhi())*(t))/(getCp()*getVolume()*getRho()))*Pup2 - getCp()*T0*VdotUsed*getRho() +
   getCp()*std::exp(((getCp()*VdotUsed*getRho() + getPhi())*(t))/(getCp()*getVolume()*getRho()))*Tf*VdotUsed*getRho() +
   getCp()*Tinlet*VdotUsed*getRho() - getCp()*std::exp(((getCp()*VdotUsed*getRho() + getPhi())*(t))/(getCp()*getVolume()*getRho()))*Tinlet*VdotUsed*
    getRho() - T0*getPhi() + Tamb*getPhi() - std::exp(((getCp()*VdotUsed*getRho() + getPhi())*(t))/(getCp()*getVolume()*getRho()))*Tamb*getPhi() +
   std::exp(((getCp()*VdotUsed*getRho() + getPhi())*(t))/(getCp()*getVolume()*getRho()))*Tf*getPhi())/
  (-1.0 + std::exp(((getCp()*VdotUsed*getRho() + getPhi())*(t))/(getCp()*getVolume()*getRho())));

}

double Tank::time(double Tf, double VdotUsed,double Pp2, double Pup2, double T0, double Tinlet, double Tamb) {

    return (Cp*volume*rho)*(log((Pp2 - phi*T0 - Cp*rho*VdotUsed*T0 + Cp*rho*VdotUsed*Tinlet + Pup2 + phi*Tamb)
                                           /(Pp2 - phi*Tf - Cp*rho*VdotUsed*Tf + Cp*rho*VdotUsed*Tinlet + Pup2 + phi*Tamb)))/(phi + Cp*rho*VdotUsed);

}

double TankPCM::temperature(double t, double VdotUsed, double Pp2, double Pup2, double T0, double Tinlet, double Tamb) {

    // calcul du crit�re de stabilit�

    double deltat = std::min( 300.0, 2.0*(getRho()*getVolume()*getCp() + mass*Cp(T0))/(getPhi() + getRho()*VdotUsed*getCp()) );

    std::cerr << "Temperature - explicit stability: " << 2.0*(getRho()*getVolume()*getCp() + mass*Cp(T0))/getPhi() << std::endl;

    if ( t < deltat ) return T0 + t*(Pp2 + Pup2 + getPhi()*(Tamb - T0) + getRho()*VdotUsed*getCp()*(Tinlet - T0))/(getRho()*getVolume()*getCp() + mass*Cp(T0));
    else {

        unsigned int steps = int(std::ceil(t/deltat));

        //std::cerr << "Steps: " << steps << std::endl;

        vector<double> Tn;

        Tn.push_back( T0 + (t/steps)*(Pp2 + Pup2 + getPhi()*(Tamb - T0) + getRho()*VdotUsed*getCp()*(Tinlet - T0))/(getRho()*getVolume()*getCp() + mass*Cp(T0)) );

        for (unsigned int i=1; i< steps; i++) {

            Tn.push_back( Tn[i-1] + (t/steps)*(Pp2 + Pup2 + getPhi()*(Tamb - Tn[i-1]) + getRho()*VdotUsed*getCp()*(Tinlet - Tn[i-1]))/(getRho()*getVolume()*getCp() + mass*Cp(Tn[i-1])) );

        }

        return std::min(Tn.back(), getTmax());

    }

}

double TankPCM::power(double t, double Tf, double VdotUsed, double Pup2, double T0, double Tinlet, double Tamb) {

    std::cerr << "rho: " << getRho() << "\tphi: " << getPhi() << std::endl;

    // calcul du crit�re de stabilit�

    double deltat = std::min( 300.0, 2.0*(getRho()*getVolume()*getCp() + mass*Cp(T0))/(getPhi()+ getRho()*VdotUsed*getCp()) );

    std::cerr << "Heating - explicit stability: " << 2.0*(getRho()*getVolume()*getCp() + mass*Cp(T0))/getPhi() << std::endl;

    if ( t < deltat ) return ((Tf-T0)/t)*(getRho()*getVolume()*getCp() + mass*Cp(T0)) - Pup2 - getPhi()*(Tamb - T0) - getRho()*VdotUsed*getCp()*(Tinlet - T0);
    else {

        unsigned int steps = int(std::ceil(t/deltat));

        //std::cerr << "Steps: " << steps << std::endl;

        vector<double> Hn;

        Hn.push_back( ((Tf-T0)/(t/steps))*(getRho()*getVolume()*getCp() + mass*Cp(T0)) - Pup2 - getPhi()*(Tamb - T0) - getRho()*VdotUsed*getCp()*(Tinlet - T0) );

        for (unsigned int i=1; i< steps; i++) {

            Hn.push_back( - Pup2 - getPhi()*(Tamb - Tf) - getRho()*VdotUsed*getCp()*(Tinlet - T0) );

        }

        return accumulate(Hn.begin(), Hn.end(), 0.0)/double(steps);

    }

}

double Tank::domesticHotWater(double t, double Tf, double VdotUsedUp, double Pp2, double Pup2, double T0, double Tinlet, double Tamb) {

    double VdotUp = VdotUsedUp;
    double VdotDown = 0.0;

    double VdotMid = 0.5*(VdotUp + VdotDown);

    do {

        if ( temperature(t, VdotMid, Pp2, Pup2, T0, Tinlet, Tamb) < Tf ) VdotUp = VdotMid;
        else VdotDown = VdotMid;

        VdotMid = 0.5*(VdotUp + VdotDown);

        //std::cerr << "VdotMid: " << VdotMid << std::endl;

    }
    while ( std::abs(temperature(t, VdotMid, Pp2, Pup2, T0, Tinlet, Tamb)-(Tf-0.01)) > 1e-2 && VdotMid > 1e-20 && VdotMid < (VdotUsedUp - 1e-20) );

    if ( VdotMid < 1e-20 ) return 0.0;
    else if ( VdotMid > (VdotUsedUp - 1e-20) ) return VdotUsedUp;
    else return VdotMid;

}


// Cognet: Start added content.
double Tank::maxSolPowerToNotExceedTcrit(double t, double VdotUsed, double power, double solPower, double T0, double Tinlet, double Tamb){
    double solPowerToUse;

    double T1 = temperature(t, VdotUsed, power, solPower, T0, Tinlet, Tamb);

    if ( T1 < Tcritical ) { // If stay below Tcritical.
        solPowerToUse = solPower; // Use all thermal power.
     }
    else { // If reaches Tcritical, only heat partially.
        double powHigh = solPower;
        double powLow = 0; //solPower*(Tcritical-T0)/(T1-T0);
        double powMid;
        int nbLoops = 0;
        bool notConverged = true;
        while ( notConverged ) {
            powMid = (powHigh+powLow)*0.5;
            T1 = temperature(t, VdotUsed, power, powMid, T0, Tinlet, Tamb);
            if ( T1<Tcritical ) { powLow = powMid; }
            else { powHigh = powMid; }
            nbLoops++;
            notConverged = ( abs(T1-Tcritical)>0.05 and nbLoops<10 );
        }
        solPowerToUse = powMid;
    }
    return solPowerToUse;
}
// Cognet: End added content


Boiler::Boiler(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr):EnergyConversionSystem(beginD, endD, pLogStr){
    //cout << "Boiler..." << endl << flush;
    boilerThermalPower = to<double>(hdl.ToElement()->Attribute("Pmax"));
    boilerThermalEfficiency = to<double>(hdl.ToElement()->Attribute("eta_th"));
    if (hdl.ToElement()->Attribute("name")){ hdl.ToElement()->QueryStringAttribute("name",&name);}
}

void Boiler::writeXML(ofstream& file, string tab=""){
    file << tab << "<Boiler name=\"" << name << "\" Pmax=\"" << boilerThermalPower << "\" eta_th=\"" << boilerThermalEfficiency << "\"/>" << endl;
}

void Boiler::writeGML(ofstream& file, string tab="") {

    file << tab << "<energy:Boiler>" << endl;
    file << tab << "\t<energy:installedNominalPower uom=\"W\">" << this->boilerThermalPower << "</energy:installedNominalPower>" << endl;
    file << tab << "\t<energy:nominalEfficiency uom=\"ratio\">" << this->boilerThermalEfficiency << "</energy:nominalEfficiency>" << endl;
    file << tab << "</energy:Boiler>" << endl;

}

CoGenerationAndBoiler::CoGenerationAndBoiler(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr):EnergyConversionSystem(beginD,endD,pLogStr){
    logStream << "Cogen + HP." << endl << flush;

    coGenThermalPower           = to<double>(hdl.ToElement()->Attribute("Pmax"));
    coGenElectricalEfficiency   = to<double>(hdl.ToElement()->Attribute("eta_el"));
    coGenThermalEfficiency      = to<double>(hdl.ToElement()->Attribute("eta_th"));
    coGenMinPartLoadCoefficient = to<double>(hdl.ToElement()->Attribute("minPartLoadCoeff"));

    boilerThermalPower      = to<double>(hdl.FirstChildElement("Boiler").ToElement()->Attribute("Pmax"));
    boilerThermalEfficiency = to<double>(hdl.FirstChildElement("Boiler").ToElement()->Attribute("eta_th"));
}

void CoGenerationAndBoiler::writeXML(ofstream& file, string tab){
    file << tab << "<CoGenerationAndBoiler Pmax=\"" << coGenThermalPower << "\" eta_el=\"" << coGenElectricalEfficiency
         << "\" eta_th=\"" << coGenThermalEfficiency << "\" minPartLoadCoeff=\"" << coGenMinPartLoadCoefficient << "\">";
    file << tab << "<Boiler Pmax=\"" << boilerThermalPower << "\" eta_th=\"" << boilerThermalEfficiency << "\"/>" << endl;
    file << "</CoGenerationAndBoiler>" << endl;
}

CoGenerationAndBoiler::CoGenerationAndBoiler(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient, double boilerThermalPower, double boilerThermalEfficiency) {

    this->coGenThermalPower           = coGenThermalPower;
    this->coGenElectricalEfficiency   = coGenElectricalEfficiency;
    this->coGenThermalEfficiency      = coGenThermalEfficiency;
    this->coGenMinPartLoadCoefficient = coGenMinPartLoadCoefficient;

    this->boilerThermalPower        = boilerThermalPower;
    this->boilerThermalEfficiency   = boilerThermalEfficiency;

}

void CoGenerationAndBoiler::getMaxThermalPower(double thermalPower1Needed, double thermalPower2Needed, double &thermalPower1Available, double &thermalPower2Available, double sourceTemp) {

    if ( (thermalPower1Needed+thermalPower2Needed) < coGenThermalPower*coGenMinPartLoadCoefficient ) {

        thermalPower1Available = 0.0;
        thermalPower2Available = 0.0;

    }
    else if ( (thermalPower1Needed+thermalPower2Needed) <= (coGenThermalPower + boilerThermalPower) ) {

        thermalPower1Available = thermalPower1Needed;
        thermalPower2Available = thermalPower2Needed;

    }
    else {

        if ( coGenThermalPower + boilerThermalPower >= thermalPower1Needed ) {
            // th1 satisfait
            thermalPower1Available = thermalPower1Needed;
            thermalPower2Available = coGenThermalPower + boilerThermalPower - thermalPower1Needed;
        }
        else {
            // th1 meme pas satisfait, donne le max
            thermalPower1Available = coGenThermalPower + boilerThermalPower;
            thermalPower2Available = 0.0;
        }
    }
}

double CoGenerationAndBoiler::getElectricConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) {

    if (thermalPower1+thermalPower2 < coGenMinPartLoadCoefficient*coGenThermalPower) {
        return 0.0;
    }
    else if (thermalPower1+thermalPower2 > coGenThermalPower) {
        return -time*(coGenThermalPower/coGenThermalEfficiency)*coGenElectricalEfficiency;
    }
    else {
        return -time*((thermalPower1+thermalPower2)/coGenThermalEfficiency)*coGenElectricalEfficiency;
    }

}

double CoGenerationAndBoiler::getCO2Production(double time, double thermalPower1, double thermalPower2, double sourceTemp, double outputTemp1, double outputTemp2) {

    if (thermalPower1+thermalPower2 < coGenMinPartLoadCoefficient*coGenThermalPower) {
        // no operation
        return 0.0;
    }
    else if (thermalPower1+thermalPower2 > coGenThermalPower) {
        return (time*coGenThermalPower/coGenThermalEfficiency)*coGenCO2EmissionCoefficient + (time*(thermalPower1+thermalPower2-coGenThermalPower)/boilerThermalEfficiency)*boilerCO2EmissionCoefficient;
    }
    else {
        return (time*(thermalPower1+thermalPower2)/coGenThermalEfficiency)*coGenCO2EmissionCoefficient;
    }
}

double CoGenerationAndBoiler::getFuelConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp, double outputTemp1, double outputTemp2) {

    if (thermalPower1+thermalPower2 < coGenMinPartLoadCoefficient*coGenThermalPower) {
        // no operation
        return 0.0;
    }
    else if (thermalPower1+thermalPower2 > coGenThermalPower) {
        return (time*coGenThermalPower/coGenThermalEfficiency);
    }
    else {
        return (time*(thermalPower1+thermalPower2)/coGenThermalEfficiency);
    }
}

void CoGenerationAndBoiler::getMaxThermalPower(vector<double> thermalPowerNeeded, vector<double> &thermalPowerAvailable) {

    // district version with many needs, the needs are the sum of the HS and DHW
    // if the one need is not satisfied, we have to give priority to HS or DHW... depending on the strategy

    // clear the vector before being filled again
    thermalPowerAvailable.clear();

    if ( accumulate(thermalPowerNeeded.begin(),thermalPowerNeeded.end(), 0.0) < coGenThermalPower*coGenMinPartLoadCoefficient ) {

        // assign everything to zero
        thermalPowerAvailable.assign(thermalPowerNeeded.size(), 0.0);

    }
    else if ( accumulate(thermalPowerNeeded.begin(),thermalPowerNeeded.end(), 0.0) <= (coGenThermalPower + boilerThermalPower) ) {

        // demande enti�rement satisfaite
        thermalPowerAvailable = thermalPowerNeeded;

    }
    else {

        double accu = 0.0;

        for (unsigned int i=0; i<thermalPowerNeeded.size(); i++) {

            accu += thermalPowerNeeded[i];
            if ( coGenThermalPower + boilerThermalPower >= accu) {
                // ce bout satisfait
                thermalPowerAvailable.push_back(thermalPowerNeeded[i]);
            }
            else {
                // ce bout pas satisfait
                thermalPowerAvailable.push_back( max(coGenThermalPower + boilerThermalPower - (accu-thermalPowerNeeded[i]), 0.0) );
            }
        }
    }

    return;
}

double CoGenerationAndBoiler::getElectricProduction(vector<double> thermalPower) {

    if ( accumulate(thermalPower.begin(),thermalPower.end(), 0.0) < coGenMinPartLoadCoefficient*coGenThermalPower) {
        return 0.0;
    }
    else if ( accumulate(thermalPower.begin(),thermalPower.end(), 0.0) > coGenThermalPower) {
        return (coGenThermalPower/coGenThermalEfficiency)*coGenElectricalEfficiency;
    }
    else {
        return (accumulate(thermalPower.begin(),thermalPower.end(), 0.0)/coGenThermalEfficiency)*coGenElectricalEfficiency;
    }

}

double CoGenerationAndBoiler::getCO2Production(double time, vector<double> thermalPower) {

    if ( accumulate(thermalPower.begin(),thermalPower.end(), 0.0) < coGenMinPartLoadCoefficient*coGenThermalPower) {
        // no operation
        return 0.0;
    }
    else if ( accumulate(thermalPower.begin(),thermalPower.end(), 0.0) > coGenThermalPower) {
        return (time*coGenThermalPower/coGenThermalEfficiency)*coGenCO2EmissionCoefficient + (time*(accumulate(thermalPower.begin(),thermalPower.end(), 0.0)-coGenThermalPower)/boilerThermalEfficiency)*boilerCO2EmissionCoefficient;
    }
    else {
        return (time*(accumulate(thermalPower.begin(),thermalPower.end(), 0.0))/coGenThermalEfficiency)*coGenCO2EmissionCoefficient;
    }
}

HeatPumpAndElectricElement::HeatPumpAndElectricElement(double heatPumpElectricPower, double heatPumpCOP, double heatPumpSrcTemp, double heatPumpOutputTemp, double electricElementPower) {

    ground = false;

    this->heatPumpElectricPower = heatPumpElectricPower;
    this->electricElementPower = electricElementPower;

    etaTech = heatPumpCOP / epsilonC(heatPumpSrcTemp, heatPumpOutputTemp);
    targetTemp = heatPumpOutputTemp;

}

void HeatPumpAndElectricElement::getMaxThermalPower(double thermalPower1Needed, double thermalPower2Needed, double &thermalPower1Available, double &thermalPower2Available, double sourceTemp) {

    double outputTemp1 = targetTemp;
    double outputTemp2 = targetTemp;

    if ( getWorkNeeded(thermalPower1Needed, sourceTemp, outputTemp1) + getWorkNeeded(thermalPower2Needed, sourceTemp, outputTemp2) <= heatPumpElectricPower ) {
        // HP only operation - everything provided
        thermalPower1Available = thermalPower1Needed;
        thermalPower2Available = thermalPower2Needed;
    }
    else {
        // meilleur temps de chauffer la maison d'abord, et en plus meilleur COP!!

        double work1, work2;
        double thermalPower1Missing, thermalPower2Missing;

        work1 = getWorkNeeded(thermalPower1Needed, sourceTemp, outputTemp1);
        if (work1 <= heatPumpElectricPower) { // thermalPower1 given...
            work2 = heatPumpElectricPower-work1;
            thermalPower1Missing = 0.0;
            thermalPower2Missing = thermalPower2Needed - getHeatProduced(work2, sourceTemp, outputTemp2);
        }
        else {
            // work1 insufficient for thermalPower1, limited to heatPumpElectricPower
            work1 = heatPumpElectricPower;
            work2 = 0.0;
            thermalPower1Missing = thermalPower1Needed - getHeatProduced(heatPumpElectricPower, sourceTemp, outputTemp1);
            thermalPower2Missing = thermalPower2Needed;
        }

        // on utilise l'element electrique pour promouvoir la chaleur dans 1 puis dans 2
        if (electricElementPower > 0.0) { // mode chauffage
            if ( electricElementPower >= thermalPower1Missing ) {

                // on a tout pour th1
                thermalPower1Available = thermalPower1Needed;

                if ( electricElementPower - thermalPower1Missing >= thermalPower2Missing ) thermalPower2Available = thermalPower2Needed;
                else thermalPower2Available = electricElementPower - thermalPower1Missing;
            }
            else {
                thermalPower1Available = thermalPower1Needed - thermalPower1Missing + electricElementPower;
                thermalPower2Available = 0.0;
            }
        }
        else {
            // seulement en mode refroidissement
            thermalPower1Available = thermalPower1Needed - thermalPower1Missing;
            thermalPower2Available = thermalPower2Needed - thermalPower2Missing;

        }
    }

}


double HeatPumpAndElectricElement::getCO2Production(double time, double thermalPower1, double thermalPower2, double sourceTemp) {

    double outputTemp1 = targetTemp;
    double outputTemp2 = targetTemp;

    if ( getWorkNeeded(thermalPower1, sourceTemp, outputTemp1) + getWorkNeeded(thermalPower2, sourceTemp, outputTemp2) <= heatPumpElectricPower ) {
        // HP only operation - everything provided
        return time*( thermalPower1/( etaTech*epsilonC(sourceTemp, outputTemp1) ) + thermalPower2/( etaTech*epsilonC(sourceTemp, outputTemp2)))*electricCO2EmissionCoefficient;
    }
    else {
        // meilleur temps de chauffer la maison d'abord, et en plus meilleur COP!!

        double work1, work2;
        double thermalPower1Missing, thermalPower2Missing;

        work1 = getWorkNeeded(thermalPower1, sourceTemp, outputTemp1);
        if (work1 <= heatPumpElectricPower) { // thermalPower1 given...
            work2 = heatPumpElectricPower-work1;
            thermalPower1Missing = 0.0;
            thermalPower2Missing = thermalPower2 - getHeatProduced(work2, sourceTemp, outputTemp2);
        }
        else {
            // work1 insufficient for thermalPower1, limited to heatPumpElectricPower
            work1 = heatPumpElectricPower;
            work2 = 0.0;
            thermalPower1Missing = thermalPower1 - getHeatProduced(heatPumpElectricPower, sourceTemp, outputTemp1);
            thermalPower2Missing = thermalPower2;
        }

        if ( electricElementPower == 0.0 ) return (time*(work1+work2))*electricCO2EmissionCoefficient;
        else return (time*(work1+work2))*electricCO2EmissionCoefficient + (time*(std::max(thermalPower1Missing + thermalPower2Missing, electricElementPower)))*electricCO2EmissionCoefficient;
    }
}

double HeatPumpAndElectricElement::getElectricConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) {

    double outputTemp1 = targetTemp;
    double outputTemp2 = targetTemp;

    if ( getWorkNeeded(thermalPower1, sourceTemp, outputTemp1) + getWorkNeeded(thermalPower2, sourceTemp, outputTemp2) <= heatPumpElectricPower ) {
        // HP only operation - everything provided
        return time*( thermalPower1/( etaTech*epsilonC(sourceTemp, outputTemp1) ) + thermalPower2/( etaTech*epsilonC(sourceTemp, outputTemp2)));
    }
    else {
        // meilleur temps de chauffer la maison d'abord, et en plus meilleur COP!!

        double work1, work2;
        double thermalPower1Missing, thermalPower2Missing;

        work1 = getWorkNeeded(thermalPower1, sourceTemp, outputTemp1);
        if (work1 <= heatPumpElectricPower) { // thermalPower1 given...
            work2 = heatPumpElectricPower-work1;
            thermalPower1Missing = 0.0;
            thermalPower2Missing = thermalPower2 - getHeatProduced(work2, sourceTemp, outputTemp2);
        }
        else {
            // work1 insufficient for thermalPower1, limited to heatPumpElectricPower
            work1 = heatPumpElectricPower;
            work2 = 0.0;
            thermalPower1Missing = thermalPower1 - getHeatProduced(heatPumpElectricPower, sourceTemp, outputTemp1);
            thermalPower2Missing = thermalPower2;
        }

        if ( electricElementPower == 0.0 ) return (time*(work1+work2));
        else return (time*(work1+work2)) + (time*(std::max(thermalPower1Missing + thermalPower2Missing, electricElementPower)));
    }
}

// version district... 80�C!!

void HeatPumpAndElectricElement::getMaxThermalPower(vector<double> thermalPowerNeeded, vector<double> &thermalPowerAvailable, double sourceTemp) {

    double outputTemp = targetTemp;

    if ( getWorkNeeded(accumulate(thermalPowerNeeded.begin(),thermalPowerNeeded.end(), 0.0), sourceTemp, outputTemp) <= heatPumpElectricPower ) {
        // HP only operation - everything provided
        thermalPowerAvailable = thermalPowerNeeded;
    }
    else {
        // meilleur temps de chauffer la maison d'abord, et en plus meilleur COP!!

        double accuMissing = 0.0, accuWork = 0.0;

        for (unsigned int i=0; i<thermalPowerNeeded.size(); i++) {

            accuWork += getWorkNeeded(thermalPowerNeeded[i], sourceTemp, outputTemp);

            if ( heatPumpElectricPower >= accuWork) {
                // ce bout satisfait
                thermalPowerAvailable.push_back(thermalPowerNeeded[i]);
            }
            else {


                accuMissing += thermalPowerNeeded[i] - getHeatProduced(max(heatPumpElectricPower-(accuWork-getWorkNeeded(thermalPowerNeeded[i], sourceTemp, outputTemp)),0.0), sourceTemp, outputTemp);

                if ( electricElementPower > 0.0 ) {
                    if ( electricElementPower >= accuMissing ) {
                        // on a tout
                        thermalPowerAvailable.push_back(thermalPowerNeeded[i]);
                    }
                    else {
                        thermalPowerAvailable.push_back( max(electricElementPower - (accuMissing-thermalPowerNeeded[i]), 0.0) );
                    }
                }
                else {

                    thermalPowerAvailable.push_back( min((accuMissing-thermalPowerNeeded[i]), 0.0) );
                }
            }
        }
    }

}

double HeatPumpAndElectricElement::getCO2Production(double time, vector<double> thermalPower, double sourceTemp) {

    double outputTemp = targetTemp;

    double workHP=0.0, heatEL=0.0;

    if ( getWorkNeeded(accumulate(thermalPower.begin(),thermalPower.end(), 0.0), sourceTemp, outputTemp) <= heatPumpElectricPower ) {
        // HP only operation - everything provided
        workHP = getWorkNeeded( accumulate(thermalPower.begin(),thermalPower.end(), 0.0), sourceTemp, outputTemp );
        heatEL = 0.0;
        cerr << "Everything provided by HP!" << endl;
    }
    else {

        double accuMissing = 0.0, accuWork = 0.0, accu = 0.0;

        for (unsigned int i=0; i<thermalPower.size(); i++) {

            accu     += thermalPower[i];
            accuWork += getWorkNeeded(thermalPower[i], sourceTemp, outputTemp);

            if ( heatPumpElectricPower >= accuWork) {
                // ce bout satisfait
                workHP += getWorkNeeded(thermalPower[i], sourceTemp, outputTemp);
            }
            else {

                workHP = heatPumpElectricPower;
                accuMissing = accu - getHeatProduced(heatPumpElectricPower, sourceTemp, outputTemp);

                if ( electricElementPower > 0.0 ) {
                    if ( electricElementPower >= accuMissing ) {
                        // on a tout
                        heatEL = accuMissing;
                    }
                    else {
                        heatEL += accuMissing - electricElementPower;
                    }
                }
            }
        }
    }

    cerr << "Work HP: " << workHP << "\tHeat Electric: " << heatEL << endl;

    if ( accumulate(thermalPower.begin(), thermalPower.end(), 0.0) < 0.0 ) return (time*(workHP))*electricCO2EmissionCoefficient;
    else return (time*(workHP))*electricCO2EmissionCoefficient + (time*(heatEL))*electricCO2EmissionCoefficient;

}

CoGenerationHeatPumpAndBoiler::CoGenerationHeatPumpAndBoiler(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient, double heatPumpCOP, double heatPumpSrcTemp, double heatPumpOutputTemp, double boilerThermalPower, double boilerThermalEfficiency) {

    this->coGenThermalPower           = coGenThermalPower;
    this->coGenElectricalEfficiency   = coGenElectricalEfficiency;
    this->coGenThermalEfficiency      = coGenThermalEfficiency;
    this->coGenMinPartLoadCoefficient = coGenMinPartLoadCoefficient;

    etaTech = heatPumpCOP / epsilonC(heatPumpSrcTemp, heatPumpOutputTemp);

    this->boilerThermalPower        = boilerThermalPower;
    this->boilerThermalEfficiency   = boilerThermalEfficiency;

}

CoGenerationHeatPumpAndBoiler::CoGenerationHeatPumpAndBoiler(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient, double heatPumpEtaTech, double boilerThermalPower, double boilerThermalEfficiency) {

    this->coGenThermalPower           = coGenThermalPower;
    this->coGenElectricalEfficiency   = coGenElectricalEfficiency;
    this->coGenThermalEfficiency      = coGenThermalEfficiency;
    this->coGenMinPartLoadCoefficient = coGenMinPartLoadCoefficient;

    this->etaTech = heatPumpEtaTech;

    this->boilerThermalPower        = boilerThermalPower;
    this->boilerThermalEfficiency   = boilerThermalEfficiency;

}

void CoGenerationHeatPumpAndBoiler::getMaxThermalPower(double thermalPower1Needed, double thermalPower2Needed, double &thermalPower1Available, double &thermalPower2Available, double sourceTemp) {

    double outputTemp1 = targetTemp;
    double outputTemp2 = targetTemp;

    double COP1 = etaTech*epsilonC(sourceTemp, outputTemp1);
    double COP2 = etaTech*epsilonC(sourceTemp, outputTemp2);

    double Pgross; // gross energy input
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency;

    // cas P1*dt >= th1 ou < th1
    Pgross = (thermalPower1Needed+thermalPower2Needed)/(coGenElectricalEfficiency*COP2+coGenThermalEfficiency);
    if ( Pgross < thermalPower1Needed ) Pgross = (thermalPower1Needed*(COP2/COP1)+thermalPower2Needed)/(coGenElectricalEfficiency*COP2+coGenThermalEfficiency*(COP2/COP1));

    // voyons les diff�rents cas de travail de notre machine
    if ( Pgross > PgrossMax ) {
        // pas assez d'�nergie, regardons ce qu'il manque pour le boiler
        if (PgrossMax*coGenThermalEfficiency >= thermalPower1Needed) { // d�passe la zone critique, les besoins1 sont satisfaits
            thermalPower1Available = thermalPower1Needed;
            thermalPower2Available = PgrossMax*(coGenElectricalEfficiency*COP2+coGenThermalEfficiency) - thermalPower1Needed;
            if ( (thermalPower2Needed - thermalPower2Available) > boilerThermalPower ) thermalPower2Available += boilerThermalPower;
            else thermalPower2Available = thermalPower2Needed;
        }
        else { // ne d�passe pas th1, formule2
            thermalPower2Available = PgrossMax*(coGenElectricalEfficiency*COP2+coGenThermalEfficiency*(COP2/COP1)) - thermalPower1Needed*(COP2/COP1);
            if ( thermalPower2Available < 0.0 ) {
                thermalPower2Available = 0.0;
                thermalPower1Available = PgrossMax*(coGenElectricalEfficiency*COP2+coGenThermalEfficiency*(COP2/COP1))*(COP1/COP2);
                if ( (thermalPower1Needed - thermalPower1Available) > boilerThermalPower ) thermalPower1Available += boilerThermalPower;
                else {
                    double boilerThermalPowerLeft = boilerThermalPower - (thermalPower1Needed-thermalPower1Available);
                    thermalPower1Available = thermalPower1Needed;
                    thermalPower2Available = std::max(thermalPower2Needed, boilerThermalPowerLeft);
                }
            }
            else { // compl�ter avec du boiler pour thermalPower2
                if ( (thermalPower2Needed - thermalPower2Available) > boilerThermalPower ) thermalPower2Available += boilerThermalPower;
                else thermalPower2Available = thermalPower2Needed;
            }
        }
    }
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) {

        // ne fonctionne pas en mode SOFC, tout doit etre fourni par le boiler
        if ( (thermalPower1Needed+thermalPower2Needed) <= boilerThermalPower ) {

            thermalPower1Available = thermalPower1Needed;
            thermalPower2Available = thermalPower2Needed;

        }
        else {

            thermalPower2Available = boilerThermalPower - thermalPower1Needed;
            if ( thermalPower2Available < 0.0 ) {
                thermalPower1Available = thermalPower1Needed + thermalPower2Available;
                thermalPower2Available = 0.0;
            }
            else {
                thermalPower1Available = thermalPower1Needed;
            }
        }

    }
    else { // bonne formule, dans la zone de travail id�ale!

        thermalPower1Available = thermalPower1Needed;
        thermalPower2Available = thermalPower2Needed;

    }

}

double CoGenerationHeatPumpAndBoiler::getCO2Production(double time, double thermalPower1, double thermalPower2, double sourceTemp, double outputTemp1, double outputTemp2) {

    double COP1 = etaTech*epsilonC(sourceTemp, outputTemp1);
    double COP2 = etaTech*epsilonC(sourceTemp, outputTemp2);

    double Pgross; // gross energy input
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency;

    // cas P1*dt >= th1 ou < th1
    Pgross = (thermalPower1+thermalPower2)/(coGenElectricalEfficiency*COP2+coGenThermalEfficiency);
    if ( Pgross < thermalPower1 ) Pgross = (thermalPower1*(COP2/COP1)+thermalPower2)/(coGenElectricalEfficiency*COP2+coGenThermalEfficiency*(COP2/COP1));

    // voyons les diff�rents cas de travail de notre machine
    if ( Pgross > PgrossMax ) {
        // pas assez d'�nergie, tourne au max
        return time*PgrossMax*coGenCO2EmissionCoefficient;
    }
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) {
        // ne fonctionne pas en mode SOFC
        return 0.0;
    }
    else { // bonne formule, dans la zone de travail id�ale!

        return time*Pgross*coGenCO2EmissionCoefficient;

    }
}

double CoGenerationHeatPumpAndBoiler::getFuelConsumption(double time, double thermalPower1, double thermalPower2, double sourceTemp) {

    double outputTemp1 = targetTemp;
    double outputTemp2 = targetTemp;

    double COP1 = etaTech*epsilonC(sourceTemp, outputTemp1);
    double COP2 = etaTech*epsilonC(sourceTemp, outputTemp2);

    double Pgross; // gross energy input
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency;

    // cas P1*dt >= th1 ou < th1
    Pgross = (thermalPower1+thermalPower2)/(coGenElectricalEfficiency*COP2+coGenThermalEfficiency);
    if ( Pgross < thermalPower1 ) Pgross = (thermalPower1*(COP2/COP1)+thermalPower2)/(coGenElectricalEfficiency*COP2+coGenThermalEfficiency*(COP2/COP1));

    // voyons les diff�rents cas de travail de notre machine
    if ( Pgross > PgrossMax ) {
        // pas assez d'�nergie, tourne au max
        return time*PgrossMax;
    }
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) {
        // ne fonctionne pas en mode SOFC
        return 0.0;
    }
    else { // bonne formule, dans la zone de travail id�ale!

        return time*Pgross;

    }
}

// cas � plusieurs demandes, pour le district...
void CoGenerationHeatPumpAndBoiler::getMaxThermalPower(vector<double> thermalPowerNeeded, vector<double> &thermalPowerAvailable, double sourceTemp) {

    double outputTemp = targetTemp;

    double COP = etaTech*epsilonC(sourceTemp, outputTemp);

    cerr << "COP: " << COP << endl;

    double Pgross; // gross energy input
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency;

    // calcul du Pgross
    Pgross = accumulate(thermalPowerNeeded.begin(), thermalPowerNeeded.end(), 0.0) / (coGenElectricalEfficiency*COP+coGenThermalEfficiency);

    cerr << "Pgross: " << Pgross << "\tPgrossMax: " << PgrossMax << endl;

    // voyons les diff�rents cas de travail de notre machine
    if ( Pgross > PgrossMax ) {

        // pas assez d'�nergie, regardons ce qu'il manque pour le boiler
        double accu = 0.0;

        for (unsigned int i=0; i<thermalPowerNeeded.size(); i++) {

            accu += thermalPowerNeeded[i];
            if ( PgrossMax*(coGenThermalEfficiency + coGenElectricalEfficiency*COP) + boilerThermalPower >= accu) {
                // ce bout satisfait
                thermalPowerAvailable.push_back(thermalPowerNeeded[i]);
            }
            else {
                // ce bout pas satisfait
                thermalPowerAvailable.push_back( max(PgrossMax*(coGenThermalEfficiency + coGenElectricalEfficiency*COP) + boilerThermalPower - (accu-thermalPowerNeeded[i]), 0.0) );
            }

        }
    }
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) {

        thermalPowerAvailable.assign(thermalPowerNeeded.size(), 0.0);

    }
    else { // bonne formule, dans la zone de travail id�ale!

        thermalPowerAvailable = thermalPowerNeeded;

    }

}

double CoGenerationHeatPumpAndBoiler::getCO2Production(double time, vector<double> thermalPower, double sourceTemp, double outputTemp) {

    double COP = etaTech*epsilonC(sourceTemp, outputTemp);

    double Pgross; // gross energy input
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency;

    // calcul de Pgross
    Pgross = accumulate(thermalPower.begin(), thermalPower.end(), 0.0) / (coGenElectricalEfficiency*COP+coGenThermalEfficiency);

    // voyons les diff�rents cas de travail de notre machine
    if ( Pgross > PgrossMax ) {
        // pas assez d'�nergie, tourne au max
        return time*PgrossMax*coGenCO2EmissionCoefficient;
    }
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) {
        // ne fonctionne pas en mode SOFC
        return 0.0;
    }
    else { // bonne formule, dans la zone de travail id�ale!

        return time*Pgross*coGenCO2EmissionCoefficient;

    }
}

CoGeneration::CoGeneration(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr):EnergyConversionSystem(beginD, endD, pLogStr) {
    coGenThermalPower = to<double>(hdl.ToElement()->Attribute("Pmax"));
    coGenElectricalEfficiency = to<double>(hdl.ToElement()->Attribute("eta_el"));
    coGenThermalEfficiency = to<double>(hdl.ToElement()->Attribute("eta_th"));
    coGenMinPartLoadCoefficient = to<double>(hdl.ToElement()->Attribute("minPartLoadCoeff"));
}

CoGeneration::CoGeneration(double coGenThermalPower,double coGenElectricalEfficiency,double coGenThermalEfficiency,double coGenMinPartLoadCoefficient,unsigned int beginDay,unsigned int endDay):
    EnergyConversionSystem(beginDay,endDay),
    coGenThermalPower(coGenThermalPower),
    coGenElectricalEfficiency(coGenElectricalEfficiency),
    coGenThermalEfficiency(coGenThermalEfficiency),
    coGenMinPartLoadCoefficient(coGenMinPartLoadCoefficient)
{}

void CoGeneration::writeXML(ofstream& file, string tab){
    file << tab << "<CoGeneration Pmax=\"" << coGenThermalPower << "\" eta_el=\"" << coGenElectricalEfficiency
         << "\" eta_th=\"" << coGenThermalEfficiency << "\" minPartLoadCoeff=\"" << coGenMinPartLoadCoefficient << "\" />" << endl;
}


//double CoGeneration::getThermalPower(double thermalPowerNeeded, double sourceTemp) { // Cognet: Deleted this.
double CoGeneration::getThermalPower(double sourceTemp) { // Cognet: Added this. thermalPowerNeeded is now a class attribute set beforehand.

    if ( thermalPowerNeeded < coGenThermalPower*coGenMinPartLoadCoefficient ) return 0.0; // under the min part load coefficient
    else if ( thermalPowerNeeded <= coGenThermalPower )  return thermalPowerNeeded;
    else return coGenThermalPower;

}

double CoGeneration::getFuelConsumption(double time, double thermalPower, double sourceTemp) {

    if ( thermalPower < coGenThermalPower*coGenMinPartLoadCoefficient ) return 0.0;
    else if ( thermalPower <= coGenThermalPower) return (time*thermalPower/coGenThermalEfficiency);
    else return (time*coGenThermalPower/coGenThermalEfficiency);

}

double CoGeneration::getElectricConsumption(double time, double thermalPower, double sourceTemp) {

    if ( thermalPower < coGenThermalPower*coGenMinPartLoadCoefficient) return 0.0;
    else if ( thermalPower <= coGenThermalPower) return -time*(thermalPower/coGenThermalEfficiency)*coGenElectricalEfficiency;
    else return -time*(coGenThermalPower/coGenThermalEfficiency)*coGenElectricalEfficiency;

}

HeatPump::HeatPump(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr):EnergyConversionSystem(beginD, endD, pLogStr) {
    logStream << "Heat Pump, ";

    if (to<double>(hdl.ToElement()->Attribute("Pmax"))) heatPumpElectricPower = to<double>(hdl.ToElement()->Attribute("Pmax"));
    etaTech = to<double>(hdl.ToElement()->Attribute("eta_tech"));
    targetTemp = to<double>(hdl.ToElement()->Attribute("Ttarget"));
    // output of the standard performance
    logStream << "COP (2C/" << targetTemp << "C): " << etaTech*fabs(epsilonC(2.,targetTemp)) << ", COP (35C/" << targetTemp << "C): " << etaTech*fabs(epsilonC(35., targetTemp)) << endl;

    if (string(hdl.ToElement()->Attribute("Tsource"))==string("ground")) {
        logStream << "Ground source." << endl << flush;
        ground = true;

        z0 = to<double>(hdl.ToElement()->Attribute("depth"));
        alpha = to<double>(hdl.ToElement()->Attribute("alpha"));

        if (hdl.ToElement()->Attribute("z1")) {   // vertical pipes in the ground
            z1 = to<double>(hdl.ToElement()->Attribute("z1"));
        }
        else {      // horizontal pipes in the ground
            z1 = to<double>(hdl.ToElement()->Attribute("depth"));
        }
    }
//    else{  // Cognet: Deleted this, to improve the if to verify that it is "air", else throw error
    else if (string(hdl.ToElement()->Attribute("Tsource"))==string("air")) { // Cognet: Added this.
        logStream << "Air source." << endl << flush;
        ground = false;
    }
    else { throw string("Error in the XML file: a HeatPump has a Tsource attribute with neither 'ground' nor 'air'."); } // Cognet: Added this.
}

void HeatPump::writeXML(ofstream& file, string tab){
    file << tab << "<HeatPump Pmax=\"" << heatPumpElectricPower << "\" eta_tech=\"" << etaTech
         << "\" Ttarget=\"" << targetTemp << "\" ";
    if(ground){
        file << "Tsource=\"ground\" depth=\""<< z0 << "\" alpha=\"" << alpha << "\" ";
        if (z1 != z0) file << "position=\"vertical\" z1=\"" << z1 << "\" ";
    }
    else file << "Tsource=\"air\"";
    file << "/>" << endl;
}

void HeatPump::writeGML(ofstream& file, string tab) {

    file << tab << "<energy:HeatPump>" << endl;
    file << tab << "\t<energy:installedNominalPower uom=\"W\">" << this->heatPumpElectricPower << "</energy:installedNominalPower>" << endl;
    file << tab << "\t<energy:nominalEfficiency uom=\"ratio\">" << this->etaTech << "</energy:nominalEfficiency>" << endl;
    if (targetTemp >= 35.) {
        file << tab << "\t<energy:carnotEfficiency>" << this->epsilonC(2.,targetTemp) << "</energy:carnotEfficiency>" << endl;
    }
    else {
        file << tab << "\t<energy:carnotEfficiency>" << fabs(this->epsilonC(35.,targetTemp)) << "</energy:carnotEfficiency>" << endl;
    }
    file << tab << "\t<energy:heatSource>";
    if (ground) {
        if (z0 == z1) file << "HorizontalGroundCollector";
        else          file << "VerticalGroundCollector";
    }
    else file << "AmbientAir";
    file << "</energy:heatSource>" << endl;
    file << tab << "</energy:HeatPump>" << endl;

}

HeatPump::HeatPump(double heatPumpElectricPower,double heatPumpCOP,double heatPumpSrcTemp,double heatPumpOutputTemp,unsigned int beginDay,unsigned int endDay):
    EnergyConversionSystem(beginDay,endDay),
    heatPumpElectricPower(heatPumpElectricPower),
    targetTemp(heatPumpOutputTemp)
{
    ground = false;
    etaTech = heatPumpCOP / epsilonC(heatPumpSrcTemp, heatPumpOutputTemp);
}

HeatPump::HeatPump(double heatPumpElectricPower,double heatPumpEtaTech,double targetTemp,unsigned int beginDay,unsigned int endDay):
    EnergyConversionSystem(beginDay,endDay),
    heatPumpElectricPower(heatPumpElectricPower),
    targetTemp(targetTemp),
    etaTech(heatPumpEtaTech)
{ ground = false; }

double HeatPump::getHeatProduced(double work, double sourceTemp, double outputTemp) {

    // as epsilonC follows the sign of the heat/cold demand, so does the heat/cold produced
    return work*etaTech*epsilonC(sourceTemp, outputTemp);

}

double HeatPump::getWorkNeeded(double thermalPower, double sourceTemp, double outputTemp) {

    // the work is always positive
    return std::abs(thermalPower/(etaTech*epsilonC(sourceTemp, outputTemp)));

}

//double HeatPump::getThermalPower(double thermalPowerNeeded, double sourceTemp) { // Cognet: Deleted this.
double HeatPump::getThermalPower(double sourceTemp) { // Cognet: Added this. thermalPowerNeeded is now a class attribute set beforehand.

    if ( getWorkNeeded(thermalPowerNeeded, sourceTemp, targetTemp) <= heatPumpElectricPower ) return thermalPowerNeeded;
    else return getHeatProduced(heatPumpElectricPower, sourceTemp, targetTemp);

}

double HeatPump::getElectricConsumption(double time, double thermalPower, double sourceTemp) {

    if ( getWorkNeeded(thermalPower, sourceTemp, targetTemp) <= heatPumpElectricPower ) return time*(thermalPower/(etaTech*epsilonC(sourceTemp, targetTemp)));
    else return time*heatPumpElectricPower;

}

CoGenerationHeatPump::CoGenerationHeatPump(TiXmlHandle hdl, unsigned int beginD, unsigned int endD, ostream* pLogStr):CoGeneration(hdl,beginD,endD,pLogStr),HeatPump(hdl.FirstChildElement("HeatPump"),beginD,endD,pLogStr) {
    logStream << "Cogen + HP." << endl << flush;
    // Ignore heatPumpElectricPower read in xml file ? Or keep line bellow only if no value in xml file ?
    //if(heatPumpElectricPower==0)
    heatPumpElectricPower = coGenThermalPower/coGenThermalEfficiency*coGenElectricalEfficiency;
}

CoGenerationHeatPump::CoGenerationHeatPump(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient,
                     double heatPumpCOP, double heatPumpSrcTemp, double heatPumpOutputTemp,unsigned int beginDay,unsigned int endDay):
    CoGeneration(coGenThermalPower, coGenElectricalEfficiency, coGenThermalEfficiency, coGenMinPartLoadCoefficient,beginDay,endDay),
    HeatPump(coGenThermalPower/coGenThermalEfficiency*coGenElectricalEfficiency, heatPumpCOP, heatPumpSrcTemp, heatPumpOutputTemp,beginDay,endDay) {}

CoGenerationHeatPump::CoGenerationHeatPump(double coGenThermalPower, double coGenElectricalEfficiency, double coGenThermalEfficiency, double coGenMinPartLoadCoefficient,
                     double heatPumpEtaTech, double targetTemp,unsigned int beginDay,unsigned int endDay) :
    CoGeneration(coGenThermalPower, coGenElectricalEfficiency, coGenThermalEfficiency, coGenMinPartLoadCoefficient,beginDay,endDay),
    HeatPump(coGenThermalPower/coGenThermalEfficiency*coGenElectricalEfficiency, heatPumpEtaTech, targetTemp,beginDay,endDay) {}


void CoGenerationHeatPump::writeXML(ofstream& file, string tab){
    file << tab << "<CHP-HP Pmax=\"" << coGenThermalPower << "\" eta_el=\"" << coGenElectricalEfficiency
         << "\" eta_th=\"" << coGenThermalEfficiency << "\" minPartLoadCoeff=\"" << coGenMinPartLoadCoefficient << "\">" << endl;
    this->HeatPump::writeXML(file, tab+"\t");
    file << tab << "</CHP-HP>" << endl;
}

//double CoGenerationHeatPump::getThermalPower(double thermalPowerNeeded, double sourceTemp) { // Cognet: Deleted this.
double CoGenerationHeatPump::getThermalPower(double sourceTemp) { // Cognet: Added this. thermalPowerNeeded is now a class attribute set beforehand.
    double COP = etaTech*epsilonC(sourceTemp, targetTemp);

    double Pgross = thermalPowerNeeded/(coGenElectricalEfficiency*COP+coGenThermalEfficiency); // gross energy needed to satisfy the needs
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency; // maximal gross energy taken by the machine

    if ( Pgross > PgrossMax ) return PgrossMax*(coGenElectricalEfficiency*COP+coGenThermalEfficiency);
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) return 0.0;
    else return thermalPowerNeeded;

}

double CoGenerationHeatPump::getFuelConsumption(double time, double thermalPower, double sourceTemp) {

    double COP = etaTech*epsilonC(sourceTemp, targetTemp);

    double Pgross = thermalPower/(coGenElectricalEfficiency*COP+coGenThermalEfficiency);
    double PgrossMax = coGenThermalPower/coGenThermalEfficiency;

    if ( Pgross > PgrossMax ) return time*PgrossMax;
    else if ( (Pgross/PgrossMax) < coGenMinPartLoadCoefficient ) return 0.0;
    else return time*Pgross;

}





Pump::Pump(TiXmlHandle hdl) {

    string name = "efficiencyPump";
    if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &efficiencyPump) ) { throw string("Error in the XML file: a Pump doesn't have attribute: '"+name+"'."); }
    if ( efficiencyPump<=0 )    { throw string("Error in the XML file: a Pump has "+name+"<=0."); }

    name = "n0";
    if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &n0) ) { throw string("Error in the XML file: a Pump doesn't have attribute: '"+name+"'."); }
    if ( n0<=0 )    { throw string("Error in the XML file: a Pump has "+name+"<=0."); }

    name = "a0";
    if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &a0) ) { throw string("Error in the XML file: a Pump doesn't have attribute: '"+name+"'."); }
    if ( a0<=0 )    { throw string("Error in the XML file: a Pump "+name+"<=0."); }

    name = "a1";
    if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &a1) ) { throw string("Error in the XML file: a Pump doesn't have attribute: '"+name+"'."); }

    name = "a2";
    if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &a2) ) { throw string("Error in the XML file: a Pump doesn't have attribute: '"+name+"'."); }
    if ( a2>0 )    { throw string("Error in the XML file: a Pump has "+name+">0."); }

    n = 0.5f*n0; // TODO improve this?
}

void Pump::computeElectricAndThermalPower(float const& pressureDiff, float const& massFlow, float const& rho, float& electricPow, float& thermalPow) {
    if (massFlow>=0.f and pressureDiff<=0.f) { // Normal conditions, mass flow in correct direction and pump increases the pressure and temperature.
        // Hydraulic power needed is pressureDiff*massFlow/rho
        electricPow = -pressureDiff*massFlow/(rho*efficiencyPump); // Electric power supplied. (minus sign due to conventiions)
        thermalPow = electricPow*(1.f-efficiencyPump); // Since pumpPower = hydraulicPower+thermalPower = hydraulicPower/efficiency
    } else { // Pump is resisting the flow (similar to a valve), dissipating heat.
        electricPow = 0.f; // For simplifications, this no electricity, only heat dissipation (even if it could technically generate electricity, like a turbine).
        thermalPow = pressureDiff*massFlow/rho; // Should always be positive.
    }
}

void Pump::computePressureDiff(float const& massFlow, float& deltaP, float& dDeltaP_dm) {
    float x = n/n0;
    if (massFlow>=0.f) {
        deltaP = -(a0*x*x + a1*x*massFlow + a2*massFlow*massFlow);
        dDeltaP_dm = -(a1*x + 2.f*a2*massFlow);
    } else {
//        deltaP = -a0*x*x;
//        dDeltaP_dm = 0.f;
        deltaP = -a0*x*x + a2*massFlow*massFlow*1000000.f; // TODO: improve this. This is to ensure that the friction increases if flow goes backwards through the pump.
        dDeltaP_dm = 2.f*a2*massFlow*1000000.f;
    }
}

float Pump::computeIdealN(float const& massFlow, float const& pressureDiff) {
    float a = -a0;
    float b = -a1*massFlow;
    float c = -a2*massFlow*massFlow - pressureDiff;
    float delta = abs(b*b-4*a*c); // SHould always be positive even without abs.
    float idealN = n0 * ( -b-sqrt(delta) )/(2.f*a);
    return idealN;
}

void Pump::updateRpms(float& sumDeltaRpm, float& sumRpm, float const& learningRate, float const& targetMassFlow, float const& targetPressureDiff) {
    float nMax = n0;
    float nMin = computeNMin();
    float nPrev = n;

    n = computeIdealN(targetMassFlow, targetPressureDiff);
    if (n<nMin) { n = nMin; }
    else if (n>nMax) { n = nMax; }

    n = nPrev*(1.f-learningRate)+n*learningRate;

    sumDeltaRpm += abs(n-nPrev);
    sumRpm += n;
}

void Pump::setNToMax(float& sumDeltaRpm, float& sumRpm, float const& learningRate) {
    float nPrev = n;
    n = nPrev*(1.f-learningRate)+n0*learningRate;
    sumDeltaRpm += abs(n-nPrev);
    sumRpm += n;
}




float Valve::deltaP0_invRho0_36002 = 1296000000.f ; // deltaP0=1e5, invRho0=1/1e3, 3600^2=12960000, multiply them.

void Valve::computePressureDiffAndDerivative(float const& m, float const& rho, float& deltaP, float& dDeltaP_dm) {
    dDeltaP_dm = deltaP0_invRho0_36002*2.f*abs(m)/(rho*kv*kv);
    deltaP = 0.5f*m*dDeltaP_dm;
}

float Valve::computeIdealKv(float const& rho, float const& massFlow, float const& pressureDiff) {
    return massFlow*sqrt(abs(deltaP0_invRho0_36002/(rho*pressureDiff)));
}

void Valve::updateKv(float const& rho, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& targetMassFlow, float const& targetPressureDiff) {
    float kvMin = computeKvMin();
    float kvPrev = kv;

//    if (massFlow<=0.f or pressureDiff<=0.f) { // If the pressure drops the wrong way, or if the supply temperatures don't match : close the valve.
//        kv = kvMin;
//    }
//    else {
    kv = computeIdealKv(rho, targetMassFlow, targetPressureDiff);
    if (kv<kvMin) { kv = kvMin; }
    else if (kv>kvMax) { kv = kvMax; }

    kv = kvPrev*(1.f-learningRate)+kv*learningRate;
//    }
    sumDeltaKv += abs(kv-kvPrev);
    sumKv += kv;
}

void Valve::setKvToMin(float& sumDeltaKv, float& sumKv, float const& learningRate) {
    float kvPrev = kv;
    kv = kvPrev*(1.f-learningRate)+computeKvMin()*learningRate;
    sumDeltaKv += abs(kv-kvPrev);
    sumKv += kv;
}







float PIDController::computeControlVariable(float const& desiredSetpoint, float const& processVariable, float const& kp, float const& ki, float const& kd) {
    float err = desiredSetpoint-processVariable;
    integralErr += err;
    float deriv = err-prevErr;
    prevErr = err;
    float controlVariable = kp*err + ki*integralErr + kd*deriv;
    return controlVariable;
}

float Carla::integrate(vector<float> const& f, float const& xmin, float const& xmax) {
    float integ = -0.5f*(f.front()+f.back());
    for (auto const& el : f) { integ += el; }
    integ *= (xmax-xmin)/(f.size()-1);
    return integ;
}

float Carla::selectAction() {
    float z = randomUniform(0,1);
    float cumulProba = 0.f;
    float probaToNextStep = 0.f;
    size_t i = 0;
    do {
        probaToNextStep = 0.5f*(probaDensity_fx[i]+probaDensity_fx[i+1])*(position_x[i+1]-position_x[i]);
        cumulProba += probaToNextStep;
        i++;
    } while (cumulProba<z and i<probaDensity_fx.size()-1);
    float newAction;
    cumulProba -= probaToNextStep;
    if (cumulProba<z) {
        i--;
        newAction = solveQuadratic(position_x[i], position_x[i+1], probaDensity_fx[i], probaDensity_fx[i+1], z-cumulProba);
    } else {
        newAction = xmax(); // Somehow the integral is bigger than one, just take xmax.
    }
    return newAction;
}

float Carla::solveQuadratic(float const& xi, float const& xip1, float const& fi, float const& fip1, float const& area) {
    float slope = (fip1-fi)/(xip1-xi);
    float a = slope*0.5f;
    float b = fi-xi*slope;
    float c = -( area + xi*b + xi*xi*a );
    float sol;
    if (a==0.f) { // Not quadratic equation.
        if (b==0.f) { // The affine function is zero everywhere (should not happen).
            sol = xi;
        } else {
            sol = -c/b;
        }
    }
    else { // Quadratic equation.
        float tmp = sqrt(b*b-4*a*c);
        sol = ( -b+tmp )/(2.f*a);
        if (sol<xi or sol>xip1) { // Check whether this root is the correct solution.
            sol = ( -b-tmp )/(2.f*a);
        }
    }
    return sol;
}

Carla::Carla(size_t const& n, float const& xmin, float const& xmax) : position_x(n), probaDensity_fx(n) {
    if (n<=1) { throw string ("Error in Carla::Carla(...) , n is <=1 but should be >1"); }
    for (size_t i=0; i<n; i++) {
        position_x[i] = xmin + (xmax - xmin)*i/(n-1);
        probaDensity_fx[i] = 1.f/(xmax - xmin);
    }
    actionChoice_r = 0.5f*(xmin+xmax);
}

float Carla::step(float const& beta) {
    float g_h = 0.3f, g_w = 0.02f;
    float a = g_h/(xmax()-xmin()), b = -0.5f/( g_w*g_w*(xmax()-xmin())*(xmax()-xmin()) );
    for (size_t i=0; i<position_x.size(); i++) {
        float dx = position_x[i]-actionChoice_r;
        float gauss_H = a*exp(b*dx*dx);
        probaDensity_fx[i] += beta*gauss_H;
    }
    float alpha = 1.f/integrate(probaDensity_fx, xmin(), xmax());
    for (auto & el : probaDensity_fx) { el *= alpha; }
    actionChoice_r = selectAction();
    return actionChoice_r;
}

PIDControllerCarla::PIDControllerCarla(size_t const& n, size_t const& r, float const& xmin, float const& xmax)
    : PIDController(), carlaKp(n, xmin, xmax), carlaKi(n, xmin, xmax), mm(r) { }
//    : PIDController(), carlaKp(n, xmin, xmax), carlaKi(n, xmin, xmax), carlaKd(n, xmin, xmax), mm(r) { }

float PIDControllerCarla::evaluatePerformanceBeta(float const& costJ) {
    float beta;
    if (mm.nbMemorized()>3) {
        float jMed = mm.findMedian(), jMin = mm.findMin();
        beta = min(  max( 0.f , (jMed-costJ)/(jMed-jMin) )  ,  1.f  );
    }
    else {
        beta = 0.f;
    }
    mm.addElement(costJ);
    return beta;
}

float PIDControllerCarla::computeControlVariable(float const& desiredSetpoint, float const& processVariable) {
    float cost = abs(desiredSetpoint-processVariable);
    float beta = evaluatePerformanceBeta(cost);
    float kp = carlaKp.step(beta);
    float ki = carlaKi.step(beta);
    //    float kd = carlaKd.step(beta);
        float kd = 0.f;
    return PIDController::computeControlVariable(desiredSetpoint, processVariable, kp, ki, kd);
}





Substation* Substation::createNewSubstation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr) {
    Substation* newSubstation;

    string type;
    if ( hdl.ToElement()->QueryStringAttribute("type", &type) )   { throw string("Error in the XML file: a Substation doesn't have attribute: 'type'."); }

    vector<string> possibleTypes = { "simple", "regulatedElement", "prosumer" };
    string possibleTypesString = possibleTypes[0];
    for (size_t i=1; i<possibleTypes.size(); i++) { possibleTypesString += ", "+possibleTypes[i]; }

    if (type==possibleTypes[0]) {
        newSubstation = new Substation(hdl, pBui, beginD, endD, pLogStr);
    } else if (type==possibleTypes[1]) {
        newSubstation = new RegulatedPressureSubstation(hdl, pBui, beginD, endD, pLogStr);
    } else if (type==possibleTypes[2]) {
        newSubstation = new ProsumerSubstation(hdl, pBui, beginD, endD, pLogStr);
    } else {
        throw string("Error in the XML file: a Substation has attribute type="+type+", which isn't valid. Possible types are : "+possibleTypesString);
    }

    return newSubstation;
}

Substation::Substation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr) : EnergyConversionSystem(beginD, endD, pLogStr), pBuilding(pBui), valve(nullptr) { // pid(200, 100, 0.f, 0.5f)

    if(hdl.ToElement() ){
        logStream <<"Substation" << endl << flush;

        // DEC it is linked to.
        unsigned int linkedNodeId;
        if ( hdl.ToElement()->QueryUnsignedAttribute("linkedNodeId", &linkedNodeId) ) { throw string("Error in the XML file: a Substation doesn't have attribute: 'linkedNodeId'."); }

        // Thermal power designed (maximal power during operation).
        if ( hdl.ToElement()->QueryFloatAttribute("designThermalPower", &designThermalPower) ) { throw string("Error in the XML file: a Substation doesn't have attribute: 'designThermalPower'."); }
        if( designThermalPower<=0. ) { throw string("Error in the XML file: a Substation has designThermalPower<=0."); }

        // Designed temperature difference on the primary network side.
        if ( hdl.ToElement()->QueryFloatAttribute("designTempDifference", &designTempDifference) ) { throw string("Error in the XML file: a Substation doesn't have attribute: 'designTempDifference'."); }
        if( designTempDifference<=0. ) { throw string("Error in the XML file: a Substation has designTempDifference<=0."); }

        // Epsilon of the heat exchanger, linked to its efficiency.
        if ( hdl.ToElement()->QueryFloatAttribute("designEpsilon", &designEpsilon) ) { throw string("Error in the XML file: a Substation doesn't have attribute: 'designEpsilon'."); }
        if( designEpsilon<=0. or designEpsilon>1. ) { throw string("Error in the XML file: a Substation has designEpsilon<=0. or designEpsilon>1."); }

        // Find and set pDEC and pNode that this Substation is linked to.
        size_t i=0;
        while( (i<pBui->getDistrict()->getnDECs()) && (not pBui->getDistrict()->getDEC(i)->addSubstationIfContainNode(this, linkedNodeId)) ) { i++; }
        if (i==pBui->getDistrict()->getnDECs()) {
            throw string("Error in the XML file: a Substation is trying to link to a NodePair that has not been found in DistrictEnergyCenter.");
        }

        valve = new Valve(this->maxKv(pDEC->getRho()));
    }
}

void Substation::setThermalPowerNeeded(double tPN) {
    EnergyConversionSystem::setThermalPowerNeeded(tPN);

    desiredMassFlow = mc( abs(thermalPowerNeeded), designThermalPower, designThermalPower/(designTempDifference*pDEC->getCp()) );

    // Set secondarySideInputTemp.
    if( thermalPowerNeeded>0. ) {
        // TODO, improve this block ! These formulas make too big approximations.
        float T_HS = getBuilding()->getHeatStockTemperature();
        if (getBuilding()->getDHWHeatStock()!=NULL) {
            float T_DHW = getBuilding()->getDHWStockT(), Pow_DHW = getBuilding()->getDHW_needs(), Pow_HS = getBuilding()->getHS_needs();
            if (Pow_HS+Pow_DHW==0.f) {
                secondarySideInputTemp = (T_HS+T_DHW)*0.5; // This may happen if imposed heat demands are used.
            } else {
                secondarySideInputTemp = (T_HS*Pow_HS+T_DHW*Pow_DHW)/(Pow_HS+Pow_DHW); // The temperature entering the heat exchanger secondary side, weighted average.
            }
        } else {
            secondarySideInputTemp = T_HS;
        }
    }
    else if ( thermalPowerNeeded<0. ) {
        secondarySideInputTemp = getBuilding()->getColdStockTemperature();
    }
    else {
        secondarySideInputTemp = pNode->getSupplyTemperature();
    }
}


double Substation::getThermalPower(double sourceTemp) {
    return thermalPowerExchanged; // Computations are done in ComputeHeatExchanged().
}


double Substation::getFuelConsumption(double time, double thermalPower, double sourceTemp) {
   return 0.;
}

double Substation::getElectricConsumption(double time, double thermalPower, double sourceTemp) {
    return 0.;
}


double Substation::mc(double thermalPowerNeeded, double designThermalPower, double m_c_n) {
    double ratio = thermalPowerNeeded / designThermalPower;
    if (ratio >= 0.75)
        return m_c_n;
    else if ( 0.5 <= ratio )
        return 0.75 * m_c_n;
    else if ( 0.25 <= ratio )
        return 0.5 * m_c_n;
    else
        return 0.25 * m_c_n;

// Fully proportional, but many changes in flows meaning longer simulation times.
//    if (ratio>1.0) { ratio = 1.0; }
//    else if (ratio<0.01) { ratio = 0.01; }
//    return ratio*m_c_n;
}

void Substation::computeHeatExchanged(float const& primarySideCp, float const& primarySideRho, float const& primarySideInputTemp, float const& primarySideMassFlow, float const& primarySidePressureDiff, Climate* pClim, unsigned int day, unsigned int hour) {

    if (primarySideMassFlow>0.f) {

        float secondarySideMassFlow = primarySideMassFlow; // TODO improve this ?
        float secondarySideCp = primarySideCp; // TODO: Improve this ? Make possibility to set in xml ?

        float sgnNeeds;
        if( thermalPowerNeeded>0. ) {
            sgnNeeds = 1.;
        }
        else if ( thermalPowerNeeded<0. ) {
            sgnNeeds = -1.;
        }
        else {
            sgnNeeds = 0.;
        }


        if (  (primarySideInputTemp-secondarySideInputTemp)*sgnNeeds<=0. ) { // If the temperature difference has opposite sign to the needs, then the heat tranfer cannot take place in the correct direction.
            thermalPowerExchanged = 0.f;

        } else { // In the normal case.
            double C_min = min(secondarySideCp * secondarySideMassFlow  , primarySideCp * primarySideMassFlow);
            double C_max = max(secondarySideCp * secondarySideMassFlow  , primarySideCp * primarySideMassFlow);
            double UA = 2 * designEpsilon * designThermalPower / ((1 - designEpsilon) * designTempDifference)
                    / (pow(designThermalPower / (designTempDifference*secondarySideCp*secondarySideMassFlow)  , 0.583) +
                       pow(designThermalPower / (designTempDifference*primarySideCp*primarySideMassFlow)  , 0.583));

            double HXEpsilon;
            if (C_min > 0.9999*C_max) { HXEpsilon = UA/C_min / (1 + UA/C_min); }
            else { HXEpsilon = (1 - exp(-(UA/C_min) * (1 - C_min/C_max))) /(1 - (C_min/C_max) * exp(-(UA/C_min) * (1 - C_min/C_max)) ); }

            thermalPowerExchanged = HXEpsilon  * C_min  * (primarySideInputTemp - secondarySideInputTemp ); // Positive for heating, negative for cooling.

            if ( thermalPowerExchanged*thermalPowerNeeded<0. ) { thermalPowerExchanged = 0.f; } // If the thermal power computed has opposite sign to the needs, set to zero (corresponds to closing the valve on secondary network).
            if ( abs(thermalPowerExchanged)>abs(designThermalPower) ) { thermalPowerExchanged = sgnNeeds*designThermalPower; } // Don't exceed design thermal power.
            if ( abs(thermalPowerExchanged)>abs(thermalPowerNeeded) ) { thermalPowerExchanged = thermalPowerNeeded; } // Don't exceed the thermal power demanded by consumer.
        }

        primarySideOutputTemp = computeOutputTemp(primarySideInputTemp, thermalPowerExchanged, primarySideMassFlow, primarySidePressureDiff, primarySideRho, primarySideCp);

    }
    else if (primarySideMassFlow<0.f) {
        thermalPowerExchanged = 0.f;
        primarySideOutputTemp = computeOutputTemp(primarySideInputTemp, thermalPowerExchanged, primarySideMassFlow, primarySidePressureDiff, primarySideRho, primarySideCp);
    }
    else { // primarySideMassFlow==0.f
        thermalPowerExchanged = 0.f;
    }
}

float Substation::maxKv(float const& rho) {
    return 3600.f*computeDesignMassFlow()*sqrt(100000.f/10000.f)/rho; // Kv in m^3/h, massFlow in kg/s. TODO improve this (this means there is at most designMassFlow for deltaP=0.1bar)
}

void Substation::updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& massFlowSupplyToReturn) {
    float currentPressureDiffControlElement = computePressureDiffControlElement(deltaP);
    float targetMassFLow, targetPressureDiff;
    if (massFlow>0.f) {
        targetMassFLow = desiredMassFlow;
        targetPressureDiff = currentPressureDiffControlElement;
        if(targetPressureDiff<=0.f) { targetPressureDiff = 1.0f; }
    }
    else {
        targetMassFLow = -desiredMassFlow*0.001f;
        targetPressureDiff = abs(deltaP)+1.f;
    } // If flow is going in the wrong direction, close the valve.
    valve->updateKv(rho, sumDeltaKv, sumKv, learningRate, targetMassFLow, targetPressureDiff);

    sumDeltaRpm += 0.f; sumRpm += 0.f;
}


void Substation::relativeErrorMassFlow(float const& massFlow, float& relErr, float& sumErr) {
    relErr = relativeErrorMassFlow(massFlow);
    sumErr += abs(massFlow-desiredMassFlow);
}


float Substation::computeDesignMassFlow() {
    return designThermalPower/(designTempDifference*pDEC->getCp());
}

void Substation::computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm) {
    float m_ = (*m);
    float deltaP_, dDeltaP_dm_;
    valve->computePressureDiffAndDerivative(m_, rho, deltaP_, dDeltaP_dm_);
    (*deltaP) = deltaP_;
    (*dDeltaP_dm) = dDeltaP_dm_;
}

float Substation::computeOutputTemp(float const& inputTemp, float const& thermalPowerExchanged, float const& massFlow, float const& pressureDiff, float const& rho, float const& cp) {
    float sgnMassFLow;
    if (massFlow>0.f) { sgnMassFLow = 1.f; }
    else if (massFlow<0.f) { sgnMassFLow = -1.f; }
    else { sgnMassFLow = 0.f; }

    float deltaTempHX = sgnMassFLow * (-thermalPowerExchanged) / (cp*massFlow); // Temperature increase due to heat exchanger [K].
    float deltaTempValve = sgnMassFLow * valve->computeTemperatureIncrease(pressureDiff, cp, rho); // Temperature increase du to valve. The pressure drop dissipates into thermal energy [K].
    return inputTemp + deltaTempHX + deltaTempValve;
}







RegulatingElement::RegulatingElement(TiXmlHandle hdl) {
    if ( hdl.ToElement()->QueryFloatAttribute("targetRegulatedPathPressureDiff", &targetRegulatedPathPressureDiff) ) { throw string("Error in the XML file: a RegulatingElement doesn't have attribute: 'targetRegulatedPathPressureDiff'."); }
    pressureDiff = targetRegulatedPathPressureDiff + randomUniform(-1,1);
}




RegulatedPressureSubstation::RegulatedPressureSubstation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr)
    : Substation(hdl, pBui, beginD, endD, pLogStr) {
    if (hdl.FirstChildElement("RegulatingElement").ToElement()) {
        regEle = new RegulatingElement(hdl.FirstChildElement("RegulatingElement"));
        pDEC->getNetwork()->addRegulatingElement(regEle);
    } else {
        throw string("Error in the XML file: a RegulatedPressureSubstation doesn't contain any RegulatingElement child element.");
    }
}

void RegulatedPressureSubstation::computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm) {
    // Case of substation with two edges : a valve and a differential pressure regulator.
    Substation::computePressureDiffAndDerivative(m, rho, deltaP, dDeltaP_dm);
    *(deltaP+1) = regEle->getPressureDiff();
    *(dDeltaP_dm+1) = 0.f;
}






void ProsumerSubstation::deleteDynAllocated() {
    if (temperatureSetpoint!=nullptr) { delete temperatureSetpoint; }
    if (pressureSetpoint!=nullptr) { delete pressureSetpoint; }
    if (pump!=nullptr) { delete pump; }
    if (pumpFlowControlValve!=nullptr) { delete pumpFlowControlValve; }
}

bool ProsumerSubstation::isHeatSource() {
    if (pBuilding->getHeatingUnit()==this) { return true; }
    else if (pBuilding->getCoolingUnit()==this) { return false; }
    else { throw string("Error in ProsumerSubstation::isHeatSource(), 'this' is neither heatingUnit nor coolingUnit."); }
}

float ProsumerSubstation::computeSolarThermalPower(float const& targetSupplyTemp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour) {
    float avgFluidTempInSolCollector = 0.5f*(targetSupplyTemp + primarySideReturnTemp); // Hypothesis that the temperature enters the panel at primarySideReturnTemp, exits at targetSupplyTemp.
    float solarThermalPower = pBuilding->getSolTherFracLeft() * Model::computeSolarThermalPower(pBuilding, pClim, day, hour, avgFluidTempInSolCollector, false); // Defined as positive to inject heat on the network.
    return solarThermalPower;
}

float ProsumerSubstation::minDesiredPumpFlow() {
//    return -0.01f*abs(pDEC->getNetwork()->computeMassFlowSupplyToReturn());
//    return -computeDesignMassFlow()*0.0001f;
    return -0.01f; // TODO: find a better value?
}

ProsumerSubstation::ProsumerSubstation(TiXmlHandle hdl,  Building* pBui, unsigned int beginD, unsigned int endD, ostream* pLogStr)
    : Substation(hdl, pBui, beginD, endD, pLogStr), temperatureSetpoint(nullptr), pressureSetpoint(nullptr), pump(nullptr), pumpFlowControlValve(nullptr), producerModeOn(false), pumpElectricPower(0.f) { // , turnOfPumpFlowControlValve(false) {

    try {
        if (hdl.FirstChildElement("TemperatureSetpoint").ToElement()) {
            temperatureSetpoint = TemperatureSetpoint::createNewTemperatureSetpoint(hdl.FirstChildElement("TemperatureSetpoint")); // Dynamical allocation.
        } else {
            throw string("Error in the XML file: a ProsumerSubstation doesn't contain any TemperatureSetpoint child element.");
        }

        if (hdl.FirstChildElement("PressureSetpoint").ToElement()) {
            pressureSetpoint = PressureSetpoint::createNewPressureSetpoint(hdl.FirstChildElement("PressureSetpoint")); // Dynamical allocation.
        } else {
            throw string("Error in the XML file: a ProsumerSubstation doesn't contain any PressureSetpoint child element.");
        }

        if (hdl.FirstChildElement("Pump").ToElement()) {
            pump = new Pump(hdl.FirstChildElement("Pump"));
        } else {
            throw string("Error in the XML file: a ProsumerSubstation doesn't contain any Pump child element.");
        }

        float kvMaxFlowControlValve;
        if ( hdl.ToElement()->QueryFloatAttribute("kvMaxFlowControlValve", &kvMaxFlowControlValve) ) { throw string("Error in the XML file: a Substation of type prosumer doesn't have attribute: 'kvMaxFlowControlValve' [m^3/h]."); }
        if (kvMaxFlowControlValve<=0) { throw string("Error in the XML file: a Substation of type prosumer has kvMaxFlowControlValve<=0, it must be positive by convention."); }
        pumpFlowControlValve = new Valve(kvMaxFlowControlValve);

    } catch (...) {
        deleteDynAllocated();
        throw;
    }

}

void ProsumerSubstation::setThermalPowerNeeded(double tPN) {
    Substation::setThermalPowerNeeded(tPN);
    producerModeOn = false;
}

double ProsumerSubstation::getThermalPower(double sourceTemp) {
    if (producerModeOn) {
        return 0.f; // thermalPowerExchanged variable contains the heat fed into the network, not consumed by substation.
    } else { // Consumer or passive mode.
        return Substation::getThermalPower(sourceTemp);
    }
}

double ProsumerSubstation::getElectricConsumption(double time, double thermalPower, double sourceTemp) {
    return time*pumpElectricPower;
}

void ProsumerSubstation::updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate, float const& massFlowSupplyToReturn) {
    if (not producerModeOn) { // Use valve (consumer or passive).
        Substation::updateControlVariable(massFlow, deltaP, rho, sumDeltaRpm, sumRpm, sumDeltaKv, sumKv, learningRate, massFlowSupplyToReturn);
    } else { // Use pump and pumpFlowControlValve.
//        if (turnOfPumpFlowControlValve) { // Update only the pump flow control valve.

            float targetMassFlow, targetPressureDiff;
            if (massFlow>=0.f) { // Flow is going the wrong direction (supply to return).
                targetMassFlow = desiredMassFlow*0.001f; // Shut the valve.
                targetPressureDiff = abs(deltaP)+1.f; // To avoid sqrt(negative number).

            } else { // Flow is going the correct direction (return to supply).
                targetMassFlow = -desiredMassFlow;
                float pressureDiffThroughPumpAtTargetFlow, dpressureDiffdmThroughPumpAtTargetFlow;
                pump->computePressureDiff(targetMassFlow, pressureDiffThroughPumpAtTargetFlow, dpressureDiffdmThroughPumpAtTargetFlow);
                targetPressureDiff = - pressureDiffThroughPumpAtTargetFlow - deltaP*targetMassFlow*targetMassFlow/(massFlow*massFlow);
                if(targetPressureDiff<=0.f) { targetPressureDiff = 1.0f; }
            }
            pumpFlowControlValve->updateKv(rho, sumDeltaKv, sumKv, learningRate, targetMassFlow, targetPressureDiff);

//        } else { // Update only the pump.


                float setpointPressureDiff = pressureSetpoint->computeTargetPressureDiff(-massFlow);

                if (massFlow>0.f) { // Flow is going the wrong direction (supply to return).
                    targetMassFlow = massFlow*0.001f; // Try to reverse the direction of the flow.
                }
                else { // Flow is going the correct direction (return to supply).
//                    targetMassFlow = -massFlow; // Keep same mass flow.
                    targetMassFlow = -massFlow*sqrt(abs(setpointPressureDiff/deltaP)); // Keep mass flow if same resistance from rest of network.
                }

                float pressureDiffThroughValveAtTargetFlow, dpressureDiffdmThroughValveAtTargetFlow;
                pumpFlowControlValve->computePressureDiffAndDerivative(targetMassFlow, rho, pressureDiffThroughValveAtTargetFlow, dpressureDiffdmThroughValveAtTargetFlow);

                targetPressureDiff = setpointPressureDiff - pressureDiffThroughValveAtTargetFlow;
                if (targetPressureDiff>0) { targetPressureDiff = -1.f; }

                pump->updateRpms(sumDeltaRpm, sumRpm, learningRate, targetMassFlow, targetPressureDiff);
//        }
//        turnOfPumpFlowControlValve = not turnOfPumpFlowControlValve;
    }
}

void ProsumerSubstation::computeHeatExchanged(float const& primarySideCp, float const& primarySideRho, float const& primarySideInputTemp, float const& primarySideMassFlow, float const& primarySidePressureDiff, Climate* pClim, unsigned int day, unsigned int hour) {

    if (not producerModeOn) { // Consumer or passive mode.
        pumpElectricPower = 0.f;
        Substation::computeHeatExchanged(primarySideCp, primarySideRho, primarySideInputTemp, primarySideMassFlow, primarySidePressureDiff, pClim, day, hour);
    }
    else { // Producer mode.

        // Compute valve or pump temperature increase.
        float controlElementThermalPower; // "cp*(mass flow tail to head)*(Temperature_head - Temperature_tail)", positive when temperature increases in the direction of the flow [W]. // TODO: Use this!
        if (desiredMassFlow>=0.f) { // Valve is used.
            controlElementThermalPower = primarySideCp*primarySideMassFlow*valve->computeTemperatureIncrease(primarySidePressureDiff, primarySideCp, primarySideRho);
            pumpElectricPower = 0.f;
        } else { // Pump is used.
            pump->computeElectricAndThermalPower(-primarySidePressureDiff, -primarySideMassFlow, primarySideRho, pumpElectricPower, controlElementThermalPower); // TODO, this formula needs to be updated to take into account pumpFlowControlValve
        }

        if (isHeatSource()) {
            if (primarySideMassFlow<0.f) { // If the effective mass flow is from return to supply.

                float targetSupplyTemp = temperatureSetpoint->computeTargetTemperature(pClim, day, hour);
                float totalThermalPowerNeeded = primarySideCp*(-primarySideMassFlow)*(targetSupplyTemp - primarySideInputTemp);
                float solarThermalPowerNeeded = totalThermalPowerNeeded;
//                float solarThermalPowerNeeded = totalThermalPowerNeeded - controlElementThermalPower; // Pump or valve preheats the fluid. // TODO : Implement this line instead of the previous one.


                if (solarThermalPowerNeeded<0.f) { // The fluid's temperature will exceed the target temperature.
                    thermalPowerExchanged = 0.f;
                    primarySideOutputTemp = primarySideInputTemp;
                } else {

                    float solarThermalPowerAvailable = computeSolarThermalPower(targetSupplyTemp, primarySideInputTemp, pClim, day, hour);

                    if (solarThermalPowerNeeded<=solarThermalPowerAvailable) { // Sufficient solar thermal to heat to the target temperature.
                        thermalPowerExchanged = -solarThermalPowerNeeded; // Only use the necessary solar thermal power, to not exceed the target temperature. (Minus sign due to conventions)
                        primarySideOutputTemp = targetSupplyTemp;

                    } else { // Not enough solar thermal to heat to the target temperature.

                        float outputTempHigh = targetSupplyTemp;
                        float outputTempLow = primarySideInputTemp + solarThermalPowerAvailable/(-primarySideMassFlow*primarySideCp);
                        float outputTemp;
                        float newOutputTemp;
                        do { // Binary search.
                            outputTemp = 0.5f*(outputTempHigh+outputTempLow);
                            solarThermalPowerAvailable = computeSolarThermalPower(outputTemp, primarySideInputTemp, pClim, day, hour);
                            newOutputTemp = primarySideInputTemp + solarThermalPowerAvailable/(-primarySideMassFlow*primarySideCp);
                            if (outputTemp>newOutputTemp) { outputTempHigh = outputTemp;  }
                            else { outputTempLow = outputTemp;  }
                        } while (abs(outputTemp-newOutputTemp)>0.05f and outputTempHigh-outputTempLow>0.05f);

                        thermalPowerExchanged = -solarThermalPowerAvailable; // Use all solar thermal available. (Minus sign due to conventions)
                        primarySideOutputTemp= newOutputTemp;
                    }
                    // TODO : Impose a maximal designThermalPower ?
                }
            }
            else { // Flow is from supply to return, which is the wrong direction, cannot feed in.
                thermalPowerExchanged = 0.f;
                primarySideOutputTemp = primarySideInputTemp;
            }
        }
        else {
            throw string("ProsumerSubstation for cooling are not yet implemented.");
        }
    }
}

void ProsumerSubstation::computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm) {
    if (not producerModeOn) { // Use the valve (consumer or passive).
        Substation::computePressureDiffAndDerivative(m, rho, deltaP, dDeltaP_dm);
    } else { // Use the pump (producer).
        float m_ = -(*m); // Convention of the pump and substation are opposite.
        float deltaP_pump, dDeltaP_dm_pump, deltaP_valve, dDeltaP_dm_valve;
        pump->computePressureDiff(m_, deltaP_pump, dDeltaP_dm_pump);
        pumpFlowControlValve->computePressureDiffAndDerivative(m_, rho, deltaP_valve, dDeltaP_dm_valve);
        *deltaP = -deltaP_pump-deltaP_valve; // Convention of the pump and substation are opposite.
        *dDeltaP_dm = dDeltaP_dm_pump+dDeltaP_dm_valve; // Has double sign negation, so stays the same.
    }
}

void ProsumerSubstation::updateDesiredMassFlow(float const& cp, float const& primarySideReturnTemp, Climate* pClim, unsigned int day, unsigned int hour) {

    if (thermalPowerNeeded==0.f) { // Not consmuming, check if can be producer (feed into the network).
        producerModeOn = true; // before debug
        if (isHeatSource()) {

            // TODO : Improve by taking into account the heating from the pump or valve.

            float targetSupplyTemp = temperatureSetpoint->computeTargetTemperature(pClim, day, hour);
            float targetFeedInPower = computeSolarThermalPower(targetSupplyTemp, primarySideReturnTemp, pClim, day, hour);
            float tempIncrease = targetSupplyTemp - primarySideReturnTemp; // If temperature increase is too small, or negative : the desired mass flow will be too negative, or impossible.

            if (0.f<targetFeedInPower and 0.05f<tempIncrease) {
                // If desiredMassFlow is too negative, put a cap? eg minus the design mass flow? Or simply let the pump reach it's maximum capacity?
                desiredMassFlow = -targetFeedInPower/(cp*tempIncrease);
            }
            else {
                // Don't try to inject heat, just ensure a small flow.
                desiredMassFlow = minDesiredPumpFlow(); // Pump only a very small flow. // before debug
            }

        }
        else {
            throw string("ProsumerSubstation for cooling are not yet implemented.");
        }
    }
    else { // Substation in consumer mode, desiredMassFlow is left unchanged.
        producerModeOn = false;
    }
}

void ProsumerSubstation::setProsumerSolarThermal() {
    if (producerModeOn) {
        float oldValue = pBuilding->getSolarThermalProduction();
        float newValue = oldValue + Model::dt * (-thermalPowerExchanged); // Minus sign because of convention.
        pBuilding->eraseSolarThermalProduction_back();
        pBuilding->setSolarThermalProduction(newValue);
    } // Else it is consumer mode and the exchanged heat is consumed not injected solar thermal.
}

void ProsumerSubstation::relativeErrorMassFlow(float const& massFlow, float& relErr, float& sumErr) {
    if (producerModeOn and (desiredMassFlow==minDesiredPumpFlow()) and (massFlow<0.f) and pump->nIsAlmostMax() and pumpFlowControlValve->kvIsAlmostMin()) { // The pump is already at minimum power, we cannot go lower. (to avoid getting stuck). // before debug
        relErr = 0.f;
        sumErr += 0.f;
    } else {
        Substation::relativeErrorMassFlow(massFlow, relErr, sumErr);
    }
}






TemperatureSetpoint* TemperatureSetpoint::createNewTemperatureSetpoint(TiXmlHandle hdl) {
    TemperatureSetpoint* newTemperatureSetpoint;

    string type;
    if ( hdl.ToElement()->QueryStringAttribute("type", &type) )   { throw string("Error in the XML file: a TemperatureSetpoint doesn't have attribute: 'type'."); }

    if (type=="constant") {
        newTemperatureSetpoint = new ConstantTemperatureSetpoint(hdl);
    } else if (type=="affine") {
        newTemperatureSetpoint = new AffineTemperatureSetpoint(hdl);
    } else if (type=="affineWinterConstantSummer") {
        newTemperatureSetpoint = new AffineWinterConstantSummerSetpoint(hdl);
    } else if (type=="imposedValuesOrConstant") {
        newTemperatureSetpoint = new ImposedValuesOrConstantSetpoint(hdl);
    }
    else {
        throw string("Error in the XML file: a TemperatureSetpoint has attribute type="+type+", which isn't valid.");
    }

    return newTemperatureSetpoint;
}


ConstantTemperatureSetpoint::ConstantTemperatureSetpoint(TiXmlHandle hdl) : TemperatureSetpoint() {
    if ( hdl.ToElement()->QueryFloatAttribute("targetSupplyTemp", &targetSupplyTemp) ) { throw string("Error in the XML file: a TemperatureSetpoint of type constant doesn't have attribute: 'targetSupplyTemp'."); }
}



AffineTemperatureSetpoint::AffineTemperatureSetpoint(TiXmlHandle hdl) : TemperatureSetpoint() {
    if ( hdl.ToElement()->QueryFloatAttribute("lowExtTemp", &lowExtTemp) ) { throw string("Error in the XML file: a TemperatureSetpoint of type affine doesn't have attribute: 'lowExtTemp'."); }
    if ( hdl.ToElement()->QueryFloatAttribute("highExtTemp", &highExtTemp) ) { throw string("Error in the XML file: a TemperatureSetpoint of type affine doesn't have attribute: 'highExtTemp'."); }
    if ( hdl.ToElement()->QueryFloatAttribute("lowExtTempSupplyTemp", &lowExtTempSupplyTemp) ) { throw string("Error in the XML file: a TemperatureSetpoint of type affine doesn't have attribute: 'lowExtTempSupplyTemp'."); }
    if ( hdl.ToElement()->QueryFloatAttribute("highExtTempSupplyTemp", &highExtTempSupplyTemp) ) { throw string("Error in the XML file: a TemperatureSetpoint of type affine doesn't have attribute: 'highExtTempSupplyTemp'."); }
}

float AffineTemperatureSetpoint::computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) {
    float extTemp = avgExtTempLast24Hours(pClimate, day, hour);
    float targetTemp;
    if (extTemp<lowExtTemp) {
        targetTemp = lowExtTempSupplyTemp;
    } else if (extTemp>highExtTemp) {
        targetTemp = highExtTempSupplyTemp;
    } else {
        targetTemp = lowExtTempSupplyTemp + (highExtTempSupplyTemp-lowExtTempSupplyTemp) * (extTemp-lowExtTemp)/(highExtTemp-lowExtTemp);
    }
    return targetTemp;
}

float AffineTemperatureSetpoint::avgExtTempLast24Hours(Climate* pClimate, unsigned int day, unsigned int hour) {
    float extTemp = pClimate->getToutCelsius(day,hour);
    for (size_t i=0 ; i<23 ; i++) {
        if (hour==1) {
            hour = 24;
            if (day==1) { day = 365; }
            else { day--; }
        }
        else { hour--; }
        extTemp += pClimate->getToutCelsius(day,hour);
    }
    extTemp *= 0.041667; // 0.041667 = 1/24
    return extTemp;
}



AffineWinterConstantSummerSetpoint::AffineWinterConstantSummerSetpoint(TiXmlHandle hdl) : AffineTemperatureSetpoint(hdl), notYetActivated(true) {
    if ( hdl.ToElement()->QueryFloatAttribute("startSummerTempThreshold", &startSummerTempThreshold) ) { throw string("Error in the XML file: a TemperatureSetpoint of type affineWinterConstantSummer doesn't have attribute: 'startSummerTempThreshold'."); }
    if ( hdl.ToElement()->QueryFloatAttribute("endSummerTempThreshold", &endSummerTempThreshold) ) { throw string("Error in the XML file: a TemperatureSetpoint of type affineWinterConstantSummer doesn't have attribute: 'endSummerTempThreshold'."); }
    if ( startSummerTempThreshold<endSummerTempThreshold ) { throw string("Error in the XML file: a TemperatureSetpoint of type affineWinterConstantSummer has startSummerTempThreshold<endSummerTempThreshold."); }
}

float AffineWinterConstantSummerSetpoint::computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) {
    if (notYetActivated) { // Only useful when function is called for the first time.
        notYetActivated = false;
        if (day>=136 and day<=288) { summerModeOn = true; }
        else { summerModeOn = false; }
    }

    float extTemp = AffineTemperatureSetpoint::avgExtTempLast24Hours(pClimate, day, hour);
    if (summerModeOn) { // Hysteresis to change between winter and summer modes.
        if (extTemp<endSummerTempThreshold) { summerModeOn = false; }
    } else {
        if (extTemp>startSummerTempThreshold) { summerModeOn = true; }
    }

    if (summerModeOn) { // Determine target temperature.
        return highExtTempSupplyTemp;
    } else {
        return AffineTemperatureSetpoint::computeTargetTemperature(pClimate, day, hour);
    }
}



bool ImposedValuesOrConstantSetpoint::hasImposedValue(unsigned int day, unsigned int hour, float& retValue) {
    string dayHour = "d"+toString(day)+"h"+toString(hour);
    if (imposedValues.count(dayHour.c_str())==0) {
        return false;
    } else {
        retValue = imposedValues[dayHour];
        return true;
    }
}

ImposedValuesOrConstantSetpoint::ImposedValuesOrConstantSetpoint(TiXmlHandle hdl) : TemperatureSetpoint() {
    if ( hdl.ToElement()->QueryFloatAttribute("constantTempIfNoImposed", &constantTempIfNoImposed) ) { throw string("Error in the XML file: a TemperatureSetpoint of type imposedValuesOrConstant doesn't have attribute: 'constantTempIfNoImposed'."); }

    if(hdl.FirstChildElement("ImposedValues").ToElement()) {
        TiXmlAttribute* attrib = hdl.FirstChildElement("ImposedValues").ToElement()->FirstAttribute();
        while ( attrib!=NULL ) {
            string dayHour = attrib->Name();
            double temp;
            if ( attrib->QueryDoubleValue(&temp) ) { throw string("Error in the XML file: a TemperatureSetpoint of type imposedValuesOrConstant has attribute "+dayHour+" which isn't a number."); }
            if (imposedValues.count(dayHour) != 0) {
                throw string("Error in the XML file: a TemperatureSetpoint of type imposedValuesOrConstant has on multiple times the attribute "+dayHour+", there must be at most one.");
            } else {
                imposedValues[dayHour] = temp;
            }
            attrib = attrib->Next();
        }
    } else {
        throw string("Error in the XML file: a TemperatureSetpoint of type imposedValuesOrConstant doesn't have child element: 'ImposedValues'.");
    }
}

float ImposedValuesOrConstantSetpoint::computeTargetTemperature(Climate* pClimate, unsigned int day, unsigned int hour) {
    float temp;
    if (not hasImposedValue(day, hour, temp)) {
        temp = constantTempIfNoImposed;
    }
    return temp;
}






PressureSetpoint* PressureSetpoint::createNewPressureSetpoint(TiXmlHandle hdl) {
    PressureSetpoint* newPressureSetpoint;

    string type;
    if ( hdl.ToElement()->QueryStringAttribute("type", &type) )   { throw string("Error in the XML file: a PressureSetpoint doesn't have attribute: 'type'."); }

    if (type=="constant") {
        newPressureSetpoint = new ConstantPressureSetpoint(hdl);
    } else if (type=="affine") {
        newPressureSetpoint = new AffinePressureSetpoint(hdl);
    } else {
        throw string("Error in the XML file: a PressureSetpoint has attribute type="+type+", which isn't valid.");
    }

    return newPressureSetpoint;
}


ConstantPressureSetpoint::ConstantPressureSetpoint(TiXmlHandle hdl) : PressureSetpoint() {
    if ( hdl.ToElement()->QueryFloatAttribute("targetPressureDiff", &targetPressureDiff) ) { throw string("Error in the XML file: a PressureSetpoint of type doesn't have attribute: 'targetPressureDiff'."); }
    if (targetPressureDiff<=0) { throw string("Error in the XML file: a PressureSetpoint of type constant has targetPressureDiff<=0, it must be positive by convention."); }
    targetPressureDiff = -targetPressureDiff;
}


AffinePressureSetpoint::AffinePressureSetpoint(TiXmlHandle hdl) : PressureSetpoint() {
    bool keepGoing = true;
    unsigned int i = 1;
    while (keepGoing) {
        string name = "massFlow"+to_string(i);
        float value;
        if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &value) ) {
            keepGoing = false;
        } else {
            massFlows.push_back(value);
            name = "pressureDiff"+to_string(i);
            if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &value) ) {
                throw string("Error in the XML file: a PressureSetpoint of type affine has attribute 'massFlow"+to_string(i)+"' but not '"+name+"'.");
            } else {
                pressureDiffs.push_back(value);
            }
        }
        i++;
    }
    if (i<=2) {
        throw string("Error in the XML file: a PressureSetpoint of type affine must have at least 2 points, so attributes: 'massFlow1', 'pressureDiff1', 'massFlow2' and 'pressureDiff2'.");
    }
    for (size_t i=1; i<massFlows.size(); i++) {
        if (massFlows[i]<=massFlows[i-1]) {
            throw string("Error in the XML file: a PressureSetpoint of type affine has 'massFlow"+to_string(i)+"'<='massFlow"+to_string(i-1)+". But massFlows should be strictly increasing.");
        }
    }

    for (auto const& m : massFlows) {
        if (m<=0.f) { throw string("Error in the XML file: a PressureSetpoint of type affine has a massFlow<=0, by convention they must be positive."); }
    }

    for (auto& dP : pressureDiffs) {
        if (dP<=0.f) { throw string("Error in the XML file: a PressureSetpoint of type affine has a pressureDiff<=0, by convention they must be positive."); }
        else { dP = -dP; } // In the xml they are positive for userfriendlyness, but in the code they are negative when the pump increases the pressure.
    }
}

float AffinePressureSetpoint::computeTargetPressureDiff(float const& massFlow) {
    float deltaP;
    if (massFlow<=massFlows[0]) { // Flat.
        deltaP = pressureDiffs[0];
    } else if (massFlow>massFlows[massFlows.size()-1]) { // Flat.
        deltaP = pressureDiffs[massFlows.size()-1];
    } else { // Affine by parts.
        size_t i = 0;
        while (massFlow>massFlows[i+1]) { i++; }
        deltaP = pressureDiffs[i]+(massFlow-massFlows[i])*(pressureDiffs[i+1]-pressureDiffs[i])/(massFlows[i+1]-massFlows[i]);
    }
    return deltaP;
}






Storage* Storage::createNewStorage(TiXmlHandle hdl) {
    Storage* newStorage;

    string type;
    if ( hdl.ToElement()->QueryStringAttribute("type", &type) )   { throw string("Error in the XML file: a Storage doesn't have attribute: 'type'."); }

    if (type=="simple") {
        newStorage = new SimpleStorage(hdl);
    } else {
        throw string("Error in the XML file: a Storage has attribute type="+type+", which isn't valid.");
    }

    return newStorage;
}






SimpleStorage::SimpleStorage(TiXmlHandle hdl) : Storage(), temperatureRecord(0) {
    if ( hdl.ToElement()->QueryFloatAttribute("initialTemperature", &temperature) ) { throw string("Error in the XML file: a Storage of type simple doesn't have attribute: 'initialTemperature' [degree C]."); }
    if ( temperature<=-273.15f ) { throw string("Error in the XML file: a Storage of type simple has initialTemperature<=-273.15, which isn't physically possible."); }
    if ( hdl.ToElement()->QueryFloatAttribute("heatCapacity", &heatCapacity) ) { throw string("Error in the XML file: a Storage of type simple doesn't have attribute: 'heatCapacity' [J/kg]."); }
    if ( heatCapacity<=0.f ) { throw string("Error in the XML file: a Storage of type simple has heatCapacity<=0, which isn't physically possible."); }
}

float SimpleStorage::computeOutputTemperature(bool storageHeatsUp, float const& inputTemp, float const& m, float const& cp) {
    // m must be strictly positive
    float outputTemp;
    float tempDiff = inputTemp - temperature;
    if ( (storageHeatsUp and tempDiff>0.f)  or  ((not storageHeatsUp) and tempDiff<0.f) ) {
        float x = m*cp*Model::dt/heatCapacity;
        outputTemp = inputTemp - tempDiff * (1.f-exp(-x)) / x;
    }
    else {
        outputTemp = inputTemp; // Don't exchange heat (eg by using a bypass), since the temperature difference isn't correct.
    }
    return outputTemp;
}

void SimpleStorage::confirmStoredHeat(float const& tempDiff, float const& m, float const& cp) {
    // m must be positive
    temperature = temperature - tempDiff * m * cp * Model::dt / heatCapacity;
}

void SimpleStorage::writeTHHeaderText(fstream& textFile, string prefix) {
//    Storage::writeTHHeaderText(textFile, string);
    textFile << prefix << ":StorageTemp(celsius)" <<"\t";
}

void SimpleStorage::writeTHResultsText(fstream& textFile, unsigned int i) {
//    Storage::writeTHResultsText(textFile, i);
    textFile << fixed << setprecision(2) << getTemperature(i) <<"\t";
}










void ThermalStation::deleteDynAllocated() {
    if (pump!=nullptr) { delete pump; }
    if (ecs!=nullptr) { delete ecs; }
    if (temperatureSetpoint!=nullptr) { delete temperatureSetpoint; }
    if (pressureSetpoint!=nullptr) { delete pressureSetpoint; }
}

ThermalStation* ThermalStation::createNewThermalStation(TiXmlHandle hdl, Network* net, ostream* pLogStr) {
    ThermalStation* newThermalStation;

    string type;
    if ( hdl.ToElement()->QueryStringAttribute("type", &type) )   { throw string("Error in the XML file: a ThermalStation doesn't have attribute: 'type'."); }

    if (type=="simple") {
        newThermalStation = new ThermalStation(hdl, net, pLogStr);
    } else if (type=="seasonalStorageHeating") {
        newThermalStation = new SeasonalStorageHeatingThermalStation(hdl, net, pLogStr);
    } else {
        throw string("Error in the XML file: a ThermalStation has attribute type="+type+", which isn't valid.");
    }

    return newThermalStation;
}

ThermalStation::ThermalStation(TiXmlHandle hdl, Network* net, ostream* pLogStr) : pump(nullptr), ecs(nullptr), temperatureSetpoint(nullptr), pressureSetpoint(nullptr), logStream(std::cout.rdbuf()){ // pid(200, 100, 0.f, 0.5f),
    try {
        logStream.rdbuf(pLogStr->rdbuf()); // Add the logStream

        unsigned int beginDay, endDay;
        if ( hdl.ToElement()->Attribute("beginDay") == NULL ) beginDay = 1;
        else beginDay = to<unsigned int>(hdl.ToElement()->Attribute("beginDay"));
        if ( hdl.ToElement()->Attribute("endDay") == NULL )   endDay   = 365;
        else endDay   = to<unsigned int>(hdl.ToElement()->Attribute("endDay"));

        unsigned int linkedNodeId;
        string name = "linkedNodeId";
        if ( hdl.ToElement()->QueryUnsignedAttribute(name.c_str(), &linkedNodeId) ) { throw string("Error in the XML file: a ThermalStation doesn't have attribute: '"+name+"'."); }
        node = net->pointerOfThermalStationNodeWithId(linkedNodeId);
        if (node==nullptr) { throw string("Error in the XML file: a ThermalStation is trying to link to a ThermalStationNodePair id="+to_string(linkedNodeId)+" that has not been found in the DistrictEnergyCenter id="+to_string(net->getDEC()->getId())+"."); }
        node->setThermalStation(this);

        if (hdl.FirstChildElement("Pump").ToElement()) {
            pump = new Pump(hdl.FirstChildElement("Pump"));
        } else {
            throw string("Error in the XML file: a ThermalStation doesn't contain any Pump child element.");
        }

        // Load different possible types of EnergyConversionSystem.
        bool ecsFound = false;
        string errorMsg = "Error in the XML file: a DistrictEnergyCenter has a ThermalStation with multiple EnergyConversionSystems, there must be only one (eg only a boiler, not boiler and heat pump).";
        if (hdl.FirstChildElement("Boiler").ToElement()) {
            if ( ecsFound ) { delete ecs; throw errorMsg; }
            ecsFound = true;
            ecs = new Boiler(hdl.FirstChildElement("Boiler"),beginDay,endDay,&(logStream));
        }
        if (hdl.FirstChildElement("HeatPump").ToElement()) {
            if ( ecsFound ) { delete ecs; throw errorMsg; }
            ecsFound = true;
            ecs = new HeatPump(hdl.FirstChildElement("HeatPump"),beginDay,endDay,&(logStream));
        }
        if (hdl.FirstChildElement("CHP").ToElement()) {
            if ( ecsFound ) { delete ecs; throw errorMsg; }
            ecsFound = true;
            ecs = new CoGeneration(hdl.FirstChildElement("CHP"),beginDay,endDay,&(logStream));
        }
        if (hdl.FirstChildElement("CHP-HP").ToElement()) {
            if ( ecsFound ) { delete ecs; throw errorMsg; }
            ecsFound = true;
            ecs = new CoGenerationHeatPump(hdl.FirstChildElement("CHP-HP"),beginDay,endDay,&(logStream));
        }
        if ( not ecsFound ) {
            throw string("Error in the XML file: in a DistrictEnergyCenter, a ThermalStation, that doesn't contain any EnergyConversionSystem.");
        }


        if (hdl.FirstChildElement("TemperatureSetpoint").ToElement()) {
            temperatureSetpoint = TemperatureSetpoint::createNewTemperatureSetpoint(hdl.FirstChildElement("TemperatureSetpoint")); // Dynamical allocation.
        } else {
            throw string("Error in the XML file: a ThermalStation doesn't contain any TemperatureSetpoint child element.");
        }


        if (hdl.FirstChildElement("PressureSetpoint").ToElement()) {
            pressureSetpoint = PressureSetpoint::createNewPressureSetpoint(hdl.FirstChildElement("PressureSetpoint")); // Dynamical allocation.
        } else {
            throw string("Error in the XML file: a ThermalStation doesn't contain any PressureSetpoint child element.");
        }

    } catch (...) {
        ThermalStation::deleteDynAllocated();
        throw;
    }
}

void ThermalStation::computeThermal(float pressureDiff, float massFlow, float rho, float cp,  float inputTemp, Climate* pClimate, unsigned int day, unsigned int hour) {

    float thermalPowerNeeded;

    float thermalPowerPump;
    pump->computeElectricAndThermalPower(pressureDiff, massFlow, rho, pumpPower, thermalPowerPump);

    if (massFlow>0.f) {
        // Energy conversion system
        float targetTemp = temperatureSetpoint->computeTargetTemperature(pClimate, day, hour);
        float totalThermalPowerNeeded = cp*massFlow*(targetTemp - inputTemp);
        thermalPowerNeeded = totalThermalPowerNeeded - thermalPowerPump; // Pump preheats the fluid. We can end up in situation where the thermal station usually heats the fluid, but the pump friction heats more than the target temperature, meaning that we ask the boiler for a negative power...
    }
    else {
        thermalPowerNeeded = 0.f;
    }

    computeThermalPowerProvided(thermalPowerNeeded, pClimate, day, hour);

    if (massFlow!=0.f) {
        float deltaT = (thermalPowerProvided+thermalPowerPump) / (cp * abs(massFlow));
        outputTemperature = inputTemp + deltaT;
    } // Else the output temperature is not changed. (this avoids division by zero)
}

void ThermalStation::computePressureDiff(float const& massFlow, float const& rho, float& deltaP, float& dDeltaP_dm) {
    pump->computePressureDiff(massFlow, deltaP, dDeltaP_dm);
}

void ThermalStation::updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate) {
    float targetMassFlow;
    if (massFlow>=0.f) { targetMassFlow = massFlow; }
    else { targetMassFlow = abs(massFlow)*0.001f; } // If negative mass flow, try to get it positive? TODO improve? How to chose target mass flow?
    pump->updateRpms(sumDeltaRpm, sumRpm, learningRate, targetMassFlow, computeTargetPressureDiff(targetMassFlow));

    sumDeltaKv += 0.f; sumKv += 0.f;
}


void ThermalStation::relativeErrorPressureDiff(float const& massFlow, float const& pressureDiff, float& relErr, float& sumErrP, float& sumDesiredP) {
    float goal = -computeTargetPressureDiff(massFlow); // Minus signs to have positive pressures.
    relErr = (-pressureDiff-goal)/goal;
    sumErrP += abs(-pressureDiff-goal);
    sumDesiredP += goal;
}


void ThermalStation::setMachinePower_FuelConsumption_ElectricConsumption(Climate* pClim, unsigned int day, unsigned int hour) {
    // TODO, improve this, make the heatpump be able to access the temperature on its own.
    // the heatPumpSrcTemp depends on the type of heatpump (air or ground vertical/horizontal)
    double heatPumpSrcTemp = pClim->getToutCelsius(day,hour); // air temperature
    float z0=0.f,z1=0.f,alpha=0.f;
    if ( ecs->getGround(z0,z1,alpha) ) {
        heatPumpSrcTemp = pClim->getTgroundCelsius(day,hour,z0,alpha,z1);
    }
    setMachinePower(getThermalPowerProvided());
    setFuelConsumption( ecs->getFuelConsumption(double (Model::dt), getThermalPowerProvided(), heatPumpSrcTemp) );
    setElectricConsumption( ecs->getElectricConsumption(double (Model::dt), getThermalPowerProvided(), heatPumpSrcTemp) + pumpPower*double (Model::dt) );
}

void ThermalStation::computeThermalPowerProvided(float const& thermalPowerNeeded, Climate* pClimate, unsigned int day, unsigned int hour) {

    ecs->setThermalPowerNeeded(thermalPowerNeeded);

    // TODO, improve this, make the heatpump be able to access the temperature on its own.
    // the heatPumpSrcTemp depends on the type of heatpump (air or ground vertical/horizontal)
    float heatPumpSrcTemp = pClimate->getToutCelsius(day,hour); // air temperature
    float z0=0.f,z1=0.f,alpha=0.f;
    if (ecs->getGround(z0,z1,alpha)) {
        heatPumpSrcTemp = pClimate->getTgroundCelsius(day,hour,z0,alpha,z1);
    }

    thermalPowerProvided = ecs->getThermalPower(heatPumpSrcTemp);
}

void ThermalStation::writeTHHeaderText(fstream& textFile, string prefix) {
    textFile << prefix <<":MachinePower(W)" <<"\t"
             << prefix <<":PumpPower(W)" <<"\t"
             << prefix <<":FuelConsumption(J)" <<"\t"
             << prefix <<":ElectricConsumption(J)" <<"\t";
}

void ThermalStation::writeTHResultsText(fstream& textFile, unsigned int i) {
    textFile << fixed << setprecision(0) << getMachinePower(i)<<"\t";
    textFile << fixed << setprecision(2) << getPumpPower(i) <<"\t";
    textFile << fixed << setprecision(0) << getFuelConsumption(i)<<"\t";
    textFile << fixed << setprecision(0) << getElectricConsumption(i)<<"\t";
}






SeasonalStorageHeatingThermalStation::SeasonalStorageHeatingThermalStation(TiXmlHandle hdl, Network* net, ostream* pLogStr) : ThermalStation(hdl, net, pLogStr), storage(nullptr), valve(nullptr), storageModeOn(false), tempDiffAroundStorage(1.f), flowToStore(0.001f) {
    try {
        string name = "kvMax";
        float kvMax;
        if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &kvMax) ) { throw string("Error in the XML file: a ThermalStation of type seasonalStorage doesn't have attribute: '"+name+"'."); }
        if ( kvMax<=0.f ) { throw string("Error in the XML file: a ThermalStation of type seasonalStorage has kvMax<=0.f, which isn't physically possible."); }

        if (hdl.FirstChildElement("Storage").ToElement()) {
            storage = Storage::createNewStorage(hdl.FirstChildElement("Storage")); // Dynamical allocation.
        } else {
            throw string("Error in the XML file: a SeasonalStorageThermalStation doesn't contain any Storage child element.");
        }

        valve = new Valve(kvMax);

    } catch (...) {
        SeasonalStorageHeatingThermalStation::deleteDynAllocated();
        throw;
    }
}


void SeasonalStorageHeatingThermalStation::computePressureDiff(float const& massFlow, float const& rho, float& deltaP, float& dDeltaP_dm) {
    if (storageModeOn) { // Use the valve.
        valve->computePressureDiffAndDerivative(-massFlow, rho, deltaP, dDeltaP_dm);
        deltaP *= -1.f; // Sign conventions are oppostite for the valve. dDeltaP_dm has double sign negation, so stays the same.
    } else { // Use the pump.
        ThermalStation::computePressureDiff(massFlow, rho, deltaP, dDeltaP_dm);
    }
}

void SeasonalStorageHeatingThermalStation::updateControlVariable(float const& massFlow, float const& deltaP, float const& rho, float& sumDeltaRpm, float& sumRpm, float& sumDeltaKv, float& sumKv, float const& learningRate) {
    if (storageModeOn) { // Use valve.
        float targetMassFLow, targetPressureDiff;
        if (massFlow<0.f) { // Correct direction (supply to return)
            targetMassFLow = -flowToStore;
            targetPressureDiff = -deltaP;
            if(targetPressureDiff<=0.f) { targetPressureDiff = 1.0f; }
        }
        else { // Wrong direction (return to supply).
            targetMassFLow = flowToStore*0.001f; // Want a small flow, close valve.
            targetPressureDiff = abs(deltaP)+1.f;
        } // If flow is going in the wrong direction, close the valve.

        valve->updateKv(rho, sumDeltaKv, sumKv, learningRate, targetMassFLow, targetPressureDiff);
//        valve->keepKvConst(sumDeltaKv, sumKv); // Keep valve always at const value.
//        valve->setKvToMax(sumDeltaKv, sumKv); // Keep valve always at max.

    } else { // Use pump.
        ThermalStation::updateControlVariable(massFlow, deltaP, rho, sumDeltaRpm, sumRpm, sumDeltaKv, sumKv, learningRate);
    }
}

void SeasonalStorageHeatingThermalStation::computeThermal(float pressureDiff, float massFlow, float rho, float cp, float inputTemp, Climate *pClimate, unsigned int day, unsigned int hour) {
    if (storageModeOn) {

        pumpPower = 0.f;
        float temperatureIncreaseValve = abs(valve->computeTemperatureIncrease(-pressureDiff, cp, rho));
        float tempAfterValve = inputTemp + temperatureIncreaseValve;

        float thermalPowerEcsNeeded = 0.f;
        computeThermalPowerProvided(thermalPowerEcsNeeded, pClimate, day, hour);

        if (massFlow>=0.f) { // Wrong direction.
            outputTemperature = tempAfterValve;
        }
        else if (massFlow<0.f) { // Good direction.
            // This only works for stations that heat (not cool) and that store heat (not cold).
            bool storageHeatsUp = true; // Here we store heat, so we heat up the storage.
            outputTemperature = storage->computeOutputTemperature(storageHeatsUp, tempAfterValve, -massFlow, cp);
        }
        tempDiffAroundStorage = outputTemperature - tempAfterValve;

        // TODO: add protection to fan off heat if there is too much? to keep at least a 5 degree differrence?

    }
    else { // Not storaging, heat the fluid to target temperature.
        float thermalPowerEcsNeeded;

        float thermalPowerPump;
        pump->computeElectricAndThermalPower(pressureDiff, massFlow, rho, pumpPower, thermalPowerPump);
        float tempAfterPump = inputTemp + thermalPowerPump/(cp*abs(massFlow));
        float tempAfterStorage;

        if (massFlow<=0.f) { // Wrong direction;
            tempAfterStorage = tempAfterPump;
            thermalPowerEcsNeeded = 0.f;
        }
        else { // Correct direction.
            float targetTemp = temperatureSetpoint->computeTargetTemperature(pClimate, day, hour);
            if (tempAfterPump>=targetTemp) { // Target temperature already exceeded.
                tempAfterStorage = tempAfterPump;
                thermalPowerEcsNeeded = 0.f;
            } else { // Still need to heat up.
                bool storageHeatsUp = false; // Here we take heat from storage, so we cool down the storage.
                tempAfterStorage = storage->computeOutputTemperature(storageHeatsUp, tempAfterPump, massFlow, cp);
                if (tempAfterStorage>targetTemp) {
                    tempAfterStorage = targetTemp;
                    thermalPowerEcsNeeded = 0.f;
                } else {
                    thermalPowerEcsNeeded = cp*massFlow*(targetTemp - tempAfterStorage);
                }
            }
        }
        tempDiffAroundStorage = tempAfterStorage - tempAfterPump;

        computeThermalPowerProvided(thermalPowerEcsNeeded, pClimate, day, hour);


        if (massFlow!=0.f) {
            outputTemperature = tempAfterStorage + thermalPowerProvided / (cp * abs(massFlow));
        } // Else don't change the temperature.
    }
}

void SeasonalStorageHeatingThermalStation::updateOperationMode(float const& sumSubstationDesiredMassFlow, size_t const& nbThermalStations) {
    // The substation (of which some are prosumers) are overall injecting more mass flow in the supply, so thermal stations cannot also inject, store the extra heat.
    if (sumSubstationDesiredMassFlow<=0.f) {
        storageModeOn = true;
        flowToStore = sumSubstationDesiredMassFlow/float(nbThermalStations);
    } else {
        storageModeOn = false;
        flowToStore = 0.f;
    }
}


void SeasonalStorageHeatingThermalStation::relativeErrorPressureDiff(float const& massFlow, float const& pressureDiff, float& relErr, float& sumErrP, float& sumDesiredP) {
    if (storageModeOn) {
        // The goal is set to the current pressure, so is always achieved.
        float goal = -flowToStore; // Minus signs to have positive mass flow.
        relErr = (-massFlow-goal)/goal;
        sumErrP += abs(-pressureDiff-goal);
        sumDesiredP += goal;
    } else {
        ThermalStation::relativeErrorPressureDiff(massFlow, pressureDiff, relErr, sumErrP, sumDesiredP);
    }
}

void SeasonalStorageHeatingThermalStation::writeTHHeaderText(fstream& textFile, string prefix) {
    ThermalStation::writeTHHeaderText(textFile, prefix);
    storage->writeTHHeaderText(textFile, prefix);
}

void SeasonalStorageHeatingThermalStation::writeTHResultsText(fstream& textFile, unsigned int i) {
    ThermalStation::writeTHResultsText(textFile, i);
    storage->writeTHResultsText(textFile, i);
}









vector<string> DistrictEnergyCenter::muPossibilities = { "0.0004", "Vogel-Fulcher-Tammann" }; // 0.0004 is a constant value in Pa*s, Vogel-Fulcher-Tammann in an empirical formula.

DistrictEnergyCenter::DistrictEnergyCenter(TiXmlHandle hdl, District* pDistrict, ostream* pLogStr) :
    pDistrict(pDistrict), thermalStations(0), pipelineNetwork(nullptr), totalThermalLossRecord(0), logStream(std::cout.rdbuf()) {

    logStream.rdbuf(pLogStr->rdbuf()); // Add the logStream

    if( hdl.ToElement() ){

        if ( hdl.ToElement()->QueryUnsignedAttribute("id", &id) )   { throw string("Error in the XML file: a DistrictEnergyCenter doesn't have attribute: 'id'."); }
        if ( hdl.ToElement()->QueryDoubleAttribute("Cp", &cp) )     { throw string("Error in the XML file: a DistrictEnergyCenter doesn't have attribute: 'Cp' [J/(kg K)]."); }
        if ( hdl.ToElement()->QueryDoubleAttribute("rho", &rho) )   { throw string("Error in the XML file: a DistrictEnergyCenter doesn't have attribute: 'rho' [kg/m^3]."); }
        if ( hdl.ToElement()->QueryStringAttribute("mu", &mu) )   { throw string("Error in the XML file: a DistrictEnergyCenter doesn't have attribute: 'mu'."); }

        if ( cp<=0 )  { throw string("Error in the XML file: DistrictEnergyCenter id="+to_string(id)+" has Cp<=0."); }
        if ( rho<=0 ) { throw string("Error in the XML file: DistrictEnergyCenter id="+to_string(id)+" has rho<=0."); }
        if ( mu!=muPossibilities[0] and mu!=muPossibilities[1] ) { throw string("Error in the XML file: DistrictEnergyCenter id="+to_string(id)+" has mu!="+muPossibilities[0]+" and mu!="+muPossibilities[1]+"."); }

        logStream << "District energy center: id=" << id << ", Cp=" << cp << ", rho=" << rho << endl;


        try {
            if(hdl.FirstChildElement("Network").ToElement() ) { pipelineNetwork = new Network(hdl.FirstChildElement("Network"), this, &(this->logStream)); }
            else { throw string("Error in the XML file: a DistrictEnergyCenter, doesn't have any 'Network' child element."); }

            string name = "ThermalStation";
            TiXmlElement* xmlEl = hdl.ToElement()->FirstChildElement(name);
            if ( not xmlEl ) { throw string("Error in the XML file: a DistrictEnergyCenter doesn't have any '"+name+"' child elements."); }
            else {
                while( xmlEl ) {
                    thermalStations.push_back(ThermalStation::createNewThermalStation(xmlEl, pipelineNetwork, &(this->logStream)));
                    xmlEl = xmlEl->NextSiblingElement(name); // Go to next thermal station.
                }
            }
        } catch (...) {
            deleteDynAllocated();
            throw;
        }
    }
    else {
        // Throw error ?
    }
}

void DistrictEnergyCenter::deleteDynAllocated() {
    while(!thermalStations.empty()) {
        delete thermalStations.back();
        thermalStations.pop_back();
    }
    if (pipelineNetwork!=nullptr) { delete pipelineNetwork; }
}

DistrictEnergyCenter::~DistrictEnergyCenter() {
    logStream << "Destructor of DistrictEnergyCenter" << endl;
    deleteDynAllocated();
}

double DistrictEnergyCenter::getMu(double temp) {
    float mu_ = 0.0004; // JK - added default value to 0.0004
    if (mu==muPossibilities[0]){ mu_ = 0.0004; }
    else if (mu==muPossibilities[1]) { mu_ = 0.00002939*exp( 507.88/(temp-(149.3-273.15)) ); } // Check for temperature range where this formula is valid?
    return mu_;
}

void DistrictEnergyCenter::recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour) {
    totalThermalLossRecord.push_back(pipelineNetwork->computeTotalThermalLoss());
    pipelineNetwork->recordTimeStep(pClim, day, hour, cp);
}

void DistrictEnergyCenter::eraseRecords(unsigned int keepValue) {
    totalThermalLossRecord.erase(totalThermalLossRecord.begin(),totalThermalLossRecord.end()-min(keepValue,(unsigned int)totalThermalLossRecord.size()));
    pipelineNetwork->eraseRecords(keepValue);
}

void DistrictEnergyCenter::eraseRecords_back() {
    totalThermalLossRecord.pop_back();
    pipelineNetwork->eraseRecords_back();
}

void DistrictEnergyCenter::convergeToEquilibrium(Climate *pClim, unsigned int day, unsigned int hour) {
    pipelineNetwork->convergeToEquilibrium(rho, cp, pClim, day, hour);
    recordTimeStep(pClim, day, hour);
}


bool DistrictEnergyCenter::addSubstationIfContainNode(Substation* pSub, unsigned int nodeId) {
    SubstationNodePair* node = pipelineNetwork->pointerOfSubstationNodeWithId(nodeId);
    if ( node==nullptr ) {
        return false;
    }
    else {
        pSub->setDEC(this);
        pSub->setNode(node);

        node->setSubstation(pSub);
//        pipelineNetwork->initializeMassFlows(node);
        return true;
    }
}







void Network::deleteNodesAndPipes() {
    while(!pipePairs.empty()) {
        delete pipePairs.back();
        pipePairs.pop_back();
    }
    while(!separateNodePairs.empty()) {
        delete separateNodePairs.back();
        separateNodePairs.pop_back();
    }
    while(!thermalStationNodePairs.empty()) {
        delete thermalStationNodePairs.back();
        thermalStationNodePairs.pop_back();
    }
    while(!substationNodePairs.empty()) {
        delete substationNodePairs.back();
        substationNodePairs.pop_back();
    }
}

Network::Network(TiXmlHandle hdl, DistrictEnergyCenter *pDEC, ostream *pLogStr) :
    pDEC(pDEC), loopMassFlows(0), loopMatrix(0), regulatedPathMatrix(0), jacobianRegulatedPathColumns(0), logStream(std::cout.rdbuf()) {

    logStream.rdbuf(pLogStr->rdbuf());
    logStream << "Constructing Network" << endl;

    // Loading the network properties.
    if(hdl.ToElement() ) {

        if ( hdl.ToElement()->QueryFloatAttribute("soilkValue", &soilkValue) )     { throw string("Error in the XML file: a Network doesn't have attribute: 'soilkValue'."); }

        if (soilkValue<=0) { throw string("Error in the XML file: Network of DistrictEnergyCenter id="+to_string(pDEC->getId())+" has soilkValue<=0."); }

        try {
            string name = "NodePair";
            TiXmlElement* xmlEl = hdl.ToElement()->FirstChildElement(name);
            while( xmlEl ) {
                separateNodePairs.push_back( new NodePair( xmlEl, pDEC->getInitTemp() ) );
                xmlEl = xmlEl->NextSiblingElement(name); // Go to next network node.
            }

            name = "SubstationNodePair";
            xmlEl = hdl.ToElement()->FirstChildElement(name);
            if ( not xmlEl ) { throw string("Error in the XML file: a Network doesn't have any '"+name+"' child elements."); }
            else {
                while( xmlEl ) {
                    substationNodePairs.push_back( new SubstationNodePair( xmlEl, pDEC->getInitTemp() ) );
                    xmlEl = xmlEl->NextSiblingElement(name); // Go to next network node.
                }
            }

            name = "ThermalStationNodePair";
            xmlEl = hdl.ToElement()->FirstChildElement(name);
            if ( not xmlEl ) { throw string("Error in the XML file: a Network doesn't have any '"+name+"' child elements."); }
            else {
                while( xmlEl ) {
                    thermalStationNodePairs.push_back( new ThermalStationNodePair( xmlEl, pDEC->getInitTemp() ) );
                    xmlEl = xmlEl->NextSiblingElement(name); // Go to next network node.
                }
            }

            name = "PipePair";
            xmlEl = hdl.ToElement()->FirstChildElement(name);
            if ( not xmlEl ) { throw string("Error in the XML file: a Network doesn't have any '"+name+"' child elements."); }
            else {
                while( xmlEl ) {
                    pipePairs.push_back( new PipePair( xmlEl, this ) );
                    xmlEl = xmlEl->NextSiblingElement(name); // Go to next pipe pair.
                }
            }

            for (auto node : separateNodePairs){ node->checkNbPipes(); }
            for (auto node : substationNodePairs){ node->checkNbPipes(); }
            for (auto node : thermalStationNodePairs){ node->checkNbPipes(); }

        } catch(...) { // To delete the unique ids created.
            deleteNodesAndPipes();
            throw;
        }
    }
    else {
        // should it throw error ?
    }
}


Network::~Network() {
    logStream << "Destructor of Network" << endl;
    deleteNodesAndPipes();
}


void Network::writeTHHeaderText(fstream& textFile) {
    unsigned int decId = pDEC->getId();

    for (auto node : thermalStationNodePairs) {
        node->writeTHHeaderText(textFile, decId);
    }

    for (auto node : substationNodePairs) {
        node->writeTHHeaderText(textFile, decId);
    }

    for (auto node : separateNodePairs) {
        node->writeTHHeaderText(textFile, decId);
    }

    for (auto pair : pipePairs) {
        unsigned int pairId = pair->getId();
        textFile << "DEC" << decId <<":PipePair" << pairId << ":SupplyMassFlow(kg/s)" <<"\t"
                 << "DEC" << decId <<":PipePair" << pairId << ":ReturnMassFlow(kg/s)" <<"\t"
                 << "DEC" << decId <<":PipePair" << pairId << ":SupplyPressureDiff(Pa)" <<"\t"
                 << "DEC" << decId <<":PipePair" << pairId << ":ReturnPressureDiff(Pa)" <<"\t"
                 << "DEC" << decId <<":PipePair" << pairId << ":SupplyThermalLoss(W)" <<"\t"
                 << "DEC" << decId <<":PipePair" << pairId << ":ReturnThermalLoss(W)" <<"\t";
    }
}

void Network::writeTHResultsText(fstream& textFile, unsigned int i) {

    for (auto node : thermalStationNodePairs) {
        node->writeTHResultsText(textFile, i);
    }

    for (auto node : substationNodePairs) {
        node->writeTHResultsText(textFile, i);
    }

    for (auto node : separateNodePairs) {
        node->writeTHResultsText(textFile, i);
    }

    for (auto pair : pipePairs) {
        textFile << fixed << setprecision(5) << pair->getSupplyPipe()->getMassFlowRecord(i) <<"\t";
        textFile << fixed << setprecision(5) << pair->getReturnPipe()->getMassFlowRecord(i) <<"\t";
        textFile << fixed << setprecision(2) << pair->getSupplyPipe()->getPressureDiffRecord(i) <<"\t";
        textFile << fixed << setprecision(2) << pair->getReturnPipe()->getPressureDiffRecord(i) <<"\t";
        textFile << fixed << setprecision(1) << pair->getSupplyPipe()->getThermalLossRecord(i) <<"\t";
        textFile << fixed << setprecision(1) << pair->getReturnPipe()->getThermalLossRecord(i) <<"\t";
    }

}

int Network::nbFloatsRecorded() { // TODO: Improve this similarly to writeTHHeaderText, so that the code adapts better to future changes.
    return thermalStationNodePairs.size()*8+substationNodePairs.size()*4+separateNodePairs.size()*2+pipePairs.size()*6;
}

void Network::initializeMassFlows() {
    for (size_t i=0; i<substationNodePairs.size(); i++) {
        SubstationNodePair* ssnp = substationNodePairs[i];
        float massFlow = ssnp->getDesignMassFlow();
        loopMassFlows[ loopMassFlows.size()-substationNodePairs.size()+i ] = massFlow; // In the loop passing through the substation, put the design mass flow.

        massFlow /= thermalStationNodePairs.size(); // Distribute equally the flow between all thermal stations.
        for (size_t j=1; j<thermalStationNodePairs.size(); j++) {
            loopMassFlows[ loopMassFlows.size()-substationNodePairs.size()-thermalStationNodePairs.size()+j ] -= massFlow; // In the loops from the zeroth thermal station to other thermal stations, add this flow so that all thermal stations contribute equally to produce this flow.
        }
    }
}

void Network::propagateNetwork(float cp, float soilTemp, bool isSupply) {
    // To keep information of which components have already been traversed.
    for (auto pipe : pipePairs) { pipe->setAlreadyTraversed(false); }
    for (auto tsnp : thermalStationNodePairs) { tsnp->setAlreadyTraversed(false); }
    for (auto ssnp : substationNodePairs) { ssnp->setAlreadyTraversed(false); }

    // Propagates from high pressures to low pressures.
    for (auto tsnp : thermalStationNodePairs) { tsnp->propagateNetwork(soilTemp, cp, isSupply); }
    for (auto ssnp : substationNodePairs) { ssnp->propagateNetwork(soilTemp, cp, isSupply); }
}

void Network::recordTimeStep(Climate* pClim, unsigned int day, unsigned int hour, float const& cp) {
    for(auto pipePair : pipePairs) { pipePair->recordTimeStep(); } // recordThermalLosses  recordMassFlows  recordPressureDiffs
    for(auto node : separateNodePairs) { node->recordTimeStep(pClim, day, hour, cp); } // recordTemperatures
    for(auto node : thermalStationNodePairs) { node->recordTimeStep(pClim, day, hour, cp); } // recordTemperatures  recordMassFlow  recordPressureDiff  setMachinePower_FuelConsumption_ElectricConsumption setPumpPowers  confirmStorage
    for(auto node : substationNodePairs) { node->recordTimeStep(pClim, day, hour, cp); } // recordTemperatures  recordMassFlow  recordPressureDiff  setProsumerSolarThermal
}

void Network::eraseRecords(unsigned int keepValue) {
    for(auto pipePair : pipePairs) { pipePair->eraseRecords(keepValue); }
    for(auto node : separateNodePairs) { node->eraseRecords(keepValue); }
    for(auto node : thermalStationNodePairs) { node->eraseRecords(keepValue); }
    for(auto node : substationNodePairs) { node->eraseRecords(keepValue); }
}

void Network::eraseRecords_back() {
    for(auto pipePair : pipePairs) { pipePair->eraseRecords_back(); }
    for(auto node : separateNodePairs) { node->eraseRecords_back(); }
    for(auto node : thermalStationNodePairs) { node->eraseRecords_back(); }
    for(auto node : substationNodePairs) { node->eraseRecords_back(); }
}

float Network::computeTotalThermalLoss() {
    float total = 0.f;
    for ( auto pipePair : pipePairs ) {
        total += pipePair->getSupplyPipe()->getThermalLoss() + pipePair->getReturnPipe()->getThermalLoss();
    }
    return total;
}

void Network::checkThermalStationsFound() {
    for (auto node : thermalStationNodePairs) {
        if (node->thermalStationNotSet()) { throw string("Error in the XML file: ThermalStationNodePair of id="+to_string(node->getId())+" has no ThermalStation pointing to it in the same DistrictEnergyCenter."); }
    }
}

void Network::substationsConstructionFinished() {

    for (auto node : substationNodePairs) {
        if (node->substationNotSet()) { throw string("Error in the XML file: SubstationNodePair of id="+to_string(node->getId())+" has no Substation pointing to it."); }
    }

    map<PipePair*,size_t> edgeIdxPipe;
    map<NodePairConnectingSupplyReturn*,size_t> edgeIdxNode;
    size_t nbEdges;
    computeEdgeIdx(edgeIdxPipe, edgeIdxNode, nbEdges);
    findLoops(loopMatrix, edgeIdxPipe, edgeIdxNode, nbEdges);
    findRegulatedPaths(regulatedPathMatrix, edgeIdxNode, nbEdges);
    computeJacobianRegulatedPathColumns(jacobianRegulatedPathColumns, loopMatrix, regulatedPathMatrix,  edgeIdxNode, nbEdges);
    for (size_t i=0; i<nbLoops(); i++) { loopMassFlows.push_back(randomUniform(0.001,0.002)); }

    if (regulatedPathMatrix.size()!=regulatingElements.size()) { throw string("Error in Network::checkSubstationsFound(), regulatedPathMatrix.size()!=regulatingElements.size()."); }

    initializeMassFlows();
}

NodePair* Network::pointerOfNodeWithId(unsigned int nodeId) {
    for(auto node : separateNodePairs) {
        if (nodeId==node->getId()) {
            return node;
        }
    }
    for(auto node : thermalStationNodePairs) {
        if (nodeId==node->getId()) {
            return node;
        }
    }
    for(auto node : substationNodePairs) {
        if (nodeId==node->getId()) {
            return node;
        }
    }
    return nullptr;
}

SubstationNodePair* Network::pointerOfSubstationNodeWithId(unsigned int nodeId) {
    for(auto node : substationNodePairs) {
        if (nodeId==node->getId()) {
            return node;
        }
    }
    return nullptr;
}

ThermalStationNodePair* Network::pointerOfThermalStationNodeWithId(unsigned int nodeId) {
    for(auto node : thermalStationNodePairs) {
        if (nodeId==node->getId()) {
            return node;
        }
    }
    return nullptr;
}

ostream& Network::print(ostream& stream) {
    stream << "Network, DECid = " << getDEC()->getId() << endl;
    return stream;
}

//size_t Network::computeEdgeIdx(PipePair* pp, bool isSupply) {
//    size_t idx = 0;
//    while (pipePairs[idx]!=pp and idx<pipePairs.size()) { idx++; }
//    if (idx==pipePairs.size()) { throw string("Error in Network::computeEdgeIdx(PipePair*,bool), the PipePair given is not in pipePairs."); }
//    idx *= 2;
//    if (not isSupply) { idx++; }
//    return idx;
//}

//size_t Network::computeEdgeIdx(ThermalStationNodePair* np) {
//    size_t idx = 0;
//    while (thermalStationNodePairs[idx]!=np and idx<thermalStationNodePairs.size()) { idx++; }
//    if (idx==thermalStationNodePairs.size()) { throw string("Error in Network::computeEdgeIdx(ThermalStationNodePair*), the node given is not in thermalStationNodePairs."); }
//    idx += pipePairs.size()*2;
//    return idx;
//}

//size_t Network::computeEdgeIdx(SubstationNodePair* np) {
//    size_t idx = 0;
//    while (substationNodePairs[idx]!=np and idx<substationNodePairs.size()) { idx++; }
//    if (idx==substationNodePairs.size()) { throw string("Error in Network::computeEdgeIdx(SubstationNodePair*), the node given is not in substationNodePairs."); }
//    idx += pipePairs.size()*2 + thermalStationNodePairs.size();
//    return idx;
//}

void Network::computeEdgeIdx(map<PipePair*,size_t>& edgeIdxPipe, map<NodePairConnectingSupplyReturn*,size_t>& edgeIdxNode, size_t& nbEdges) {
    nbEdges = 0;
    for(auto const& pp : pipePairs) {
        edgeIdxPipe[pp] = nbEdges; // Idx for the supply pipe.
        nbEdges += 2; // Skip over the idx for the return pipe.
    }
    for(auto const& tsnp : thermalStationNodePairs) {
        edgeIdxNode[tsnp] = nbEdges;
        nbEdges++;
    }
    for(auto const& ssnp : substationNodePairs) {
        edgeIdxNode[ssnp] = nbEdges; // Idx corresponds to the first edge.
        nbEdges += ssnp->nbEdges(); // Skip over the idxs of other edges.
    }
}


void Network::findLoops(vector<vector<float>>& matB, map<PipePair*,size_t> const& edgeIdxPipe, map<NodePairConnectingSupplyReturn*,size_t> const& edgeIdxNode, size_t const& nbEdges) {

    matB = vector<vector<float>> (0);

    map<NodePair*,PipePair*> previousPipe;

    map<NodePair*,bool> nodePairVisited;
    for (auto nodePair : separateNodePairs) {       nodePairVisited[nodePair] = false;  previousPipe[nodePair] = nullptr; }
    for (auto nodePair : thermalStationNodePairs) { nodePairVisited[nodePair] = false;  previousPipe[nodePair] = nullptr; }
    for (auto nodePair : substationNodePairs) {     nodePairVisited[nodePair] = false;  previousPipe[nodePair] = nullptr; }
    map<PipePair*,bool> pipePairVisited;
    for (auto pipePair : pipePairs) { pipePairVisited[pipePair] = false; }

    ThermalStationNodePair* startNode = thermalStationNodePairs[0];
    queue<NodePair*> q;
    q.push(startNode);
    nodePairVisited[startNode] = true;


    while (not q.empty()) { // Breadth first search to get the spanning tree.
        NodePair* currentNode = q.front();
        q.pop();
        for (auto pipePair : currentNode->getPipes()) { // Loop on pipes connected to the current node.
            NodePair* otherNode = pipePair->getOtherNode(currentNode);
            if (not nodePairVisited[otherNode]) {
                pipePairVisited[pipePair] = true;
                nodePairVisited[otherNode] = true;
                previousPipe[otherNode] = pipePair;
                q.push(otherNode);
            }
        }
    }

    for (auto el : nodePairVisited) { // Check that all nodes were reached.
        bool nodeVisited = el.second;
        if (not nodeVisited) { throw string("Error in the XML file: the NodePairs id="+to_string(startNode->getId())+" and id="+to_string(el.first->getId())+" cannot be connected via the pipe network."); }
    }

    // Add loops whose pipes are only present in supply network (and the corresponding loop in return networks).
    for (auto el : pipePairVisited) { // Check for pipePairs that were not visited, meaning thay are not in the spanning tree, each one corresponds to a loop.
        bool pipeVisited = el.second;
        if (not pipeVisited) {
            PipePair* pipe;
            NodePair* node;
            stack<PipePair*> predecessorEdgesTailSide, predecessorEdgesHeadSide;
            stack<NodePair*> predecessorNodesTailSide, predecessorNodesHeadSide;

            node = el.first->getTailNode();
            pipe = previousPipe[node];
            while (pipe!=nullptr) { // Find path from the unvisited edge, all the way to the starting node, on the tail side.
                predecessorEdgesTailSide.push(pipe); // Push in stack.
                node = pipe->getOtherNode(node); // Get previous node.
                predecessorNodesTailSide.push(node); // Push in stack.
                pipe = previousPipe[node]; // Get previous pipe.
            }

            node = el.first->getHeadNode();
            pipe = previousPipe[node];
            while (pipe!=nullptr) { // Find path from the unvisited edge, all the way to the starting node, on the head side.
                predecessorEdgesHeadSide.push(pipe); // Push in stack.
                node = pipe->getOtherNode(node); // Get previous node.
                predecessorNodesHeadSide.push(node); // Push in stack.
                pipe = previousPipe[node]; // Get previous pipe.
            }

            while ( (not predecessorEdgesTailSide.empty())  and  (not predecessorEdgesHeadSide.empty())  and  predecessorEdgesTailSide.top()==predecessorEdgesHeadSide.top() ) { // Pop the common elements, that are not part of the loop.
                predecessorEdgesTailSide.pop();
                predecessorNodesTailSide.pop();
                predecessorEdgesHeadSide.pop();
                predecessorNodesHeadSide.pop();
            }


            vector<float> rowToAddSupply (nbEdges, 0.f);
            vector<float> rowToAddReturn (nbEdges, 0.f);

            float direction = 1.f; // Add the unvisited edge.
//            size_t idx = computeEdgeIdx(el.first, true);
//            rowToAddSupply[idx] = direction;
//            idx = computeEdgeIdx(el.first, false);
//            rowToAddReturn[idx] = direction;
            size_t idx = edgeIdxPipe.at(el.first);
            rowToAddSupply[idx] = direction;
            rowToAddReturn[idx+1] = direction;

            while (not predecessorEdgesTailSide.empty()) { // Add tail side of the loop.
                pipe = predecessorEdgesTailSide.top(); predecessorEdgesTailSide.pop();
                node = predecessorNodesTailSide.top(); predecessorNodesTailSide.pop();
                if (pipe->getTailNode()==node) { direction = 1.f; }
                else { direction = -1.f; }
//                idx = computeEdgeIdx(pipe, true);
//                rowToAddSupply[idx] = direction;
//                idx = computeEdgeIdx(pipe, false);
//                rowToAddReturn[idx] = direction;
                idx = edgeIdxPipe.at(pipe);
                rowToAddSupply[idx] = direction;
                rowToAddReturn[idx+1] = direction;
            }

            while (not predecessorEdgesHeadSide.empty()) { // Add head side of the loop.
                pipe = predecessorEdgesHeadSide.top(); predecessorEdgesHeadSide.pop();
                node = predecessorNodesHeadSide.top(); predecessorNodesHeadSide.pop();
                float direction;
                if (pipe->getHeadNode()==node) { direction = 1.f; }
                else { direction = -1.f; }
//                idx = computeEdgeIdx(pipe, true);
//                rowToAddSupply[idx] = direction;
//                idx = computeEdgeIdx(pipe, false);
//                rowToAddReturn[idx] = direction;
                idx = edgeIdxPipe.at(pipe);
                rowToAddSupply[idx] = direction;
                rowToAddReturn[idx+1] = direction;
            }
            matB.push_back(rowToAddSupply);
            matB.push_back(rowToAddReturn);
        }
    }


    // Add loops that happen at thermal stations other than the startNode.
    for (size_t i=1 ; i<thermalStationNodePairs.size() ; i++) {
        ThermalStationNodePair* tsnode = thermalStationNodePairs[i];
        vector<float> rowToAdd(nbEdges, 0.f);

        float direction = -1.f; // Add the other thermal station edge (direction is from supply to return).
//        size_t idx = computeEdgeIdx(tsnode);
        size_t idx = edgeIdxNode.at(tsnode);
        rowToAdd[idx] = direction;

        direction = 1.f;
//        idx = computeEdgeIdx(startNode);
        idx = edgeIdxNode.at(startNode);
        rowToAdd[idx] = direction; // Add the thermal station edge (direction is from return to supply).


        NodePair* node = tsnode;
        PipePair* pipe = previousPipe[node];
        while (pipe!=nullptr) {

            if (pipe->getHeadNode()==node) { direction = 1.f; } // Direction for the supply.
            else { direction = -1.f; }
//            idx = computeEdgeIdx(pipe, true);
//            rowToAdd[idx] = direction;
//            idx = computeEdgeIdx(pipe, false);
//            rowToAdd[idx] = -1.f*direction; // Return has opposite direction.
            idx = edgeIdxPipe.at(pipe);
            rowToAdd[idx] = direction; // Supply.
            rowToAdd[idx+1] = -direction; // Return has opposite direction.

            node = pipe->getOtherNode(node); // Get previous node.
            pipe = previousPipe[node]; // Get previous pipe.
        }
        matB.push_back(rowToAdd);
    }


    // Add loops that happen at every substation and the other thermal stations.
    for (auto subsnode : substationNodePairs) {
        vector<float> rowToAdd(nbEdges, 0.f);

        float direction = 1.f; // Add the substation edge (direction is from supply to return).
//        size_t idx = computeEdgeIdx(subsnode);
//        rowToAdd[idx] = direction;
        size_t idx = edgeIdxNode.at(subsnode);
        for(size_t i=0; int(i)<subsnode->nbEdges(); i++) { rowToAdd[idx+i] = direction; } // A substation may contain multiple edges (valve, differential pressure regulator, etc).

//        idx = computeEdgeIdx(startNode);
        idx = edgeIdxNode.at(startNode);
        rowToAdd[idx] = direction; // Add the thermal station edge (direction is from return to supply).


        NodePair* node = subsnode;
        PipePair* pipe = previousPipe[node];
        while (pipe!=nullptr) {

            if (pipe->getHeadNode()==node) { direction = 1.f; } // Direction for the supply.
            else { direction = -1.f; }
//            idx = computeEdgeIdx(pipe, true);
//            rowToAdd[idx] = direction;
//            idx = computeEdgeIdx(pipe, false);
//            rowToAdd[idx] = -1.f*direction; // Return has opposite direction.
            idx = edgeIdxPipe.at(pipe);
            rowToAdd[idx] = direction; // Supply.
            rowToAdd[idx+1] = -direction; // Return has opposite direction.

            node = pipe->getOtherNode(node); // Get previous node.
            pipe = previousPipe[node]; // Get previous pipe.
        }

        matB.push_back(rowToAdd);
    }

//    logStream << "matB" << endl;
//    for (auto const& row : matB) {
//        for (auto const& el : row) {
//            logStream << el << " ";
//        }
//        logStream << endl;
//    }

}

void Network::findRegulatedPaths(vector<vector<float>>& matR, map<NodePairConnectingSupplyReturn*,size_t> const& edgeIdxNode, size_t const& nbEdges) {
    matR = vector<vector<float>> (0);
    for (auto subsnode : substationNodePairs) {

        if (subsnode->hasRegEle()) {
            vector<float> rowToAdd(nbEdges, 0.f);
            size_t idx = edgeIdxNode.at(subsnode);
            float direction = 1.f; // Add the substation edge(s) (direction is from supply to return).
            for(size_t i=0; (int)i+1<subsnode->nbEdges(); i++) { rowToAdd[idx+i] = direction; } // The regulated path contains all edges except the last, which is the regulating element.
            matR.push_back(rowToAdd);
        }
    }
    logStream << "matR" << endl;
    for (auto const& row : matR) {
        for (auto const& el : row) {
            logStream << el << " ";
        }
        logStream << endl;
    }
}

void Network::computeJacobianRegulatedPathColumns(vector<vector<float>>& matJRPC, vector<vector<float>> const& matB, vector<vector<float>> const& matR,  map<NodePairConnectingSupplyReturn*,size_t> const& edgeIdxNode, size_t const& nbEdges) {
    vector<vector<float>> matRegEleTransposed (0); // nb rows = nb reg paths, nb cols = nb edges.
    for (auto subsnode : substationNodePairs) {
        if (subsnode->hasRegEle()) {
            vector<float> rowToAdd(nbEdges, 0.f);
            size_t idx = edgeIdxNode.at(subsnode);
            rowToAdd[idx+subsnode->nbEdges()-1] = 1.f; // The regulating element is the last edge.
            matRegEleTransposed.push_back(rowToAdd);
        }
    }

    matJRPC = vector<vector<float>> (0);
    if (matRegEleTransposed.size()>0) { // This condition is added so that if there are zero regulating elements, then matJRPC=vector<vector<float>>(0) not vector<vector<float>>(n,vector<vector<float>>(0)).
        for (auto const& rowB : matB) {
            vector<float> rowToAdd(matR.size(), 0.f);
            for (size_t j=0; j<matRegEleTransposed.size(); j++) {
                vector<float> colRegEle = matRegEleTransposed[j];
                for (size_t k=0; k<rowB.size(); k++) {
                    rowToAdd[j] += rowB[k]*colRegEle[k];
                }
            }
            matJRPC.push_back(rowToAdd);
        }

        for (auto const& rowR : matR) { // Under current implementation of regulating paths, this matrix should be full of zeros.
            vector<float> rowToAdd(matR.size(), 0.f);
            for (size_t j=0; j<matRegEleTransposed.size(); j++) {
                vector<float> colRegEle = matRegEleTransposed[j];
                for (size_t k=0; k<rowR.size(); k++) {
                    rowToAdd[j] += rowR[k]*colRegEle[k];
                }
            }
            matJRPC.push_back(rowToAdd);
        }
    }
    logStream << "matJRPC" << endl;
    for (auto const& row : matJRPC) {
        for (auto const& el : row) {
            logStream << el << " ";
        }
        logStream << endl;
    }
}


void Network::matTransposedVecMult(vector<vector<float>> const& matTransposed, vector<float> const& vec, vector<float>& result) {
    for (size_t i=0 ; i<result.size() ; i++) {
        float sum = 0.f;
        for (size_t j=0 ; j<vec.size() ; j++) {
            sum += matTransposed[j][i] * vec[j];
        }
        result[i] = sum;
    }
}

void Network::matVecMult(vector<vector<float>> const& mat, vector<float> const& vec, vector<double>& result, size_t const& iResStart) {
    for (int i=0 ; i<int(mat.size()) ; i++) {
        double sum = 0.;
        for (int j=0 ; j<int(vec.size()) ; j++) {
            sum += mat[i][j] * vec[j];
        }
        result[iResStart+i] = sum;
    }
}

void Network::matDiagMatTransposedMult(vector<vector<float>> const& mat, vector<float> const& diag, vector<vector<float>> const& matT, double *resultMat, int resultMat_n, size_t const& iResStart, size_t const& jResStart) {
    for (int i=0; i<(int)mat.size() ; i++) { // Loop on rows.

        vector<float> row = mat[i]; // Save the row.
        for (size_t k=0 ; k<row.size() ; k++) { // Multiply by the diagonal matrix only once.
            row[k] *= diag[k];
        }

        for (int j=0 ; j<(int)matT.size() ; j++) { // Loop on columns.

            double sum = 0.;
            for (size_t k=0 ; k<row.size() ; k++) {
                sum += row[k] * matT[j][k];
            }
            resultMat[(iResStart+i)+resultMat_n*(jResStart+j)] = sum;
        }
    }
}

float Network::computeDynamicViscosityWater(float const& temp) {
    return 0.00002939f*exp(507.88f/(temp+123.85f)); // Semi-empirical Vogel-Fulcher-Tammann relation, giving dynamic viscosity in Pa*s.
}

float Network::computeDarcyFrictionFactor(float const& massFlow, float const& radius, float const& temp) {
//    float mu = 0.0004; // Dynamic viscocity of water [Pa*s] TODO improve this formula.
    float mu = pDEC->getMu(temp);
    float reynolds (0.63662*massFlow/(radius*mu)); // 2*massflow/(pi*radius*mu) [ ]
    float eps_diam = 0.000045f/(2.f*radius); // Relative roughness = absolute roughness/diameter [ ]

    float darcyFricFact;

    if (reynolds<=2320.f) { // Laminar flow.
        darcyFricFact = 64.f/reynolds;
    }
    else if (reynolds>=4000.f) { // Turbulent flow.
        darcyFricFact = pow( 1.8*log10( pow(eps_diam*0.27027f, 1.11)+6.9f/reynolds ), -2); // Haaland equation.
    }
    else { // Transition region.
        float darcy4000 = pow( 1.8*log10( pow(eps_diam*0.27027f, 1.11)+0.0017250), -2); // Haaland equation with reynolds=4000.
        float darcy2320 = 0.027586; // Is 64/2320.
        darcyFricFact = darcy2320 + (darcy4000-darcy2320)*(reynolds-2320.f)*0.00059524f; // Linear interpolation, since 0.00059524=1/(4000-2320)

    }
    return darcyFricFact;
}

void Network::computePressureDiff(vector<float> const& massFlows, float const& rho, vector<float>& deltaP, vector<float>& dDeltaP_dm) {
    float cst (19.739f*rho); // 2*pi*pi*rho
    for (size_t i=0 ; i<pipePairs.size() ; i++) {
        float r ( pipePairs[i]->getInnerRadius() );
        float tmp ( pipePairs[i]->getLength()/(cst*r*r*r*r*r) );
        float altitudePressureLoss = pipePairs[i]->computeAltitudePressureLoss(rho);

        size_t idx = 2*i;
        float massFlow ( massFlows[idx] );
        dDeltaP_dm[idx] = abs(massFlow)*computeDarcyFrictionFactor(abs(massFlow), r, pipePairs[i]->getSupplyPipe()->getDownstreamTemp())*tmp; // Usually the fluid stays longer close to the downstream temperature then the upstream temperature. This could be improve.
        deltaP[idx] = 0.5f*massFlow*dDeltaP_dm[idx] + altitudePressureLoss;

        idx = 2*i+1;
        massFlow = massFlows[idx];
        dDeltaP_dm[idx] = abs(massFlow)*computeDarcyFrictionFactor(abs(massFlow), r, pipePairs[i]->getReturnPipe()->getDownstreamTemp())*tmp; // Usually the fluid stays longer close to the downstream temperature then the upstream temperature. This could be improve.
        deltaP[idx] = 0.5f*massFlow*dDeltaP_dm[idx] + altitudePressureLoss;
    }

    for (size_t i=0 ; i<thermalStationNodePairs.size() ; i++) {
        size_t idx = pipePairs.size()*2+i;
        thermalStationNodePairs[i]->getThermalStation()->computePressureDiff(massFlows[idx], rho, deltaP[idx], dDeltaP_dm[idx]);
    }

    int iBeg = pipePairs.size()*2+thermalStationNodePairs.size();
    for (size_t i=0 ; i<substationNodePairs.size() ; i++) {
        substationNodePairs[i]->computePressureDiffAndDerivative(massFlows.begin()+iBeg, rho, deltaP.begin()+iBeg, dDeltaP_dm.begin()+iBeg);
        iBeg += substationNodePairs[i]->nbEdges();
    }
}

void Network::convergeThermal(float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour) {

    float soilTemp = pClim->getTgroundCelsius(day, hour, 0.f); // TODO : Improve this, with a better formula, as this may be wrong in some climate files.
//    float soilTemp = 8.f; // debug

    size_t nNodes = separateNodePairs.size()+thermalStationNodePairs.size()+substationNodePairs.size();
    vector<float> nodeSupplyTemperature_prev (nNodes);
    vector<float> nodeReturnTemperature_prev (nNodes);
    size_t nLoops = 0;
    bool hasConverged;

    do {
        // Remember values of previous step (to check for convergence).
        for (size_t i=0 ; i<separateNodePairs.size() ; i++) {
            nodeSupplyTemperature_prev[i] = separateNodePairs[i]->getSupplyTemperature();
            nodeReturnTemperature_prev[i] = separateNodePairs[i]->getReturnTemperature();
        }
        for (size_t i=0 ; i<thermalStationNodePairs.size() ; i++) {
            nodeSupplyTemperature_prev[separateNodePairs.size()+i] = thermalStationNodePairs[i]->getSupplyTemperature();
            nodeReturnTemperature_prev[separateNodePairs.size()+i] = thermalStationNodePairs[i]->getReturnTemperature();
        }
        for (size_t i=0 ; i<substationNodePairs.size() ; i++) {
            nodeSupplyTemperature_prev[separateNodePairs.size()+thermalStationNodePairs.size()+i] = substationNodePairs[i]->getSupplyTemperature();
            nodeReturnTemperature_prev[separateNodePairs.size()+thermalStationNodePairs.size()+i] = substationNodePairs[i]->getReturnTemperature();
        }

        // From supply to return.
        for(auto ssnp : substationNodePairs) { ssnp->computeThermal(cp, rho, pClim, day, hour, true); }
        for(auto tsnp : thermalStationNodePairs) { tsnp->computeThermal(rho, cp, pClim, day, hour, true); }
        // Through return.
        propagateNetwork(cp, soilTemp, false);
        // From return to supply.
        for(auto ssnp : substationNodePairs) { ssnp->computeThermal(cp, rho, pClim, day, hour, false); }
        for(auto tsnp : thermalStationNodePairs) { tsnp->computeThermal(rho, cp, pClim, day, hour, false); }
        // Through supply.
        propagateNetwork(cp, soilTemp, true);

        nLoops++;
        // Check for convergence (if values didn't move too much since previous step)
        hasConverged = true;
        for ( size_t i=0 ; hasConverged and (i<separateNodePairs.size()) ; i++ ) {
            hasConverged = hasConverged and ( abs(nodeSupplyTemperature_prev[i] - separateNodePairs[i]->getSupplyTemperature())<0.2f );
            hasConverged = hasConverged and ( abs(nodeReturnTemperature_prev[i] - separateNodePairs[i]->getReturnTemperature())<0.2f );
        }
        for ( size_t i=0 ; hasConverged and (i<thermalStationNodePairs.size()) ; i++ ) {
            hasConverged = hasConverged and ( abs(nodeSupplyTemperature_prev[separateNodePairs.size()+i] - thermalStationNodePairs[i]->getSupplyTemperature())<0.2f );
            hasConverged = hasConverged and ( abs(nodeReturnTemperature_prev[separateNodePairs.size()+i] - thermalStationNodePairs[i]->getReturnTemperature())<0.2f );
        }
        for ( size_t i=0 ; hasConverged and (i<substationNodePairs.size()) ; i++ ) {
            hasConverged = hasConverged and ( abs(nodeSupplyTemperature_prev[separateNodePairs.size()+thermalStationNodePairs.size()+i] - substationNodePairs[i]->getSupplyTemperature())<0.2f );
            hasConverged = hasConverged and ( abs(nodeReturnTemperature_prev[separateNodePairs.size()+thermalStationNodePairs.size()+i] - substationNodePairs[i]->getReturnTemperature())<0.2f );
        }

    } while ( nLoops<70 and (not hasConverged) );
}

void Network::computeResiduals(vector<float>& edgeMassFlows, vector<double>& residuals, vector<float>& deltaP, vector<float>& dDeltaP_dm, float const& rho) {
    // Compute edgeMassFlows from loopMassFlows.
    matTransposedVecMult(loopMatrix, loopMassFlows, edgeMassFlows);

    // Compute pressure difference and the derivative from the edgeMassFlows.
    computePressureDiff(edgeMassFlows, rho, deltaP, dDeltaP_dm);

    // Compute pressure difference.
    matVecMult(loopMatrix, deltaP, residuals, 0); // In loops.
    matVecMult(regulatedPathMatrix, deltaP, residuals, nbLoops()); // In regulated paths.

    // Compute residual. In loops, the pressure difference must be zero (so subtract zero).
    for (size_t i=0; i<nbRegPaths(); i++) { residuals[nbLoops()+i] -= regulatingElements[i]->getTargetRegulatedPathPressureDiff(); } // In regulated paths, get residual by subtracting the target pressure difference.
}

double Network::computeResidualNorm(vector<double> const& residuals) {
    double residualNorm = 0.0;
    for (size_t i=0; i<residuals.size(); i++) { residualNorm += residuals[i]*residuals[i]; } // TODO make a scaled/normalized version ?
    residualNorm = sqrt(residualNorm);
    return residualNorm;
}

double Network::computeMuPrime(vector<double> const& residuals_n, vector<double> const& residuals_np1, double const& lambda_n, double const& residualNorm_n) {
    double denominator = 0.0;
    for (size_t i=0; i<residuals_np1.size(); i++) {
        double ele = residuals_np1[i]-(1-lambda_n)*residuals_n[i];
        denominator += ele*ele;
    }
    denominator = sqrt(denominator);
    return 0.5*residualNorm_n*lambda_n*lambda_n/denominator;
}

bool Network::hasConverged(vector<double> const& residuals) {
    for (auto const& residual : residuals) {
        if ( abs(residual)>100.f ) { return false; } // All residuals must be less than 100Pa=0.01bar to converge.
    }
    return true;
}

void Network::convergeHydraulic(float rho) {
    vector<float> edgeMassFlows (nbEdges());
    vector<float> deltaP (nbEdges());
    vector<float> dDeltaP_dm (nbEdges());

    vector<double> residuals_n(nbLoops()+nbRegPaths()), residuals_np1(nbLoops()+nbRegPaths());

    int jacobian_n = nbLoops()+nbRegPaths();
    double deltaX_n[jacobian_n];
    double jacobian[jacobian_n*jacobian_n];

    // NLEQ-RES algorithm (based on Newton's algorithm)
    double lambda_n = 1.0, lambdaPrime_n, mu_n, muPrime_n, muPrime_nm1=0., residualNorm_nm1=0., residualNorm_n, residualNorm_np1, theta_n;
    vector<float> x_n (nbLoops()+nbRegPaths());

    // Initial step.
    size_t nLoops = 0;
    size_t debugNumberResidualsComputed = 1;
    bool regularityTestFailed = false;
    for (size_t i=0; i<nbLoops(); i++) { x_n[i] = loopMassFlows[i]; }
    for (size_t i=0; i<nbRegPaths(); i++) { x_n[nbLoops()+i] = regulatingElements[i]->getPressureDiff(); }
    computeResiduals(edgeMassFlows, residuals_n, deltaP, dDeltaP_dm, rho);
    residualNorm_n = computeResidualNorm(residuals_n);

    // Loop until converged.
    while ( nLoops<70 and (not hasConverged(residuals_n)) and (not regularityTestFailed) ) {

        // Compute Jacobian.
        matDiagMatTransposedMult(loopMatrix, dDeltaP_dm, loopMatrix, jacobian, jacobian_n, 0, 0);
        matDiagMatTransposedMult(regulatedPathMatrix, dDeltaP_dm, loopMatrix, jacobian, jacobian_n, nbLoops(), 0);
        for (size_t i=0; i<nbLoops()+nbRegPaths(); i++) { for (size_t j=0; j<nbRegPaths(); j++) { jacobian[i+jacobian_n*(nbLoops()+j)] = jacobianRegulatedPathColumns[i][j]; } } // These columns are always the same (where precomputed).

        for (size_t i=0; i<residuals_n.size(); i++) { deltaX_n[i] = -residuals_n[i]; }

        solve_Ax_equal_b(jacobian,deltaX_n,jacobian_n);

        if (nLoops>0) {
            mu_n = muPrime_nm1*residualNorm_nm1/residualNorm_n;
            lambda_n = min(1., mu_n);
        }

        size_t nbGoToRegularityTest = 0, nbGoToTrialIterate = 0; // To avoid staying stuck in loop.

        // Regularity test
        regularityTest:
        if (lambda_n<0.00001f) {
            regularityTestFailed = true;
        } else {

            // Compute trial iterate x_np1 = x_n + lambda_n * Delatx_n
            trialIterate:
            for (size_t i=0; i<nbLoops(); i++) { loopMassFlows[i] = x_n[i] + lambda_n*deltaX_n[i]; }
            for (size_t i=0; i<nbRegPaths(); i++) { regulatingElements[i]->setPressureDiff(x_n[nbLoops()+i] + lambda_n*deltaX_n[nbLoops()+i]); }

            // Compute monitoring quantities.
            debugNumberResidualsComputed++; // debug
            computeResiduals(edgeMassFlows, residuals_np1, deltaP, dDeltaP_dm, rho);
            residualNorm_np1 = computeResidualNorm(residuals_np1);
            theta_n = residualNorm_np1/residualNorm_n;
            muPrime_n = computeMuPrime(residuals_n, residuals_np1, lambda_n, residualNorm_n);

            if (theta_n>=1.0 or theta_n>1.0-lambda_n*0.25) {
                lambdaPrime_n = min(muPrime_n, 0.5*lambda_n);
                lambda_n = lambdaPrime_n;
                if (nbGoToRegularityTest<70) {
                    nbGoToRegularityTest++;
                    goto regularityTest;
                } // Else continue and accept the step, instead of staying stuck in a loop.
            } else {
                lambdaPrime_n = min(1.0, muPrime_n);
                if (lambdaPrime_n>=4.0*lambda_n) {
                    lambda_n = lambdaPrime_n;
                    if (nbGoToTrialIterate<70) {
                        nbGoToTrialIterate++;
                        goto trialIterate;
                    } // Else continue and accept the step, instead of staying stuck in a loop.
                }
            }

            // Accept x_np1 as next iterate.
            for (size_t i=0; i<nbLoops(); i++) { x_n[i] = loopMassFlows[i]; }
            for (size_t i=0; i<nbRegPaths(); i++) { x_n[nbLoops()+i] = regulatingElements[i]->getPressureDiff(); }
            residuals_n = residuals_np1;
            residualNorm_nm1 = residualNorm_n;
            residualNorm_n = residualNorm_np1;
            muPrime_nm1 = muPrime_n;

            nLoops++;
        }
    }


    // Set the results
    for (size_t i=0 ; i<pipePairs.size() ; i++) {
        pipePairs[i]->hydraulicConverged(edgeMassFlows[2*i], deltaP[2*i], edgeMassFlows[2*i+1], deltaP[2*i+1], rho);
    }
    for (size_t i=0 ; i<thermalStationNodePairs.size() ; i++) {
        thermalStationNodePairs[i]->setMassFlow(edgeMassFlows[pipePairs.size()*2+i]);
        thermalStationNodePairs[i]->setPressureDiff(deltaP[pipePairs.size()*2+i]);
    }
    int iBeg = pipePairs.size()*2+thermalStationNodePairs.size();
    for (size_t i=0 ; i<substationNodePairs.size() ; i++) {
        substationNodePairs[i]->hydraulicConverged(edgeMassFlows.begin()+iBeg, deltaP.begin()+iBeg);
        iBeg += substationNodePairs[i]->nbEdges();
    }
}

void Network::convergeToEquilibrium(float rho, float cp, Climate *pClim, unsigned int day, unsigned int hour) {

    float sumErrM, sumDesiredM = 0.f/*, relErrM = 1.f, derivRelErrM = -1.f*/, sumDeltaKv, sumKv/*, relFluctuationKv = 1.f, derivRelErrP = -1.f*/;
    vector<float> prevRelErrM (0), prevRelErrP (0);
    float sumErrP, sumDesiredP/*, relErrP = 1.f*/, sumDeltaRpm, sumRpm/*, relFluctuationRpm = 1.f*/;
    for (auto& ssnp : substationNodePairs) { sumDesiredM += ssnp->getDesiredMassFlow(); } // Substations update valve positions.

    updateDesiredMassFlow(cp, pClim, day, hour); // Inform prosumers, so that they can determine the supply target temperature.
    updateOperationMode(); // Inform thermal stations whether they must switch to storage mode.

    size_t nLoops = 0;
    bool hasConverged;

    vector<float> relErr (substationNodePairs.size()+thermalStationNodePairs.size()); // debug

    do {
        float learningRate = pow(1.f+nLoops, -0.75); // Similar to a step-size or learning rate in gradient descent. This fulfills the Robbins-Monroe condition. //        float learningRate = pow(1.f+nLoops, -0.75); // Similar to a step-size or learning rate in gradient descent. This fulfills the Robbins-Monroe condition.

        sumDeltaKv = 0.f; sumKv = 0.f;
        sumDeltaRpm = 0.f; sumRpm = 0.f;
        sumErrM = 0.f;
        sumErrP = 0.f; sumDesiredP = 0.f;

        float massFlowSupplyToReturn = computeMassFlowSupplyToReturn(); // debug
        for (auto& ssnp : substationNodePairs) { ssnp->updateControlVariable(rho, sumDeltaRpm, sumRpm, sumDeltaKv, sumKv, learningRate, massFlowSupplyToReturn); } // Substations update the control variable (Mostly valve position. Can be pump rotation if prosumer feeding in heat).
        for (auto& tsnp : thermalStationNodePairs) { tsnp->updateControlVariable(rho, sumDeltaRpm, sumRpm, sumDeltaKv, sumKv, learningRate); } // Thermal stations update the control variable (Mostly pump rotational speed. Can be valve position if in storage mode).
        convergeHydraulic(rho);
        convergeThermal(rho, cp, pClim, day, hour);


        updateDesiredMassFlow(cp, pClim, day, hour);
        updateOperationMode();
        for (size_t i=0; i<substationNodePairs.size(); i++) { substationNodePairs[i]->relativeErrorMassFlow(relErr[i], sumErrM); }
        for (size_t i=0; i<thermalStationNodePairs.size(); i++) { thermalStationNodePairs[i]->relativeErrorPressureDiff(relErr[substationNodePairs.size()+i], sumErrP, sumDesiredP); }

        nLoops++;
        // Check for convergence (if substation massflows, and thermal station pressures, are within two percent of the setpoint values)
        hasConverged = true;
        for (size_t i=0; hasConverged and (i<relErr.size()); i++) {
            hasConverged = hasConverged and ( abs(relErr[i])<0.02f );
        }

    } while ( nLoops<150 and (not hasConverged) );
}

void Network::computeDerivative(float& deriv, vector<float>& prev, float const& curr, size_t const& idx) {
    if (idx==0) {
        prev.push_back(curr);
    }
    else if (idx==1) {
        prev.push_back(curr);
        deriv = prev[1]-prev[0];
    }
    else if (idx==2) {
        prev.push_back(curr);
        deriv = 0.5f*(prev[2]-prev[0]);
//        deriv = 0.5f*(prev[0]-4.f*prev[1]+3.f*prev[2]);
    } else {
        size_t idx2 = idx%3;
        size_t idx0 = (idx2+1)%3;
//        size_t idx1 = (idx0+1)%3;
        prev[idx2] = curr;
        deriv = 0.5f*(prev[idx2]-prev[idx0]);
//        deriv = 0.5f*(prev[idx0]-4.f*prev[idx1]+3.f*prev[idx2]);
    }
}

void Network::updateDesiredMassFlow(float const& cp, Climate* pClim, unsigned int day, unsigned int hour) {
    // Inform prosumers of the temperature at return node, the climate, day, and hour ; so that the producer substations can determine the desired mass flow.
    for (auto& ssnp : substationNodePairs) {
        ssnp->updateDesiredMassFlow(cp, pClim, day, hour);
    }
}

void Network::updateOperationMode() {
    // If the substation (of which some are prosumers) are overall injecting more mass flow in the supply, then the thermal stations cannot also inject, they more consume and store the extra heat.
    float sumSubstationDesiredMassFlow = 0.f;
    for (auto ssnp : substationNodePairs) {
        sumSubstationDesiredMassFlow += ssnp->getDesiredMassFlow(); // Sum all the desired mass flows.
    }
    for (auto tsnp : thermalStationNodePairs) {
        tsnp->getThermalStation()->updateOperationMode(sumSubstationDesiredMassFlow, thermalStationNodePairs.size()); // Inform the thermal stations.
    }
}

float Network::computeMassFlowSupplyToReturn() {
    float sum = 0.f;
    for (auto tsnp : thermalStationNodePairs) {
        float m = tsnp->getMassFlow();
        if (m<0.f) {
            sum -= m;
        }
    }
    for (auto ssnp : substationNodePairs) {
        float m = ssnp->getMassFlow();
        if (m>0.f) {
            sum += m;
        }
    }
    return sum;
}

float Network::avgPressureDiffSupplyReturn(Substation* apartFromMe) {
    float sum = 0.f, nb = 0.f;
    for (auto tsnp : thermalStationNodePairs) {
        sum -= tsnp->getPressureDiff();
        nb += 1;
    }
    for (auto ssnp : substationNodePairs) {
        if (not ssnp->hasSubstation(apartFromMe)) {
            sum += ssnp->getPressureDiff();
            nb += 1;
        }
    }
    return sum/nb;
}










ostream& operator<<(ostream& stream, vector<unsigned int> const& v) {
    for(size_t i=0 ; i+1<v.size() ; i++ ) { // Print all except last element.
        stream << v[i] << " ";
    }
    if ( v.size()>0 ) { stream << v.back(); } // Print last element.
    return stream;
}

ostream& operator<<(ostream& stream, vector<double> const& v) {
    for(size_t i=0 ; i+1<v.size() ; i++ ) { // Print all except last element.
        stream << v[i] << " ";
    }
    if ( v.size()>0 ) { stream << v.back(); } // Print last element.
    return stream;
}

ostream& operator<<(ostream& stream, Network* net) {
    return net->print(stream);
}




vector<unsigned int> NodePair::ids = {};

NodePair::NodePair(TiXmlHandle hdl, float initTemp) :
    connectedPipes(0), supplyTemperature(initTemp), returnTemperature(initTemp), supplyTemperatureRecord(0), returnTemperatureRecord(0) {

    if(hdl.ToElement() ) {

        // Get parameters.
        string name = "z";
        if ( hdl.ToElement()->QueryFloatAttribute(name.c_str(), &z) ) { throw string("Error in the XML file: a NodePair doesn't have attribute: '"+name+"'."); }

        name = "id";
        if ( hdl.ToElement()->QueryUnsignedAttribute(name.c_str(), &id) )   { throw string("Error in the XML file: a NodePair doesn't have attribute: '"+name+"'."); }

        // Check that id is unique.
        if ( std::find(ids.begin(), ids.end(), id) != ids.end() ) {
            throw string("Error in the XML file: multiple NodePairs have the same id="+to_string(id)+".");
        }
        ids.push_back(id);
    }
    else {
        // throw error ?
    }
}


NodePair::~NodePair() {
    ids.erase( find( ids.begin(), ids.end(), id ) );
}

void NodePair::addConnectedPipe(PipePair *pipe) {
    connectedPipes.push_back(pipe);
}

bool NodePair::areAllUpstreamsComputed(bool isSupply) {

    for (auto pipe : connectedPipes) { // Loop on inflowing pipes.
        if (pipe->massFlowsToMe(this, isSupply)) { // Flow is incoming.
            if (not pipe->getAlreadyTraversed()) { // If temperature and heat loss have not yet been computed.
                return false;
            }
        }
    }
    return true;
}

void NodePair::computeAndSetTemperature(float initSumTM, float initSumM, bool isSupply) {
    float sumTM = initSumTM;
    float sumM = initSumM;
    for (auto pipe : connectedPipes) { // Loop on inflowing pipes (compute their temperatures).
        if (pipe->massFlowsToMe(this, isSupply)) {
            float massFlow = abs(pipe->getPipe(isSupply)->getMassFlow());
            float incomingTemp = pipe->getPipe(isSupply)->getDownstreamTemp();
            sumTM += massFlow * incomingTemp;
            sumM += massFlow;
        }
    }
    if (abs(sumM)>0.f) { // To avoid division by zero, this should not happen.
        setTemperature(sumTM/sumM, isSupply);
    } // Else do nothing, so keep the old temperature.
}

void NodePair::propagateDownstream(float soilTemp, float cp, bool isSupply) {
    // Compute heat losses in ouflowing pipes.
    for (auto pipe : connectedPipes) { // Loop on outflowing pipes (compute their thermal losses).
        if (not pipe->massFlowsToMe(this, isSupply)) { // Flow is outgoing.
            pipe->computeThermalLoss(soilTemp, cp, isSupply);
            pipe->setAlreadyTraversed(true); // It is important that all these flags are set before doing pipe->getOtherNode(this)->propagateNetwork
        }
    }

    // Propagate to the ouflowing nodes.
    for (auto pipe : connectedPipes) {
        if (not pipe->massFlowsToMe(this, isSupply)) { //Flow is outgoing.
            pipe->getOtherNode(this)->propagateNetwork(soilTemp, cp, isSupply);
        }
    }
}

void NodePair::writeTHHeaderText(fstream& textFile, unsigned int decId) {
    string prefix = getPrefix(decId);
    textFile << prefix << ":SupplyTemp(celsius)" <<"\t"
             << prefix << ":ReturnTemp(celsius)" <<"\t";
}

void NodePair::writeTHResultsText(fstream& textFile, unsigned int i) {
    textFile << fixed << setprecision(1) << getSupplyTemperatureRecord(i) <<"\t";
    textFile << fixed << setprecision(1) << getReturnTemperatureRecord(i) <<"\t";
}

bool NodePair::propagateNetwork(float soilTemp, float cp, bool isSupply) {
    bool isTraversed;
    if (areAllUpstreamsComputed(isSupply)) {
        float initSumTM = 0.f;
        float initSumM = 0.f;
        computeAndSetTemperature(initSumTM, initSumM, isSupply);

        propagateDownstream(soilTemp, cp, isSupply);
        isTraversed = true;
    }
    else { isTraversed = false; }
    return isTraversed;
}

void NodePair::checkNbPipes() {
    if (connectedPipes.size()<2) {
        throw string("Error in the XML file: NodePair of id="+to_string(id)+" is only connected to less than 2 PipePair, there must be at least 2 (or the node should be a SubstationNodePair or ThermalStationNodePair).");
    }
}






void NodePairConnectingSupplyReturn::writeTHHeaderText(fstream& textFile, unsigned int decId) {
    NodePair::writeTHHeaderText(textFile, decId);
    string prefix = getPrefix(decId);
    textFile << prefix << ":MassFlow(kg/s)" <<"\t"
             << prefix << ":PressureDiff(Pa)" <<"\t";
}

void NodePairConnectingSupplyReturn::writeTHResultsText(fstream& textFile, unsigned int i) {
    NodePair::writeTHResultsText(textFile, i);
    textFile << fixed << setprecision(5) << getMassFlowRecord(i) <<"\t";
    textFile << fixed << setprecision(0) << getPressureDiffRecord(i) <<"\t";
}






void SubstationNodePair::setSubstation(Substation* sub) {
    if (substation!=nullptr)  { throw string("Error in the XML file: the NodePair of id="+to_string(getId())+" is being referred by multiple Substations, there can only be one (if necessary, add small pipes to the network)."); }
    substation = sub;
}

void SubstationNodePair::computePressureDiffAndDerivative(vector<float>::const_iterator m, float const& rho, vector<float>::iterator deltaP, vector<float>::iterator dDeltaP_dm) {
    substation->computePressureDiffAndDerivative(m, rho, deltaP, dDeltaP_dm);
}

void SubstationNodePair::hydraulicConverged(vector<float>::const_iterator m, vector<float>::const_iterator deltaP) {
    for (size_t i=1; int(i)<nbEdges(); i++) {
        if ( abs((*m)-(*(m+i))) > 0.0001f*(*m) ) { throw string("Error in SubstationNodePair::hydraulicConverged, not all massflows are equal."); }
    }
    massFlow = *m;
    pressureDiff = 0.f;
    for (int i=0; i<nbEdges(); i++) {
        pressureDiff += *deltaP;
        deltaP++;
    }
}

void SubstationNodePair::computeThermal(float const& cp, float const& rho, Climate* pClim, unsigned int day, unsigned int hour, bool supplyToReturn) {
    if (  (supplyToReturn and massFlow>0.f)  or  ((not supplyToReturn) and massFlow<0.f)  ) {
        float inputTemp;
        if (massFlow>=0.f) { inputTemp = getSupplyTemperature(); }
        else { inputTemp = getReturnTemperature(); }

        substation->computeHeatExchanged(cp, rho, inputTemp, massFlow, pressureDiff, pClim, day, hour);
    }
}


bool SubstationNodePair::propagateNetwork(float soilTemp, float cp, bool isSupply) {
    bool isTraversed = false;
    if (not alreadyTraversed) {
        if ( (isSupply and getMassFlow()<0.f) or ((not isSupply) and getMassFlow()>0.f) ) { // TODO: could make a unified version for NodePairConnectingSupplyReturn
            if (areAllUpstreamsComputed(isSupply)) {
                float initSumTM = abs(getMassFlow()) * substation->getPrimarySideOutputTemp();
                float initSumM = abs(getMassFlow());
                computeAndSetTemperature(initSumTM, initSumM, isSupply);

                alreadyTraversed = true;
                isTraversed = true;

                propagateDownstream(soilTemp, cp, isSupply);
            }
        }
        else {
            isTraversed = NodePair::propagateNetwork(soilTemp, cp, isSupply);
            alreadyTraversed = isTraversed;
        }
    }
    return isTraversed;
}

void SubstationNodePair::checkNbPipes() {
    if (getNbPipes()<1) {
        throw string("Error in the XML file: SubstationNodePair of id="+to_string(getId())+" is connected to "+to_string(getNbPipes())+" PipePair, it must be connected to at least 1 PipePair.");
    }
}

void SubstationNodePair::writeTHHeaderText(fstream& textFile, unsigned int decId) {
    string prefix = getPrefix(decId);
    NodePairConnectingSupplyReturn::writeTHHeaderText(textFile, decId);
    substation->writeTHHeaderText(textFile, prefix);
}

void SubstationNodePair::writeTHResultsText(fstream& textFile, unsigned int i) {
    NodePairConnectingSupplyReturn::writeTHResultsText(textFile, i);
    substation->writeTHResultsText(textFile, i);
}







void ThermalStationNodePair::setThermalStation(ThermalStation* st) {
    if (thermalStation!=nullptr)  { throw string("Error in the XML file: a NodePair is being referred by multiple thermal stations, there can only be one."); }
    thermalStation = st;
}

void ThermalStationNodePair::computeThermal(float rho, float cp, Climate* pClim, unsigned int day, unsigned int hour, bool supplyToReturn) {
    if (  (supplyToReturn and massFlow<0.f)  or  ((not supplyToReturn) and massFlow>0.f)  ) {
        float inputTemp;
        if (massFlow>=0.f) { inputTemp = getReturnTemperature(); }
        else { inputTemp = getSupplyTemperature(); }

        thermalStation->computeThermal(pressureDiff, massFlow, rho, cp, inputTemp, pClim, day, hour);
    }
}


bool ThermalStationNodePair::propagateNetwork(float soilTemp, float cp, bool isSupply) {
    bool isTraversed = false;
    if (not alreadyTraversed) {
        if ( (isSupply and getMassFlow()>0.f) or ((not isSupply) and getMassFlow()<0.f) ) {
            if (areAllUpstreamsComputed(isSupply)) {
                float initSumTM = abs(getMassFlow()) * thermalStation->getOutputTemperature();
                float initSumM = abs(getMassFlow());
                computeAndSetTemperature(initSumTM, initSumM, isSupply);

                alreadyTraversed = true;
                isTraversed = true;

                propagateDownstream(soilTemp, cp, isSupply);
            }
        }
        else {
            isTraversed = NodePair::propagateNetwork(soilTemp, cp, isSupply);
            alreadyTraversed = isTraversed;
        }
    }
    return isTraversed;
}

void ThermalStationNodePair::checkNbPipes() {
    if (getNbPipes()<1) {
        throw string("Error in the XML file: ThermalStationNodePair of id="+to_string(getId())+" is connected to "+to_string(getNbPipes())+" PipePair, it must be connected to at least 1 PipePair.");
    }
}

void ThermalStationNodePair::writeTHHeaderText(fstream& textFile, unsigned int decId) {
    string prefix = getPrefix(decId);
    NodePairConnectingSupplyReturn::writeTHHeaderText(textFile, decId);
    thermalStation->writeTHHeaderText(textFile, prefix);
}

void ThermalStationNodePair::writeTHResultsText(fstream& textFile, unsigned int i) {
    NodePairConnectingSupplyReturn::writeTHResultsText(textFile, i);
    thermalStation->writeTHResultsText(textFile, i);
}







vector<unsigned int> PipePair::ids = {};

PipePair::PipePair(TiXmlHandle hdl, Network* net) : alreadyTraversed(false) {

    if(hdl.ToElement() ) {

        unsigned int node1Id, node2Id;

        // Get parameters.
        if ( hdl.ToElement()->QueryUnsignedAttribute("id", &id) )                               { throw string("Error in the XML file: a PipePair doesn't have attribute: 'id'."); }
        if ( hdl.ToElement()->QueryUnsignedAttribute("node1", &node1Id) )             { throw string("Error in the XML file: a PipePair doesn't have attribute: 'node1'."); }
        if ( hdl.ToElement()->QueryUnsignedAttribute("node2", &node2Id) )           { throw string("Error in the XML file: a PipePair doesn't have attribute: 'node2'."); }
        if ( hdl.ToElement()->QueryFloatAttribute("length", &length) )                          { throw string("Error in the XML file: a PipePair doesn't have attribute: 'length' in meters."); }
        if ( hdl.ToElement()->QueryFloatAttribute("innerRadius", &innerRadius) )                { throw string("Error in the XML file: a PipePair doesn't have attribute: 'innerRadius' in meters."); }
        if ( hdl.ToElement()->QueryFloatAttribute("interPipeDistance", &interPipeDistance) )    { throw string("Error in the XML file: a PipePair doesn't have attribute: 'interPipeDistance' in meters."); }

        connectedNodes[0] = net->pointerOfNodeWithId(node1Id);
        if (connectedNodes[0]==nullptr) { throw string("Error in the XML file: in Network of PipePair id="+to_string(id)+" the node1="+to_string(node1Id)+" is referenced, but isn't present in the Network of DistrictEnergyCenter id="+to_string(net->getDEC()->getId())+"."); }
        connectedNodes[1] = net->pointerOfNodeWithId(node2Id);
        if (connectedNodes[1]==nullptr) { throw string("Error in the XML file: in Network of PipePair id="+to_string(id)+" the node2="+to_string(node2Id)+" is referenced, but isn't present in the Network of DistrictEnergyCenter id="+to_string(net->getDEC()->getId())+"."); }

        for (auto pipe : connectedNodes[0]->getPipes()) { // Check for duplicated pipes.
            if (pipe->getOtherNode(connectedNodes[0])==connectedNodes[1]) {
                throw string("Error in the XML file: PipePair id="+to_string(id)+" connects the same two nodes as PipePair id="+to_string(pipe->getId())+".");
            }
        }
        connectedNodes[0]->addConnectedPipe(this);
        connectedNodes[1]->addConnectedPipe(this);

        // Sanity check.
        if ( length<=0 )            { throw string("Error in the XML file: PipePair id="+to_string(id)+" has length<=0."); }
        if ( innerRadius<=0 )       { throw string("Error in the XML file: PipePair id="+to_string(id)+" has innerRadius<=0."); }
        if ( interPipeDistance<=0 ) { throw string("Error in the XML file: PipePair id="+to_string(id)+" has interPipeDistance<=0."); }

        float initDownstreamTemp = 0.5f*(getHeadNode()->getSupplyTemperature()+getTailNode()->getSupplyTemperature()); // Arbitrarily set downstream temperature to average between the two node temperatures.
        TiXmlElement* supplyPipeHdl = hdl.ToElement()->FirstChildElement( "SupplyPipe" );
        if ( not supplyPipeHdl ) { throw string("Error in the XML file: a PipePair doesn't have any 'SupplyPipe' child element."); }
        else { supplyPipe = new Pipe( supplyPipeHdl, this, net->getSoilkValue(), initDownstreamTemp ); }

        initDownstreamTemp = 0.5f*(getHeadNode()->getReturnTemperature()+getTailNode()->getReturnTemperature()); // Arbitrarily set downstream temperature to average between the two node temperatures.
        TiXmlElement* returnPipeHdl = hdl.ToElement()->FirstChildElement( "ReturnPipe" );
        if ( not returnPipeHdl ) { throw string("Error in the XML file: a PipePair doesn't have any 'ReturnPipe' child element."); }
        else { returnPipe = new Pipe( returnPipeHdl, this, net->getSoilkValue(), initDownstreamTemp ); }

        interPipeThermalResistance = supplyPipe->computeInsulationResistance(innerRadius) + returnPipe->computeInsulationResistance(innerRadius) + computeInterPipeSoilThermalResistance(supplyPipe, returnPipe, interPipeDistance, innerRadius, net->getSoilkValue());

        // Check that id is unique. This is done last so that the id is added to ids only if the constructor finishes successfully.
        if ( std::find(ids.begin(), ids.end(), id) != ids.end() ) { throw string("Error in the XML file: multiple PipePairs have the same id="+to_string(id)+"."); }
        ids.push_back(id);
    }
    else {
        // throw error ?
    }
}

PipePair::~PipePair() {
    ids.erase( find( ids.begin(), ids.end(), id ) );
}


float PipePair::computeInterPipeSoilThermalResistance(Pipe* p1, Pipe* p2, float interPipeDistance, float innerRadius, float soilConductivity) {
    float OuterRad1 = innerRadius + p1->getInsulationThick();
    float OuterRad2 = innerRadius + p2->getInsulationThick();
    float x = (interPipeDistance*interPipeDistance - OuterRad1*OuterRad1 - OuterRad2*OuterRad2)/(2*OuterRad1*OuterRad2); // Or maybe if the interPipeDistance is the horizontal distance, then we should also take into account the difference in buried depth, so real distance between the pipes is sqrt(interpipehorizontaldist^2+diffBuriedDepth^2)
    return log(x + pow(x*x-1, 0.5)) / (2*3.14*soilConductivity);
}


NodePair* PipePair::getOtherNode(NodePair* np) {
    if (np==connectedNodes[0]) { return connectedNodes[1]; }
    if (np==connectedNodes[1]) { return connectedNodes[0]; }
    else { throw string("Error in PipePair::getOtherNode, the node given as input is not connected to this PipePair."); }
}


void PipePair::computeThermalLoss(float soilTemp, float cp, bool isSupply) {
    Pipe* pipeToCompute = getPipe(isSupply);
    Pipe* pipeTwin = getPipe(not isSupply);

    NodePair* inflowingNode;
    if (pipeToCompute->getMassFlow()>0) { inflowingNode = getTailNode(); }
    else { inflowingNode = getHeadNode(); }

    float inputTemp = inflowingNode->getTemperature(isSupply);
    float twinNearInputTemp, twinNearOutputTemp;

    if (pipeToCompute->getMassFlow()*pipeTwin->getMassFlow()>0) {
        twinNearInputTemp = inflowingNode->getTemperature(not isSupply);
        twinNearOutputTemp = pipeTwin->getDownstreamTemp();
    } else {
        twinNearInputTemp =  pipeTwin->getDownstreamTemp();
        twinNearOutputTemp = getOtherNode(inflowingNode)->getTemperature(not isSupply);
    }

    pipeToCompute->computeThermalLoss(inputTemp, twinNearInputTemp, twinNearOutputTemp, soilTemp, cp, length, interPipeThermalResistance);
}

bool PipePair::massFlowsToMe(NodePair* n, bool isSupply) {
    // Or simply check whether pressure of opposite node is higher ?
    bool nodeIsHead;
    if (n==connectedNodes[0]) { nodeIsHead = false; }
    else if (n==connectedNodes[1]) { nodeIsHead = true; }
    else { throw string("Error in PipePair::massFlowsToMe, the node given as input is not connected to this PipePair."); }

    Pipe* pipe;
    if (isSupply) { pipe = supplyPipe; }
    else { pipe = returnPipe; }
    bool flowsToHead = (pipe->getMassFlow()>0);

    return (nodeIsHead and flowsToHead) or ((not nodeIsHead) and (not flowsToHead));
}

float PipePair::computeAltitudePressureLoss(float const& rho) {
    return rho*9.81f*(getHeadNode()->getZ() - getTailNode()->getZ());
}

void PipePair::hydraulicConverged(float const& supplyMassFlow, float const& supplyDeltaP, float const& returnMassFlow, float const& returnDeltaP, float const& rho) {
    float altitudePressureLoss = computeAltitudePressureLoss(rho);
    supplyPipe->hydraulicConverged(supplyMassFlow, supplyDeltaP, rho, length, altitudePressureLoss);
    returnPipe->hydraulicConverged(returnMassFlow, returnDeltaP, rho, length, altitudePressureLoss);
}








float Pipe::computeSoilResistance(float pipeRadius, float soilThermalConductivity) {
    float x = buriedDepth/pipeRadius;
    if (x<=3) {
        return log(x + pow(x*x-1, 0.5)) / (2*3.14*soilThermalConductivity);
    } else {
        return log(2*x) / (2*3.14*soilThermalConductivity);
    }
}

float Pipe::computeInsulationResistance(float pipeRadius) {
    return log((pipeRadius+insulationThick)/pipeRadius) / (2*3.14*insulationkValue);
}

Pipe::Pipe(TiXmlHandle hdl, PipePair* parent, const float& soilThermalConductivity, const float& downstreamTemp) :
    massFlow(0.f), pressureDiff(0.f), thermalLoss(0.f), downstreamTemp(downstreamTemp), pressureLossHeat(0.f), massFlowRecord(0), pressureDiffRecord(0), thermalLossRecord(0) {
    if(hdl.ToElement() ) {
        if ( hdl.ToElement()->QueryFloatAttribute("insulationThick", &insulationThick) )        { throw string("Error in the XML file: a Pipe doesn't have attribute: 'insulationThick' in meters."); }
        if ( hdl.ToElement()->QueryFloatAttribute("insulationkValue", &insulationkValue) )      { throw string("Error in the XML file: a Pipe doesn't have attribute: 'insulationkValue' in W/(K*m)."); }
        if ( hdl.ToElement()->QueryFloatAttribute("buriedDepth", &buriedDepth) )                { throw string("Error in the XML file: a Pipe doesn't have attribute: 'buriedDepth' in meters."); }

        // Sanity check.
        unsigned int id = parent->getId();
        if ( insulationThick<=0 )   { throw string("Error in the XML file: PipePair id="+to_string(id)+" has a Pipe with insulationThick<=0."); }
        if ( insulationkValue<=0 )  { throw string("Error in the XML file: PipePair id="+to_string(id)+" has a Pipe with insulationkValue<=0."); }
        if ( buriedDepth<=0 )       { throw string("Error in the XML file: PipePair id="+to_string(id)+" has a Pipe with buriedDepth<=0."); }
        if ( buriedDepth <= parent->getInnerRadius() + insulationThick ) { throw string("Error in the XML file: a PipePair id="+to_string(id)+" has: buriedDepth<=innerRadius+insulationThick"); }

        thermalResistance = computeSoilResistance(parent->getInnerRadius()+insulationThick, soilThermalConductivity) + computeInsulationResistance(parent->getInnerRadius());
    }
    else {
        // throw error ?
    }

}

void Pipe::computeThermalLoss(float inputTemp, float twinNearInputTemp, float twinNearOutputTemp, float soilTemp, float cp, float length, float interPipeThermalResistance) {
    if (massFlow!=0.f) { // To avoid division by zero.
        float mcpInv = 1.f/(abs(massFlow)*cp);
        float c1 = - mcpInv*(1.f/thermalResistance+1.f/interPipeThermalResistance);
        float c1InvSq = 1.f/(c1*c1);
        float c2 = mcpInv*(soilTemp/thermalResistance + twinNearInputTemp/interPipeThermalResistance + pressureLossHeat);
        float c3 = mcpInv/(interPipeThermalResistance*length) * (twinNearOutputTemp-twinNearInputTemp);
        downstreamTemp = - ( c3 + c1*(c2 + c3*length) )*c1InvSq + exp(c1*length)*( inputTemp + (c3 + c1*c2)*c1InvSq );

        thermalLoss = (inputTemp-downstreamTemp)*abs(massFlow)*cp;
    }
    else { // If no massflow, suppose no heat losses during the overall time step, as no liquid is passing through.
        thermalLoss = 0.f;
    }
}

void Pipe::hydraulicConverged(float const& massFlow_, float const& deltaP, float const& rho, float const& length, float const& altitudePressureLoss) {
    massFlow = massFlow_;
    pressureDiff = deltaP;
    pressureLossHeat = abs( massFlow*(deltaP-altitudePressureLoss)/(rho*length) );
}






#pragma GCC diagnostic warning "-Wunused-parameter"
