#include "models.h"

#include "climate.h"
#include "scene.h"
#include "building.h"
#include "occupants.h"

// *** Model class, CitySim  *** //
// *** jerome.kaempf@epfl.ch *** //

float temperaturePrecision = 0.1f; // temperature precision for the models

bool Model::thermalExplicit = false; // is the thermal implicit model used

int Model::ThermalWarmUpTime(Building *pBuilding) {

    // liens
    unsigned int NP = pBuilding->getnNodes();

    // initialise the Atilde matrix
    double Atilde[NP*NP]; // new for dgesv
    for (unsigned int i=0;i<NP*NP;++i) Atilde[i]=0.;

    // initialise calculation matrix
    for (unsigned int i=0; i<pBuilding->getnZones(); ++i) {
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnNodes(); ++j) {
            for (unsigned int k=0; k<pBuilding->getZonenNodes(i); ++k) {
                // put the values in the G matrix
                Atilde[(Thermal_getMatrixPosition(pBuilding,i)+j)+NP*(Thermal_getMatrixPosition(pBuilding,i)+k)] = pBuilding->getZone(i)->getVariableMatrixElement(j,k);
            }
        }
    }

    //print_matrix("Atilde - 1",NP,NP,Atilde,NP);

    // Atilde will contain C^1 * G, with G being the total conductances (fixed and variable)
    for (unsigned int i=0;i<NP;++i) {
        for (unsigned int j=0;j<NP;++j) {
            // we sum the (timely) variable conductances and the (timely) fixed conductances of the building, and divide by the capacitance
            Atilde[i+NP*j] = (Atilde[i+NP*j]+pBuilding->getG1(i,j)) / pBuilding->getC(i,i);
        }
    }

    //print_matrix("Atilde - 2",NP,NP,Atilde,NP);

    // we solve the problem Atilde * tau = (-1, ..., -1)
    double btilde[NP];
    for (unsigned int i=0; i<NP; ++i) btilde[i] = -1.;
    solve_Ax_equal_b(Atilde,btilde,NP);

    // takes the biggest warm up time of all the nodes
    int warmUpDays = 0;
    for (unsigned int i=0;i<NP;++i) {
        //cout << "tau " << i << " = " << ceil(btilde[i]/3600./24.) << " days" << endl;
        warmUpDays = max(warmUpDays, static_cast<int>(ceil(btilde[i]/3600./24.)));
    }
    //cout << "warmUpDays: " << warmUpDays << endl;
    return warmUpDays;

}

void Model::ThermalStepImplicit(Building *pBuilding, Climate *pClimate, unsigned int day, unsigned int hour) {

    //outside temperature. In this model it doesn't matter if the temperature is in ∞C or in K, it solve an equadiff with difference.
    float Tout = pClimate->getToutCelsius(day,hour);
    float Tground = pClimate->getTgroundCelsius(day,hour);
    //cerr << "Tout: " << Tout << "\tTground: " << Tground << endl;

    // Evapotranspiration process on the roofs, sets the climate dependent variables
    for (size_t i=0;i<pBuilding->getnZones();++i) {
        for (size_t j=0;j<pBuilding->getZone(i)->getnRoofs();++j) pBuilding->getZone(i)->getRoof(j)->set_YX(pClimate,day,hour);
    }

    // liens
    unsigned int NP = pBuilding->getnNodes();

    // propre au problËme Ax=b
    double G[NP][NP];
    for (unsigned int i=0;i<NP;++i) { for (unsigned int j=0;j<NP;++j) G[i][j]=0.; }

    double b[NP], T[NP];

    /// propre au systËme ‡ rÈsoudre

    int np = pBuilding->getnZones();

    int mp = NP - np;

    double TnInf[np], TnSup[np];
    double heating[np], cooling[np];

    double Asecond[mp*mp];
    double phi[mp];

    // initializes the building electricity consumption and production (prior to the occupancy model)
    pBuilding->setElectricConsumption(0.f);
    pBuilding->setSolarPVProduction(0.f);

    for (unsigned int i=0;i<pBuilding->getnZones();++i) {

        // *** computing the occupants + activity heat gains *** //
        pBuilding->getZone(i)->setOccupantsCountAndActivity(day,hour);
        // JK - 10.07.2015: Note to be added to Lc: the gains due to lights

        //cerr << "Building id= " << pBuilding->getId() << "\t";
        //cerr << "Persons: " << pBuilding->getZone(i)->getOccupantsNumber() << "\tFraction present: " << pBuilding->getZone(i)->getOccupantsFraction(day,hour) << "\tLc: " << Lc << endl;

        // *** blinds procedure for the walls - deterministic fashion *** //
        if (!Model::thermalExplicit) pBuilding->deterministicShadingAction(/*day*/);

        // saving the internal gains in a vector (except the Qs)
        pBuilding->getZone(i)->setQi(pBuilding->getZone(i)->getQsun2() + pBuilding->getZone(i)->getConvectiveInternalHeatGains() + pBuilding->getZone(i)->getRadiativeInternalHeatGains());

        for (unsigned int j=0;j<pBuilding->getZone(i)->getnNodes();j++) {

            if (pBuilding->getZone(i)->getnT(j) == 0) T[Thermal_getMatrixPosition(pBuilding,i)+j]=15.0*pBuilding->getZoneC(i,j);
            else T[Thermal_getMatrixPosition(pBuilding,i)+j]=pBuilding->getZoneT(i,j)*pBuilding->getZoneC(i,j);

        }

        //cerr << "After definition of T." << endl;

        TnInf[i]=pBuilding->getZone(i)->getTmin()*pBuilding->getZone(i)->getC(0);
 //        cerr<<"TnInf " << TnInf[i]<<endl;

        TnSup[i]=pBuilding->getZone(i)->getTmax()*pBuilding->getZone(i)->getC(0);
//        cerr<<"TnSup " << TnSup[i]<<endl;

        // convectiveKe calculated with WindSpeed and WindDirection

//        cerr << "WindClarke: " << Thermal_KeClarke(pClimate->getWindSpeed(step*dt),pClimate->getWindDirection(step*dt), Tout, pBuilding->getZone(i)->getAzimuthWall(), pBuilding->getZone(i)->getSurfaceWall(), pBuilding->getZone(i)->getSwa()) << endl;
//        cerr << "Wind:       " << Thermal_Ke(pClimate->getWindSpeed(step*dt),pClimate->getWindDirection(step*dt), Tout, pBuilding->getZone(i)->getAzimuthWall(), pBuilding->getZone(i)->getSurfaceWall() ) << endl;

        pBuilding->getZone(i)->setKe(Thermal_KeWalls(pClimate,pBuilding->getZone(i),day,hour)); // convective
        pBuilding->getZone(i)->setHr(Thermal_HrWalls(pBuilding->getZone(i))); // radiative
        //pBuilding->getZone(i)->setKe(Thermal_KeClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour), azimuthWall, surfaceWall));
        //pBuilding->getZone(i)->setKe(Thermal_KeClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour), Tout, azimuthWall, surfaceWall, pBuilding->getZone(i)->getSwa() ));
        //pBuilding->getZone(i)->setKe(Thermal_KeClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour), Tout, pBuilding->getZone(i)->getWalls(), pBuilding->getZone(i)->getSwa() ));
        //pBuilding->getZone(i)->setKi(Thermal_KiFixed( pBuilding->getZone(i)->getSwa() )); - as fixed, put in the constructor of Zone JK - 18.04.2009
//        cerr << "Kappatot: " << pBuilding->getZone(i)->getKappatot() << endl;

        // no ventilation in the building for the prevision of the heating
        pBuilding->getZone(i)->setNvent(0.f);

        // calculation of UA
        //vector<float> azimuthWindow = pBuilding->getZone(i)->getAzimuthWindow();
        //vector<float> surfaceWindow = pBuilding->getZone(i)->getSurfaceWindow();
        // JK - 16.03.2014 - removed for BESTEST compatibility
        //pBuilding->getZone(i)->setKwindow(Thermal_Kwindows(pClimate,pBuilding->getZone(i),day,hour));
        pBuilding->getZone(i)->setKroof(Thermal_Kroofs(pClimate,pBuilding->getZone(i),day,hour));
//        cerr << "UA : "<<pBuilding->getZone(i)->getUA()<<endl;

     // show Kappa1 et Kappa2
//        cerr << "Kappa1: " << pBuilding->getZone(i)->getKappa1() << "\tKappa2: " << pBuilding->getZone(i)->getKappa2() << endl;

        // source term b
        for (unsigned int j=0;j<pBuilding->getZone(i)->getnNodes();++j) {
             b[Thermal_getMatrixPosition(pBuilding,i)+j] = pBuilding->getZone(i)->getSourceTerm(j,Tout,Tground);
        }

        //cerr << "After definition of the source term b." << endl;

    }

    // shows the source vector
    //cout << "source vector:" << b << endl;

    // initialise calculation matrix
    for (unsigned int i=0; i<pBuilding->getnZones(); ++i) {
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnNodes(); ++j) {
            for (unsigned int k=0; k<pBuilding->getZonenNodes(i); ++k) {
                // put the values in the G matrix
                G[Thermal_getMatrixPosition(pBuilding,i)+j][Thermal_getMatrixPosition(pBuilding,i)+k] = pBuilding->getZone(i)->getVariableMatrixElement(j,k);
            }
        }
    }

    // shows the conductance matrix
    //cout << "variable conductance matrix: " << G << endl;

    for (unsigned int i=0;i<NP;i++) {
        for (unsigned int j=0;j<NP;j++) {
            // G1 contains the (timely) fixed conductances of the building, that are added to the variable conductances initialised in G
            G[i][j]+=pBuilding->getG1(i,j); // previously G[i][j]=G1[i][j]+G2[i][j]
        }
    }

    // shows the conductance matrix
    //cout << "variable + fixed conductance matrix: " << G << endl;

    // Calculation of heating needs to reach TsetInf
	for (int i=0;i<np;i++) {
	  heating[i]=TnInf[i]-T[Thermal_getMatrixPosition(pBuilding,i)]-dt*b[Thermal_getMatrixPosition(pBuilding,i)];
	  for (int j=0;j<np;j++) {
	    heating[i]-=dt*G[Thermal_getMatrixPosition(pBuilding,i)][Thermal_getMatrixPosition(pBuilding,j)]*pBuilding->getZone(j)->getTmin();
	  }
	}

    for (int i=0;i<np;i++) {
        for (unsigned int j=1; j< pBuilding->getZonenNodes(i); j++) {

            phi[Thermal_getSubMatrixPosition(pBuilding ,i)+j-1]=T[Thermal_getMatrixPosition(pBuilding,i)+j] + dt*b[Thermal_getMatrixPosition(pBuilding,i)+j];
            for (int k=0;k<np;k++) {
                phi[Thermal_getSubMatrixPosition(pBuilding ,i)+j-1]+=dt*G[Thermal_getMatrixPosition(pBuilding,i)+j][Thermal_getMatrixPosition(pBuilding,k)]*pBuilding->getZone(k)->getTmin();
                for (unsigned int l=1; l< pBuilding->getZonenNodes(k); l++) {
                    Asecond[Thermal_getSubMatrixPosition(pBuilding ,i)+j-1 + mp*(Thermal_getSubMatrixPosition(pBuilding,k)+l-1)]=dt*G[Thermal_getMatrixPosition(pBuilding,i)+j][Thermal_getMatrixPosition(pBuilding,k)+l]-pBuilding->getC(Thermal_getMatrixPosition(pBuilding,i)+j,Thermal_getMatrixPosition(pBuilding,k)+l);
                }
            }
        }
    }

    // inverse with LAPACK
    inverse_square_matrix(Asecond,mp);
    // end of inverse with LAPACK

    double AsecondPhi[mp];
    // Product of Asecond^{-1} with phi
	for (int i=0;i<mp;i++) {
      AsecondPhi[i]=0.0;
	  for (int j=0;j<mp;j++) {
        AsecondPhi[i]+=Asecond[i+mp*j]*phi[j];
	  }
	}

    // Product of AsecondPhi with A matrix
	for (int i=0;i<np;i++) {
        for (unsigned int j=1; j< pBuilding->getZonenNodes(i); j++) {
	        heating[i]+=dt*G[Thermal_getMatrixPosition(pBuilding,i)][Thermal_getMatrixPosition(pBuilding,i)+j]*AsecondPhi[Thermal_getSubMatrixPosition(pBuilding ,i)+j-1];
        }
	}
    // End of calculation of heating needs to reach TsetInf

    // Calculation of cooling needs to reach TsetSup
	for (int i=0;i<np;i++) {
	  cooling[i]=TnSup[i]-T[Thermal_getMatrixPosition(pBuilding,i)]-dt*b[Thermal_getMatrixPosition(pBuilding,i)];
	  for (int j=0;j<np;++j) {
	    cooling[i]-=dt*G[Thermal_getMatrixPosition(pBuilding,i)][Thermal_getMatrixPosition(pBuilding,j)]*pBuilding->getZone(j)->getTmax();
	  }
	}

    for (int i=0;i<np;i++) {
        for (unsigned int j=1; j< pBuilding->getZonenNodes(i); j++) {

            phi[Thermal_getSubMatrixPosition(pBuilding ,i)+j-1]=T[Thermal_getMatrixPosition(pBuilding,i)+j] + dt*b[Thermal_getMatrixPosition(pBuilding,i)+j];

            for (int k=0;k<np;k++) {
                phi[Thermal_getSubMatrixPosition(pBuilding ,i)+j-1]+=dt*G[Thermal_getMatrixPosition(pBuilding,i)+j][Thermal_getMatrixPosition(pBuilding,k)]*pBuilding->getZone(k)->getTmax();
            }
        }
    }

    // Product of Asecond^{-1} with phi
	for (int i=0;i<mp;i++) {
      AsecondPhi[i]=0.0;
	  for (int j=0;j<mp;j++) {
        AsecondPhi[i]+=Asecond[i+mp*j]*phi[j];
	  }
	}

    // Product of AsecondPhi with A matrix
	for (int i=0;i<np;i++) {
        for (unsigned int j=1; j< pBuilding->getZonenNodes(i); j++) {
	        cooling[i]+=dt*G[Thermal_getMatrixPosition(pBuilding,i)][Thermal_getMatrixPosition(pBuilding,i)+j]*AsecondPhi[Thermal_getSubMatrixPosition(pBuilding,i)+j-1];
        }
	}
    // End of calculation of cooling needs to reach TsetSup

    //	 Choosing the heating or cooling load
	for (int i=0;i<np;i++) {
        // heating necessary if positive
        pBuilding->getZone(i)->setHeating(max(heating[i]/dt,0.));
	    //cerr << "Heating air: " << max(heating[i]/dt,0.) << " W constant during the time step "<< (24*(day-1)+hour) << endl;
	    b[Thermal_getMatrixPosition(pBuilding,i)]+=max(heating[i]/dt,0.);

        // cooling necessary if negative
	    pBuilding->getZone(i)->setCooling(min(cooling[i]/dt,0.));
	    //cerr << "Cooling air: " << min(cooling[i]/dt,0.) << " W constant during the time step "<< (24*(day-1)+hour) << endl;
        b[Thermal_getMatrixPosition(pBuilding,i)]+=min(cooling[i]/dt,0.);
	}

    // sauvegarde de la tempÈrature prÈvue (qui sera la bonne si tout va bien)
    double Aprime[NP*NP];
    double bprime[NP];
    for (unsigned int i=0;i<NP;i++) {
      for (unsigned int j=0;j<NP;j++) {
        Aprime[i+NP*j]=pBuilding->getC(i,j)-dt*G[i][j];
      }
      bprime[i]=T[i]+dt*b[i];
    }

    solve_Ax_equal_b(Aprime,bprime,NP);

    for (unsigned int i=0; i<pBuilding->getnZones(); i++) {
            pBuilding->getZone(i)->setTaForeseen(bprime[Thermal_getMatrixPosition(pBuilding,i)]);
            //cerr << "Ta Foreseen: " << bprime[Thermal_getMatrixPosition(pBuilding,i)] << endl;
    }

//    // shows the heating vector
//    cout << "heating vector:" << heating << endl;
//    cout << "cooling vector:" << cooling << endl;
//    cout << "T foreseen vector:" << bprime << endl;

    return;

}

void Model::ThermalStepImplicitTemperature(Building *pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {

    // outside air temperature
    float Tout = pClimate->getToutCelsius(day,hour);
    float Tground = pClimate->getTgroundCelsius(day,hour);
    //cerr << "Tout: " << Tout << "\tTground: " << Tground << endl;

    // liens
    unsigned int NP = pBuilding->getnNodes();

	// propre au problËme Ax=b
    double G[NP][NP];
    for (unsigned int i=0;i<NP;++i) { for (unsigned int j=0;j<NP;++j) G[i][j]=0.; }
    double b[NP], T[NP];

    // sensible heat provided to each thermal zone
    float Qs = 0.f;

    // loop on the thermal zones
    for (unsigned int i=0;i<pBuilding->getnZones();i++) {

        Qs = pBuilding->getZone(i)->getQs();
//        cout << "Qs: " << Qs*3600.f << endl;

        for (unsigned int j=0;j<pBuilding->getZone(i)->getnNodes();j++) {

            if (pBuilding->getZone(i)->getnT(j) == 0) T[Thermal_getMatrixPosition(pBuilding,i)+j]=15.0*pBuilding->getZoneC(i,j);
            else T[Thermal_getMatrixPosition(pBuilding,i)+j]=pBuilding->getZoneT(i,j)*pBuilding->getZoneC(i,j);

        }

        // convectiveKe calculated with WindSpeed and WindDirection
//        cerr << "WindClarke: " << Thermal_KeClarke(pClimate->getWindSpeed(step*dt),pClimate->getWindDirection(step*dt), pClimate->getToutCelsius(step*dt), pBuilding->getZone(i)->getAzimuthWall(), pBuilding->getZone(i)->getSurfaceWall(), pBuilding->getZone(i)->getSwa()) << endl;
//        cerr << "Wind:       " << Thermal_Ke(pClimate->getWindSpeed(step*dt),pClimate->getWindDirection(step*dt), pClimate->getToutCelsius(step*dt), pBuilding->getZone(i)->getAzimuthWall(), pBuilding->getZone(i)->getSurfaceWall() ) << endl;

        // JK - 11.04.2013: remove the creation of the two vectors used by setKe
        // JK - 11.04.2013: removed setKe as we are keeping the same from the ThermalStepImplicit (in a vector format)
        //pBuilding->getZone(i)->setKe(Thermal_KeClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour), Tout, azimuthWall, surfaceWall, pBuilding->getZone(i)->getSwa() ));
        //pBuilding->getZone(i)->setKi(Thermal_KiFixed( pBuilding->getZone(i)->getSwa() ));
//        cerr << "Kappatot: " << pBuilding->getZone(i)->getKappatot() << endl;

// write the nvent in a file!
//    fstream textFile("nvent.dat",ios::out|ios::binary|ios::app);
//    textFile << pBuilding->getZone(i)->getNvent() << endl;
//    textFile.close();

        // *** windows opening model in a deterministic fashion *** //
        if (pBuilding->getZone(i)->getTaForeseen() > (pBuilding->getZone(i)->getTmin()+1.0) && (pBuilding->getZone(i)->getTaForeseen()-Tout) > 1.0) {
            // people present && temperature in the non-controlled zone => open the windows
            //cerr << "hour: " << hour << "\tTi: " << pBuilding->getZone(i)->getTa() << "\tTout: " << Tout << "\tNvent: " << Model::deterministicWindowsNvent(2.0,pBuilding->getZone(i)->getTa(),Tout) << endl;
            //pBuilding->getZone(i)->setNvent(Model::deterministicWindowsNvent(2.0,pBuilding->getZone(i)->getTa(),Tout));

            float VdotVent = Model::ventilationFlowRate(pBuilding->getZone(i)->getTaForeseen()+273.15f, pClimate->getToutCelsius(day,hour)+273.15f, pClimate->getWindSpeed(day,hour), pBuilding->getZone(i)->getSwiO());
            //cerr << "timeStep: " << ((day-1)*24 + hour -1) << ", day: " << day << ", hour: " << hour << "\tTi: " << pBuilding->getZone(i)->getTaForeseen() << "\tTout: " << Tout << "\tVdotVent: " << VdotVent << endl;

            float depletionTime = 0.f;
            if (VdotVent > 0.f) {
                // linear approximation
                //depletionTime = (pBuilding->getZone(i)->getVi()/VdotVent)*((pBuilding->getZone(i)->getTaForeseen()-max(pBuilding->getTmin(),Tout+temperaturePrecision))/(pBuilding->getZone(i)->getTaForeseen()-Tout));
                // log approximation
                depletionTime = (pBuilding->getZone(i)->getVi()/VdotVent)*log((Tout-pBuilding->getZone(i)->getTaForeseen())/(Tout-max(pBuilding->getZone(i)->getTmin(),Tout+temperaturePrecision)));
                //cerr << "depletion time: " << depletionTime << endl;
                depletionTime = min(depletionTime,(float)(Model::dt));
            }
            //cerr << "time (s): " << depletionTime << "\tratio: " << depletionTime/(float)(Model::dt) << "\tVdotVent\': " << (depletionTime/(float)(Model::dt))*VdotVent << endl;
            pBuilding->getZone(i)->setVdotVent((depletionTime/(float)(Model::dt))*VdotVent);
        }
        else pBuilding->getZone(i)->setVdotVent(0.f);

//        cout << "VdotVent: " << 3600.f*pBuilding->getZone(i)->getVdotVent() << endl;

        // *** determination of the lights state and electric consumption *** //
        // DP : commented for now, create trouble when DayLight is not computed (see doDayLightSim parameter of XmlScene::simulateTimeStep) // Cognet: Dapeng has removed this, should it be put back?
        //Model::lightAction_Threshold(pBuilding->getZone(i));
        //pBuilding->addElectricConsumption(Model::lightsElectricConsumption(pBuilding->getZone(i)));

        // *** calculation of UA, which depends on the VdotVent *** -> in getUA() //

        // show Kappa1 et Kappa2
//        cerr << "Kappa1: " << pBuilding->getZone(i)->getKappa1() << "\tKappa2: " << pBuilding->getZone(i)->getKappa2() << endl;

        // source terms b, here we find the same source as for the prevision PLUS the term Qs
        for (unsigned int j=0;j<pBuilding->getZone(i)->getnNodes();++j) {
             b[Thermal_getMatrixPosition(pBuilding,i)+j] = pBuilding->getZone(i)->getSourceTerm(j,Tout,Tground,Qs);
        }

    }

    // prints source terms
    //cout << "source term: " << b << endl;

    // initialise calculation matrix
    for (unsigned int i=0; i<pBuilding->getnZones(); ++i) {
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnNodes(); ++j) {
            for (unsigned int k=0; k<pBuilding->getZonenNodes(i); ++k) {
                // put the values in the G matrix
                G[Thermal_getMatrixPosition(pBuilding,i)+j][Thermal_getMatrixPosition(pBuilding,i)+k] = pBuilding->getZone(i)->getVariableMatrixElement(j,k);
            }
        }
    }

    for (unsigned int i=0;i<NP;i++) {
        for (unsigned int j=0;j<NP;j++) {
            // G1 contains the (timely) fixed conductances of the building, that are added to the variable conductances initialised in G
            G[i][j]+=pBuilding->getG1(i,j); // previously G[i][j]=G1[i][j]+G2[i][j]
        }
    }

    // evaluation of the temperatures
    double Aprime[NP*NP];
    double bprime[NP];
    for (unsigned int i=0;i<NP;i++) {
      for (unsigned int j=0;j<NP;j++) {
        Aprime[i+NP*j]=pBuilding->getC(i,j)-dt*G[i][j];
      }
      bprime[i]=T[i]+dt*b[i];
    }

    solve_Ax_equal_b(Aprime,bprime,NP);

    for (unsigned int i=0; i<pBuilding->getnZones(); i++) {
        // saves the results in the zone (Ta,Tw, eventually Tw2)
        for (unsigned int j=0; j<pBuilding->getZonenNodes(i); j++) {
            //if (pBuilding->getZone(i)->getTaForeseen() > (pBuilding->getTmin()+1.0) && (pBuilding->getZone(i)->getTaForeseen()-Tout) > 1.0)
            //    cerr << "Temperature: " << bprime[Thermal_getMatrixPosition(pBuilding,i)+j] << endl;
            pBuilding->getZone(i)->setT(j,bprime[Thermal_getMatrixPosition(pBuilding,i)+j]);
//            cerr << "Zone: " << i << "\tNode: " << j << "\tT imp2: " << bprime[Thermal_getMatrixPosition(pBuilding,i)+j] << "\t Text: " << pClimate->getToutCelsius(step*dt) << endl;
        }

        // saves the outside surface temperature in the surface corresponding to the zone (which is calculated from Tw - wall and Ta - roofs)
        pBuilding->getZone(i)->setTos(Tout);

        // Evapotranspiration process on the roofs, sets the water evaporated
        for (size_t j=0;j<pBuilding->getZone(i)->getnRoofs();++j) pBuilding->getZone(i)->getRoof(j)->setWaterEvapotranspiration();

    // end of the loop on Zones
    }

    return;

}

// JHK - modified for ground model temperature 07.10.2013
void Model::ThermalStepImplicitTemperature(Ground *pGround, Climate* pClimate, unsigned int day, unsigned int hour) {

	// set hc for the ground
    Thermal_Kgrounds(pClimate, pGround, day, hour);

    // checks if the type of Ground is defined
    if (pGround->getComposite()==NULL) {
        // returns the surface temperature from the climate file and do not set the layer temperature as they don't exist
        pGround->setTemperature(pClimate->getTgroundCelsius(day,hour));
        // saves the water consumed
        pGround->setWaterEvapotranspiration(0.f); // 694.5 is the latent heat of water, considered as constant with water temperature
        return;
    }

    // outside air temperature
    float Tout = pClimate->getToutCelsius(day,hour);
    // gets the ground temperature according to Darren's model
    float Tground = pClimate->getTgroundCelsius(day,hour,pGround->getComposite()->getDepth(),pGround->getComposite()->getDiffusivity(),pGround->getComposite()->getDepth());

    // Evapotranspiration process on the grounds, sets the climate dependent variables
    //pGround->set_YX(pClimate,day,hour);

    // define Psi and Lambda, Psi + Lambda*g1*(theta_1 - theta_s)
    float deltaT = ((17.502*240.97)/pow(Tout+240.97,2))*pClimate->getSaturatedVapourPressure(Tout);
    float psychrometricConstant = pClimate->getPatm(day,hour)/1000.*0.665e-3;
    float localWindSpeed = ((pClimate->getWindSpeed(day,hour)>0)?pClimate->getWindSpeed(day,hour):0.01);
    float term2 = ((pClimate->getIdh(day,hour)>0.f)?((1.+0.34*localWindSpeed)*psychrometricConstant + deltaT)
                                                   :((1.+1.70*localWindSpeed)*psychrometricConstant + deltaT) );
    float latentHeat = (2501.-2.37*Tout)/3.6; // in Wh/kg
    float Psi = (((psychrometricConstant*(37./(Tout+273.15))*localWindSpeed)*(pClimate->getSaturatedVapourPressure(Tout)-pClimate->getVapourPressure(day,hour)))/term2)*latentHeat;
    float Lambda = deltaT/term2;

    // get the number of layers in the ground
    int NP = pGround->getComposite()->getnLayers();

    // source term of the equation (first and last are non-zero)
    double b[NP];
//  b[0] = (pGround->getk(0)/(pGround->get_hc()+pGround->get_hr()))                   *(pGround->get_hc()*Tout+(1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance()+pGround->get_hr()*pGround->getEnvironmentalTemperature());
    b[0] = (pGround->getG(0)/( (1.+Lambda)*pGround->getG(0)+pGround->get_hr()+pGround->get_hc()))*(pGround->get_hc()*Tout+(1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance()+pGround->get_hr()*pGround->getEnvironmentalTemperature()+Psi);
    for (int i=1;i<NP-1;++i) b[i]=0.;
    b[NP-1] = pGround->getk(NP)*Tground;

    // evaluation of the temperatures
    double Aprime[NP*NP];
    double bprime[NP];
    for (int i=0;i<NP;++i) {
      for (int j=0;j<NP;++j) {
        Aprime[i+NP*j] = pGround->getCapacitance(i,j) - dt*pGround->getConductance(i,j);
      }
      bprime[i] = pGround->getLayerTemperature(i)*pGround->getCapacitance(i,i) + dt*b[i];
    }

    // evapotranspiration (ET) modification term by GU
    Aprime[0] = pGround->getCapacitance(0,0) - dt*( -pGround->getk(1) - pGround->getG(0)*(pGround->get_hr()+pGround->get_hc())/((1.+Lambda)*pGround->getG(0)+pGround->get_hr()+pGround->get_hc()) );

    solve_Ax_equal_b(Aprime,bprime,NP);

    // saves the results in the ground
    for (int i=0; i<NP; ++i) {
        pGround->setLayerTemperature(i,bprime[i]);
    }

    // saves the outside surface temperature in the surface corresponding to the ground, bprime[0] is the first layer temperature
    //float surfaceTemperature = (bprime[0]*pGround->getG(0) + Tout*pGround->get_hc() + (1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance() + pGround->get_hr()*pGround->getEnvironmentalTemperature())/(pGround->get_hc()+pGround->getG(0)+pGround->get_hr());
    float surfaceTemperature = ((1.+Lambda)*pGround->getG(0)*bprime[0] + Tout*pGround->get_hc() + (1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance() + pGround->get_hr()*pGround->getEnvironmentalTemperature() + Psi)/((1.+Lambda)*pGround->getG(0)+pGround->get_hr()+pGround->get_hc());
    pGround->setTemperature(surfaceTemperature);

    // saves the water consumed
    pGround->setWaterEvapotranspiration((Psi+Lambda*pGround->getG(0)*(bprime[0]-surfaceTemperature))/latentHeat); // latent heat of water in Wh/kg -> W/m2 * Wh/kg = kg/(m2 h) = mm/h

    return;

}

// JHK - modified for a simplified version of the ground model temperature 08.04.2014
void Model::ThermalStepImplicitTemperature_simplified(Ground *pGround, Climate* pClimate, unsigned int day, unsigned int hour) {

	// set hc for the ground
    Thermal_Kgrounds(pClimate, pGround, day, hour);

    // checks if the type of Ground is defined
    if (pGround->getComposite()==NULL) {
        // returns the surface temperature from the climate file and do not set the layer temperature as they don't exist
        pGround->setTemperature(pClimate->getTgroundCelsius(day,hour));
        // saves the water consumed -> 0
        pGround->setWaterEvapotranspiration(0.f);
        return;
    }

    // outside air temperature
    float Tout = pClimate->getToutCelsius(day,hour);
    // gets the ground temperature according to Darren's model
    float Tground = pClimate->getTgroundCelsius(day,hour,pGround->getComposite()->getDepth(),pGround->getComposite()->getDiffusivity(),pGround->getComposite()->getDepth());

    // Evapotranspiration process on the grounds, sets the climate dependent variables
    pGround->set_YX(pClimate,day,hour);

	// element of the matrix A
	//float alpha1 = -pGround->getG1()*(pGround->get_hc()+pGround->get_hr())/(pGround->get_hc()+pGround->getG1()+pGround->get_hr()) - pGround->getG2();
    float alpha1 = -pGround->getG1() + (pGround->getG1()*pGround->getG1()/(pGround->get_hc()+pGround->getG1()+pGround->get_hr()+pGround->get_X())) - pGround->getG2();

    // source term of the equation
    //float beta1 = pGround->getG2()*Tground + pGround->getG1()*(pGround->get_hc()*Tout+(1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance()+pGround->get_hr()*pGround->getEnvironmentalTemperature())/(pGround->get_hc()+pGround->getG1()+pGround->get_hr());
    float beta1 = pGround->getG2()*Tground + pGround->getG1()*(pGround->get_hc()*Tout+(1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance()+pGround->get_hr()*pGround->getEnvironmentalTemperature()-pGround->get_Y())/(pGround->get_hc()+pGround->getG1()+pGround->get_hr()+pGround->get_X());
    // evaluation of the temperature
    float theta1 = (pGround->getC1()*pGround->getLayerTemperature(0) + dt*beta1)/(pGround->getC1() - dt*alpha1);

    // saves the results in the ground
    pGround->setLayerTemperature(0,theta1);

    // saves the outside surface temperature in the surface corresponding to the ground
    //float surfaceTemperature = (theta1*pGround->getG1() + Tout*pGround->get_hc() + (1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance() + pGround->get_hr()*pGround->getEnvironmentalTemperature())/(pGround->get_hc()+pGround->getG1()+pGround->get_hr());
    float surfaceTemperature = (theta1*pGround->getG1() + Tout*pGround->get_hc() + (1.f-pGround->getShortWaveReflectance())*pGround->getShortWaveIrradiance() + pGround->get_hr()*pGround->getEnvironmentalTemperature()-pGround->get_Y())/(pGround->get_hc()+pGround->getG1()+pGround->get_hr()+pGround->get_X());
    pGround->setTemperature(surfaceTemperature);

    // saves the water consumed
    pGround->setWaterEvapotranspiration((pGround->get_Y() +  pGround->get_X()*surfaceTemperature)/694.5); // 694.5 is the latent heat of water, considered as constant with water temperature

    // for debugging purposes, writes all temperatures in the soil
    #ifdef DEBUG
    ostringstream saveTemp;
    saveTemp << day << "\t" << hour << "\t" << pGround->getId() << "\t" << surfaceTemperature;
    saveTemp << "\t" << theta1;
    saveTemp << "\t" << Tground << endl;
    save(string("groundTemp.dat"),saveTemp,false);
    #endif

    return;

}

// JHK - thermal model for the Trees by SC
void Model::ThermalStepTree(Tree *pTree, Climate* pClimate, unsigned int day, unsigned int hour) {

    // outside air temperature
    float Ta = pClimate->getToutCelsius(day,hour);
    float Cp_a = 29.3f; // 29.3 is the Cp of air in J/(mol K)

    // computes averages on the top layer
    float totalArea = 0.f, SWa = 0.f, LWa = 0.f, LWe = 0.f, g_r = 0.f;
    for (size_t i=0;i<pTree->getLeaves()->size();++i) {
        totalArea+=pTree->getLeaves()->at(i)->getArea();
        // SW absorbed, note: 0.2 is the Transmittance of the leaves by default from Oke
        SWa+=(1.f-pTree->getLeaves()->at(i)->getShortWaveReflectance()-0.2f)*pTree->getLeaves()->at(i)->getShortWaveIrradiance()*pTree->getLeaves()->at(i)->getArea();
        // LW absorbed
        LWa+=pTree->getLeaves()->at(i)->getLongWaveAbsorbed()*pTree->getLeaves()->at(i)->getArea();
        // LW emitted at the air temperature
        LWe+=(pTree->getLeaves()->at(i)->getLongWaveEmissivity()*5.670373e-8*pow(Ta+273.15,4))*pTree->getLeaves()->at(i)->getArea();
        // LW radiative conductance
        g_r+=(4.*1.*5.670373e-8*pow(Ta+273.15,3)/Cp_a)*pTree->getLeaves()->at(i)->getArea();
    }
    SWa /= totalArea; // average absorbed shortwave irradiance on the top layer
    LWa /= totalArea; // average absorbed longwave irradiance on the top layer
    LWe /= totalArea; // average emitted longwave irradiance on the top layer
    g_r /= totalArea; // average radiative conductance ot the top layer

    // compute the conductances
    float g_Hr = 1.4*0.135*sqrt(max(pClimate->getWindSpeed(day,hour),0.01f)/(0.72*pTree->getLeafWidth())) + g_r;
    float g_va = 1.4*0.147*sqrt(max(pClimate->getWindSpeed(day,hour),0.01f)/(0.72*pTree->getLeafWidth()));
    float g_vs_ab = ((pClimate->getIdh(day,hour)>0.f)||(pClimate->getRelativeHumidity(day,hour)>0.8f))?0.3f:0.01f;
    float g_vs_ad = g_vs_ab; // equal conductances
    float g_v = (0.5*g_vs_ab*g_va)/(g_vs_ab+g_va)+(0.5*g_vs_ad*g_va)/(g_vs_ad+g_va);
    // gets the apparent psychometric constant (°C⁻¹)
    float gammaStar = 6.67e-4 * g_Hr / g_v;

    // compute the slope of saturation mole fraction function s (°C⁻¹)
    float s = pClimate->getSaturatedVapourPressureDerivative(Ta)/(pClimate->getPatm(day,hour)/1000.);

    // compute the vapour deficit
    float D = pClimate->getSaturatedVapourPressure(pClimate->getToutCelsius(day,hour))-pClimate->getSaturatedVapourPressure(pClimate->getTd(day,hour))+6.66e-4*(pClimate->getPatm(day,hour)/1000.)*(Ta-pClimate->getTd(day,hour));

    // compute the surface temperature
    float surfaceTemperature = Ta + gammaStar/(s+gammaStar)*((SWa+LWa-LWe)/(g_Hr*Cp_a)-D/(pClimate->getPatm(day,hour)/1000.*gammaStar));

    // save the temperature in all the layers, sublayers and trunc
    for (size_t i=0;i<pTree->getnSurfaces();++i) pTree->getSurface(i)->setTemperature(surfaceTemperature);

    #ifdef DEBUG
    ostringstream debugFile;
    if (day==1 && hour==1)
        debugFile << "day\thour\tid\tSWa\tTenv\tLWa\tLWe\tg_r\tg_Hr\tg_va\tg_vs_ab\tg_vs_ad\tg_v\tgammaStar\ts\tD\tsurfaceTemperature" << endl;
    debugFile << day << "\t" << hour << "\t" << pTree->getId() << "\t"
              << SWa << "\t" << pTree->getLeaves()->at(0)->getEnvironmentalTemperature() << "\t" << LWa << "\t" << LWe << "\t" << g_r << "\t" << g_Hr << "\t" << g_va << "\t" << g_vs_ab << "\t"
              << g_vs_ad << "\t" << g_v << "\t" << gammaStar << "\t" << s << "\t" << D << "\t" << surfaceTemperature << endl;
    save(string("tree.dat"),debugFile,false);
    #endif

    return;

}

void Model::ThermalExplicitStability(Building *pBuilding) {

    // liens
    int NP = pBuilding->getnNodes();

    // initialise the Atilde matrix Atilde = C^1 * G
    double Atilde[NP*NP];

    // gets the FixedMatrixComponents for the conductances
    for (int i=0;i<NP;i++) {
        for (int j=0;j<NP;j++) {
            // G1 contains the (timely) fixed conductances of the building and C the conductances (here C^-1 * K)
            Atilde[i+NP*j] = pBuilding->getG1(i,j) / pBuilding->getC(i,i);
        }
    }

    // output vectors
    double d[NP];

    // compute the eigenvalues using LAPACK
    eigenvalues_A(Atilde,d,NP);

    // shows the eigenvalues
    for (int i=0;i<NP;i++) {
        cout << "eigenvalue " << i << " = " << d[i] << "\ttime = " << (-2./d[i]) << " s" << endl;
    }

    return;
}

void Model::ThermalStepExplicitTemperature(Building *pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {

    // outside air temperature
    float Tout = pClimate->getToutCelsius(day,hour);
    float Tground = pClimate->getTgroundCelsius(day,hour);
    //cerr << "Tout: " << Tout << "\tTground: " << Tground << endl;

    // liens
    unsigned int NP = pBuilding->getnNodes();

    // propre au problËme
    double G[NP][NP];
    for (unsigned int i=0;i<NP;++i) { for (unsigned int j=0;j<NP;++j) G[i][j]=0.; }
    double b[NP], Tc[NP], T[NP];

    // propre ‡ la rÈsolution explicite
    double Texpl[NP];

    // loop on the smaller time steps for temperature determination
    for (unsigned int step2=0; step2 < dt/dt2; ++step2) {

        // loop on the zones
        for (unsigned int i=0;i<pBuilding->getnZones();++i) {

            // the occupants presence -> heat released
            /// commented: if needed again, this must be uncommented and adapted
            //Lr = pBuilding->getZone(i)->getOccupantsSensibleHeatRadiative()  * pBuilding->getZone(i)->getOccupants()->getOccupantsFraction(day,hour,step2);
            //Lc = pBuilding->getZone(i)->getOccupantsSensibleHeatConvective() * pBuilding->getZone(i)->getOccupants()->getOccupantsFraction(day,hour,step2);
            //pBuilding->getZone(i)->setRadiativeInternalHeatGains(Lr);
            //pBuilding->getZone(i)->setConvectiveInternalHeatGains(Lc);

            // *** the blind procedure *** //
            for (unsigned int wallIndex=0; wallIndex<pBuilding->getZone(i)->getnWalls(); ++wallIndex) {
                // loop on the external walls
                Wall* thisWall = pBuilding->getZone(i)->getWall(wallIndex);
                // determination of the blinds state which changes every step2
                Model::lowerShadingAction(pBuilding->getZone(i), thisWall, pClimate, day, hour, step2);
                /// TODO: see what we can do with the upperShadingAction
                //Model::upperShadingAction(pBuilding->getZone(i), thisWall, pClimate, day, hour, step2);
                //cerr << "LowerShadingState: " << pBuilding->getZone(i)->getLowerShadingState() << "\tVdotVent: " << pBuilding->getZone(i)->getVdotVent() << endl;
            }
            float Qs = pBuilding->getZone(i)->getQs(); // the heat/cold provided by the system

            //cerr << "Zone: " << i << "\tQsun1: " << pBuilding->getZone(i)->getQsun1() << "\tQsun2: " << pBuilding->getZone(i)->getQsun2() << "\tQs: " << Qs << endl;

            // preparation of the vector of temperatures
            for (unsigned int j=0;j<pBuilding->getZonenNodes(i);j++) {
                Tc[Thermal_getMatrixPosition(pBuilding,i)+j] = pBuilding->getZone(i)->getTExpl(j)*pBuilding->getZoneC(i,j);
                T[Thermal_getMatrixPosition(pBuilding,i)+j]  = pBuilding->getZone(i)->getTExpl(j);
            }

            // JK - 10.04.2013: removed the setKe as we are keeping the same from the ThermalStepImplicit (vector format

            // *** stochastic determination of the window state *** //
            string windowModel = pBuilding->getZone(i)->getOccupants()->getStochasticWindowParameters()->getWindowModel();
            if      (windowModel == "Markov")
                    Model::windowAction_Markov(pBuilding->getZone(i), pClimate, day, hour, step2);
            else if (windowModel == "Hybrid")
                    Model::windowAction_Hybrid(pBuilding->getZone(i), pClimate, day, hour, step2);
            else if (windowModel == "Bernoulli")
                    Model::windowAction_Bernoulli(pBuilding->getZone(i), pClimate, day, hour, step2);
            else if (windowModel == "Humphreys")
                    Model::windowAction_Humphreys(pBuilding->getZone(i), pClimate, day, hour, step2);
            else throw(string("No valid Stochastic Window Model, candidates are: Markov, Hybrid, Bernoulli, Humphreys"));

            // *** according to the window state defined previously, compute the ventilation flow rate *** //
            if ( pBuilding->getZone(i)->getWindowState() > 0.f ) {
                pBuilding->getZone(i)->setVdotVent( Model::ventilationFlowRate(pBuilding->getZone(i)->getTaExpl()+273.15f, pClimate->getToutCelsius(day,hour)+273.15f, pClimate->getWindSpeed(day,hour), pBuilding->getZone(i)->getSwiO()) );
            }
            else pBuilding->getZone(i)->setVdotVent(0.f);

            // *** determination of the lights state and electric consumption *** //
            // DP : commented for now, create trouble when DayLight is not computed (see doDayLightSim parameter of XmlScene::simulateTimeStep)
            //Model::lightAction_Lightswitch2002(pBuilding->getZone(i), day, hour, step2);
            //pBuilding->addElectricConsumption(Model::lightsElectricConsumption(pBuilding->getZone(i)));

            // *** calculation of UAm which depends on VdotVent *** -> in getUA() /

            // source terms b
            for (unsigned int j=0;j<pBuilding->getZone(i)->getnNodes();++j) {
                 b[Thermal_getMatrixPosition(pBuilding,i)+j] = pBuilding->getZone(i)->getSourceTerm(j,Tout,Tground,Qs);
            }

        }

        // initialise calculation matrix
        for (unsigned int i=0; i<pBuilding->getnZones(); ++i) {
            for (unsigned int j=0; j<pBuilding->getZone(i)->getnNodes(); ++j) {
                for (unsigned int k=0; k<pBuilding->getZonenNodes(i); ++k) {
                    // put the values in the G matrix
                    G[Thermal_getMatrixPosition(pBuilding,i)+j][Thermal_getMatrixPosition(pBuilding,i)+k] = pBuilding->getZone(i)->getVariableMatrixElement(j,k);
                }
            }
        }

        for (unsigned int i=0;i<NP;i++) {
            for (unsigned int j=0;j<NP;j++) {
                // G1 contains the (timely) fixed conductances of the building, that are added to the variable values initialised in G
                G[i][j]+=pBuilding->getG1(i,j); // previously G[i][j]=G1[i][j]+G2[i][j]
            }
        }

        // EXPLICIT Scheme
        for (unsigned int i=0;i<NP;++i) {
          Texpl[i]=Tc[i]+dt2*b[i];
              for (unsigned int j=0;j<NP;++j) {
                Texpl[i]+=dt2*G[i][j]*T[j];
              }
        }

        // save the output values
        for (unsigned int i=0; i<pBuilding->getnZones(); ++i) {
            for (unsigned int j=0; j<pBuilding->getZonenNodes(i); ++j) {
                Texpl[Thermal_getMatrixPosition(pBuilding,i)+j]/=pBuilding->getZoneC(i,j);
                pBuilding->getZone(i)->setTExpl(j,Texpl[Thermal_getMatrixPosition(pBuilding,i)+j]);
                //cerr << "Step: " << 24*(day-1)+hour << "\tStep2: " << step2 << "\tTemperature: " << Texpl[Thermal_getMatrixPosition(pBuilding,i)+j] << "\tsource term: " << b[Thermal_getMatrixPosition(pBuilding,i)+j] << endl;
            }
        }

    } // ends the loop on the internal time steps @ 5 min

    // saves the final temperature (after 1 hour)
    for (unsigned int i=0; i<pBuilding->getnZones(); ++i) {
        // saves the results in the zone (Ta,Tw, eventually Tw2)
        for (unsigned int j=0; j<pBuilding->getZonenNodes(i); ++j) {
            pBuilding->getZone(i)->setT(j,Texpl[Thermal_getMatrixPosition(pBuilding,i)+j]);
        }
        // saves the outside surface temperature in the surface corresponding to the zone (which is calculated from Tw - wall and Ta - roofs)
        pBuilding->getZone(i)->setTos(Tout);
    }

    return;

}

int Model::Thermal_getMatrixPosition(Building *pBuilding, int zoneNumber) {
    int sum = 0;
    for (int j=0; j<zoneNumber; j++) {
        sum += pBuilding->getZonenNodes(j);
    }
    return sum;
}

int Model::Thermal_getSubMatrixPosition(Building *pBuilding, int zoneNumber) {

  int sum = 0;
  for (int j=0; j<zoneNumber; j++) {
        sum += pBuilding->getZonenNodes(j)-1;
  }
  return sum;

}

double Model::Thermal_KiFixed(const double& surfaceWalls) {
    return 3.0*surfaceWalls;
}

double Model::Thermal_KeClarke(const double& WindSpeed, const double& WindDir/*, const double& Tout*/, const vector<float> &azimuth, const vector<float> &surface, const double& totalSurfaceWalls) {

    double hc=0., Ke=0., verticalSurface=0.;
    vector<double> dr(azimuth.size(), 0.), W(azimuth.size(), 0.);

    for (unsigned int i=0;i<azimuth.size();i++) {
        dr[i]=abs(azimuth[i]-WindDir);
        if ( dr[i] > 180 ) dr[i]=360-dr[i];

        if (dr[i] < 10) {
            W[i]=0.5*WindSpeed;
            if ( W[i] > 0.5 ) W[i]=0.5;
            if ( WindSpeed > 2 ) W[i]=0.25*WindSpeed;
        }

        else if (dr[i] > 90) W[i]=WindSpeed*0.25*abs(sin(dr[i]/180*M_PI));

        else W[i]=WindSpeed*sin(dr[i]/180*M_PI);

        hc = 2.8 + 3.0*W[i];
        Ke += hc*surface[i];

        verticalSurface+=surface[i];
    }

    hc = 2.8 + 3.0*WindSpeed;
    Ke += (totalSurfaceWalls - verticalSurface) * hc; // horizontalSurface * hc

    return Ke;
}

float Model::Thermal_hcClarke(const float& WindSpeed, const float& WindDir, const float& azimuth) {

    // McAdams in ESP-r, see Applied Thermal Engineering 56 (2013) pp.134-151
    float W=0.f;
    float dr=abs(azimuth-WindDir);
    if (dr > 180.f) dr=360.f-dr; // between 0 and 180∞

    if (dr < 10.f) { // Windward surface with 0∞ < dr <= 10∞
        W=0.5f*WindSpeed;
        if ( W > 0.5f ) W=0.5f;
        if ( WindSpeed > 2.f ) W=0.25f*WindSpeed;
    }
    else if (dr > 90.f) W=WindSpeed*0.25f*abs(sin(dr/180.f*M_PI)); // Leeward surface (90∞ < dr <= 180∞)
    else W=WindSpeed*sin(dr/180.f*M_PI); // Windward surface with 10∞ < dr <= 90∞

    return (2.8f + 3.0f*W);
}

float Model::Thermal_hcLiuHarris(const float& WindSpeed, const float& WindDir, const float& azimuth) {

    // Liu and Harris, see Applied Thermal Engineering 56 (2013) pp.134-151
    float dr=abs(azimuth-WindDir);
    if (dr > 180.f) dr=360.f-dr; // between 0 and 180∞

    if (dr <= 90.f) return 1.53*WindSpeed + 1.43; // Windward surface with 0∞ < dr <= 90∞
    else return 0.90*WindSpeed + 3.28; // Leeward surface (90∞ < dr <= 180∞)

}

float Model::Thermal_KeClarke(const float& WindSpeed, const float& WindDir, const vector<float> &azimuth, const vector<float> &surface) {

    float Ke=0.f;
    // make the sum on all surfaces
    for (unsigned int i=0;i<azimuth.size();i++) {
        Ke += Thermal_hcClarke(WindSpeed,WindDir,azimuth[i])*surface[i];
    }
    return Ke;
}

float Model::Thermal_hcCli2(Climate* pClimate, unsigned int surfaceId, unsigned int day, unsigned int hour) {

    return (pClimate->getKeCoeff1(surfaceId,day,hour)*pow(pClimate->getWindSpeed(day,hour),pClimate->getKeCoeff2(surfaceId,day,hour)) + pClimate->getKeCoeff3(surfaceId,day,hour));

}

float Model::Thermal_KeWalls(Climate* pClimate, Zone* pZone, unsigned int day, unsigned int hour) {

    // sets the hc for the Surfaces (needed by the PVT)
    for (size_t i=0; i<pZone->getnSurfaces(); ++i) {
        // test if we have the cli2 data
        if (pClimate->isEmptyCli2())
            pZone->getSurface(i)->set_hc(Thermal_hcClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour),pZone->getSurface(i)->getAzimuth()));
        else
            pZone->getSurface(i)->set_hc(Thermal_hcCli2(pClimate,pZone->getSurface(i)->getId(),day,hour));
    }

    float Ke = 0.f;
    // loop on the walls to extract the Ke
    for (size_t i=0; i<pZone->getnWalls(); ++i) {
        // test if we have the cli2 data
        if (pClimate->isEmptyCli2())
            pZone->getWall(i)->set_hc(Thermal_hcClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour),pZone->getWall(i)->getAzimuth()));
        else
            pZone->getWall(i)->set_hc(Thermal_hcCli2(pClimate,pZone->getWall(i)->getId(),day,hour));
        // sums up for Ke
        Ke += pZone->getWall(i)->get_hc()*pZone->getWall(i)->getWallArea();
    }
    return Ke;
}

float Model::Thermal_HrWalls(Zone* pZone) {

    float Hr = 0.f;
    // loop on the walls to extract the Ke
    for (size_t i=0; i<pZone->getnWalls(); ++i) {
        Hr += pZone->getWall(i)->get_hr()*pZone->getWall(i)->getWallArea();
    }
    return Hr;
}

float Model::Thermal_Uprime(Surface* pSurface, float hc) {
    // gets the modified U-value thanks to the variable hc
    return 1.f/( 1.f/pSurface->getGlazingUvalue() - 1.f/25.f + 1.f/hc + 1.f/pSurface->get_hr() - 1.f/8.f + 1.f/3.f); // hc_int = 3 from CIBSE guide
}

float Model::Thermal_Kwindows(Climate* pClimate, Zone* pZone, unsigned int day, unsigned int hour) {

    float Kwindow = 0.f;
    // loop on the walls to extract the Kwindow
    for (size_t i=0; i<pZone->getnWalls(); ++i) {
        // test if we have the cli2 data
        if (pClimate->isEmptyCli2()) {
            Kwindow += Thermal_Uprime(pZone->getWall(i), Thermal_hcClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour),pZone->getWall(i)->getAzimuth()))
                       *pZone->getWall(i)->getGlazingArea();
        }
        else {
            Kwindow += Thermal_Uprime(pZone->getWall(i), Thermal_hcCli2(pClimate,pZone->getWall(i)->getId(),day,hour))
                       *pZone->getWall(i)->getGlazingArea();
        }
    }
    // loop on the roofs to extract the Kwindow
    for (size_t i=0; i<pZone->getnRoofs(); ++i) {
        // test if we have the cli2 data
        if (pClimate->isEmptyCli2())  {
            Kwindow += Thermal_Uprime(pZone->getRoof(i), Thermal_hcClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour),pZone->getRoof(i)->getAzimuth()))
                       *pZone->getRoof(i)->getGlazingArea();
        }
        else {
            Kwindow += Thermal_Uprime(pZone->getRoof(i), Thermal_hcCli2(pClimate,pZone->getRoof(i)->getId(),day,hour))
                       *pZone->getRoof(i)->getGlazingArea();
        }
    }
    return Kwindow;
}

float Model::Thermal_Kroofs(Climate* pClimate, Zone* pZone, unsigned int day, unsigned int hour) {

    float Kroof = 0.f;
    // loop on the roofs to extract the Ke
    for (size_t i=0; i<pZone->getnRoofs(); ++i) {
        // test if we have the cli2 data, and save the hc in the surface itself
        if (pClimate->isEmptyCli2())
            pZone->getRoof(i)->set_hc(Thermal_hcClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour),pZone->getRoof(i)->getAzimuth()));
        else
            pZone->getRoof(i)->set_hc(Thermal_hcCli2(pClimate,pZone->getRoof(i)->getId(),day,hour));
        // adds the roof in the Kroof (global)
        Kroof += pZone->getRoof(i)->getKappa()*pZone->getRoof(i)->get_hc()*pZone->getRoof(i)->getRoofArea();
    }
    return Kroof;
}

void Model::Thermal_Kgrounds(Climate* pClimate, Ground* pGround, unsigned int day, unsigned int hour) {

    // test if we have the cli2 data, and save the hc in the surface itself
    if (pClimate->isEmptyCli2())
        pGround->set_hc(Thermal_hcClarke(pClimate->getWindSpeed(day,hour),pClimate->getWindDirection(day,hour), pGround->getAzimuth()));
    else
        pGround->set_hc(Thermal_hcCli2(pClimate, pGround->getId(),day,hour));

}

    double Model::Thermal_Ke(const double& WindSpeed, const double& WindDir, const double& Tout, const vector<double> &azimuth, const vector<double> &surface) {
        double hc=0, Ke=0;
        vector<double> dr(azimuth.size(), 0.), W(azimuth.size(), 0.);
        for (unsigned int i=0;i<azimuth.size();i++) {
            dr[i]=azimuth[i]-WindDir;
            if (abs(dr[i])>180) dr[i]=360-abs(dr[i]);
            if (abs(dr[i])<90) {
                if (WindSpeed > 2) W[i]=0.25*WindSpeed;
                else W[i]=0.5; }
            else W[i]=0.3+0.05*WindSpeed;

            W[i]=294.26*W[i]/(273.16+Tout);
            if ( W[i] < 4.88 ) hc=5.678*(0.99+0.21*W[i]/0.3048);
            else hc=5.678*0.5*pow(W[i]/0.3048, 0.78);

            Ke += hc*surface[i];
        }
        return Ke;
    }


void Model::ThermalAllAvailable(Building *pBuilding) {

    for (unsigned int i=0;i<pBuilding->getnZones();i++) {

        pBuilding->getZone(i)->setQs(pBuilding->getZone(i)->getHeating() + pBuilding->getZone(i)->getCooling());

    }

    pBuilding->setMachinePower(0.0);
    pBuilding->setFuelConsumption(0.0);
    pBuilding->addElectricConsumption(0.f);
    /// TODO: check here the consistency with addElectricConsumption and SolarPVProduction
    pBuilding->setSolarThermalProduction(0.f);

}

void Model::HVAC_Needs(Building *pBuilding,Climate* pClimate,unsigned int day,unsigned int hour) {

    // gets the needed information for this timestep in the climate file and in the building
    double t1 =  pClimate->getToutCelsius(day,hour);
    double hr =  pClimate->getRelativeHumidity(day,hour);
    double patm = pClimate->getPatm(day,hour);

    double t3heat = pBuilding->getTmaxSupply();
    double t3cool = pBuilding->getTminSupply();
    double deltat = pBuilding->getDeltaT();

    bool evaporativeCooling = pBuilding->getEvaporativeCooling();

    for (unsigned int i=0;i<pBuilding->getnZones();++i){

        double t5 = pBuilding->getZone(i)->getTaForeseen();
        bool mControl = true;
        if (t5 > pBuilding->getZone(i)->getTmax() || t5 < pBuilding->getZone(i)->getTmin()) mControl = false;

        double t5prev = pBuilding->getZone(i)->getTa();

        double qs = pBuilding->getZone(i)->getHeating() + pBuilding->getZone(i)->getCooling(); //Sensible load

        double np = pBuilding->getZone(i)->getOccupantsCount(); // occupants number * occupants fraction
        double ql = np*pBuilding->getZone(i)->getOccupantsLatentHeat(); // Latent load

        // moisture content w1
        double w1 = HVAC_moistureContent(t1, hr, patm);

        // set point temperature
        //if (t5 > pBuilding->getTmax()) t5 = pBuilding->getTmax();
        //else if (t5 < pBuilding->getTmin()) t5 = pBuilding->getTmin();

        // off-coil temperature
        double t2 = pBuilding->getCoilEfficiency()*(t5prev-t1)+t1;
        double w2 = min(w1, HVAC_moistureContent(t2, 1.0, patm)); // condensation dans le coil
        // hydroThermalWheel
        if (pBuilding->getCoilHydroThermalWheelEfficiency() > 0. && pBuilding->getZone(i)->getnMoistureContent() != 0)
            w2 = min(pBuilding->getCoilHydroThermalWheelEfficiency()*(pBuilding->getZone(i)->getMoistureContent()-w1)+w1, HVAC_moistureContent(t2, 1.0, patm));

        // temperature already in the OK range at the exit of the heat exchanger
        //if ( qs == 0. && ( t2 < pBuilding->getTmax() && t2 > pBuilding->getTmin() ) ) t5 = t2;

        // supply air temperature
        double t3;
        if ( qs > 0.0 ) t3 = t3heat;
        else if ( qs < 0.0 ) t3 = t3cool;
        else t3 = t5;

        // flow rate - minimal fresh air requirements
        float m1dot = HVAC_massFlowRate(t2, patm, w2, qs, t3, t5, np);
        if (m1dot == 0.f) {
            // system switched off
            pBuilding->getZone(i)->setHVACHeat(0.);
            pBuilding->getZone(i)->setHVACCool(0.);
            pBuilding->getZone(i)->setHVACReheat(0.);
            pBuilding->getZone(i)->setHVACHumidification(0.);
            pBuilding->getZone(i)->setHVACEvaporation(0.);
            pBuilding->getZone(i)->setHVACMassFlowRate(0.);
            return;
        }

        // moisture content w5
        double w5prime = HVAC_totalMoistureContent(w2, ql, m1dot);

        // moisture control w3 - temperature control t3
        double deltatHVAC = 0.0;
        double evl = 0.0;
        double w3;
        if (mControl) w3 = HVAC_moistureControl(pBuilding->getZone(i)->getTmax(), t2, t3, deltatHVAC, HVAC_temperatureChange(deltat, t2, t3), t5, w2, w5prime, patm, evl, evaporativeCooling);
        else w3 = w2;

        double tSupplyHVAC = t3 + deltatHVAC + HVAC_temperatureChange(deltat, t2, t3);
        double ws = std::min(w3, HVAC_moistureContent(tSupplyHVAC, 1.0, patm));
        w3 = ws;

        // correction dans m1dot pour tenir compte d'un Èventuel changement d'humiditÈ du systËme
        m1dot = HVAC_massFlowRate(t3,patm,ws,qs,t3,t5,np);

        // recalculate final state (5) for VERIFICATION
        t5 = tSupplyHVAC - HVAC_temperatureChange(deltat, t2, t3) - qs/(m1dot*1000.0*(1.007+ws*1.84));
        //double w5 = HVAC_totalMoistureContent(ws, ql, m1dot);

        // heating & cooling routines
        double hel = 0.0;
        double hul = 0.0;
        double col = 0.0;
        double rel = 0.0;
        HVAC_heat(t2, w2, tSupplyHVAC, ws, /*patm,*/ hel, hul);
        HVAC_cool(t2, w2, tSupplyHVAC, ws, /*t5, w5,*/ patm, col, rel, hul);

        // saving the results into the zone
        pBuilding->getZone(i)->setHVACHeat(hel*m1dot);
        pBuilding->getZone(i)->setHVACCool(col*m1dot);
        pBuilding->getZone(i)->setHVACReheat(rel*m1dot);
        pBuilding->getZone(i)->setHVACHumidification(hul*m1dot);
        pBuilding->getZone(i)->setHVACEvaporation(evl*m1dot);
        pBuilding->getZone(i)->setHVACMassFlowRate(m1dot);

    }
}

void Model::HVAC_Control(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {

    double heatingNeeds = 0.;
    double coolingNeeds = 0.;
    // computation of the total demand in heat and cold for the building
    for (unsigned int i=0; i<pBuilding->getnZones(); i++) heatingNeeds += pBuilding->getZone(i)->getHVACHeat();
    for (unsigned int i=0; i<pBuilding->getnZones(); i++) coolingNeeds += pBuilding->getZone(i)->getHVACCool();

    // heating with heatTank
    if (heatingNeeds > 0.) {

        double HS_Pup = 0.; // this will be provided by the solar panel
        float Tamb = 18.f; // 18∞C in the cave

        // lets start with the heat stock, we compute the temperature after the time step
        double HS_T1 = pBuilding->getHeatStock()->temperature(dt, 0.0, 0.0, -heatingNeeds, pBuilding->getHeatStockTemperature(), 0.0, Tamb);

        // compute the HS heating needs
        double HS_needs = 0.;
        if ( HS_T1 < pBuilding->getHeatStock()->getTmin() )
            HS_needs = pBuilding->getHeatStock()->power(dt, pBuilding->getHeatStock()->getTmax(), 0.0, HS_Pup-heatingNeeds, pBuilding->getHeatStockTemperature(), 0.0, Tamb);

        // note: the domestic water consumption is not taken into account now
        //float DHW_needs = 0.f;
        //float DHW_Tmax = 0.f;

        // the heatPumpSrcTemp depends on the type of heatpump (air or ground vertical/horizontal)
        double heatPumpSrcTemp = pClimate->getToutCelsius(day,hour); // air temperature
        float z0=0.f,z1=0.f,alpha=0.f;
        if ( pBuilding->getHeatingUnit() != NULL ) {
            if ( pBuilding->getHeatingUnit()->getGround(z0,z1,alpha) ) {
                heatPumpSrcTemp = pClimate->getTgroundCelsius(day,hour,z0,alpha,z1);
                //cerr << "Heat Pump Src Temp: " << heatPumpSrcTemp << endl;
            }
        }

        double HS_Pp = 0.;//, DHW_Pp = 0.;
        // see how the energy unit can provide, HS_Pp is the max heating power that can give the MACHINE (if working)
        if (HS_needs > 0. && pBuilding->getHeatingUnit()!=NULL && pBuilding->getHeatingUnit()->isWorking(day)) { /// TODO: use here a more complex version for also taking into account domestic water needs
            pBuilding->getHeatingUnit()->setThermalPowerNeeded(HS_needs); // Cognet: Added this.
            HS_Pp = pBuilding->getHeatingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//            HS_Pp = pBuilding->getHeatingUnit()->getThermalPower(HS_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
        }

        // total energy available in the tank with HS_Pp
        double heatingAvailable = max(-pBuilding->getHeatStock()->power(dt, pBuilding->getHeatStock()->getTmin(),
                                                                        0., HS_Pup+HS_Pp, pBuilding->getHeatStockTemperature(),
                                                                        0., Tamb),0.);

        // the energy that can be used, the minimum between the needs and the available
        double heatingUsed = min(heatingNeeds, heatingAvailable);
        // take that heatingUsed from the stock with HS_Pp provided, and calculate new temperature
        HS_T1 = pBuilding->getHeatStock()->temperature(dt, 0., HS_Pp, -heatingUsed, pBuilding->getHeatStockTemperature(), 0., Tamb);
        // save the new stock temperature
        pBuilding->setHeatStockTemperature(HS_T1);
        pBuilding->setColdStockTemperature(pBuilding->getColdStock()->temperature(dt, 0., 0., 0., pBuilding->getColdStockTemperature(), 0., Tamb));

        // Energy consumption in primary resource to provide HS_Pp
        pBuilding->setMachinePower(HS_Pp);
        if (pBuilding->getHeatingUnit() != NULL) {
            pBuilding->setFuelConsumption(pBuilding->getHeatingUnit()->getFuelConsumption(double(dt),HS_Pp,heatPumpSrcTemp));
            pBuilding->addElectricConsumption(pBuilding->getHeatingUnit()->getElectricConsumption(double(dt),HS_Pp,heatPumpSrcTemp));
        }
        else {
            pBuilding->setFuelConsumption(0.);
            pBuilding->addElectricConsumption(0.);
        }

        // repartition of the energy in the zones
        for (unsigned int i=0; i<pBuilding->getnZones(); i++) {

            /// TODO: depending on the type of systems present, improve this bit
            // for the HEATING
            if ( heatingUsed > pBuilding->getZone(i)->getHVACHeat() ) { // demand satisfied
                pBuilding->getZone(i)->setHVACHeatAvailable(pBuilding->getZone(i)->getHVACHeat());
                heatingUsed -= pBuilding->getZone(i)->getHVACHeat();
            }
            else {
                pBuilding->getZone(i)->setHVACHeatAvailable(heatingUsed);
                heatingUsed -= heatingUsed;
            }

            // for the connex elements
            pBuilding->getZone(i)->setHVACCoolAvailable(0.);
            pBuilding->getZone(i)->setHVACReheatAvailable(0.);
            pBuilding->getZone(i)->setHVACHumidificationAvailable(0.);
            pBuilding->getZone(i)->setHVACEvaporationAvailable(0.);

            // for the mass flow rate (according to the fan power) - TO BE CHANGED
            pBuilding->getZone(i)->setHVACMassFlowRateAvailable(pBuilding->getZone(i)->getHVACMassFlowRate());

        }

        if (heatingUsed > 1.) throw(string("Error in the repartition of the heating energy - in the control model (HVAC): " + toString(heatingUsed)));
    }
    else if (coolingNeeds < 0.) { // cooling with coldTank

        double CS_Pup = 0.; // energy removed by the solar panel
        float Tamb = 18.f; // 18∞C in the cave
        double heatPumpSrcTemp = pClimate->getToutCelsius(day,hour);

        // calculation of the new stock temperature
        double CS_T1 = pBuilding->getColdStock()->temperature(dt, 0., 0., -coolingNeeds, pBuilding->getColdStockTemperature(), 0., Tamb);

        // calculation of the needs to reach the lowest tank temperature of upper limit overtaken
        double CS_needs = 0.;
        if ( CS_T1 >= pBuilding->getColdStock()->getTmax() )
            CS_needs = pBuilding->getColdStock()->power(dt, pBuilding->getColdStock()->getTmin(), 0., CS_Pup-coolingNeeds, pBuilding->getColdStockTemperature(), 0., Tamb);

        // calculation of the refreshing available
        double CS_Pp = 0.;
        if (CS_needs < 0. && pBuilding->getCoolingUnit()!=NULL && pBuilding->getCoolingUnit()->isWorking(day)) {
            pBuilding->getCoolingUnit()->setThermalPowerNeeded(CS_needs); // Cognet: Added this.
            CS_Pp = pBuilding->getCoolingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//            CS_Pp = pBuilding->getCoolingUnit()->getThermalPower(CS_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
        }

        // total energy available in the tank with CS_Pp provided
        double coolingAvailable = min(-pBuilding->getColdStock()->power(dt, pBuilding->getColdStock()->getTmax(),
                                                                     0., CS_Pup+CS_Pp, pBuilding->getColdStockTemperature(),
                                                                     0., Tamb),0.);

        double coolingUsed = max(coolingNeeds, coolingAvailable);
        // compute the new stock temperature with CS_Pp provided and the coolingUsed
        CS_T1 = pBuilding->getColdStock()->temperature(dt, 0., CS_Pp, -coolingUsed, pBuilding->getColdStockTemperature(), 0., Tamb);
        // save the new stock temperature
        pBuilding->setColdStockTemperature(CS_T1);
        pBuilding->setHeatStockTemperature(pBuilding->getHeatStock()->temperature(dt, 0., 0., 0., pBuilding->getHeatStockTemperature(), 0., Tamb));

        // Energy consumption to provide CS_Pp
        pBuilding->setMachinePower(CS_Pp);
        if (pBuilding->getCoolingUnit() != NULL) {
            pBuilding->setFuelConsumption(pBuilding->getCoolingUnit()->getFuelConsumption(double(dt),CS_Pp,heatPumpSrcTemp));
            pBuilding->addElectricConsumption(pBuilding->getCoolingUnit()->getElectricConsumption(double(dt),CS_Pp,heatPumpSrcTemp));
        }
        else {
            pBuilding->setFuelConsumption(0.);
            pBuilding->addElectricConsumption(0.);
        }

        // repartition of the energy in the zones
        for (unsigned int i=0; i<pBuilding->getnZones(); i++) {

            // for the COOLING
            if ( coolingUsed < pBuilding->getZone(i)->getHVACCool() ) { // demand satisfied
                pBuilding->getZone(i)->setHVACCoolAvailable(pBuilding->getZone(i)->getHVACCool());
                coolingUsed -= pBuilding->getZone(i)->getHVACCool();
            }
            else {
                pBuilding->getZone(i)->setHVACCoolAvailable(coolingUsed);
                coolingUsed -= coolingUsed;
            }

            // for the connex elements
            pBuilding->getZone(i)->setHVACHeatAvailable(0.);
            pBuilding->getZone(i)->setHVACReheatAvailable(0.);
            pBuilding->getZone(i)->setHVACHumidificationAvailable(0.);
            pBuilding->getZone(i)->setHVACEvaporationAvailable(0.);

            // for the mass flow rate (according to the fan power) - TO BE CHANGED
            pBuilding->getZone(i)->setHVACMassFlowRateAvailable(pBuilding->getZone(i)->getHVACMassFlowRate());

        }

        if (coolingUsed < -1.) throw(string("Error in the repartition of the cooling energy - in the control model (HVAC): " + toString(coolingUsed)));
    }
    else {

        float Tamb = 18.f; // 18∞C in the cave

        pBuilding->setMachinePower(0.);
        pBuilding->setFuelConsumption(0.);
        pBuilding->addElectricConsumption(0.);

        pBuilding->setHeatStockTemperature(pBuilding->getHeatStock()->temperature(dt, 0., 0., 0., pBuilding->getHeatStockTemperature(), 0., Tamb));
        pBuilding->setColdStockTemperature(pBuilding->getColdStock()->temperature(dt, 0., 0., 0., pBuilding->getColdStockTemperature(), 0., Tamb));

        // repartition of the energy in the zones
        for (unsigned int i=0; i<pBuilding->getnZones(); i++) {

            // for the heating
            pBuilding->getZone(i)->setHVACHeatAvailable(0.);

            // for the COOLING
            pBuilding->getZone(i)->setHVACCoolAvailable(0.);
            pBuilding->getZone(i)->setHVACReheatAvailable(0.);

            // for the connex elements
            pBuilding->getZone(i)->setHVACHumidificationAvailable(0.);
            pBuilding->getZone(i)->setHVACEvaporationAvailable(0.);

            // for the mass flow rate (according to the fan power) - TO BE CHANGED
            pBuilding->getZone(i)->setHVACMassFlowRateAvailable(pBuilding->getZone(i)->getHVACMassFlowRate());
        }
    }

    // compute the electricity produced by the PV
    for (unsigned int i=0; i<pBuilding->getnZones(); i++) {
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnRoofs(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
            pBuilding->addElectricConsumption(-float(dt)*pBuilding->getZone(i)->getRoof(j)->getWindElectricPower(pClimate->getWindSpeed(day,hour),pClimate->getAltitude()) ); // in joules
//            if ( pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tRoof: " << j << "\tIrradiation: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getRoof(j)->getArea() <<
//            "\tIrradiation Tot: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnWalls(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getWall(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
//            if ( pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tRoof: " << j << "\tIrradiation: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getRoof(j)->getArea() <<
//            "\tIrradiation Tot: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnSurfaces(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getSurface(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
//            if ( pBuilding->getZone(i)->getSurface(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tSurface: " << j << "\tIrradiation (W/m2): " << pBuilding->getZone(i)->getSurface(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getSurface(j)->getArea() <<
//            "\tIrradiation (W): " << pBuilding->getZone(i)->getSurface(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getSurface(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
    }

}

void Model::HVAC_Available(Building *pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {

    double t1 =  pClimate->getToutCelsius(day,hour);
    double hr =  pClimate->getRelativeHumidity(day,hour);
    double patm = pClimate->getPatm(day,hour);

    double deltat = pBuilding->getDeltaT();

    for(unsigned int i=0;i<pBuilding->getnZones();i++){

        double t5 = pBuilding->getZone(i)->getTaForeseen();
        double t5prev = pBuilding->getZone(i)->getTa();

        // set point temperature
        //if (t5 > pBuilding->getTmax()) t5 = pBuilding->getTmax();
        //else if (t5 < pBuilding->getTmin()) t5 = pBuilding->getTmin();

        double m1dot = pBuilding->getZone(i)->getHVACMassFlowRateAvailable();
        double hel = pBuilding->getZone(i)->getHVACHeatAvailable();
        double col = pBuilding->getZone(i)->getHVACCoolAvailable();
        double rel = pBuilding->getZone(i)->getHVACReheatAvailable();
        double hul = pBuilding->getZone(i)->getHVACHumidificationAvailable();
        double evl = pBuilding->getZone(i)->getHVACEvaporationAvailable();

        // moisture content w1
        double w1 = HVAC_moistureContent(t1, hr, patm);

        // off-coil temperature
        double t2 = pBuilding->getCoilEfficiency()*(t5prev-t1)+t1;
        double w2 = min(w1, HVAC_moistureContent(t2, 1.0, patm)); // condensation dans le coil
        // presence of the hydroThermalWheel
        if (pBuilding->getCoilHydroThermalWheelEfficiency() > 0. && pBuilding->getZone(i)->getnMoistureContent() != 0)
            w2 = min(pBuilding->getCoilHydroThermalWheelEfficiency()*(pBuilding->getZone(i)->getMoistureContent()-w1)+w1, HVAC_moistureContent(t2, 1.0, patm));

        // on part de t2,w2 pour arriver sur ts,ws
        // calcul de t2prime, w2prime
        double H2 = HVAC_enthalpyDryAir(t2) + w2*HVAC_enthalpyWaterVapour(t2); // enthalpy of Point 2
        double w2prime;
        if (m1dot > 0.) w2prime = w2 + evl/m1dot; else w2prime = w2; // check if mass flow rate exists
        double t2prime = (H2 + 0.026 - w2prime*2501.0) / (1.007 + w2prime*1.84);

        double ts, ws;

        if (hel > 0.) { // heating
            double Cpa = 1007. + w2prime*1840.;
            ts = hel/(m1dot*Cpa) + t2prime;
            double Cpw = 2.501e6 + 1840.*t2prime; // changed ts in t2prime, JK - 7/04/09
            ws = hul/(m1dot*Cpw) + w2prime;
        }
        else if (col < 0.) { // cooling
            if (hul > 0.) { // humidification
                double Cpa = 1007. + w2prime*1840.;
                ts = col/(m1dot*Cpa) + t2prime;
                double Cpw = 2.501e6 + 1840.*t2prime; // changed ts in t2prime, JK - 7/04/09
                ws = hul/(m1dot*Cpw) + w2prime;
            }
            else { // cooling, dehumidification and reheat
                double Cpa = 1007. + w2prime*1840.;
                if ( col/m1dot > Cpa*(HVAC_bulbTemperature(w2prime,100.,patm)-t2prime) ) { // cooling seulement
                    ts = col/(m1dot*Cpa) + t2prime ;
                    ws = w2prime;
                }
                else { // dehumidification
                   ts = HVAC_bulbTemperature(col/m1dot,t2prime,w2prime,patm);
                    ws = HVAC_moistureContent(ts,1.0,patm);
                    ts += (rel/m1dot)*(1007.+ws*1840.);
                }
            }
        }
        else { // ni heating, ni cooling
            ts = t2prime;
            ws = w2prime;
        }

        double np = pBuilding->getZone(i)->getOccupantsCount(); // number of people * presence fraction (or probability)
        double ql = np*pBuilding->getZone(i)->getOccupantsLatentHeat(); // Latent load

        double w5 = HVAC_totalMoistureContent(ws, ql, m1dot);
        // saves the moisture content in the zone
        pBuilding->getZone(i)->setMoistureContent(w5);

        double qs = m1dot*(1007.+ws*1840.)*(ts-deltat-t5);
        // saves the Qs given to the system
        pBuilding->getZone(i)->setQs(qs);

    }
}

void Model::noHVAC_Control(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {

    // double precision required in this method for a precise repartition of the energy
    double heatingNeeds = 0.;
    double coolingNeeds = 0.;

    float Tamb = 18.f; // 18°C in the cave
    double HS_Pp = 0., DHW_Pp = 0.;
    double HS_Pup = 0., CS_Pup = 0.; // this will be provided by the solar thermal panel (heating and cooling)

    // initialize heat stock temperature
    const double HS_T0 = pBuilding->getHeatStockTemperature();
    double DHW_T0 = numeric_limits<double>::signaling_NaN();
    if (pBuilding->getDHWHeatStock()!= NULL) DHW_T0 = pBuilding->getDHWStockT();
    const double CS_T0 = pBuilding->getColdStockTemperature();

    // initialize the machine power to zero
    pBuilding->setMachinePower(0.);
    pBuilding->setFuelConsumption(0.);

    // Heat stock temperature after one time step
    double HS_T1, DHW_T1=DHW_T0, CS_T1;

    // computation of the total demand in heat and cold for the building
    for (size_t i=0; i<pBuilding->getnZones(); ++i) {
        heatingNeeds += pBuilding->getZone(i)->getHeating();
        coolingNeeds += pBuilding->getZone(i)->getCooling();
    }

    // Solar Thermal gains for heating -> HS_Pup
    for (unsigned int i=0; i<pBuilding->getnZones(); i++) { // Cognet: This is computed for the hot storage (since use HS_T0), but then later under certain conditions, it is given to the DHW tank. This may cause errors.
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnRoofs(); ++j)
            HS_Pup += pBuilding->getZone(i)->getRoof(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour), pClimate->getWindSpeed(day,hour), HS_T0);
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnWalls(); ++j)
            HS_Pup += pBuilding->getZone(i)->getWall(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour), pClimate->getWindSpeed(day,hour), HS_T0);
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnSurfaces(); ++j)
            HS_Pup += pBuilding->getZone(i)->getSurface(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour), pClimate->getWindSpeed(day,hour), HS_T0);
    }
    pBuilding->setSolarThermalProduction(float(dt)*HS_Pup);

    // Space Heating Energy used
    double heatingUsed = 0.;

    // Domestic Hot Water (DHW) usage in the zones
    double VdotUsed = 0.; // in m³/s
    for (size_t i=0; i<pBuilding->getnZones(); ++i) VdotUsed += pBuilding->getZone(i)->getDHWConsumption(day,hour)/(1000.*3600.); // DHW consumption in l/h -> m³/s

    // heating with heatTank

    // the heatPumpSrcTemp depends on the type of heatpump (air or ground vertical/horizontal)
    double heatPumpSrcTemp = pClimate->getToutCelsius(day,hour); // air temperature
    float z0=0.f,z1=0.f,alpha=0.f;
    if ( pBuilding->getHeatingUnit() != NULL ) {
        if ( pBuilding->getHeatingUnit()->getGround(z0,z1,alpha) ) {
            heatPumpSrcTemp = pClimate->getTgroundCelsius(day,hour,z0,alpha,z1);
            //cerr << "Heat Pump Src Temp: " << heatPumpSrcTemp << endl;
        }
    }

    // lets start with the heat stock, we compute the temperature after the time step, no losses of fluid from the tank, solar energy provided, heating needs extracted
    HS_T1 = pBuilding->getHeatStock()->temperature(dt, 0., HS_Pup, -heatingNeeds, HS_T0, 0.0, Tamb);

    // compute the HS heating needs
    double HS_needs = 0., DHW_needs = 0.;
    double machineThermalPower = 0.;
    if (HS_T1 < pBuilding->getHeatStock()->getTmin() && heatingNeeds > 0.) {
        // if the temperature is below the minimum, compute the energy needed to reach the maximum with added solar power and extraction of heating needs
        HS_needs = pBuilding->getHeatStock()->power(dt, pBuilding->getHeatStock()->getTmax(), 0., HS_Pup-heatingNeeds, HS_T0, 0.0, Tamb);

        // if the DHWtank exists, compute its needs
        if (pBuilding->getDHWHeatStock() != NULL)
        {
            DHW_T1 = pBuilding->getDHWHeatStock()->temperature(dt, VdotUsed, 0., 0., DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
            if (DHW_T1 < pBuilding->getDHWHeatStock()->getTmin())
                DHW_needs = pBuilding->getDHWHeatStock()->power(dt, pBuilding->getDHWHeatStock()->getTmax(), VdotUsed, 0.0, DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
        }

        // see how the energy unit can provide, HS_Pp is the max heating power that can give the MACHINE hopefully equal to the needs
        if (pBuilding->getHeatingUnit()!=NULL) {
            if (pBuilding->getHeatingUnit()->isWorking(day)) {
                pBuilding->getHeatingUnit()->setThermalPowerNeeded(HS_needs+DHW_needs); // Cognet: Added this.
                machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//                machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(HS_needs+DHW_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
            }
            else {
                // Cognet: So in this case the solar thermal energy is not used ? because it was used in the computation HS_needs ? Possible error
                pBuilding->getHeatingUnit()->setThermalPowerNeeded(DHW_needs); // Cognet: Added this.
                machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//                machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(DHW_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
            }
        }

        // heat provided to the heat stocks, priority to the DHW to avoid legionella
        DHW_Pp = fmin(DHW_needs, machineThermalPower);
        HS_Pp = fmin(HS_needs, machineThermalPower - DHW_Pp);

        // the energy that can be used, the minimum between the needs and the available in the space heating tank with HS_Pp to reach Tmin
        heatingUsed = fmin(heatingNeeds, fmax(-pBuilding->getHeatStock()->power(dt, pBuilding->getHeatStock()->getTmin(),0., HS_Pup+HS_Pp, HS_T0, 0., Tamb), 0.));
        // take that heatingUsed from the stock with HS_Pp + HS_Pup provided, and calculate new temperature
        HS_T1 = pBuilding->getHeatStock()->temperature(dt, 0., HS_Pp + HS_Pup, -heatingUsed, HS_T0, 0., Tamb);

        if(pBuilding->getDHWHeatStock() != NULL) //check if DHWtank exists
        {
            // take the VdotUsed from the stock with DHW_Pp provided, and calculate new temperature
            DHW_T1 = pBuilding->getDHWHeatStock()->temperature(dt, VdotUsed, DHW_Pp, 0.0, DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb); // Cognet: If the heatingUnit is not working, then the solar power is not given to the DHW! Possible error.
        }
    }
    else {
        // the heating used is equal to the needs
        heatingUsed = heatingNeeds;
        // checks if DHWtank exists
        if(pBuilding->getDHWHeatStock() != NULL)
        {
            // first priority DHW, compute the new DHW stock temperature with VdotUsed and solar energy provided HS_Pup
            DHW_T1 = pBuilding->getDHWHeatStock()->temperature(dt, VdotUsed, 0., HS_Pup, DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
            HS_T1 = pBuilding->getHeatStock()->temperature(dt, 0., 0., -heatingNeeds, HS_T0, 0.0, Tamb); // Cognet: But here we don't recheck if the temperature is below Tmin ! Because of this we don't have the thermal solar. It is only done if (DHW_T1<Tmin).
            if (DHW_T1 < pBuilding->getDHWHeatStock()->getTmin()) {
                // compute the needs to reach Tmax
                DHW_needs = pBuilding->getDHWHeatStock()->power(dt, pBuilding->getDHWHeatStock()->getTmax(), VdotUsed, HS_Pup, DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
                if (HS_T1 < pBuilding->getHeatStock()->getTmin())
                    HS_needs = pBuilding->getHeatStock()->power(dt, pBuilding->getHeatStock()->getTmax(), 0., -heatingNeeds, HS_T0, 0.0, Tamb);
                // see how the energy unit can provide, HS_Pp is the max heating power that can give the MACHINE hopefully equal to the needs
                if (pBuilding->getHeatingUnit()!=NULL) {
                    if (pBuilding->getHeatingUnit()->isWorking(day)) {
                        pBuilding->getHeatingUnit()->setThermalPowerNeeded(HS_needs+DHW_needs); // Cognet: Added this.
                        machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//                        machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(HS_needs+DHW_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
                    }
                    else {
                        pBuilding->getHeatingUnit()->setThermalPowerNeeded(DHW_needs); // Cognet: Added this.
                        machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//                        machineThermalPower = pBuilding->getHeatingUnit()->getThermalPower(DHW_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
                    }
                }
                // heat provided to the heat stocks, priority to the DHW to avoid legionella
                DHW_Pp = fmin(DHW_needs, machineThermalPower);
                HS_Pp = fmin(HS_needs, machineThermalPower - DHW_Pp);

                // the energy that can be used, the minimum between the needs and the available in the space heating tank with HS_Pp to reach Tmin
                heatingUsed = fmin(heatingNeeds, fmax(-pBuilding->getHeatStock()->power(dt, pBuilding->getHeatStock()->getTmin(),0., HS_Pp, HS_T0, 0., Tamb), 0.));
                // take that heatingUsed from the stock with HS_Pp provided, and calculate new temperature
                HS_T1 = pBuilding->getHeatStock()->temperature(dt, 0., HS_Pp, -heatingUsed, HS_T0, 0., Tamb);
                // take that VdotUsed from the stock with DHW_Pp provided, and calculate new temperature
                DHW_T1 = pBuilding->getDHWHeatStock()->temperature(dt, VdotUsed, DHW_Pp, HS_Pup, DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
            }
            else if (DHW_T1 > pBuilding->getDHWHeatStock()->getTcritical()) {
                // compute the remaining time for the providing HS_Pup to the DHW until Tcritical and the new stock temperature after dt-DHW_time for the decay
                double DHW_time = pBuilding->getDHWHeatStock()->time(pBuilding->getDHWHeatStock()->getTcritical(), VdotUsed, 0., HS_Pup, DHW_T0, pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
                DHW_T1 = pBuilding->getDHWHeatStock()->temperature(dt-DHW_time, VdotUsed, 0., 0., pBuilding->getDHWHeatStock()->getTcritical(), pBuilding->getDHWHeatStock()->getTinlet(), Tamb);
                // compute the new stock temperature achieved after dt-DHW_time considering the initial temperature
                HS_T1 = pBuilding->getHeatStock()->temperature(DHW_time, 0., 0.,-heatingUsed, HS_T0, pBuilding->getHeatStock()->getTinlet(), Tamb);
                HS_T1 = pBuilding->getHeatStock()->temperature(dt-DHW_time, 0., HS_Pup,-heatingUsed, HS_T1, pBuilding->getHeatStock()->getTinlet(), Tamb);
            }
        }
        // checks that the temperature of the main tank does not exceed Tcritical
        HS_T1 = fmin(HS_T1,pBuilding->getHeatStock()->getTcritical());
    }

    // save the new stock temperature
    pBuilding->setHeatStockTemperature(HS_T1);
    if (pBuilding->getDHWHeatStock() != NULL) pBuilding->setDHWStockT(DHW_T1);

    // Energy consumption in primary resource to provide HS_Pp
    pBuilding->addMachinePower(HS_Pp+DHW_Pp);
    if (pBuilding->getHeatingUnit() != NULL) {
        pBuilding->addFuelConsumption(pBuilding->getHeatingUnit()->getFuelConsumption(double(dt),HS_Pp+DHW_Pp,heatPumpSrcTemp));
        pBuilding->addElectricConsumption(pBuilding->getHeatingUnit()->getElectricConsumption(double(dt),HS_Pp+DHW_Pp,heatPumpSrcTemp));
    }

    // deal with the cool tank

    // Solar Thermal Calculation for cooling
//        for (size_t i=0; i<pBuilding->getnZones(); ++i)
//        for (unsigned int j=0; j<pBuilding->getZone(i)->getnRoofs(); ++j)
//        CS_Pup +=(pBuilding->getZone(i)->getRoof(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour),50));
//        pBuilding->setSolarThermalProduction(double(dt)*CS_Pup);

    // we only simulate an air heat pump for the cooling system // Cognet: Why is it necessarily an air heat pump ?
    heatPumpSrcTemp = pClimate->getToutCelsius(day,hour);
    // calculation of the new stock temperature
    CS_T1 = pBuilding->getColdStock()->temperature(dt, 0., 0., -coolingNeeds, CS_T0, 0., Tamb); // Cognet: Should we add CS_Pup ?

    // calculation of the needs to reach the lowest tank temperature of upper limit overtaken
    double CS_needs = 0.;
    if ( CS_T1 >= pBuilding->getColdStock()->getTmax() )
        CS_needs = pBuilding->getColdStock()->power(dt, pBuilding->getColdStock()->getTmin(), 0., CS_Pup-coolingNeeds, CS_T0, 0., Tamb);

    // calculation of the refreshing available from the MACHINE
    double CS_Pp = 0.;
    if (CS_needs < 0. && pBuilding->getCoolingUnit()!=NULL && pBuilding->getCoolingUnit()->isWorking(day)) {
        pBuilding->getCoolingUnit()->setThermalPowerNeeded(CS_needs); // Cognet: Added this.
        CS_Pp = pBuilding->getCoolingUnit()->getThermalPower(heatPumpSrcTemp); // Cognet: Added this.
//        CS_Pp = pBuilding->getCoolingUnit()->getThermalPower(CS_needs, heatPumpSrcTemp); // Cognet: Deleted this, the process is now done in two steps.
    }
    // the energy that can be used, the maximum between the needs and the available (as negative) in the tank with CS_Pp provided
    double coolingUsed = max(coolingNeeds, min(-pBuilding->getColdStock()->power(dt, pBuilding->getColdStock()->getTmax(),0., CS_Pup+CS_Pp, CS_T0,0., Tamb),0.));
    // take that coolingUsed from the stock with CS_Pp provided, and calculate new temperature
    CS_T1 = pBuilding->getColdStock()->temperature(dt, 0., CS_Pup+CS_Pp, -coolingUsed, pBuilding->getColdStockTemperature(), 0., Tamb);
    // save the new stock temperature
    pBuilding->setColdStockTemperature(CS_T1);

    // Energy consumption to provide CS_Pp, the machine power is always positive, unlike the CS_Pp which is always negative
    pBuilding->addMachinePower(-CS_Pp);
    if (pBuilding->getCoolingUnit() != NULL) {
        pBuilding->addFuelConsumption(pBuilding->getCoolingUnit()->getFuelConsumption(double(dt),CS_Pp,heatPumpSrcTemp));
        pBuilding->addElectricConsumption(pBuilding->getCoolingUnit()->getElectricConsumption(double(dt),CS_Pp,heatPumpSrcTemp));
    }

   // repartition of the energy in the zones
    for (size_t i=0; i<pBuilding->getnZones(); ++i) {
        // for the HEATING
        if (pBuilding->getZone(i)->getHeating() > 0) {
            if ( heatingUsed > pBuilding->getZone(i)->getHeating() ) { // demand satisfied
                pBuilding->getZone(i)->setQs(pBuilding->getZone(i)->getHeating());
                heatingUsed -= pBuilding->getZone(i)->getHeating();
            }
            else {
                pBuilding->getZone(i)->setQs(heatingUsed);
                heatingUsed -= heatingUsed;
            }
        }
        else if (pBuilding->getZone(i)->getCooling() <0) {
            // for the COOLING
            if ( coolingUsed < pBuilding->getZone(i)->getCooling() ) { // demand satisfied
                pBuilding->getZone(i)->setQs(pBuilding->getZone(i)->getCooling());
                coolingUsed -= pBuilding->getZone(i)->getCooling();
            }
            else {
                pBuilding->getZone(i)->setQs(coolingUsed);
                coolingUsed -= coolingUsed;
            }
        }
        else pBuilding->getZone(i)->setQs(0.); // nothing is provided
    }
    if (heatingUsed > 1.) throw(string("Error in the repartition of the heating energy - in the control model (noHVAC): ")+toString(heatingUsed));
    if (coolingUsed < -1.) throw(string("Error in the repartition of the cooling energy - in the control model (noHVAC): " + toString(coolingUsed)));

    // compute the electricity produced by the PV
    for (unsigned int i=0; i<pBuilding->getnZones(); i++) {
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnRoofs(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
            pBuilding->addElectricConsumption(-float(dt)*pBuilding->getZone(i)->getRoof(j)->getWindElectricPower(pClimate->getWindSpeed(day,hour),pClimate->getAltitude()) ); // in joules
//            if ( pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tRoof: " << j << "\tIrradiation: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getRoof(j)->getArea() <<
//            "\tIrradiation Tot: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnWalls(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getWall(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
//            if ( pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tRoof: " << j << "\tIrradiation: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getRoof(j)->getArea() <<
//            "\tIrradiation Tot: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnSurfaces(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getSurface(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
        }
    }

}

/**
 * Main control loop. Supposes the buildings' needs are already computed. Simulates the heat and cool tanks, the demands that the heating/cooling units can satisfy, their primary energy consumption. The solar and wind energy production...
 * @param pDistrict Pointer to district to be compute.
 * @param pClimate Pointer to climate used.
 * @param day Day of the year.
 * @param hour Hour of the day.
 */
void Model::noHVAC_Control_EnergyHub(District* pDistrict, Climate* pClimate, unsigned int day, unsigned int hour) {

    // compute the Qs, electricalConsumption, fuelConsumption, tanks temperatures for building that are not connected

    bool useOldAlgo = false; // Cognet: Added this, harcoded, to chose whether to use the old algorithm or not (which doesn't work for DECs).
    if ( useOldAlgo ) { // Cognet: Added this.

        for (size_t i=0; i<pDistrict->getnBuildings(); ++i) {
// Cognet: Start deleted content not used anymore.
//            //beginning of contents added by Dapeng
//            if( (pDistrict->getBuilding(i)->getHeatingUnit()!=NULL && pDistrict->getBuilding(i)->getHeatingUnit()->getLabel() == "Substation") || // Cognet: Changed from small s to capital S in substation
//                (pDistrict->getBuilding(i)->getCoolingUnit()!=NULL && pDistrict->getBuilding(i)->getCoolingUnit()->getLabel() == "Substation") ) {
//    //                    logStream <<"substation is simulating."<<endl << flush;

//            }
//            else {
//            //ending of contents added by Dapeng
// Cognet: End deleted content not used anymore.
            if (pDistrict->getBuilding(i)->getHeatStock() == NULL && pDistrict->getBuilding(i)->getColdStock() == NULL) {
                // no heat/cold tank, only the demands are calculated
                Model::ThermalAllAvailable(pDistrict->getBuilding(i)); // Cognet: Maybe this was for an older version of the code. It seems now Buildings must have tanks when initialized.
            }
            else {
                // controls the machines to produce the demand
                if (pDistrict->getBuilding(i)->getHVACpresence()) {
                    Model::HVAC_Needs(pDistrict->getBuilding(i),pClimate,day,hour);
                    Model::HVAC_Control(pDistrict->getBuilding(i),pClimate,day,hour);
                    Model::HVAC_Available(pDistrict->getBuilding(i),pClimate,day,hour);
                }
                else Model::noHVAC_Control(pDistrict->getBuilding(i),pClimate,day,hour);
            }
        }
    }

// Cognet: Start deleted content, not used anymore.
//        //beginning of contents added by Dapeng
//        //calculate distric energy center
//        for(unsigned int i=0; i<pDistrict->getnDECs(); ++i) {
//            pDistrict->getDEC(i)->setSubInletTemp(day);
//            for(unsigned int j=0; j<pDistrict->getDEC(i)->getnBuildingsLinkedDEC(); ++j) {
//                if (pDistrict->getBuilding(j)->getHVACpresence() &&
//                    ( pDistrict->getDEC(i)->getHeatingUnit()->isWorking(day) ||
//                      pDistrict->getDEC(i)->getCoolingUnit()->isWorking(day) ) ) {
//                    Model::HVAC_Needs(pDistrict->getBuilding(j),pClimate,day,hour);
//                    Model::HVAC_Control_ForDEC(pDistrict->getDEC(i), pDistrict->getBuilding(j), pClimate,day,hour);
//                    Model::HVAC_Available(pDistrict->getBuilding(j),pClimate,day,hour);
//                }
//                else Model::noHVAC_Control_ForDEC(pDistrict->getDEC(i), pDistrict->getDEC(i)->getBuildingLinkedDEC(j), pClimate,day,hour);
//            }
//        }
//        //ending of contents added by Dapeng
// Cognet: End deleted content, not used anymore.

// Cognet: Start added content.
    else { // Use new algorithm

        for (size_t i=0; i<pDistrict->getnBuildings(); ++i) {
            Building* pBui = pDistrict->getBuilding(i);

            if (pBui->getHVACpresence()) {
                throw string("Model::noHVAC_Control_EnergyHub, a building has HVACpresence, but this was not implemented.");
            }
            else if (pBui->getHeatStock()==NULL && pBui->getColdStock()==NULL) {  // No heat/cold tank, only the demands are calculated. Although it seems Buildings must have tanks when initialized.
                Model::ThermalAllAvailable(pBui);
            }
            else {
                Model::noHVAC_Control_Init(pBui,pClimate,day,hour); // Initialization, computes solar and wind electric productions and more.
                Model::noHVAC_Control_Needs(pBui, pClimate, day, hour); // Computes heating and cooling needs for the buildings tanks.
            }
        }
        // All buildings must have their needs computed before determining the thermal power that can be provided by the heating/cooling units (simulates the DistrictEnergyCenters).
        Model::noHVAC_Control_ThermalPower(pDistrict,pClimate,day,hour);

        for (size_t i=0; i<pDistrict->getnBuildings(); ++i) {
            Building* pBui = pDistrict->getBuilding(i);
            if (pBui->getHeatStock()!=NULL or pBui->getColdStock()!=NULL) {
                Model::noHVAC_Control_Finish(pBui,pClimate,day,hour); // Computes and sets new tank temperatures, consumptions, etc.
            }
        }
    }
// Cognet: End added content.


    // compute Qs for your network
    /*
    code here

    1) get the values of the tempratures
    2) choose a temperature for each tank of the buildings connected
    3) get electricity consumption
    4) compute Qvalues
    5) decide actions for the next time step

    */

}

// Cognet: start of added content
double Model::computeSolarThermalPower(Building* pBui, Climate* pClimate, unsigned int day, unsigned int hour, double tankTemp, bool saveValue) {
    // Solar thermal panels.
    double heatingPower = 0.; // Also called "xx_Pup".
    // double coolingPower CS_Pup = 0.; // TODO to implement cooling solar panels.
    for (unsigned int i=0; i<pBui->getnZones(); i++) {
        for (unsigned int j=0; j<pBui->getZone(i)->getnRoofs(); ++j) {
            heatingPower += pBui->getZone(i)->getRoof(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour), pClimate->getWindSpeed(day,hour), tankTemp, saveValue);
        }
        for (unsigned int j=0; j<pBui->getZone(i)->getnWalls(); ++j) {
            heatingPower += pBui->getZone(i)->getWall(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour), pClimate->getWindSpeed(day,hour), tankTemp, saveValue);
        }
        for (unsigned int j=0; j<pBui->getZone(i)->getnSurfaces(); ++j) {
            heatingPower += pBui->getZone(i)->getSurface(j)->getSolarThermalProduction(pClimate->getToutCelsius(day,hour), pClimate->getWindSpeed(day,hour), tankTemp, saveValue);
        }
    }
//    pBui->setSolarThermalProduction(heatingPower); // Cognet: The solar power effectively used is set later.
    return heatingPower;
}

double Model::computeHeatPumpSrcTemp(EnergyConversionSystem* pECS, Climate* pClimate, unsigned int day, unsigned int hour) {
    // Temperatures for heatpumps, the heatPumpSrcTemp depends on the type of heatpump (air or ground vertical/horizontal).
    double heatPumpSrcTemp;
    float z0=0.f,z1=0.f,alpha=0.f;
    if ( (pECS != NULL)  &&  pECS->getGround(z0,z1,alpha) ) {
        heatPumpSrcTemp = pClimate->getTgroundCelsius(day,hour,z0,alpha,z1); // Ground temperature.
    } else {
        heatPumpSrcTemp = pClimate->getToutCelsius(day,hour); // Air temperature.
    }
    return heatPumpSrcTemp;
}

void Model::computeAndAddPhotovoltaicAndWindElectricProduction(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {
    // Compute the electricity produced by the PV and wind turbines.
    for (unsigned int i=0; i<pBuilding->getnZones(); i++) {
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnRoofs(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
            pBuilding->addElectricConsumption(-float(dt)*pBuilding->getZone(i)->getRoof(j)->getWindElectricPower(pClimate->getWindSpeed(day,hour),pClimate->getAltitude()) ); // in joules // Cognet: Should we also add this line to walls and surfaces ? If they can also have wind turbines.
//            if ( pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tRoof: " << j << "\tIrradiation: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getRoof(j)->getArea() <<
//            "\tIrradiation Tot: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnWalls(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getWall(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
//            if ( pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) > 0. )
//            cerr << "Building: " << pBuilding->getId() << "\tRoof: " << j << "\tIrradiation: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation()/pBuilding->getZone(i)->getRoof(j)->getArea() <<
//            "\tIrradiation Tot: " << pBuilding->getZone(i)->getRoof(j)->getShortWaveIrradiation() <<
//            "\tPV electric production: " << pBuilding->getZone(i)->getRoof(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) << endl;
        }
        for (unsigned int j=0; j<pBuilding->getZone(i)->getnSurfaces(); ++j) {
            pBuilding->addSolarPVProduction(float(dt)*pBuilding->getZone(i)->getSurface(j)->getPVElectricProduction(pClimate->getToutCelsius(day,hour)) ); // in joules
        }
    }
}

void Model::computeAndSetBuildingHeatingCoolingNeeds(Building* pBuilding) {
    // Computation of the total demand in heat and cold for the building (double precision required for a precise repartition of the energy)
    double heatingNeeds = 0.;
    double coolingNeeds = 0.;
    for (size_t i=0; i<pBuilding->getnZones(); ++i) {
        heatingNeeds += pBuilding->getZone(i)->getHeating();
        coolingNeeds += pBuilding->getZone(i)->getCooling();
    }
    pBuilding->setHeatingNeeds(heatingNeeds);
    pBuilding->setCoolingNeeds(coolingNeeds);
}

void Model::computeAndSetDHWUsage(Building* pBuilding, unsigned int day, unsigned int hour) {
    // Domestic hot water (DHW) usage in the zones
    double VdotUsed = 0.; // in m^3/s
    for (size_t i=0; i<pBuilding->getnZones(); ++i) {
        VdotUsed += pBuilding->getZone(i)->getDHWConsumption(day,hour)/(1000.*3600.); // DHW consumption in l/h -> m^3/s
    }
    pBuilding->setVdotUsed(VdotUsed);
}

void Model::noHVAC_Control_Init(Building* pBuilding, Climate* pClimate, unsigned int day, unsigned int hour) {
    pBuilding->setTamb(18.f); // 18°C in the cave // Cognet: can we improve this ? using the ground temperature? or climate info ?

    // Initialize the machine power to zero
    pBuilding->setMachinePower(0.);
    pBuilding->setFuelConsumption(0.);

    // Sum up needs over the zones.
    computeAndSetBuildingHeatingCoolingNeeds(pBuilding);

    // Solar and wind electricity.
    computeAndAddPhotovoltaicAndWindElectricProduction(pBuilding, pClimate, day, hour);

    // Sum up domestic hot water (DHW) usage over the zones.
    computeAndSetDHWUsage(pBuilding, day, hour);
}

void Model::noHVAC_Control_Needs(Building* pBui, Climate* pClim, unsigned int day, unsigned int hour) {
    // TODO Implement cooling solar panels, and the critical temperature cases.

    double HS_needs, DHW_needs, CS_needs; // Powers needed by the tanks (not counting solar), that will be asked the the heat/cooling units.
    double HS_SolPp, DHW_SolPp, CS_SolPp; // Solar power provided to the tanks.
    double HS_T1, DHW_T1, CS_T1; // Tank temperatures at next time step.

    // Local variables to simplify code reading.
    float Tamb = pBui->getTamb(); // Temperature of room where tanks are.
    // Power needed to heat the zones.
    double heatingNeeds = pBui->getHeatingNeeds();
    double coolingNeeds = pBui->getCoolingNeeds();

    // Hot storage tank.
    double HS_T0 = pBui->getHeatStockTemperature();
    double HS_Tmin = pBui->getHeatStock()->getTmin();
    double HS_Tmax = pBui->getHeatStock()->getTmax();
    double HS_Pup = computeSolarThermalPower(pBui, pClim, day, hour, HS_T0, true); // Solar thermal power available. Approximation, uses the initial temperature T0 during the whole time step.
    // Domestic hot water tank.
    double DHW_T0 = numeric_limits<double>::signaling_NaN();
    double DHW_Tmin = numeric_limits<double>::signaling_NaN();
    double DHW_Tmax = numeric_limits<double>::signaling_NaN();
    double DHW_Tinlet = numeric_limits<double>::signaling_NaN();
    double DHW_Pup = numeric_limits<double>::signaling_NaN();
    double VdotUsed = pBui->getVdotUsed();
    if (pBui->getDHWHeatStock()) {
        DHW_T0 = pBui->getDHWStockT();
        DHW_Tmin = pBui->getDHWHeatStock()->getTmin();
        DHW_Tmax = pBui->getDHWHeatStock()->getTmax();
        DHW_Tinlet = pBui->getDHWHeatStock()->getTinlet();
        DHW_Pup = computeSolarThermalPower(pBui, pClim, day, hour, DHW_T0, false); // Solar thermal power available. Approximation, uses the initial temperature T0 during the whole time step.
    }
    pBui->setSolTherFracLeft(1.); // Initially, 100% of the solar thermal is available.
    // Cold storage tank.
    double CS_T0 = pBui->getColdStockTemperature();
    double CS_Tmin = pBui->getColdStock()->getTmin();
    double CS_Tmax = pBui->getColdStock()->getTmax();
    double CS_Pup = 0.; // TODO, improve this.

    // Allocate solar thermal power to Heat Stock (HS), Domestic Hot Water (DHW) and Cold Stock (CS).
    if ( ! pBui->getDHWHeatStock() ) { // Without DHW tank.
        // Give all solar thermal to HS.
        HS_SolPp = pBui->getHeatStock()->maxSolPowerToNotExceedTcrit(dt, 0., -heatingNeeds, HS_Pup, HS_T0, 0., Tamb); // Check for Tcritical in HS.
        DHW_SolPp = 0.;
        if (HS_Pup!=0.) { pBui->multSolTherFracLeft(1-HS_SolPp/HS_Pup); } // If not all solar thermal power is used (eg 60%), save what fraction is left (in this case 40%).
    }
    else { // With DHW tank.
        // Heat HS and DHW tanks. (try to give HS priority on solar thermal because lower temperature more efficient)
        bool prioritizeDHWorHS;
        if (  heatingNeeds<=0.  or  ( pBui->getHeatingUnit()!=NULL and (!pBui->getHeatingUnit()->isWorking(day)) )  ) { // If no space heating needed, or if it is not a working day for heating unit.
            prioritizeDHWorHS = true;
        }
        else { // If working day and space heating needs are > 0.
            // Compute temperatures without solar power, to decide which tank to prioritize.
            HS_T1 = pBui->getHeatStock()->temperature(dt, 0., 0., -heatingNeeds, HS_T0, 0., Tamb);
            DHW_T1 = pBui->getDHWHeatStock()->temperature(dt, VdotUsed, 0., 0., DHW_T0, DHW_Tinlet, Tamb);

            if ( HS_T1>HS_Tmin and DHW_T1<DHW_Tmin ) { // If HS stays above Tmin and DHW goes below Tmin, give solar thermal in priority to DHW.
                prioritizeDHWorHS = true;
            }
            else { // If HS goes below Tmin or if both HS and DHW stay above Tmin, give solar thermal in priority to HS.
                prioritizeDHWorHS = false;
            }
        }

        // Give solar thermal according to priority.
        if ( prioritizeDHWorHS ) {
            DHW_SolPp = pBui->getDHWHeatStock()->maxSolPowerToNotExceedTcrit(dt, VdotUsed, 0., DHW_Pup, DHW_T0, DHW_Tinlet, Tamb); // DHW takes max solar thermal to stay below Tcritical.
            // Approximation, the solar power left is proportional to what was not used.
            if (DHW_Pup!=0.) { pBui->multSolTherFracLeft(1-DHW_SolPp/DHW_Pup); }
            double HS_solPowerLeft = HS_Pup*pBui->getSolTherFracLeft(); // Eg: if DHW consumed 80% of what was available to him, then HS gets 20% of what is available.
            HS_SolPp = pBui->getHeatStock()->maxSolPowerToNotExceedTcrit(dt, 0., -heatingNeeds, HS_solPowerLeft, HS_T0, 0., Tamb); // HS takes max solar thermal to stay below Tcritical.
            if (HS_Pup!=0.) { pBui->multSolTherFracLeft(1-HS_SolPp/HS_Pup); }
        }
        else {
            HS_SolPp = pBui->getHeatStock()->maxSolPowerToNotExceedTcrit(dt, 0., -heatingNeeds, HS_Pup, HS_T0, 0., Tamb); // HS takes max solar thermal to stay below Tcritical.
            // Approximation, the solar power left is proportional to what was not used.
            if (HS_Pup!=0.) { pBui->multSolTherFracLeft(1-HS_SolPp/HS_Pup); }
            double DHW_solPowerLeft = DHW_Pup*pBui->getSolTherFracLeft();
            DHW_SolPp = pBui->getDHWHeatStock()->maxSolPowerToNotExceedTcrit(dt, VdotUsed, 0., DHW_solPowerLeft, DHW_T0, DHW_Tinlet, Tamb); // DHW takes max solar thermal to stay below Tcritical.
            if (DHW_Pup!=0.) { pBui->multSolTherFracLeft(1-DHW_SolPp/DHW_Pup); }
        }

    }



//    WARNING, the Tcritical default value is 90 C, so the "CS_T1>CS_Tcritical" may cause problems if no value of Tcritical given
//    if (  (pBui->getCoolingUnit()!=NULL and pBui->getCoolingUnit()->isWorking(day))  or  pBui->getCoolingUnit()==NULL  ) { // If working day or no cooling unit, use solar thermal.
//        // Compute CS temperature with solar thermal.
//        CS_T1 = pBui->getColdStock()->temperature(dt, 0., 0., CS_Pup-coolingNeeds, CS_T0, 0., Tamb);
//        if (CS_T1 > CS_Tcritical) { // If CS stays above Tcritical.
//            CS_SolPp = CS_Pup;
//        }
//        else { // If DHW reaches Tcritical.
//            double timeToTcrit = pBui->getColdStock()->time(CS_Tcritical, 0., 0., CS_Pup-coolingNeeds, CS_T0, 0., Tamb);
//            // Approximation : instead of CS receiving 'CS_Pup' for 'timeToTcrit', then 'zero' for 'dt-timeToTcrit', it receives 'CS_Pup*timeToTcrit/dt' for 'dt' time.
//            CS_SolPp = CS_Pup*timeToTcrit/dt;
//        }
//    }
//    else {
        CS_SolPp = 0.; // Use this for now.
//    }
    pBui->setHS_SolPp( HS_SolPp );
    pBui->setDHW_SolPp( DHW_SolPp );
    pBui->setCS_SolPp( CS_SolPp );
    pBui->setSolarThermalProduction(float(dt)*(HS_SolPp+DHW_SolPp+CS_SolPp));


    // Power needs, that are asked to the heating/cooling units.
    HS_needs = 0.;
    DHW_needs = 0.;
    if ( pBui->getHeatingUnit() ) { // If heating unit exists.
        if ( pBui->getHeatingUnit()->isWorking(day) and heatingNeeds>0. ) {  // If it is a working day and space heating is needed, heat HS.
            HS_T1 = pBui->getHeatStock()->temperature(dt, 0., HS_SolPp, -heatingNeeds, HS_T0, 0., Tamb);
            if (HS_T1 < HS_Tmin) { // If HS goes below Tmin, heat to Tmax.
                HS_needs = pBui->getHeatStock()->power(dt, HS_Tmax, 0., HS_SolPp-heatingNeeds, HS_T0, 0.0, Tamb);
            }
        } // Otherwise HS_needs=0

        if ( pBui->getDHWHeatStock() ) { // Heat DHW.
            DHW_T1 = pBui->getDHWHeatStock()->temperature(dt, VdotUsed, 0., DHW_SolPp, DHW_T0, DHW_Tinlet, Tamb);
            if (DHW_T1 < DHW_Tmin) { // If DHW temperature goes below Tmin, heat to Tmax.
                DHW_needs = pBui->getDHWHeatStock()->power(dt, DHW_Tmax, VdotUsed, DHW_SolPp, DHW_T0, DHW_Tinlet, Tamb);
            }
        } // Otherwise DHW_needs=0
    }

    CS_needs = 0.;
    if ( pBui->getCoolingUnit()!=NULL and pBui->getCoolingUnit()->isWorking(day) and coolingNeeds<0  ) { // If cooling unit is working and space cooling is needed.
        CS_T1 = pBui->getColdStock()->temperature(dt, 0., 0., CS_SolPp-coolingNeeds, CS_T0, 0., Tamb);
        if (CS_T1 > CS_Tmax) { // If CS temperature goes above Tmax, cool to Tmin.
            CS_needs = pBui->getColdStock()->power(dt, CS_Tmin, 0., CS_Pup-coolingNeeds, CS_T0, 0., Tamb);
        }
    }
    pBui->setHS_needs( HS_needs );
    pBui->setDHW_needs( DHW_needs );
    pBui->setCS_needs( CS_needs );
}

void Model::noHVAC_Control_ThermalPower(District* pDis, Climate* pClim, unsigned int day, unsigned int hour) {

    // Give heating/cooling units the thermal power needs information.
    for (size_t i=0 ; i<pDis->getnBuildings() ; ++i) {
        Building* pBui = pDis->getBuilding(i);
        if (pBui->getHeatingUnit()) {
            float heatDemanded;
            float heatDemandedHS;// Added by Max
            float heatDemandedDHW;// Added by Max
            if ( not pBui->hasImposedHeatDemand(day, hour, heatDemanded) ) { // If an imposed heat demand value exists, use the imposed value. (To get rid of imposed value code, remove this line.)
                heatDemandedHS = pBui->getHS_needs();// Added by Max
                heatDemandedDHW = pBui->getDHW_needs();// Added by Max
                heatDemanded = heatDemandedHS + heatDemandedDHW;
            }
            else{
                heatDemandedHS = 0.8*heatDemanded; // Added by Max. 0.8 and 0.2 are arbitrary
                heatDemandedDHW = 0.2*heatDemanded; // Added by Max. Needed to distinguish DHW from HS from a single value of imposedHeatDemand in the XML file
            } // (To get rid of imposed value code, remove this line.)
            pBui->getHeatingUnit()->setThermalPowerNeeded( heatDemanded );
            pBui->getHeatingUnit()->setThermalPowerNeededHS(heatDemandedHS); // Added by Max
            pBui->getHeatingUnit()->setThermalPowerNeededDHW(heatDemandedDHW); // Added by Max
        }
        if (pBui->getCoolingUnit()) {
            pBui->getCoolingUnit()->setThermalPowerNeeded( pBui->getCS_needs() );
        }
    }

    // Makes the District Energy Centers converge to the solution, to determine how much thermal power each substation can provide.
    for(size_t i=0 ; i<pDis->getnDECs() ; i++) {
        pDis->getDEC(i)->convergeToEquilibrium(pClim, day, hour);
    }

    // Compute the thermal power that heating/cooling units can provide.
    for (size_t i=0 ; i<pDis->getnBuildings() ; ++i) {
        Building* pBui = pDis->getBuilding(i);
        double HS_Pp, DHW_Pp, CS_Pp;
        // Heating Unit
        if (pBui->getHeatingUnit()) {
            double heatPumpSrcTemp = computeHeatPumpSrcTemp(pBui->getHeatingUnit(), pClim, day, hour); // TODO, make it so that the heat pump can themselves go get the temperature (of air or ground)

            double machineThermalPower = pBui->getHeatingUnit()->getThermalPower(heatPumpSrcTemp);
            DHW_Pp = min(pBui->getDHW_needs(), machineThermalPower); // Priority to the DHW in machineThermalPower not enough.
            HS_Pp = min(pBui->getHS_needs(), machineThermalPower-DHW_Pp); // Give what's left of machineThermalPower to HS.
        } else {
            HS_Pp = 0.;
            DHW_Pp = 0.;
        }
        pBui->setHS_Pp( HS_Pp );
        pBui->setDHW_Pp( DHW_Pp );

        // Cooling Unit
        if (pBui->getCoolingUnit()) {
            double heatPumpSrcTemp = computeHeatPumpSrcTemp(pBui->getCoolingUnit(), pClim, day, hour); // TODO, make it so that the heat pump can themselves go get the temperature (of air or ground)
            CS_Pp = pBui->getCoolingUnit()->getThermalPower(heatPumpSrcTemp);
        } else {
            CS_Pp = 0.;
        }
        pBui->setCS_Pp( CS_Pp );
    }
}

void Model::noHVAC_Control_Finish(Building* pBui, Climate* pClim, unsigned int day, unsigned int hour) {

    // Temporary variables to simplify code reading.
    double heatingNeeds = pBui->getHeatingNeeds(); // Power needs to heat zones of the building.
    double coolingNeeds = pBui->getCoolingNeeds(); // Power needs to cool zones of the building.
    double HS_Tmin = pBui->getHeatStock()->getTmin();
    double CS_Tmax = pBui->getColdStock()->getTmax();
    double HS_Pp = pBui->getHS_Pp();
    double CS_Pp = pBui->getCS_Pp();
    double HS_SolPp = pBui->getHS_SolPp();
    double CS_SolPp = pBui->getCS_SolPp();
    double HS_T0 = pBui->getHeatStockTemperature();
    double CS_T0 = pBui->getColdStockTemperature();
    float Tamb = pBui->getTamb();

    // Thermal power that is effectively be used for space heating/cooling.
    double heatingUsed = min(heatingNeeds, max(-pBui->getHeatStock()->power(dt, HS_Tmin, 0., HS_SolPp+HS_Pp, HS_T0, 0., Tamb), 0.));
    double coolingUsed = max(coolingNeeds, min(-pBui->getColdStock()->power(dt, CS_Tmax, 0., CS_SolPp+CS_Pp, CS_T0, 0., Tamb), 0.));

    // Compute new temperatures.
    double HS_T1 = pBui->getHeatStock()->temperature(dt, 0., HS_SolPp+HS_Pp, -heatingUsed, HS_T0, 0., Tamb);
    double CS_T1 = pBui->getColdStock()->temperature(dt, 0., CS_SolPp+CS_Pp, -coolingUsed, CS_T0, 0., Tamb);

// The Tcritical check is done in noHVAC_Control_Needs
//    // checks that the temperature of the main tank does not exceed Tcritical
//    HS_T1 = min(HS_T1,pBuilding->getHeatStock()->getTcritical());

    pBui->setHeatStockTemperature(HS_T1);
    pBui->setColdStockTemperature(CS_T1);


    double DHW_Pp = pBui->getDHW_Pp();
    if (pBui->getDHWHeatStock()) {
        double DHW_Tinlet = pBui->getDHWHeatStock()->getTinlet();
        double VdotUsed = pBui->getVdotUsed();
        double DHW_SolPp = pBui->getDHW_SolPp();
        double DHW_T0 = pBui->getDHWStockT();
        double DHW_T1 = pBui->getDHWHeatStock()->temperature(dt, VdotUsed, DHW_SolPp+DHW_Pp, 0., DHW_T0, DHW_Tinlet, Tamb);
        pBui->setDHWStockT(DHW_T1);
    }

    // Energy consumption in primary resource.
    float heatdemanded; // Added by Max
    if (pBui->hasImposedHeatDemand(day,hour,heatdemanded)){ // (To get rid of imposed value code, remove this line.)
        double heatPumpSrcTemp = computeHeatPumpSrcTemp(pBui->getHeatingUnit(), pClim, day, hour); // (To get rid of imposed value code, remove this line.)
        double machineThermalPower = pBui->getHeatingUnit()->getThermalPower(heatPumpSrcTemp); // (To get rid of imposed value code, remove this line.)
        pBui->addMachinePower(machineThermalPower); // (To get rid of imposed value code, remove this line.)
        // Added by Max. In order to impose the heat demand to the Fuel and electricConsumption if it's the case.
        pBui->addFuelConsumption(pBui->getHeatingUnit()->getFuelConsumption(double(dt), heatdemanded, heatPumpSrcTemp));
        pBui->addElectricConsumption(pBui->getHeatingUnit()->getElectricConsumption(double(dt), heatdemanded, heatPumpSrcTemp)); //TODO, what if we impose the heat demand? It doesn't seem to get the electricity consumption with the imposed heat demand (Max).
    } else { // (To get rid of imposed value code, remove this line.)
        pBui->addMachinePower(HS_Pp+DHW_Pp);
        if (pBui->getHeatingUnit()) {
            double heatPumpSrcTemp = computeHeatPumpSrcTemp(pBui->getHeatingUnit(), pClim, day, hour); // TODO, make it so that the heat pump can themselves go get the temperature (of air or ground)
            pBui->addFuelConsumption(pBui->getHeatingUnit()->getFuelConsumption(double(dt), HS_Pp+DHW_Pp, heatPumpSrcTemp));
            pBui->addElectricConsumption(pBui->getHeatingUnit()->getElectricConsumption(double(dt), HS_Pp+DHW_Pp, heatPumpSrcTemp)); //TODO, what if we impose the heat demand? It doesn't seem to get the electricity consumption with the imposed heat demand (Max).
        }
    } // (To get rid of imposed value code, remove this line.)

    // Energy consumption in primary resource (machine power is always positive, unlike CS_Pp which is always negative).
    pBui->addMachinePower(-CS_Pp);
    if (pBui->getCoolingUnit()) {
        double heatPumpSrcTemp = computeHeatPumpSrcTemp(pBui->getCoolingUnit(), pClim, day, hour); // TODO, make it so that the heat pump can themselves go get the temperature (of air or ground)
        pBui->addFuelConsumption(pBui->getCoolingUnit()->getFuelConsumption(double(dt),CS_Pp,heatPumpSrcTemp));
        pBui->addElectricConsumption(pBui->getCoolingUnit()->getElectricConsumption(double(dt),CS_Pp,heatPumpSrcTemp));
    }

    // Distribution of the energy in the zones.
    for (size_t i=0 ; i<pBui->getnZones() ; ++i) {
        double zoneHeatingNeeds = pBui->getZone(i)->getHeating();
        double zoneCoolingNeeds = pBui->getZone(i)->getCooling();
        double heatingUnsatisfied(0.0); //Added by Max
        double coolingUnsatisfied(0.0); //Added by Max
        // Heating
        if (zoneHeatingNeeds > 0) {
            if ( heatingUsed > zoneHeatingNeeds ) { // Demand satisfied.
                pBui->getZone(i)->setQs(zoneHeatingNeeds);
                heatingUsed -= zoneHeatingNeeds;
            }
            else { // Not enough heat available.
                pBui->getZone(i)->setQs(heatingUsed);
                heatingUsed -= heatingUsed;
                heatingUnsatisfied += zoneHeatingNeeds - heatingUsed; //Added by Max
            }
        }
        // Cooling
        else if (zoneCoolingNeeds < 0) {
            if ( coolingUsed < zoneCoolingNeeds ) { // Demand satisfied.
                pBui->getZone(i)->setQs(zoneCoolingNeeds);
                coolingUsed -= zoneCoolingNeeds;
            }
            else { // Not enough cold available.
                pBui->getZone(i)->setQs(coolingUsed);
                coolingUsed -= coolingUsed;
                coolingUnsatisfied += zoneCoolingNeeds - coolingUsed; // Added by Max
            }
        }
        else {
            pBui->getZone(i)->setQs(0.); // No heating or cooling needed.
        }
        pBui->setHeatingDemandUnsatisfied(heatingUnsatisfied);
        pBui->setCoolingDemandUnsatisfied(coolingUnsatisfied);
     }
    // Check that all the energy has been distributed.
     if (heatingUsed > 1.) { throw(string("Error in the repartition of the heating energy - in the control model (noHVAC): ") + toString(heatingUsed)); }
     if (coolingUsed < -1.) { throw(string("Error in the repartition of the cooling energy - in the control model (noHVAC): " + toString(coolingUsed))); }

}
// Cognet: End of added content.

//end of contents added by Dapeng


double Model::HVAC_saturatedSteamPressure(double T) {
    // transformation in kelvin
    T += 273.15;
    // saturation pressure from Murphy and Koop (2005) in Pa
    // assuming liquid above 273.15K
    if   ( T >= 332.0 ) throw ("Temperature outside allowed range (>332K)");
    else if ( T > 273.15 ) return std::exp( 54.842763 - 6763.22/T - 4.210*std::log(T) + 0.000367*T + std::tanh(0.0415*(T-218.8)) * (53.878 - 1331.22/T - 9.44523*std::log(T) + 0.014025*T) );
    else if ( T > 110.0 )  return std::exp( 9.550426 - 5723.265/T + 3.53068*std::log(T) - 0.00728332*T );
    else throw ("Temperature outside allowed range (<110K)");
}

double Model::HVAC_moistureContent(double T, double RH, double Patm) {

    double pss = HVAC_saturatedSteamPressure(T);

    return 0.62197 * (RH*pss)/(Patm - RH*pss);

}

double Model::HVAC_relativeHumidity(double T, double w, double Patm) {

    double pss = HVAC_saturatedSteamPressure(T);

    return 100.0*(Patm*w/(0.62197+w))/pss;

}

    double Model::HVAC_supplyAirTemperature(double Qs, double T) {

        if ( Qs > 0.0 ) return 35.0;
        else if ( Qs < 0.0 ) return 15.0;
        else return T;

    }

    double Model::HVAC_massFlowRate(double T2, double Patm, double w2, double Qs, double &T3, double T5, double np) {

    // mass flow rate for fresh air requirements
    double Mair = 28.97; // from p.6 Jones
    double R    = 8314.41; // from p.11 Jones
    double rhoAir = Patm*Mair/(R*(T2+273.15));
    double m1dotmin = rhoAir*np*8.e-3; // 8 litres per person per second

    // humid specific heat
    double Cp = 1000.0*(1.007+w2*1.84);
    double m1dot;
    if (T3!=T5) m1dot = Qs / (Cp*(T3-T5));
    else m1dot = m1dotmin;

    if ( m1dot < m1dotmin ) {

//        std::cerr << "Minimal fresh air requirements (" << m1dotmin << ") not met: " << m1dot << std::endl;
        T3 = (Qs+m1dotmin*Cp*T5)/(m1dotmin*Cp);
        m1dot = m1dotmin;
//        std::cerr << "Adjusting supply temperature to: " << T3 << std::endl;

    }

    return m1dot;

}

double Model::HVAC_totalMoistureContent(double w, double Ql, double mdot) { if (mdot > 0.) return w + Ql / (2.26e6*mdot); else return w; }

    double Model::HVAC_temperatureChange(double deltaT, double T2, double T3) {

        if ( T3 > T2 ) { // heating case
            return +deltaT;
        }
        else if ( T3 < T2 ) { // cooling case
            return -deltaT;
        }
        else return 0.0;

    }

    double Model::HVAC_moistureControl(double Tmax, double &T2, double T3, double &deltaTHVAC, double deltaT, double T5, double &w2, double w5prime, double Patm, double &evaporation, bool evaporativeCooling) {

        // comfort zone known from 20 to 26∞C normally, here between Tmin and Tmax 70%

        if ( T3 > T2 ) { // heating case

            if ( w5prime < HVAC_moistureContent(T5, 0.3, Patm) ) return (HVAC_moistureContent(T5, 0.3, Patm) - w5prime + w2); // traitement de l'air (humiditÈ)
            else if ( w5prime > HVAC_moistureContent(T5, 0.7, Patm) && w5prime < HVAC_moistureContent(Tmax, 0.7, Patm) ) {
                // heating to avoid dehumidifying
                deltaTHVAC = (HVAC_bulbTemperature(w5prime, 0.7, Patm)-T5);
//                std::cerr << "Heating to avoid dehumidifying" << std::endl;
                return w2;
            }
            else if ( w5prime > HVAC_moistureContent(Tmax, 0.7, Patm) ) {
                // dÈshumidification
                return std::max(std::min(std::min(HVAC_moistureContent(Tmax, 0.7, Patm) - w5prime + w2, HVAC_moistureContent(T3, 1.0, Patm)), HVAC_moistureContent(T5, 0.7, Patm) - w5prime + w2), 0.0);

            }
            else return w2;                            // pas de traitement de l'air

        }
        else if ( T3 < T2 ) { // cooling case

            if (evaporativeCooling) {

                double H2 = HVAC_enthalpyDryAir(T2) + w2*HVAC_enthalpyWaterVapour(T2); // enthalpy of Point 2
                double w2prime = std::min( HVAC_moistureContent(T5,0.7,Patm) - w5prime + w2, HVAC_moistureContent(T3+deltaT,1.0,Patm) ); // min humiditÈ
                double T2prime = (H2 + 0.026 - w2prime*2501.0) / (1.007 + w2prime*1.84);

                if (T2prime < T3) { // cas spÈcial
                    T2prime = T3;
                    w2prime = (H2 - 1.007*T2prime + 0.026) / (2501.0 + 1.84*T2prime);
                }
//                std::cerr << "T2prime: " << T2prime << "\tw2prime: " << w2prime << std::endl;

                // redefine T2,w2
                if (w2prime > w2) { // si Èvaporation!
                    evaporation = w2prime - w2; // kg/kg dry air
                    T2 = T2prime;
                    w5prime = w2prime + w5prime - w2;
                    w2 = w2prime;
                }
            }

            if ( w5prime > HVAC_moistureContent(T5, 0.7, Patm) ) { // deshumidification
                return std::max(std::min(HVAC_moistureContent(T5, 0.7, Patm) - w5prime + w2, HVAC_moistureContent(T3, 1.0, Patm)), 0.0);
            }
            else if ( w5prime < HVAC_moistureContent(T5, 0.3, Patm)) { // humidification
                return (HVAC_moistureContent(T5, 0.3, Patm) - w5prime + w2);
            }
            else return std::min(w2, HVAC_moistureContent(T3, 1.0, Patm));

        }
        else return w2;

    }

    double Model::HVAC_enthalpyDryAir(double T) { // kJ/kg

        if ( T >= 0.0 ) return 1.007*T - 0.026;
        else return 1.005*T;

    }

    double Model::HVAC_enthalpyWaterVapour(double T) { // kJ/kg

        return 2501 + 1.84*T;

    }

    double Model::HVAC_humidify(/*double T2, */double w2, double T3, double w3) {

        if ( w3 > w2 ) {

            return 1000*(w3-w2)*HVAC_enthalpyWaterVapour(T3); // steam injection in J/kg

        }
        else return 0.0;

    }

    double Model::HVAC_reheat(double T3, double w3, double Patm) {

        if ( std::abs(HVAC_moistureContent(T3, 1.0, Patm)-w3) > 1e-5 ) {
            // reheat bit
            double deltaHa = HVAC_enthalpyDryAir(T3) - HVAC_enthalpyDryAir(HVAC_bulbTemperature(w3, 1.0, Patm));
            double deltaHg = w3*(HVAC_enthalpyWaterVapour(T3) - HVAC_enthalpyWaterVapour(HVAC_bulbTemperature(w3, 1.0, Patm)));

            return 1000*(deltaHa + deltaHg);
        }
        else return 0.0;

    }

    double Model::HVAC_bulbTemperature(double w, double RH, double Patm) {

        double Trefup = 58.84;
        double Trefdown = -163.14;

        // dÈbut de l'algorithme

        double Tup = Trefup;
        double Tdown = Trefdown;
        double Tmid = 0.5*(Tup + Tdown);

        do {

            if ( w < HVAC_moistureContent(Tmid, RH, Patm) ) Tup = Tmid;
            else Tdown = Tmid;

            Tmid = 0.5*(Tup + Tdown);

       //     std::cerr << "Tmid: " << Tmid << std::endl;

        }
        while ( std::abs(HVAC_moistureContent(Tmid, RH, Patm)-w) > 1e-5 && (Tmid-Trefdown) > 1e-2 && (Trefup-Tmid) > 1e-2 );

        return Tmid;

    }

    double Model::HVAC_bulbTemperature(double cooling, double t2prime, double w2prime, double Patm) {

    double Trefup = 58.84;
    double Trefdown = -163.14;

    // dÈbut de l'algorithme

    double Tup = Trefup;
    double Tdown = Trefdown;
    double Tmid = 0.5*(Tup + Tdown);

    double Ha,Hg;

    do {

        //std::cerr << "cooling: " << cooling << std::endl;
        //std::cerr << "Tmid: " << Tmid << std::endl;
        Ha = HVAC_enthalpyDryAir(Tmid) - HVAC_enthalpyDryAir(t2prime);
        Hg = HVAC_moistureContent(Tmid,1.0,Patm)*HVAC_enthalpyWaterVapour(Tmid) - w2prime*HVAC_enthalpyWaterVapour(t2prime);
        //std::cerr << "enthalpy difference: " << 1000*(Ha+Hg) << std::endl;

        if ( cooling < 1000*(Ha+Hg) ) Tup = Tmid;
        else Tdown = Tmid;

        Tmid = 0.5*(Tup + Tdown);

    }
    while ( std::abs(cooling - 1000*(Ha+Hg)) > 1e-5 && (Tmid-Trefdown) > 1e-2 && (Trefup-Tmid) > 1e-2 );

    return Tmid;

    }

    void  Model::HVAC_heat(double T2, double w2, double T3, double w3, /*double Patm,*/ double &heating, double &humidification) {

        if ( T3 > T2 && w3 >= w2 ) { // alors le heating a un sens

            // heating bit
            double deltaHa = HVAC_enthalpyDryAir(T3) - HVAC_enthalpyDryAir(T2);
            double deltaHg = w2*(HVAC_enthalpyWaterVapour(T3) - HVAC_enthalpyWaterVapour(T2));
            heating = 1000*(deltaHa + deltaHg); // heating the air and the vapour content in J/kg
            // humidification bit
            humidification = HVAC_humidify(/*T2,*/w2,T3,w3);

        }

    }

    void  Model::HVAC_cool(double T2, double w2, double T3, double &w3, /*double T5, double w5,*/ double Patm, double &cooling, double& reheating, double &humidification) {

        if ( T3 < T2 ) { // alors le cooling a un sens

            if ( w2 <= w3 ) {

                    // cas 1: cooling puis humidification if necessary

                    double deltaHa = HVAC_enthalpyDryAir(T3) - HVAC_enthalpyDryAir(T2);
                    double deltaHg = w2*(HVAC_enthalpyWaterVapour(T3) - HVAC_enthalpyWaterVapour(T2));
                    cooling = 1000*(deltaHa + deltaHg); // cool the air and the vapour content in J/kg

                    humidification = HVAC_humidify(/*T2,*/w2,T3,w3);

            }
            else { // cooling et dÈshumidification (reheat)

                if ( std::abs(HVAC_moistureContent(T3, 1.0, Patm)-w3) > 1e-5 ) { // cas avec reheat

                    reheating = HVAC_reheat(T3, w3, Patm);

                    // cooling bit
                    double deltaHa = HVAC_enthalpyDryAir(HVAC_bulbTemperature(w3, 1.0, Patm)) - HVAC_enthalpyDryAir(T2);
                    double deltaHg = w3*HVAC_enthalpyWaterVapour(HVAC_bulbTemperature(w3, 1.0, Patm)) - w2*HVAC_enthalpyWaterVapour(T2);

                    cooling = 1000*(deltaHa + deltaHg); // heating the air and the vapour content in J/kg
                }
                else { // cas sans reheat

                    double deltaHa = HVAC_enthalpyDryAir(T3) - HVAC_enthalpyDryAir(T2);
                    double deltaHg = w3*HVAC_enthalpyWaterVapour(T3) - w2*HVAC_enthalpyWaterVapour(T2);

                    cooling = 1000*(deltaHa + deltaHg); // heating the air and the vapour content in J/kg
                }
            }
        }
        else if ( w3 < w2 ) {         // cooling pour dÈshumidifier

            reheating = HVAC_reheat(T3, w3, Patm);

            // cooling bit
            double deltaHa = HVAC_enthalpyDryAir(HVAC_bulbTemperature(w3, 1.0, Patm)) - HVAC_enthalpyDryAir(T2);
            double deltaHg = w3*HVAC_enthalpyWaterVapour(HVAC_bulbTemperature(w3, 1.0, Patm)) - w2*HVAC_enthalpyWaterVapour(T2);

            cooling = 1000*(deltaHa + deltaHg); // heating the air and the vapour content in J/kg

        }

    }

double Model::deterministicWindowsNvent(double NventMax, double Tin, double Tout) { return NventMax*exp(1.459 + 0.1448*Tout - 0.1814*Tin)/(1 + exp(1.459 + 0.1448*Tout - 0.1814*Tin)); }

float Model::ventilationFlowRate(float Ti, float Te, float V, float A) {

    // FrÈdÈric Haldi, 9.02.2010

    // Ti, Te: internal and external temperatures (K)
    // Qtot: total air flow (m3/s)
    // Qwind: volume of airflow induced by wind (m3/s)
    // Qstack: volume of airflow induced by stack effect (m3/s)
    // K: coefficient of effectiveness (ranges from about 0.4 for wind hitting an opening at a 45∞
    //    angle of incidence to 0.8 for wind hitting directly at a 90∞ angle.
    // A: area of window (m2)
    // V: outdoor wind speed (m/s)
    // Cd: discharge coefficient
    // H: window height (m)

    float K  = 0.05f;
    float Cd = 0.65f;
    float g  = 9.81f;

    // intermediate variables
    float H = sqrt(A); // approximation of the window height (square window)

    // effect of the wind
    float Qwind = K*A*V;
    // Case of a single opening
    float Qstack = Cd*A/3.f*sqrt(g*H*abs(Ti-Te)/(0.5f*(Ti+Te)))*((Ti-Te)>=0.f ? 1.f : -1.f);

    //cerr << "VdotVent: " << sqrt(pow(Qwind,2) + pow(Qstack,2)) << "\tTe: " << Te-273.15f << "\tV: " << V << "\tTi: " << Ti-273.15f << "\tA: " << A << "\tH: " << H << "\tZone Volume: " << pZone->getVi() << endl;
    return sqrt(pow(Qwind,2) + pow(Qstack,2));

}

float Model::deterministicShadingAction(float irradiance, float irradianceCutOff, float lambda) {

    return 1.f-1.f/(1.f+exp(-lambda*(irradiance-irradianceCutOff))); /* 1 - sigmoid curve */

}

void Model::lowerShadingAction(Zone* pZone, Wall* pWall, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour) {

    cout << "Uh-oh: lowerShadingAction" << endl << flush;
    // This function simulates the occurence of actions on lower shading devices (covering the
    // vision window). If an action is simulated, the chosen unshaded fraction is then predicted.
    // Input parameters are: indoor workplane illuminance (Lumint) and outdoor global horizontal
    // illuminance (Evg).

    // Lumint: Indoor illuminance next to window (lux)
    float Lumint = pWall->getTotalInternalIlluminance0()*pWall->getLowerShadingState();
    // we take the previous timestep shading state to compute Lumint

   // Evg: Outdoor illuminance in the horizontal plane without obstructions (lux)
   float Evg = pClimate->getIgh_vis(day,hour);

   // --- Constants for lower blinds (covering the vision window) --------------------------------
   // XML ENTRIES TO BE ADDED

   // Probability of lowering on arrival
   float a01arr_LB    = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_LB(0);
   float b01inarr_LB  = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_LB(1);
   float b01sarr_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_LB(2);
   // a01arr_LB <- -7.411; b01inarr_LB <-  1.035E-03; b01sarr_LB <-  2.166  # Baisse arrivÈe

   // Probability of raising on arrival
   float a10arr_LB   = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_LB(0);
   float b10inarr_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_LB(1);
   float b10sarr_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_LB(2);
   // a10arr_LB <- -1.520; b10inarr_LB <- -6.535E-04; b10sarr_LB <- -3.139  # Lever arrivÈe

   // Probability of lowering during presence
   float a01int_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_LB(0);
   float b01inint_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_LB(1);
   float b01sint_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_LB(2);
   // a01int_LB <- -8.013; b01inint_LB <-  8.405E-04; b01sint_LB <-  1.270  # Baisse prÈsence

   // Probability of raising during presence
   float a10int_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_LB(0);
   float b10inint_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_LB(1);
   float b10sint_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_LB(2);
   // a10int_LB <- -3.568; b10inint_LB <- -1.567E-04; b10sint_LB <- -2.678  # Lever prÈsence

   // Probability of full lowering
   float afulllower_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfulllower_LB(0);
   float boutfulllower_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfulllower_LB(1);
   float bsfulllower_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfulllower_LB(2);
   // afulllower_LB <- -0.280; boutfulllower_LB <- 2.108e-06 ; bsfulllower_LB <- -2.239

   // Probability of full raising
   float afullraise_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfullraise_LB(0);
   float boutfullraise_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfullraise_LB(1);
   float bsfullraise_LB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfullraise_LB(2);
   // afullraise_LB <- 0.03885; boutfullraise_LB <- -2.513e-05; bsfullraise_LB <- 2.44

   // Choice of new unshaded fraction
   float aSFlower = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getDistFrac_LB(0);
   float bSFlower = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getDistFrac_LB(1);
   float shapelower = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getDistFrac_LB(2);
   // aSFlower<- -2.294; bSFlower <- 1.522; shapelower <- 1.708


   // --- Simulation of lower blinds -------------------------------------------------------------

    float problower,probraise,ptotlow,ptotraise;
    int ActLow;


   // === Case of absence ========================================================================
   if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {
       pWall->setLowerShadingState(pWall->getLowerShadingState());
   }

   // === Case of arrival ========================================================================
   else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {   problower = exp(a01arr_LB+b01inarr_LB*Lumint+b01sarr_LB*(pWall->getLowerShadingState()))/
        (1.f+exp(a01arr_LB+b01inarr_LB*Lumint+b01sarr_LB*(pWall->getLowerShadingState())));
       probraise = exp(a10arr_LB+b10inarr_LB*Lumint+b10sarr_LB*(pWall->getLowerShadingState()))/
        (1.f+exp(a10arr_LB+b10inarr_LB*Lumint+b10sarr_LB*(pWall->getLowerShadingState())));
       if ((pWall->getLowerShadingState()) == 1.f) { probraise = 0.f; }
       if ((pWall->getLowerShadingState()) == 0.f) { problower = 0.f; }

       // --- Case of problower >= probraise -----------------------------------------------------
       if (problower >= probraise)
       {   if (randomUniform(0.f,1.f) < problower) { ActLow = -1; }
           else if (randomUniform(0.f,1.f) < probraise) { ActLow = +1; }
           else { ActLow = 0; }
       }

       // --- Case of problower < probraise ------------------------------------------------------
       else
       {   if (randomUniform(0.f,1.f) < probraise) { ActLow = +1; }
           else if (randomUniform(0.f,1.f) < problower) { ActLow = -1; }
           else { ActLow = 0; }
       }

       // --- Choice of the new unshaded fraction ------------------------------------------------
       if (ActLow == -1) // If a lowering action was predicted
       {   ptotlow = exp(afulllower_LB+boutfulllower_LB*Evg+bsfulllower_LB*(pWall->getLowerShadingState()))/
           (1.f+exp(afulllower_LB+boutfulllower_LB*Evg+bsfulllower_LB*(pWall->getLowerShadingState())));
           if (randomUniform(0.f,1.f) < ptotlow) { pWall->setLowerShadingState(0.f); }
           else
           {   float Reduction = randomWeibull(exp(aSFlower + bSFlower*(pWall->getLowerShadingState())), shapelower);
               pWall->setLowerShadingState( 0.01f*round(100.f*max((pWall->getLowerShadingState())-Reduction,0.01f)) );
           }
       }

       else if (ActLow == +1) // If a raising action was predicted
       {   ptotraise = exp(afullraise_LB+boutfullraise_LB*Evg+bsfullraise_LB*(pWall->getLowerShadingState()))/
           (1.f+exp(afullraise_LB+boutfullraise_LB*Evg+bsfullraise_LB*(pWall->getLowerShadingState())));
           if (randomUniform(0.f,1.f) < ptotraise) { pWall->setLowerShadingState(1.f); }
           else pWall->setLowerShadingState( 0.01f*round(100.f*randomUniform((pWall->getLowerShadingState()),1.f)) );
       }

       else // If no action was predicted
       {   pWall->setLowerShadingState(pWall->getLowerShadingState()); }
   }

   // === Case during presence and at departure =================================================
   else
   {    problower = exp(a01int_LB+b01inint_LB*Lumint+b01sint_LB*(pWall->getLowerShadingState()))/
        (1.f+exp(a01int_LB+b01inint_LB*Lumint+b01sint_LB*(pWall->getLowerShadingState())));
        probraise = exp(a10int_LB+b10inint_LB*Lumint+b10sint_LB*(pWall->getLowerShadingState()))/
        (1.f+exp(a10int_LB+b10inint_LB*Lumint+b10sint_LB*(pWall->getLowerShadingState())));
       if (pWall->getLowerShadingState() == 1.f) { probraise = 0.f; }
       if (pWall->getLowerShadingState() == 0.f) { problower = 0.f; }

   // --- Case of problower >= probraise ---------------------------------------------------------
       if (problower >= probraise)
       {   if (randomUniform(0.f,1.f) < problower) { ActLow = -1; }
           else if (randomUniform(0.f,1.f) < probraise) { ActLow = +1; }
           else { ActLow = 0; }
       }

   // --- Case of problower < probraise ------------------------------------------------------
       else
       {   if (randomUniform(0.f,1.f) < probraise) { ActLow = +1; }
           else if (randomUniform(0.f,1.f) < problower) { ActLow = -1; }
           else { ActLow = 0; }
       }        // --- Choice of the new unshaded fraction ------------------------------------------------
       if (ActLow == -1) // If a lowering action was predicted
       {   ptotlow = exp(afulllower_LB+boutfulllower_LB*Evg+bsfulllower_LB*(pWall->getLowerShadingState()))/
           (1.f+exp(afulllower_LB+boutfulllower_LB*Evg+bsfulllower_LB*(pWall->getLowerShadingState())));
           if (randomUniform(0.f,1.f) < ptotlow) { pWall->setLowerShadingState(0.f); }
           else
           {   float Reduction = randomWeibull(exp(aSFlower + bSFlower*(pWall->getLowerShadingState())), shapelower);
               pWall->setLowerShadingState( 0.01f*round(100.f*max((pWall->getLowerShadingState())-Reduction,0.01f)) );
           }
       }

       else if (ActLow == +1) // If a raising action was predicted
       {   ptotraise = exp(afullraise_LB+boutfullraise_LB*Evg+bsfullraise_LB*(pWall->getLowerShadingState()))/
           (1.f+exp(afullraise_LB+boutfullraise_LB*Evg+bsfullraise_LB*(pWall->getLowerShadingState())));
           if (randomUniform(0.f,1.f) < ptotraise) { pWall->setLowerShadingState(1.f); }
           else pWall->setLowerShadingState( 0.01f*round(100.f*randomUniform((pWall->getLowerShadingState()),1.f)) );
       }

       else // If no action was predicted
       {    pWall->setLowerShadingState(pWall->getLowerShadingState()); }

    }

    return;

}

void Model::upperShadingAction(Zone* pZone, Wall* pWall, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour) {
   // This function simulates the occurence of actions on upper shading devices (covering
   // anidolic systems). If an action is simulated, the chosen unshaded fraction is then predicted.
   // Input parameters are: indoor workplane illuminance (Lumint) and outdoor global horizontal
   // illuminance (Evg).

   // Lumint: Indoor illuminance next to window (lux)
   float Lumint = pWall->getTotalInternalIlluminance0()*pWall->getUpperShadingState();

   // Evg: Outdoor illuminance (lux)
   float Evg = pClimate->getIgh_vis();

   // --- Constants for upper blinds (covering the anidolic system) ----------------------------------
   // XML ENTRIES TO BE ADDED

   // Probability of lowering on arrival
   float a01arr_UB    = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_UB(0);
   float b01inarr_UB  = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_UB(1);
   float b01outarr_UB  = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_UB(2);
   float b01sarr_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerarr_UB(3);
   // a01arr_UB <- -7.294; b01inarr_UB <-  9.481E-04; b01outarr_UB <-  6.66E-06; b01sarr_UB <-   2.176

   // Probability of raising on arrival
   float a10arr_UB   = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_UB(0);
   float b10inarr_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_UB(1);
   float b10outarr_UB  = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_UB(2);
   float b10sarr_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraisearr_UB(3);
   // a10arr_UB <- -1.699; b10inarr_UB <- -5.236E-04; b10outarr_UB <- -2.182E-05; b10sarr_UB <- -3.916

   // Probability of lowering during presence
   float a01int_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_UB(0);
   float b01inint_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_UB(1);
   float b01outint_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_UB(2);
   float b01sint_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPlowerint_UB(3);
   // a01int_UB <- -8.211; b01inint_UB <-  8.343E-04; b01outint_UB <-  5.692E-06; b01sint_UB <-  1.533

   // Probability of raising during presence
   float a10int_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_UB(0);
   float b10inint_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_UB(1);
   float b10outint_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_UB(2);
   float b10sint_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPraiseint_UB(3);
   // a10int_UB <- -3.629; b10inint_UB <- -2.899E-04; b10outint_UB <- -1.686E-05; b10sint_UB <- -3.365

   // Probability of full lowering
   float afulllower_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfulllower_UB(0);
   float boutfulllower_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfulllower_UB(1);
   float bsfulllower_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfulllower_UB(2);
   // afulllower_UB <- -0.4798; boutfulllower_UB <- -3.854e-07 ; bsfulllower_UB <- -1.848

   // Probability of full raising
   float afullraise_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfullraise_UB(0);
   float boutfullraise_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfullraise_UB(1);
   float bsfullraise_UB = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticBlindsParameters()->getPfullraise_UB(2);
   // afullraise_UB <- 1.032; boutfullraise_UB <- -1.324e-05; bsfullraise_UB <- -0.2356


   // --- Simulation of upper blinds -------------------------------------------------------------

    float problower,probraise,ptotlow,ptotraise;
    int ActUpp;

   // === Case of absence ========================================================================
   if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {
       pWall->setUpperShadingState(pWall->getUpperShadingState());
   }

   // === Case of arrival ========================================================================
   else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {    problower = exp(a01arr_UB+b01inarr_UB*Lumint+b01outarr_UB*Evg+b01sarr_UB*(pWall->getUpperShadingState()))/
        (1.f+exp(a01arr_UB+b01inarr_UB*Lumint+b01outarr_UB*Evg+b01sarr_UB*(pWall->getUpperShadingState())));
        probraise = exp(a10arr_UB+b10inarr_UB*Lumint+b10outarr_UB*Evg+b10sarr_UB*(pWall->getUpperShadingState()))/
        (1.f+exp(a10arr_UB+b10inarr_UB*Lumint+b10outarr_UB*Evg+b10sarr_UB*(pWall->getUpperShadingState())));
       if (pWall->getUpperShadingState() == 1.f) { probraise = 0; }
       if (pWall->getUpperShadingState() == 0.f) { problower = 0; }

       // --- Case of problower >= probraise -----------------------------------------------------
       if (problower >= probraise)
       {   if (randomUniform(0.f,1.f) < problower) { ActUpp = -1; }
           else if (randomUniform(0.f,1.f) < probraise) { ActUpp = +1; }
           else { ActUpp = 0; }
       }

       // --- Case of problower < probraise ------------------------------------------------------
       else
       {   if (randomUniform(0.f,1.f) < probraise) { ActUpp = +1; }
           else if (randomUniform(0.f,1.f) < problower) { ActUpp = -1; }
           else { ActUpp = 0; }
       }

       // --- Choice of the new unshaded fraction ------------------------------------------------
       if (ActUpp == -1) // If a lowering action was predicted
       {   ptotlow = exp(afulllower_UB+boutfulllower_UB*Evg+bsfulllower_UB*(pWall->getUpperShadingState()))/
           (1.f+exp(afulllower_UB+boutfulllower_UB*Evg+bsfulllower_UB*(pWall->getUpperShadingState())));
           if (randomUniform(0.f,1.f) < ptotlow) { pWall->setUpperShadingState(0.f); }
           else pWall->setUpperShadingState( 0.01f*round(100.f*randomUniform(0.f, (pWall->getUpperShadingState()))) );
       }

       else if (ActUpp == +1) // If a raising action was predicted
       {    ptotraise = exp(afullraise_UB+boutfullraise_UB*Evg+bsfullraise_UB*(pWall->getUpperShadingState()))/
            (1.f+exp(afullraise_UB+boutfullraise_UB*Evg+bsfullraise_UB*(pWall->getUpperShadingState())));
            if (randomUniform(0.f,1.f) < ptotraise) { pWall->setUpperShadingState(1.f); }
            else pWall->setUpperShadingState( 0.01f*round(100.f*randomUniform((pWall->getUpperShadingState()),1.f)) );
       }

       else // If no action was predicted
       {    pWall->setUpperShadingState(pWall->getUpperShadingState()); }
   }

   // === Case during presence and at departure =================================================
   else
   {    problower = exp(a01int_UB+b01inint_UB*Lumint+b01outint_UB*Evg+b01sint_UB*(pWall->getUpperShadingState()))/
        (1.f+exp(a01int_UB+b01inint_UB*Lumint+b01outint_UB*Evg+b01sint_UB*(pWall->getUpperShadingState())));
        probraise = exp(a10int_UB+b10inint_UB*Lumint+b10outint_UB*Evg+b10sint_UB*(pWall->getUpperShadingState()))/
        (1.f+exp(a10int_UB+b10inint_UB*Lumint+b10outint_UB*Evg+b10sint_UB*(pWall->getUpperShadingState())));
        if (pWall->getUpperShadingState() == 1.f) { probraise = 0.f; };
        if (pWall->getUpperShadingState() == 0.f) { problower = 0.f; };

   // --- Case of problower >= probraise ---------------------------------------------------------
       if (problower >= probraise)
       {   if (randomUniform(0.f,1.f) < problower) { ActUpp = -1; }
           else if (randomUniform(0.f,1.f) < probraise) { ActUpp = +1; }
           else { ActUpp = 0; }
       }

   // --- Case of problower < probraise ------------------------------------------------------
       else
       {   if (randomUniform(0.f,1.f) < probraise) { ActUpp = +1; }
           else if (randomUniform(0.f,1.f) < problower) { ActUpp = -1; }
           else { ActUpp = 0; }
       }
   // --- Choice of the new unshaded fraction ------------------------------------------------
       if (ActUpp == -1) // If a lowering action was predicted
       {   ptotlow = exp(afulllower_UB+boutfulllower_UB*Evg+bsfulllower_UB*(pWall->getUpperShadingState()))/
           (1.f+exp(afulllower_UB+boutfulllower_UB*Evg+bsfulllower_UB*(pWall->getUpperShadingState())));
           if (randomUniform(0.f,1.f) < ptotlow) { pWall->setUpperShadingState(0.f); }
           else pWall->setUpperShadingState( 0.01f*round(100.f*randomUniform(0.f,(pWall->getUpperShadingState()))) );
       }

       else if (ActUpp == +1) // If a raising action was predicted
       {   ptotraise = exp(afullraise_UB+boutfullraise_UB*Evg+bsfullraise_UB*(pWall->getUpperShadingState()))/
           (1.f+exp(afullraise_UB+boutfullraise_UB*Evg+bsfullraise_UB*(pWall->getUpperShadingState())));
           if (randomUniform(0.f,1.f) < ptotraise) { pWall->setUpperShadingState(1.f); }
           else pWall->setUpperShadingState( 0.01f*round(100.f*randomUniform((pWall->getUpperShadingState()),1.f)) );
       }

       else // If no action was predicted
       {    pWall->setUpperShadingState(pWall->getUpperShadingState()); }

   }
   return;
}

/// TODO: check these new models

float Model::randomWeibull(float scale, float shape) {
   // Draws a random number from a Weibull distribution,
   // proceeding by inversion of the Weibull cdf, defined as
   // F(x) = 1 - exp(-(x/scale)^shape)

   float x = randomUniform(0.f,1.f);
   float y = scale*pow(-log(1.f-x),1.f/shape);

   return y;

}

void Model::windowAction_Markov(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour) {

    // ======================================================================================
    //    SIMULATION FENETRES AVEC MODELE MARKOVIEN (Toutes variables)
    // ======================================================================================

    // Input parameters:
    // Ti: Indoor temperature (∞C) (double)
    float Ti = pZone->getTaExpl();
    // Te: Outdoor temperature (∞C) (double)
    float Te = pClimate->getToutCelsius(day,hour);
    // Tdm: Daily mean outdoor temperature (∞C) (double)
    float Tdm = pClimate->getDailyMeanTemperature(day,hour);
    // Rain: Rain presence (bool)
    bool Rain = pClimate->getRainPresence(day,hour);
    // Occdur: Current occupancy (or absence) duration (min) (double)

    // occtime: Vector of occupied times from origin (seconds) (double)

    // floor: Floor of the simulated zone (double)
    bool groundFloor = pZone->getGroundFloor();

    // --- Regression parameters -----------------------------------------------------------------------
    // Pour l'instant, fixÈs. A terme, ces paramËtres devront Ítre pris comme inputs dans cette
    // fonction WindowAction

    // P01arr
    float a01arr    = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(0);
    float b01inarr  = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(1);
    float b01outarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(2);
    float b01absprevarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(3);
    float b01rnarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(4);

    // P10arr
    float a10arr   = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10arr(0);
    float b10inarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10arr(1);
    float b10outarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10arr(2);

    // P01int
    float a01int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(0);
    float b01inint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(1);
    float b01outint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(2);
    float b01presint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(3);
    float b01rnint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(4);

    // P10int
    float a10int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(0);
    float b10inint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(1);
    float b10outint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(2);
    float b10presint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(3);

    // P01dep
    float a01dep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(0);
    float b01outdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(1);
    float b01absdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(2);
    float b01gddep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(3);

    // P10dep
    float a10dep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(0);
    float b10indep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(1);
    float b10outdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(2);
    float b10absdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(3);
    float b10gddep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(4);

    // --- Simulation ----------------------------------------------------------------------------------

    // === Cas Absent ==========================================================================
    if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
    {
        pZone->setWindowState(pZone->getWindowState());
    }

    // === Cas ArrivÈe =========================================================================
    else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
      // --- Si dÈj‡ fermÈ ‡ l'arrivÈe ---------------------------------------------------------
    {   if (pZone->getWindowState() == 0.f) // winSim[i-1] == 0
        {   float prob01arr = exp(a01arr+b01inarr*Ti+b01outarr*Te+b01rnarr*(Rain?1.f:0.f)+b01absprevarr*( (pZone->getOccupants()->getCurrentDuration(day,hour,static_cast<int>(fracHour)-1) > 8.f*60.f*60.f) ? 1.f : 0.f ) ) /
            (1.f+exp(a01arr+b01inarr*Ti+b01outarr*Te+b01rnarr*(Rain?1.f:0.f)+b01absprevarr*( (pZone->getOccupants()->getCurrentDuration(day,hour,static_cast<int>(fracHour)-1) > 8.f*60.f*60.f) ? 1.f : 0.f )));
            if (randomUniform(0.f,1.f) < prob01arr) {pZone->setWindowState(1.f); /*WinSim[i] = 1.;*/} else {pZone->setWindowState(0.f); /*[i] = 0.;*/}
        }
      // --- Si dÈj‡ ouvert ‡ l'arrivÈe --------------------------------------------------------
       else
        {   float prob10arr = exp(a10arr+b10inarr*Ti+b10outarr*Te)/(1.f+exp(a10arr+b10inarr*Ti+b10outarr*Te));
            if (randomUniform(0.f,1.f) < prob10arr) {pZone->setWindowState(0.f); /*WinSim[i] = 0.;*/} else {pZone->setWindowState(1.f); /*WinSim[i] = 1.;*/}
        }
    }

    // === Cas intermÈdiaire =====================================================================
    else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 1.f)
       // --- Si fermÈ au pas de temps prÈcÈdent ------------------------------------------------
    {
        if (pZone->getWindowState() == 0.f) // winSim[i-1] == 0
        {   float prob01int = exp(a01int+b01inint*Ti+b01outint*Te+b01presint*pZone->getOccupants()->getCurrentDuration(day,hour,fracHour)+b01rnint*(Rain?1.f:0.f))/
            (1.f+exp(a01int+b01inint*Ti+b01outint*Te+b01presint*pZone->getOccupants()->getCurrentDuration(day,hour,fracHour)+b01rnint*(Rain?1.f:0.f)));
            if (randomUniform(0.f,1.f) < prob01int) {pZone->setWindowState(1.f); /*WinSim[i] = 1.;*/} else {pZone->setWindowState(0.f); /*WinSim[i] = 0.;*/}
        }
       // --- Si ouvert au pas de temps prÈcÈdent -----------------------------------------------
        else
        {   float prob10int = exp(a10int+b10inint*Ti+b10outint*Te+b10presint*pZone->getOccupants()->getCurrentDuration(day,hour,fracHour))/
            (1.f+exp(a10int+b10inint*Ti+b10outint*Te+b10presint*pZone->getOccupants()->getCurrentDuration(day,hour,fracHour)));
            if (randomUniform(0.f,1.f) < prob10int) {pZone->setWindowState(0.f); /*WinSim[i] = 0.;*/} else {pZone->setWindowState(1.f); /*WinSim[i] = 1.;*/}
        }
    }

     // === Cas dÈpart ==========================================================================
    else
     // --- Si courte absence -----------------------------------------------------------------
    {
        if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 8.f*60.f*60.f )
        {   if (pZone->getWindowState() == 0.f) // winSim[i-1] == 0  // Si fermÈ au pas de temps prÈcÈdent
            {   float prob01depshort = exp(a01dep+b01outdep*Tdm+b01gddep*(groundFloor ? 0.f : 1.f))/
                (1.f+exp(a01dep+b01outdep*Tdm+b01gddep*(groundFloor ? 0.f : 1.f))); // as.numeric(floor>0)
                if (randomUniform(0.f,1.f) < prob01depshort) {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/} else {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/}
            }
            else // Si ouvert au pas de temps prÈcÈdent
            {   float prob10depshort = exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10gddep*(groundFloor ? 0.f : 1.f))/
                (1.f+exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10gddep*(groundFloor ? 0.f : 1.f))); // as.numeric(floor>0)
                if (randomUniform(0.f,1.f) < prob10depshort) {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/} else {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/}
            }
        }
     // --- Si longue absence -----------------------------------------------------------------
        else // If occupant leaves for more than 8 hours
        {   if (pZone->getWindowState() == 0.f) // WinSim[i-1]==0  // Si fermÈ au pas de temps prÈcÈdent
            {   float prob01deplong = exp(a01dep+b01outdep*Tdm+b01absdep+b01gddep*(groundFloor ? 0.f : 1.f))/
                (1.f+exp(a01dep+b01outdep*Tdm+b01absdep+b01gddep*(groundFloor ? 0.f : 1.f)));
                if (randomUniform(0.f,1.f) < prob01deplong) {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/} else {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/}
            }
            else // Si ouvert au pas de temps prÈcÈdent
            {   float prob10deplong = exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10absdep+b10gddep*(groundFloor ? 0.f : 1.f))/
                (1.f+exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10absdep+b10gddep*(groundFloor ? 0.f : 1.f)));
                if (randomUniform(0.f,1.f) < prob10deplong) {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/} else {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/}
            }

        }
    }

    //cerr << "Day: " << day << "\tHour: " << hour << "\tfracHour: " << fracHour << "\twindowState: " << pZone->getWindowState() << endl;

    return;

}

void Model::windowAction_Hybrid(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour) {

   // Function to call when the model choice is "Hybrid"

   // ======================================================================================
   //    SIMULATION FENETRES AVEC HYBRIDE (Toutes variables)
   // ======================================================================================

   float timestep = 5.f; // -------------------- NEW --------------------------------------
   // Input parameters:
   // Ti: Indoor temperature (∞C) (double)
   float Ti = pZone->getTaExpl();
   // Te: Outdoor temperature (∞C) (double)
   float Te = pClimate->getToutCelsius(day,hour);
   // Tdm: Daily mean outdoor temperature (∞C) (double)
   float Tdm = pClimate->getDailyMeanTemperature(day,hour);
   // Rain: Rain presence (bool)
   bool Rain = pClimate->getRainPresence(day,hour);
   // Occdur: Current occupancy (or absence) duration (min) (double)

   // occtime: Vector of occupied times from origin (seconds) (double)

   // floor: Floor of the simulated zone (double)
   bool groundFloor = pZone->getGroundFloor();

   // --- Regression parameters -----------------------------------------------------------------------

   // P01arr
   float a01arr    = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(0);
   float b01inarr  = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(1);
   float b01outarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(2);
   float b01absprevarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(3);
   float b01rnarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01arr(4);

   // P10arr
   //float a10arr   = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10arr(0);
   //float b10inarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10arr(1);
   //float b10outarr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10arr(2);

   // P01int
   float a01int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(0);
   float b01inint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(1);
   float b01outint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(2);
   float b01presint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(3);
   float b01rnint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01int(4);

   // P10int
   //float a10int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(0);
   //float b10inint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(1);
   //float b10outint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(2);
   //float b10presint = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10int(3);

   // P01dep
   float a01dep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(0);
   float b01outdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(1);
   float b01absdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(2);
   float b01gddep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP01dep(3);

   // P10dep
   float a10dep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(0);
   float b10indep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(1);
   float b10outdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(2);
   float b10absdep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(3);
   float b10gddep = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getP10dep(4);


// -------------------- NEW --------------------------------------
   // duropen
   float aop = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getDuropen(0);
   float bopout = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getDuropen(1);
   float shapeop = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticWindowParameters()->getDuropen(2);
   // aop="2.151" bopout="0.172"  shapeop="2.39"
   float duropen = 0.f; //!< Opening duration (minutes)

   // --- Simulation ----------------------------------------------------------------------------------

   // === Cas Absent ==========================================================================
   if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {
       pZone->setWindowState(pZone->getWindowState());
   }

   // === Cas ArrivÈe =========================================================================
   else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
     // --- Si dÈj‡ fermÈ ‡ l'arrivÈe ---------------------------------------------------------
   {   if (pZone->getWindowState() == 0.f) // winSim[i-1] == 0
       {   float prob01arr = exp(a01arr+b01inarr*Ti+b01outarr*Te+b01rnarr*(Rain?1.f:0.f)+b01absprevarr*( (pZone->getOccupants()->getCurrentDuration(day,hour,static_cast<int>(fracHour)-1) > 8.f*60.f*60.f) ? 1.f : 0.f ) ) /
           (1.f+exp(a01arr+b01inarr*Ti+b01outarr*Te+b01rnarr*(Rain?1.f:0.f)+b01absprevarr*( (pZone->getOccupants()->getCurrentDuration(day,hour,static_cast<int>(fracHour)-1) > 8.f*60.f*60.f) ? 1.f : 0.f )));
           if (randomUniform(0.f,1.f) < prob01arr) { duropen = randomWeibull(exp(aop+bopout*Te), shapeop); pZone->setWindowState(1.f); } else { pZone->setWindowState(0.f); }
       }
     // --- Si dÈj‡ ouvert ‡ l'arrivÈe --------------------------------------------------------
       else // ------------------- NEW ----------------------------------------
       {   duropen = randomWeibull(exp(aop+bopout*Te), shapeop);
           if (duropen<timestep) { pZone->setWindowState(0.f); duropen=0.f; }
           else { pZone->setWindowState(1.f); duropen = duropen-timestep; }
       }
   }


   // === Cas intermÈdiaire =====================================================================
   else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 1.f)
      // --- Si fermÈ au pas de temps prÈcÈdent ------------------------------------------------
   {
       if (pZone->getWindowState() == 0.f) // winSim[i-1] == 0
       {   float prob01int = exp(a01int+b01inint*Ti+b01outint*Te+b01presint*pZone->getOccupants()->getCurrentDuration(day,hour,fracHour)+b01rnint*(Rain?1.f:0.f))/
           (1.f+exp(a01int+b01inint*Ti+b01outint*Te+b01presint*pZone->getOccupants()->getCurrentDuration(day,hour,fracHour)+b01rnint*(Rain?1.f:0.f)));
           if (randomUniform(0.f,1.f) < prob01int) { duropen = randomWeibull(exp(aop+bopout*Te), shapeop); pZone->setWindowState(1.f); } else { pZone->setWindowState(0.f); }
       }
      // --- Si ouvert au pas de temps prÈcÈdent -----------------------------------------------
       else
       {
           if (duropen<timestep) { pZone->setWindowState(0.f); duropen=0.f; }
           else { pZone->setWindowState(1.f); duropen = duropen-timestep; }
       }
   }

    // === Cas dÈpart ==========================================================================
   else
    // --- Si courte absence -----------------------------------------------------------------
   {
       if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 8.f*60.f*60.f )
       {   if (pZone->getWindowState() == 0.f) // winSim[i-1] == 0  // Si fermÈ au pas de temps prÈcÈdent
           {   float prob01depshort = exp(a01dep+b01outdep*Tdm+b01gddep*(groundFloor ? 0.f : 1.f))/
               (1.f+exp(a01dep+b01outdep*Tdm+b01gddep*(groundFloor ? 0.f : 1.f))); // as.numeric(floor>0)
               if (randomUniform(0.f,1.f) < prob01depshort) {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/} else {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/}
           }
           else // Si ouvert au pas de temps prÈcÈdent
           {   float prob10depshort = exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10gddep*(groundFloor ? 0.f : 1.f))/
               (1.f+exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10gddep*(groundFloor ? 0.f : 1.f))); // as.numeric(floor>0)
               if (randomUniform(0.f,1.f) < prob10depshort) {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/} else {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/}
           }
       }
    // --- Si longue absence -----------------------------------------------------------------
       else // If occupant leaves for more than 8 hours
       {   if (pZone->getWindowState() == 0.f) // WinSim[i-1]==0  // Si fermÈ au pas de temps prÈcÈdent
           {   float prob01deplong = exp(a01dep+b01outdep*Tdm+b01absdep+b01gddep*(groundFloor ? 0.f : 1.f))/
               (1.f+exp(a01dep+b01outdep*Tdm+b01absdep+b01gddep*(groundFloor ? 0.f : 1.f)));
               if (randomUniform(0.f,1.f) < prob01deplong) {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/} else {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/}
           }
           else // Si ouvert au pas de temps prÈcÈdent
           {   float prob10deplong = exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10absdep+b10gddep*(groundFloor ? 0.f : 1.f))/
               (1.f+exp(a10dep+b10indep*Ti+b10outdep*Tdm+b10absdep+b10gddep*(groundFloor ? 0.f : 1.f)));
               if (randomUniform(0.f,1.f) < prob10deplong) {pZone->setWindowState(0.f); /*winSim[i] = 0.;*/} else {pZone->setWindowState(1.f); /*winSim[i] = 1.;*/}
           }

       }
   }

   //cerr << "Day: " << day << "\tHour: " << hour << "\tfracHour: " << fracHour << "\twindowState: " << pZone->getWindowState() << endl;

   return;

}


void Model::windowAction_Bernoulli(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour) {

   // Function to call when the model choice is "Bernoulli"

   // ======================================================================================
   //    SIMULATION FENETRES AVEC METHODE BERNOULLI
   // ======================================================================================

   // Ti: Indoor temperature (∞C) (double)
   float Ti = pZone->getTaExpl();
   // Te: Outdoor temperature (∞C) (double)
   float Te = pClimate->getToutCelsius(day,hour);

   // Parameters for the probability to observe a window open
   float a   = 0.794f;
   float bTi = -0.1541f;
   float bTe = 0.1476f;

   // --- Simulation ----------------------------------------------------------------------------------

   // === Cas Absent ==========================================================================
   if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {
       pZone->setWindowState(pZone->getWindowState());
   }

   // === Cas PrÈsent =========================================================================
   else
   {
       float Prob = exp(a+bTe*Te+bTi*Ti)/(1.f+exp(a+bTe*Te+bTi*Ti));
       if (randomUniform(0.f,1.f) < Prob) { pZone->setWindowState(1.f); }
   }
   return;
}


void Model::windowAction_Humphreys(Zone* pZone, Climate* pClimate, unsigned int day, unsigned int hour, unsigned int fracHour) {

    // Function to call when the model choice is "Humphreys"
    // The algorithm is described in detail in H. B. Rijal et al, Energy and Buildings 39 (2007) 823-836
    // It was originally developed to be run at a time step of 1 hour.

    // ======================================================================================
    //    SIMULATION FENETRES AVEC ALGORITHME DE HUMPHREYS
    // ======================================================================================

    // Ti: Indoor temperature (∞C) (double)
    float Ti = pZone->getTaExpl();
    // Te: Outdoor temperature (∞C) (double)
    float Te = pClimate->getToutCelsius(day,hour);
    // Trm: Exponentially weighted running mean outdoor temperature (∞C) (double)
    float Trm = (pClimate->getDailyMeanTemperature(((day-1-1)%365)+1,hour)
                + 0.8*pClimate->getDailyMeanTemperature(((day-2-1)%365)+1,hour)
                + 0.6*pClimate->getDailyMeanTemperature(((day-3-1)%365)+1,hour)
                + 0.5*pClimate->getDailyMeanTemperature(((day-4-1)%365)+1,hour)
                + 0.4*pClimate->getDailyMeanTemperature(((day-5-1)%365)+1,hour)
                + 0.3*pClimate->getDailyMeanTemperature(((day-6-1)%365)+1,hour)
                + 0.2*pClimate->getDailyMeanTemperature(((day-7-1)%365)+1,hour))/3.8;
   // Tcomf: Indoor comfort temperature according to the CEN standard
   float Tcomf;
   if (Trm > 10.f) { Tcomf = 0.33f*Trm+18.8f; }
   else { Tcomf = 0.09f*Trm+22.6f; }

   // Parameters for the probability to observe a window open (from SCATs data)
   float a   = -6.43f;
   float bTi = 0.171f;
   float bTe = 0.166f;

   // --- Simulation ----------------------------------------------------------------------------------

   // === Cas Absent ==========================================================================
   if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
   {
       pZone->setWindowState(pZone->getWindowState());
   }

   // === Cas PrÈsent =========================================================================
   else
   {   // If the window is closed and the occupant is too hot (Tin > Tcomf + 2∞C)
       if ( (pZone->getWindowState() == 0.f) & (Ti-Tcomf > 2.f) )
       {   float prob = exp(a+bTi*Ti+bTe*Te)/(1.f+exp(a+bTi*Ti+bTe*Te));
           if (randomUniform(0.f,1.f) < prob) { pZone->setWindowState(1.f); }
       }
       // If the window is open and the occupant is too cold  (Tin < Tcomf - 2∞C)
       else if ( (pZone->getWindowState() == 1.f) & (Ti < Tcomf-2.f) )
       {   float prob = exp(a+bTi*Ti+bTe*Te)/(1.f+exp(a+bTi*Ti+bTe*Te));
           if (randomUniform(0.f,1.f) > prob) { pZone->setWindowState(0.f); }
       }
       else { pZone->setWindowState(pZone->getWindowState()); }
   }
   return;
}

void Model::lightAction_Lightswitch2002(Zone* pZone, unsigned int day, unsigned int hour, unsigned int fracHour) {

    // SIMULATION OF ACTIONS ON ELECTRICAL LIGHTING ACCORDING TO LIGHTSWITCH-2002 MODEL
    // Christoph F. Reinhart, Lightswitch-2002: a model for manual and automated control
    // of electric lighting and blinds, Solar Energy 77 (2004) 15-28

    // Input parameters:
    // Lumint: Indoor illuminance next to window (lux)

    // *** computing the inside illuminance according as a sum of all walls' contributions *** //
    float Lumint = pZone->getTotalInternalIlluminance0();

    // =================================================================================================
    //    LOADING CALIBRATION PARAMETERS FROM THE XML FILE
    // =================================================================================================

    // --- Probability of switch-on (arrival) -----------
    float a01arr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonarr(0);
    float b01arr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonarr(1);
    float c01arr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonarr(2);
    float m01arr = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonarr(3);
    // Default values (Reinhart, 2004): a = -0.0175, b = -4.0835, c = 1.0361, m = 1.8223

    // --- Probability of switch-on (during presence) ---
    float a01int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonint(0);
    float b01int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonint(1);
    float c01int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonint(2);
    float d01int = ((StochasticOccupantsPresence*)pZone->getOccupants())->getStochasticLightsParameters()->getPonint(3);
    // Default values (Reinhart, 2004): a = 0.0027, b=0.017, c=-64.19, d=2.41


    // =================================================================================================
    //    SIMULATION PROCESS
    // =================================================================================================

    // === 1. Case of absence ==========================================================================

    if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 0.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
    {
        pZone->setLightsState(pZone->getLightsState());
    }

    // === 2. Case of arrival =========================================================================

    else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 0.f)
      // --- 2a. If the light is already off when occupant arrives -------------------------------------------------
    {   if (pZone->getLightsState() == 0.f)
        {
            float probonarr;
            if (log10(Lumint) <= 0.843f) probonarr = 1.f;
            else if (log10(Lumint) >= 2.818f) probonarr = 0.f;
            else probonarr = a01arr + c01arr/(1.f + exp(-b01arr*(log10(Lumint - m01arr))));
            // usage of the probability to define the lights state
            if (randomUniform(0.f,1.f) < probonarr) pZone->setLightsState(1.f);
            else pZone->setLightsState(0.f);
        }
      // --- 2b. If the light is already on when occupant arrives --------------------------------------------------------
       else
        {
            pZone->setLightsState(pZone->getLightsState());
        }
    }

    // === 3. Case during presence =====================================================================

    else if ( pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour) == 1.f && pZone->getOccupants()->getOccupantsFraction(day,hour,fracHour-1) == 1.f)

       // --- 3a. If the light is already off at the previous time step ------------------------------------------------
    {
        if (pZone->getLightsState() == 0.f)
        {
            float probonint;
            if (Lumint == 0.f) probonint = 1.f;
            else probonint = a01int+b01int/(c01int*(log10(Lumint-d01int)));
            // usage of the probability to define the lights state
            if (randomUniform(0.f,1.f) < probonint) pZone->setLightsState(1.f);
            else pZone->setLightsState(0.f);
        }
       // --- 3a. If the light is already on at the previous time step -----------------------------------------------
        else
        {
            pZone->setLightsState(pZone->getLightsState());
        }
    }

     // === 4. Case of departure ==========================================================================

    else
    {
        // --- 4a. If the light is already off when occupant leaves --------------------------------------
        if (pZone->getLightsState() == 0.f) // Si dÈj‡ Èteint au pas de temps prÈcÈdent
        {
            pZone->setLightsState(pZone->getLightsState());
        }

        // --- 4b. If the light is already on when occupant leaves --------------------------------------
        else
        {
            float proboffdep;
            if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 30.f*60.f )              proboffdep = 0.086f;
            else if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 60.f*60.f )         proboffdep = 0.314f;
            else if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 2.f*60.f*60.f )     proboffdep = 0.380f;
            else if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 4.f*60.f*60.f )     proboffdep = 0.600f;
            else if ( pZone->getOccupants()->getFutureDuration(day,hour,fracHour) < 12.f*60.f*60.f )    proboffdep = 0.960f;
            else proboffdep = 1.f;
            // NOTE: These probabilities are based on Pigg, Eilers, Reed, Behavioral Aspects of Lighting
            //  and Occupancy Sensors in Private Offices: A Case Study of a University Office Building
            // usage of the probability to define the lights state
            if (randomUniform(0.f,1.f) < proboffdep) pZone->setLightsState(0.f);
            else pZone->setLightsState(1.f);
        }
    }

    //cerr << "Day: " << day << "\tHour: " << hour << "\tfracHour: " << fracHour << "\twindowState: " << pZone->getWindowState() << endl;

    return;

}

void Model::lightAction_Threshold(Zone* pZone) {

    // if under threshold then lights ON otherwise lights OFF
    if (pZone->getTotalInternalIlluminance0() < pZone->getLightsThreshold()) pZone->setLightsState(1.f);
    else pZone->setLightsState(0.f);

    return;

}

float Model::lightsElectricConsumption(Zone* pZone) {

    // suppose a perfect dimming system
    if (pZone->getLightsState() > 0.f) {
        // lights ON
        // suppose a linear power dimming for the lights
        return (1.f - max(pZone->getTotalInternalIlluminance0()/pZone->getLightsThreshold(),1.f))*pZone->getLightsPowerDensity()*pZone->getFloorArea();
    }
    else return 0.f;

/*
    // Model based on the CIBSE Code for Lighting 2009, Table 2.5

    // The table needs the lamp type, the task illuminance and gives the average installed power density
    float illum = pZone->getWall(0)->getTotalInternalIlluminance0(); // Indoor illuminance next to window (lux)
    float taskIllum = 300.f; //pZone->getTaskIlluminance(); /// to be defined
    float addIllum = max(taskIllum-illum,0.f);
    float powerDensity = 2.27; //pZone->getLightsPowerDensity(); // average installed power density W/(m≤ 100 lux)

    return (addIllum/100.f)*powerDensity*pZone->getFloorArea();
*/

}

void Model::computeCMIndices(Climate* pClimate, Building* pBuilding, unsigned int i, unsigned int preTimeStepsSimulated, float &MRT, float &COMFA, float &ITS) {

    // compute the MRT for the building pBuilding and the timestep i
    // variables used with the estimation of the MRT
    float totalArea = 0.f;
    float LWabsorbed = 0.f;
    float SWabsorbed = 0.f;
    for (unsigned int zone=0; zone<pBuilding->getnZones();++zone) {
        // loop on the walls
        for (unsigned int k=0; k<pBuilding->getZone(zone)->getnWalls(); ++k) {
            Wall* thisWall = pBuilding->getZone(zone)->getWall(k);
            totalArea += thisWall->getArea();
            LWabsorbed += thisWall->getLongWaveAbsorbed(i) * thisWall->getArea();
            SWabsorbed += thisWall->getShortWaveIrradiance(i) * (1.f - thisWall->getShortWaveReflectance()) * thisWall->getArea();
        }
        // loop on the roofs
        for (unsigned int k=0; k<pBuilding->getZone(zone)->getnRoofs(); ++k) {
            Roof* thisRoof = pBuilding->getZone(zone)->getRoof(k);
            totalArea += thisRoof->getArea();
            LWabsorbed += thisRoof->getLongWaveAbsorbed(i) * thisRoof->getArea();
            SWabsorbed += thisRoof->getShortWaveIrradiance(i) * (1.f - thisRoof->getShortWaveReflectance()) * thisRoof->getArea();
        }
    } // end the loop on zones
    LWabsorbed /= totalArea;
    SWabsorbed /= totalArea;
    MRT = pow((LWabsorbed + SWabsorbed)/(pBuilding->getMRT_Epsilon()*5.670373e-8), 0.25f)-273.15f;
    // Renolds Number and co.
    float reynoldsNumber = 0.17*((pClimate->getWindSpeed(i-preTimeStepsSimulated)>0)?pClimate->getWindSpeed(i-preTimeStepsSimulated):0.01)/1.5e-5;
    float re_A = (reynoldsNumber<4000.)?0.683:((reynoldsNumber<40000.)?0.193:0.0266);
    float re_n = (reynoldsNumber<4000.)?0.466:((reynoldsNumber<40000.)?0.618:0.805);
    // output of the COMFA* (in W/m^2)

    // Metabolic
    float metabolicActivity = 70.f;
    float walkingSpeed =0.0052*(metabolicActivity-58);
    float M = (1.-(0.15-0.0173*pClimate->getVapourPressure(i-preTimeStepsSimulated)-0.0014*pClimate->getToutCelsius(i-preTimeStepsSimulated)))*metabolicActivity; // activité métabolique 70 W/m^2

    // Convection
    float coreTemperature = 36.5 + 0.0043*M;
    float airDensity = pClimate->getPatm(i-preTimeStepsSimulated)/(287.04*(pClimate->getToutCelsius(i-preTimeStepsSimulated)+273.15));
    float resistanceBodyTissue = 1000.*airDensity/((0.13*(0.42*(M-58.)))+15.);
    float boundaryLayerResistance = 0.17/(re_A*pow(reynoldsNumber,re_n)*pow(0.71,0.33)*22.e-6);
    float intrinsicClothingInsulation = (pClimate->getToutCelsius(i-preTimeStepsSimulated)>=27.f)?0.31:((1.372-(0.01866*pClimate->getToutCelsius(i-preTimeStepsSimulated))-(0.0004849*pow(pClimate->getToutCelsius(i-preTimeStepsSimulated),2))-(0.000009333*pow(pClimate->getToutCelsius(i-preTimeStepsSimulated),3)))*0.1555);
    float staticClothingResistance = 1000.*airDensity * intrinsicClothingInsulation;
    float clothingResistance = staticClothingResistance*(-0.37*(1.-exp(-walkingSpeed/0.72))+1.);
    float skinTemperature = ((((coreTemperature-pClimate->getToutCelsius(i-preTimeStepsSimulated))<0.)?0.001:(coreTemperature-pClimate->getToutCelsius(i-preTimeStepsSimulated))) / (resistanceBodyTissue+clothingResistance+boundaryLayerResistance)) * (boundaryLayerResistance+clothingResistance) + pClimate->getToutCelsius(i-preTimeStepsSimulated);
    float C = 1000.*airDensity * (skinTemperature-pClimate->getToutCelsius(i-preTimeStepsSimulated))/(clothingResistance+boundaryLayerResistance);

    // Evaporation
    float evaporativeHeatLossThroughPerspiration =0.42*(M-58.);
    float latentHeatOfVapourization =(2501.-(2.37*pClimate->getToutCelsius(i-preTimeStepsSimulated)))*1000.;
    float specificHumidityAir =0.622*(pClimate->getVapourPressure(i-preTimeStepsSimulated)/(pClimate->getPatm(i-preTimeStepsSimulated)/1000.-pClimate->getVapourPressure(i-preTimeStepsSimulated)));
    float saturatedVapourPressureSkin = pClimate->getSaturatedVapourPressure(skinTemperature);
    float vapourPressureSkin = saturatedVapourPressureSkin*pClimate->getRelativeHumidity(i-preTimeStepsSimulated);
    float specificHumiditySkin =0.622*(vapourPressureSkin/(pClimate->getPatm(i-preTimeStepsSimulated)/1000.-vapourPressureSkin));
    float airResistanceRav =0.92*boundaryLayerResistance;
    float staticClothingVapourResistance = 0.622*latentHeatOfVapourization*airDensity*intrinsicClothingInsulation*0.18/(pClimate->getPatm(i-preTimeStepsSimulated)/1000.-pClimate->getVapourPressure(i-preTimeStepsSimulated));
    float effectiveAirVelocity = pow(pow(pClimate->getWindSpeed(i-preTimeStepsSimulated),2)+pow(walkingSpeed,2),0.5);
    float resistanceVapourTransferClothing = staticClothingVapourResistance*(-0.8*(1.-exp((-effectiveAirVelocity/1.095)+1)));
    float resistanceAirClothingSkinTissue = resistanceVapourTransferClothing+airResistanceRav+7.7e+3;
    float evaporativeLossThroughSkinDiffusion = airDensity*latentHeatOfVapourization*( (((specificHumiditySkin-specificHumidityAir)/(resistanceVapourTransferClothing+airResistanceRav))>0.f)?
                                                                                       (min( (specificHumiditySkin-specificHumidityAir)/resistanceAirClothingSkinTissue
                                                                                           ,(specificHumiditySkin-specificHumidityAir)/(resistanceVapourTransferClothing+airResistanceRav))):
                                                                                       ((specificHumiditySkin-specificHumidityAir)/resistanceAirClothingSkinTissue) );
    float E = evaporativeLossThroughSkinDiffusion+evaporativeHeatLossThroughPerspiration;

    //LongwaveEmitted
    float surfaceTemperatureIndividual = max((skinTemperature-pClimate->getToutCelsius(i-preTimeStepsSimulated))/(boundaryLayerResistance+clothingResistance)*boundaryLayerResistance+pClimate->getToutCelsius(i-preTimeStepsSimulated),7.f);
    float longWaveEmitted = 0.95*5.67e-8*pow((surfaceTemperatureIndividual+273.15),4);

    // output of the COMFA (in W/m2)
    COMFA = M+ LWabsorbed + SWabsorbed - C - E - longWaveEmitted;

    // calculation of the ITS (in W)

    // Metabolic
    float metabolicRate = (metabolicActivity - (0.2 * (metabolicActivity - 80.)));

    // Convection
    float deltaT = pClimate->getToutCelsius(i-preTimeStepsSimulated) - 35.;
    float reynoldsNumberITS = 0.17*((pClimate->getWindSpeed(i-preTimeStepsSimulated)>0)?pClimate->getWindSpeed(i-preTimeStepsSimulated):0.01)/1.6e-5;
    float A_ITS = (reynoldsNumberITS>4000.)?0.17:0.62;
    float B_ITS = (reynoldsNumberITS<4000.)?0.62:0.47;
    float convectionITS = (deltaT * pClimate->getAirDensity()*1005.*2.e-5 * A_ITS * pow(reynoldsNumberITS,B_ITS) / 0.17);
    // simplified version
    //float convectionITS = (deltaT * 8.3 * pow(((pClimate->getWindSpeed(i-preTimeStepsSimulated)>0)?pClimate->getWindSpeed(i-preTimeStepsSimulated):0.01),0.6));

    // Radiation Balance
    float LWEmittedITS = 0.95*5.67e-8*pow((35.+273.15),4);
    float radiationITS = (LWabsorbed + SWabsorbed - LWEmittedITS);

    // Cooling Rate
    float coolingRate = (metabolicRate + radiationITS + convectionITS) * totalArea;

    // CoolingEfficiency
    float saturationVapourPressure = exp(16.6536-(4030.183/(pClimate->getToutCelsius(i-preTimeStepsSimulated)+235.)));
    float vapourPressureAir = 7.52*pClimate->getRelativeHumidity(i-preTimeStepsSimulated)*saturationVapourPressure;
    float evaporativeCapacityAir = 1.163*20.5*pow(((pClimate->getWindSpeed(i-preTimeStepsSimulated)>0.)?pClimate->getWindSpeed(i-preTimeStepsSimulated):0.01),0.3) * (42.-vapourPressureAir);
    float coolingEfficiencySweating = max(exp(0.6*((coolingRate/evaporativeCapacityAir)-0.12)),1.);

    // output of the ITS (in W)
    ITS = coolingRate * coolingEfficiencySweating;

    return;
}
