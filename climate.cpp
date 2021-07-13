#include "climate.h"
//#include <iomanip>
//#include <fstream>
#include "SKYSun.h"
#include "util.h"

#include <algorithm>
#include <iterator>

static const float NA=-999.0;

// *** Climate class, CitySim *** //
// *** jerome.kaempf@epfl.ch  *** //

Climate::Climate(string filename, ostream* pLogFileStream):logStream(std::cout.rdbuf()) {

    //logStream << "Climate constructor" << endl;
    // logStream is directed by default to the "cout" streambuf
    if(pLogFileStream!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogFileStream->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    // If no filename to load -> go back directly
    if (filename.empty()) {
        logStream << "WARNING: Climate file name empty." << endl << flush;
        return;

    }

    // variables used to store the elements
    char cbuffer[200];
    string buffer;

    // variables used in case only the global irradiance in the horizontal plane is given
	bool G_h_present = false;
	bool G_Dh_present = false;
	bool Tg_present = false;

    // Climate filename opening
    //throw(string("Test error"));
    logStream << "Climate file opening" << endl;
    //try{
        fstream input;

        input.open(filename.c_str(), ios::in | ios::binary);

        if (!(input.good())){
            logStream << "Bad climate file input" << endl;
            throw(string("Error opening climate file: " + filename));
        }
        //cout << "Climate: input is good" << endl;
        //logStream << "Loading: " << filename << endl << flush;

        // start the reading
        input.getline(cbuffer, 200, '\n');
        location = cbuffer;
        //cerr<<"Climate: "<<cbuffer<<endl << flush;
        input.getline(cbuffer, 200, ',');
        latitudeN = atof(cbuffer);
        //cerr<<"Latitude: "<<latitude<<"\t";
        if(latitudeN > 90.f || latitudeN < -90.f) throw(string("Climate file : wrong format latitude (degrees North)"));

        input.getline(cbuffer, 200, ',');
        longitudeE = atof(cbuffer);
        //logStream<<"Longitude: "<<longitudeE<<"\t";
        if(longitudeE > 180.f || longitudeE < -180.f) throw(string("Climate file : wrong format longitude (degrees East)"));

        input.getline(cbuffer, 200, ',');
        altitude = atof(cbuffer);
        //cerr<<"Altitude: "<<altitude<<endl << flush;
        if(altitude < 0.f) throw(string("Climate file : wrong format altitude (m above sea level)"));

        input.getline(cbuffer, 200, '\n');
        meridianE = atoi(cbuffer);
        //logStream<<"Meridian: "<<meridianE<<endl << flush;
        if(meridianE > 12.f || meridianE < -12.f) throw(string("Climate file : wrong format meridian (degrees East)"));

        // reads the header of the file
        input >> buffer;
        if(buffer != "dm") throw(string("Climate file : wrong format dm"));     //day

        input >> buffer;
        if(buffer != "m") throw(string("Climate file : wrong format m"));     //month

        input >> buffer;
        if(buffer != "h") throw(string("Climate file : wrong format h"));     //jour

        // find irradiance data type
        input >> buffer;
        if(buffer == "G_h") {
            G_h_present = true;
            G_Dh_present = false;
            logStream << "G_h only, using Maxwell model to convert to G_Bn and G_Dh." << endl << flush;
        }
        else if (buffer == "G_Dh") {
            G_Dh_present = true;
            input >> buffer;
            if (buffer == "G_h") { G_h_present = true; logStream << "G_Dh and G_h, using sun position to convert to G_Dh and G_Bn." << endl << flush; }
            else if (buffer == "G_Bn") { G_h_present = false; logStream << "G_Dh and G_Bn provided." << endl << flush; }
            else throw(string("Climate file : wrong format, missing G_h or G_Bn after G_Dh"));
        }
        else throw(string("Climate file : wrong format, missing G_h or G_Dh"));

        input >> buffer;
        if(buffer != "Ta") throw(string("Climate file : wrong format, missing Ta"));     // Tout

        input >> buffer;
        if (buffer == "Tg" || buffer == "Ts") {
            Tg_present = true;
            input >> buffer;
            if(buffer != "FF") throw(string("Climate file : wrong format FF"));     // windSpeed
        }
        else if(buffer != "FF") throw(string("Climate file : wrong format FF"));     // windSpeed

        input >> buffer;
        if(buffer != "DD") throw(string("Climate file : wrong format DD"));     // windDirection

        input >> buffer;
        if(buffer != "RH") throw(string("Climate file : wrong format RH"));     // relativeHumidity

        input >> buffer;
        if(buffer != "RR") throw(string("Climate file : wrong format RR"));     // Prec

        input >> buffer;
        if(buffer != "N") throw(string("Climate file : wrong format N"));     // Cloudiness

        logStream << "reading of the data" << endl;
        // reading of the data
        while ( !input.eof() ) {

            input >> buffer;
            //int i = atoi(buffer.c_str());         //day

            input >> buffer;
            //int j = atoi(buffer.c_str());         //month

            input >> buffer;
            //int k = atoi(buffer.c_str());         //hour

            if (G_Dh_present) {
                input >> buffer;
                Idh.push_back(atof(buffer.c_str()));
                if (G_h_present) {
                    input >> buffer;
                    Igh.push_back(atof(buffer.c_str()));
                }
                else {
                    input >> buffer;
                    Ibn.push_back(atof(buffer.c_str()));
                }
            }
            else {
                input >> buffer;
                Igh.push_back(atof(buffer.c_str()));
            }

            input >> buffer;
            Tout.push_back(atof(buffer.c_str()));

            if (Tg_present) {
                input >> buffer;
                Tground.push_back(atof(buffer.c_str()));
            }

            input >> buffer;
            windSpeed.push_back(atof(buffer.c_str()));

            input >> buffer;
            windDirection.push_back(atof(buffer.c_str()));

            input >> buffer;
            relativeHumidity.push_back(atof(buffer.c_str()));

            input >> buffer;
            Prec.push_back(atof(buffer.c_str()));

            input >> buffer;
            cloudiness.push_back(atof(buffer.c_str()));

        }

        //logStream << "last line: " << endl;
        //logStream << "Igh: " << Igh.back() << "\tTout: " << Tout.back() << "\twindSpeed: " << windSpeed.back() << "\twindDirection: " << windDirection.back()
        //     << "\tRH: " << relativeHumidity.back() << "\tprec: " << Prec.back() << "\tcloudiness: " << cloudiness.back() << endl;

        // calculation of several variables that we use in the radiation model
        meanAnnualTemperature = accumulate(Tout.begin(),Tout.end(),0.f)/float(Tout.size());
        //calculate the temperature for a day and set the coolest day in the year and a vector of daily mean temperature
        meanDailyTemperature.assign(365,0.f);
        for(unsigned int i=0;i<365;++i) {
            meanDailyTemperature[i] = accumulate(Tout.begin()+24*i,Tout.begin()+24*(i+1),0.f)/24.f;
            // note: accumulate does not take into account the last element pointed by end: [first,last)
            //logStream << *(Tout.begin()+24*i) << "\t" << *(Tout.begin()+24*(i+1)) << endl << flush;
        }
        hotDay = distance(meanDailyTemperature.begin(), max_element(meanDailyTemperature.begin(), meanDailyTemperature.end()));
        coolDay = distance(meanDailyTemperature.begin(), min_element(meanDailyTemperature.begin(), meanDailyTemperature.end()));

        // determine dew point temperature
        Td = new float[Tout.size()];
        for(unsigned int i=0;i<Tout.size();++i){
            //        double a = 17.27;
            //        double b = 237.7;
            //        double gamma = a * Tout[i] /(b+Tout[i])+ log(Hr[i]/100.0);
            //        Td[i] = b*gamma /(a-gamma);
            //        if(Td[i] < 0 || Tout[i] < 0)
            //            Td[i] = Tout[i];
            // Clausius-Clapeyron equation with L/Rw = (1/1.85e-4)
            Td[i] = 1.f / ( 1.f/(Tout[i] + 273.15f) - 1.85e-4f * log(relativeHumidity[i]/100.f) ) - 273.15f;
        }

        // case in which we have G_h and D_h, we need to compute Ibn
        if (G_Dh_present && G_h_present) {

#ifdef DEBUG
            ofstream file((filename.substr(0,filename.size()-4)+"_SW.out").c_str());
            file << "#Ibn\tIdh\tIgh" << endl << flush;
            logStream << "G_Dh and G_h present, calculating G_Bn from the solar position." << endl << flush;
#endif

            // creates a sun for the current location
            SKYSiteLocation location(GENAngle::Degrees(latitudeN), GENAngle::Degrees(longitudeE), GENAngle::Degrees(meridianE*360.f/24.f), GENAngle::Degrees(0.f));
            SKYSun sunClimate(location);
            int dayOfYear;

            // loop on all values of G_h
            for (size_t j=0; j<Igh.size(); ++j) {
                // sets the calendar conditions
                dayOfYear = (j+1)/24+1;
                sunClimate.SetDay(dayOfYear);
                sunClimate.SetClockTime(j%24);
                //cout << "ClockTime: " << j%24+0.5f << endl << flush;
                // sets the sun altitude
                float solarAltitude = sunClimate.GetPosition().Altitude().radians(); // we get the altitude of the sun
                if (solarAltitude/M_PI*180.f < 6) solarAltitude = 0.f; // we neglect the altitudes below 6°
                if (solarAltitude > 0.f)
                    Ibn.push_back((Igh[j]-Idh[j])/sin(solarAltitude));
                else {
                    Ibn.push_back(0.f);
                    Idh[j]=Igh[j];
                }

#ifdef DEBUG
                file << Ibn[j] << "\t" << Idh[j] << "\t" << Igh[j] << endl;
#endif
            }

#ifdef DEBUG
            // close the file
            file.close();
#endif

            // clear the useless Igh vector
            Igh.clear();

        }
        // If only global irradiance available, use model in dirint.c to compute Idh and Ibn
        else if (!G_Dh_present && G_h_present) {

#ifdef DEBUG
            ofstream file((filename.substr(0,filename.size()-4)+"_SW.out").c_str());

            file << "#sunAzimuth\tsunAltitude\tIbn\tIdh\tIgh" << endl;
            logStream << "G_h only present, calculating G_Dh and G_Bn." << endl << flush;
#endif

            // creates a sun for the current location
            SKYSiteLocation location(GENAngle::Degrees(latitudeN), GENAngle::Degrees(longitudeE), GENAngle::Degrees(meridianE*360.f/24.f), GENAngle::Degrees(0.f));
            SKYSun sunClimate(location);
            int dayOfYear = 1;

            float Ib = 0.f;
            float z[3] = {NA, NA, NA};  // previous, current and next data
            float Ig[3] = {NA, NA, NA}; // previous, current and next data

            for (size_t j=0; j<Igh.size()-1; ++j) {

                Ig[1] = Igh[j];
                //cerr << "Igh[" << j << "]: " << Igh[j] << endl;
                Ig[2] = Igh[j+1];

                dayOfYear = j/24+1;
                sunClimate.SetDay(dayOfYear);
                sunClimate.SetClockTime(j%24);
#ifdef DEBUG
                file << fmod(sunClimate.GetPosition().Azimuth().degrees()+360.f,360.f) << "\t" << sunClimate.GetPosition().Altitude().degrees() << "\t";
#endif
                z[1]  = M_PI/2. - sunClimate.GetPosition().Altitude().radians(); // angle entre le zénith et le soleil (RADIANS), valeur courante
                sunClimate.SetDay((j+1)/24+1);
                sunClimate.SetClockTime((j+1)%24);
                z[2]  = M_PI/2. - sunClimate.GetPosition().Altitude().radians();	// angle entre le zénith et le soleil (RADIANS), valeur suivante

                /*
             ** float dirint(float *g, float *z, float *td, int *doy, float *alt)
             ** 'g' is the address of someplace that has 3 floats representing
             ** 	global irradiance in watts. The 3 refer to the previous, current
             **	and next observed data
             ** 'z' is the address of someplace with 3 floats representing
             **	the zenith angle in RADIANS of the three observations.
             ** 'td' points to the dewpoint temperature in degrees C.
             ** 'doy' points to the day of year (from 1).
             ** 'alt' points to the altitude of the site, in meters.
             **
             ** The returned value is the modelled direct normal irradiance in watts.
            */

                // old Maxwell method
                //Ib = dirint(&Ig[0], &z[0], &Td[j], &dayOfYear, &altitude);
                // beam normal extraterrestrial irradiance is taken here as a fixed value: 1367 W/m^2
                float i0 = (cos(static_cast<double>(dayOfYear)*(2.*M_PI/365.))*0.033 + 1.)*1367.;
                Ib = dirint_(&Ig[0], &z[0], &Td[j], &altitude, &i0);
                Ibn.push_back(Ib);
                Idh.push_back((Ig[1] - Ib*cos(z[1])));
#ifdef DEBUG
                file << Ibn.back() << "\t" << Idh.back() << "\t" << Igh[j] << endl << flush;
#endif

                /*logStream << "G_h : " << setprecision(3) << setw(5) << Ig[1]
                << ", G_bn : " << setprecision(3) << setw(5) << Ib
                << ", G_dh : " << setprecision(3) << setw(5) << Ig[1] - Ib*cos(z[1])
                << ", dy : " << setprecision(3) << setw(3) << dayOfYear
                << ", h : " << setprecision(3) << setw(3) << (j+1)%24
                << ", dT : " << setprecision(3) << setw(6) <<  Td[j]
                << ", z : " << setprecision(3) << setw(5) << z[0] << " " << setprecision(3) << setw(5) << z[1] << " " << setprecision(3) << setw(5) << z[2]
                << ", Ig : " << setprecision(3) << setw(5) << Ig[0] << " " << setprecision(3) << setw(5) << Ig[1] << " " << setprecision(3) << setw(5) << Ig[2] << endl << flush;
*/
                // pass to the next time step (move forward the data)
                Ig[0] = Ig[1];
                z[0] = z[1];
            }

            //cerr << "Igh[back]" << Ig[2] << endl;

            Ig[1] = Ig[2];
            Ig[2] = NA;
            z[1] = z[2];
            z[2] = NA;
            // old Maxwell method
            //Ib = dirint(&Ig[0], &z[0], &Td[Igh.size()-1], &dayOfYear, &altitude);
            // beam normal extraterrestrial irradiance is taken here as a fixed value: 1367 W/m^2
            dayOfYear = (Igh.size()-1)/24+1;
            float i0 = (cos(static_cast<double>(dayOfYear)*(2.*M_PI/365.))*0.033 + 1.)*1367.;
            Ib = dirint_(&Ig[0], &z[0], &Td[Igh.size()-1], &altitude, &i0);
            //cerr << "Ig[0]: " << Ig[0] << "\tIg[1]: " << Ig[1] << "\tIg[2]: " << Ig[2] << "\tTd: " << Td[Igh.size()-1] << endl;
            //cerr << "Ib: " << Ib << endl;
            Ibn.push_back(Ib);
            Idh.push_back(Ig[1] - Ib*cos(z[1]));

#ifdef DEBUG
            file << fmod(sunClimate.GetPosition().Azimuth().degrees()+360.f,360.f) << "\t" << sunClimate.GetPosition().Altitude().degrees() << "\t";
            file << Ibn.back() << "\t" << Idh.back() << "\t" << Igh.back() << endl;
            // close the file
            file.close();
#endif

            // clear the useless Igh vector
            Igh.clear();

            /// TODO: copier dans méthode pour MEU

        }
        else {

#ifdef DEBUG
            ofstream file((filename.substr(0,filename.size()-4)+"_SW.out").c_str());
            file << "#dayOfYear\thour\tsunAzimuth\tsunAltitude\tIbn\tIdh\tIgh" << endl;
            logStream << "G_Dh and G_Bn present, calculating G_h using the solar position." << endl << flush;
#endif

            // creates a sun for the current location
            SKYSiteLocation location(GENAngle::Degrees(latitudeN), GENAngle::Degrees(longitudeE), GENAngle::Degrees(meridianE*360.f/24.f), GENAngle::Degrees(0.f));
            SKYSun sunClimate(location);
            int dayOfYear, clockTime;

            // loop on all values of G_h
            for (size_t j=0; j<Idh.size(); ++j) {
                // sets the calendar conditions
                dayOfYear = j/24+1;
                sunClimate.SetDay(dayOfYear);
                clockTime = j%24+1;
                sunClimate.SetClockTime1(clockTime);
                // sets the sun altitude
                Igh.push_back(Idh[j]+Ibn[j]*sin(max(sunClimate.GetPosition().Altitude().radians(),0.f)));

#ifdef DEBUG
                file << dayOfYear << "\t" << clockTime << "\t";
                file << fmod(sunClimate.GetPosition().Azimuth().degrees()+360.f,360.f) << "\t" << sunClimate.GetPosition().Altitude().degrees() << "\t";
                file << Ibn[j] << "\t" << Idh[j] << "\t" << Igh[j] << endl;
#endif
            }

#ifdef DEBUG
            // close the file
            file.close();
#endif

            // clear the useless Igh vector
            Igh.clear();

        }

    /*}
    catch(string msg) { logStream << "\n(Caught) " << msg << endl; return; }
    catch(exception &e) { logStream << "\n(Caught Standard Exception) " << e.what() << endl; return; }
    catch(...) { logStream << "\nUnhandled Exception" << endl; return; }
    */
}

float Climate::dirint(float *g, float *z, float *td, int *doy, float *alt)
{

    // ************************************************************************
    // ******** Maxwell model to compute direct normal beam irradiance ********
    //
    // Code retrieved from the file dirint.c furnished by David Lindelöf
    // ************************************************************************

	/*
	 ** Maxwell's Beam model.
	 ** Modified to interpolate between bins when kt is high.
	 **
	 ** This function was translated from FORTRAN. The maxwell model
	 ** was converted using 'f2c', and the resulting code was hand-
	 ** massaged to introduce c-isms.
	 ** The interpolation code was simply transliterated line-by-line.
	 **
	 ** float dirint(float *g, float *z, float *td, int *doy, float *alt)
	 ** 'g' is the address of someplace that has 3 floats representing
	 ** 	global irradiance in watts. The 3 refer to the previous, current
	 **	and next observed data
	 ** 'z' is the address of someplace with 3 floats representing
	 **	the zenith angle in RADIANS of the three observations.
	 ** 'td' points to the dewpoint temperature in degrees C.
	 ** 'doy' points to the day of year (from 1).
	 ** 'alt' points to the altitude of the site, in meters.
	 **
	 ** The returned value is the modelled direct normal irradiance in watts.
	 **
	 **
	 ** Three observations are used to take into consideration the time-
	 ** varying nature of the sky. The parameterization was created
	 ** with the assumption that the time series is nominally 1 hour.
	 ** However, this function may be used satifactorily with time steps
	 ** as little as 1 minute or as large as 90 minutes.
	 ** If any of the the required values are unavailable, simply send the
	 ** magic cookie "-999.0", and the function will fall back to less
	 ** sophisticated modeling (*doy, *z, and g[1] excepted, of course).
	 **
	 ** Oh, you'll need an ANSI compiler for this unless you want to
	 ** get rid of the 'const' and prototypes and stuff.
	 **
	 ** 2-12-92 rob@dinner.asrc.albany.edu
	 ** Initial translation and shakedown. v1.0
	 */

	//static const char	*version = "1.0";

static const float	cm[6][6][7][5] =
{{{{0.38523, 0.38523, 0.38523, 0.46288, 0.31744},
    {0.33839, 0.33839, 0.22127, 0.31673, 0.50365},
    {0.23568, 0.23568, 0.24128, 0.15783, 0.26944},
    {0.83013, 0.83013, 0.17197, 0.84107, 0.45737},
    {0.54801, 0.54801, 0.47800, 0.96688, 1.03637},
    {0.54801, 0.54801, 1.00000, 3.01237, 1.97654},
    {0.58269, 0.58269, 0.22972, 0.89271, 0.56995}},
    {{0.13128, 0.13128, 0.38546, 0.51107, 0.12794},
        {0.22371, 0.22371, 0.19356, 0.30456, 0.19394},
        {0.22997, 0.22997, 0.27502, 0.31273, 0.24461},
        {0.09010, 0.18458, 0.26050, 0.68748, 0.57944},
        {0.13153, 0.13153, 0.37019, 1.38035, 1.05227},
        {1.11625, 1.11625, 0.92803, 3.52549, 2.31692},
        {0.09010, 0.23700, 0.30004, 0.81247, 0.66497}},
    {{0.58751, 0.13000, 0.40000, 0.53721, 0.83249},
        {0.30621, 0.12983, 0.20446, 0.50000, 0.68164},
        {0.22402, 0.26062, 0.33408, 0.50104, 0.35047},
        {0.42154, 0.75397, 0.75066, 3.70684, 0.98379},
        {0.70668, 0.37353, 1.24567, 0.86486, 1.99263},
        {4.86440, 0.11739, 0.26518, 0.35918, 3.31082},
        {0.39208, 0.49329, 0.65156, 1.93278, 0.89873}},
    {{0.12697, 0.12697, 0.12697, 0.12697, 0.12697},
        {0.81082, 0.81082, 0.81082, 0.81082, 0.81082},
        {3.24168, 2.50000, 2.29144, 2.29144, 2.29144},
        {4.00000, 3.00000, 2.00000, 0.97543, 1.96557},
        {12.49417, 12.49417, 8.00000, 5.08352, 8.79239},
        {21.74424, 21.74424, 21.74424, 21.74424, 21.74424},
        {3.24168, 12.49417, 1.62076, 1.37525, 2.33162}},
    {{0.12697, 0.12697, 0.12697, 0.12697, 0.12697},
        {0.81082, 0.81082, 0.81082, 0.81082, 0.81082},
        {3.24168, 2.50000, 2.29144, 2.29144, 2.29144},
        {4.00000, 3.00000, 2.00000, 0.97543, 1.96557},
        {12.49417, 12.49417, 8.00000, 5.08352, 8.79239},
        {21.74424, 21.74424, 21.74424, 21.74424, 21.74424},
        {3.24168, 12.49417, 1.62076, 1.37525, 2.33162}},
    {{0.12697, 0.12697, 0.12697, 0.12697, 0.12697},
        {0.81082, 0.81082, 0.81082, 0.81082, 0.81082},
        {3.24168, 2.50000, 2.29144, 2.29144, 2.29144},
        {4.00000, 3.00000, 2.00000, 0.97543, 1.96557},
        {12.49417, 12.49417, 8.00000, 5.08352, 8.79239},
        {21.74424, 21.74424, 21.74424, 21.74424, 21.74424},
        {3.24168, 12.49417, 1.62076, 1.37525, 2.33162}}},
    {{{0.33744, 0.33744, 0.96911, 1.09719, 1.11608},
        {0.33744, 0.33744, 0.96911, 1.11603, 0.62390},
        {0.33744, 0.33744, 1.53059, 1.02442, 0.90848},
        {0.58404, 0.58404, 0.84725, 0.91494, 1.28930},
        {0.33744, 0.33744, 0.31024, 1.43502, 1.85283},
        {0.33744, 0.33744, 1.01501, 1.09719, 2.11723},
        {0.33744, 0.33744, 0.96911, 1.14573, 1.47640}},
        {{0.30000, 0.30000, 0.70000, 1.10000, 0.79694},
            {0.21987, 0.21987, 0.52653, 0.80961, 0.64930},
            {0.38665, 0.38665, 0.11932, 0.57612, 0.68546},
            {0.74673, 0.39983, 0.47097, 0.98653, 0.78537},
            {0.57542, 0.93670, 1.64920, 1.49584, 1.33559},
            {1.31967, 4.00257, 1.27639, 2.64455, 2.51867},
            {0.66519, 0.67891, 1.01236, 1.19994, 0.98658}},
        {{0.37887, 0.97406, 0.50000, 0.49188, 0.66529},
            {0.10521, 0.26347, 0.40704, 0.55346, 0.58259},
            {0.31290, 0.34524, 1.14418, 0.85479, 0.61228},
            {0.11907, 0.36512, 0.56052, 0.79372, 0.80260},
            {0.78161, 0.83739, 1.27042, 1.53798, 1.29295},
            {1.15229, 1.15229, 1.49208, 1.24537, 2.17710},
            {0.42466, 0.52955, 0.96691, 1.03346, 0.95873}},
        {{0.31059, 0.71441, 0.25245, 0.50000, 0.60760},
            {0.97519, 0.36342, 0.50000, 0.40000, 0.50280},
            {0.17558, 0.19625, 0.47636, 1.07247, 0.49051},
            {0.71928, 0.69862, 0.65777, 1.19084, 0.68111},
            {0.42624, 1.46484, 0.67855, 1.15773, 0.97843},
            {2.50112, 1.78913, 1.38709, 2.39418, 2.39418},
            {0.49164, 0.67761, 0.68561, 1.08240, 0.73541}},
        {{0.59700, 0.50000, 0.30000, 0.31005, 0.41351},
            {0.31479, 0.33631, 0.40000, 0.40000, 0.44246},
            {0.16651, 0.46044, 0.55257, 1.00000, 0.46161},
            {0.40102, 0.55911, 0.40363, 1.01671, 0.67149},
            {0.40036, 0.75083, 0.84264, 1.80260, 1.02383},
            {3.31530, 1.51038, 2.44365, 1.63882, 2.13399},
            {0.53079, 0.74585, 0.69305, 1.45804, 0.80450}},
        {{0.59700, 0.50000, 0.30000, 0.31005, 0.80092},
            {0.31479, 0.33631, 0.40000, 0.40000, 0.23704},
            {0.16651, 0.46044, 0.55257, 1.00000, 0.58199},
            {0.40102, 0.55911, 0.40363, 1.01671, 0.89857},
            {0.40036, 0.75083, 0.84264, 1.80260, 3.40039},
            {3.31530, 1.51038, 2.44365, 1.63882, 2.50878},
            {0.20434, 1.15774, 2.00308, 2.62208, 1.40938}}},
    {{{1.24221, 1.24221, 1.24221, 1.24221, 1.24221},
        {0.05698, 0.05698, 0.65699, 0.65699, 0.92516},
        {0.08909, 0.08909, 1.04043, 1.23248, 1.20530},
        {1.05385, 1.05385, 1.39969, 1.08464, 1.23334},
        {1.15154, 1.15154, 1.11829, 1.53164, 1.41184},
        {1.49498, 1.49498, 1.70000, 1.80081, 1.67160},
        {1.01845, 1.01845, 1.15360, 1.32189, 1.29467}},
        {{0.70000, 0.70000, 1.02346, 0.70000, 0.94583},
            {0.88630, 0.88630, 1.33362, 0.80000, 1.06662},
            {0.90218, 0.90218, 0.95433, 1.12669, 1.09731},
            {1.09530, 1.07506, 1.17649, 1.13947, 1.09611},
            {1.20166, 1.20166, 1.43820, 1.25628, 1.19806},
            {1.52585, 1.52585, 1.86916, 1.98541, 1.91159},
            {1.28822, 1.08281, 1.28637, 1.16617, 1.11933}},
        {{0.60000, 1.02991, 0.85989, 0.55000, 0.81360},
            {0.60445, 1.02991, 0.85989, 0.65670, 0.92884},
            {0.45585, 0.75058, 0.80493, 0.82300, 0.91100},
            {0.52658, 0.93231, 0.90862, 0.98352, 0.98809},
            {1.03611, 1.10069, 0.84838, 1.03527, 1.04238},
            {1.04844, 1.65272, 0.90000, 2.35041, 1.08295},
            {0.81741, 0.97616, 0.86130, 0.97478, 1.00458}},
        {{0.78211, 0.56428, 0.60000, 0.60000, 0.66574},
            {0.89448, 0.68073, 0.54199, 0.80000, 0.66914},
            {0.48746, 0.81895, 0.84183, 0.87254, 0.70904},
            {0.70931, 0.87278, 0.90848, 0.95329, 0.84435},
            {0.86392, 0.94777, 0.87622, 1.07875, 0.93691},
            {1.28035, 0.86672, 0.76979, 1.07875, 0.97513},
            {0.72542, 0.86997, 0.86881, 0.95119, 0.82922}},
        {{0.79175, 0.65404, 0.48317, 0.40900, 0.59718},
            {0.56614, 0.94899, 0.97182, 0.65357, 0.71855},
            {0.64871, 0.63773, 0.87051, 0.86060, 0.69430},
            {0.63763, 0.76761, 0.92567, 0.99031, 0.84767},
            {0.73638, 0.94606, 1.11759, 1.02934, 0.94702},
            {1.18097, 0.85000, 1.05000, 0.95000, 0.88858},
            {0.70056, 0.80144, 0.96197, 0.90614, 0.82388}},
        {{0.50000, 0.50000, 0.58677, 0.47055, 0.62979},
            {0.50000, 0.50000, 1.05622, 1.26014, 0.65814},
            {0.50000, 0.50000, 0.63183, 0.84262, 0.58278},
            {0.55471, 0.73473, 0.98582, 0.91564, 0.89826},
            {0.71251, 1.20599, 0.90951, 1.07826, 0.88561},
            {1.89926, 1.55971, 1.00000, 1.15000, 1.12039},
            {0.65388, 0.79312, 0.90332, 0.94407, 0.79613}}},
    {{{1.00000, 1.00000, 1.05000, 1.17038, 1.17809},
        {0.96058, 0.96058, 1.05953, 1.17903, 1.13169},
        {0.87147, 0.87147, 0.99586, 1.14191, 1.11460},
        {1.20159, 1.20159, 0.99361, 1.10938, 1.12632},
        {1.06501, 1.06501, 0.82866, 0.93997, 1.01793},
        {1.06501, 1.06501, 0.62369, 1.11962, 1.13226},
        {1.07157, 1.07157, 0.95807, 1.11413, 1.12711}},
        {{0.95000, 0.97339, 0.85252, 1.09220, 1.09659},
            {0.80412, 0.91387, 0.98099, 1.09458, 1.04242},
            {0.73754, 0.93597, 0.99994, 1.05649, 1.05006},
            {1.03298, 1.03454, 0.96846, 1.03208, 1.01578},
            {0.90000, 0.97721, 0.94596, 1.00884, 0.96996},
            {0.60000, 0.75000, 0.75000, 0.84471, 0.89910},
            {0.92680, 0.96503, 0.96852, 1.04491, 1.03231}},
        {{0.85000, 1.02971, 0.96110, 1.05567, 1.00970},
            {0.81853, 0.96001, 0.99645, 1.08197, 1.03647},
            {0.76538, 0.95350, 0.94826, 1.05211, 1.00014},
            {0.77561, 0.90961, 0.92780, 0.98780, 0.95210},
            {1.00099, 0.88188, 0.87595, 0.94910, 0.89369},
            {0.90237, 0.87596, 0.80799, 0.94241, 0.91792},
            {0.85658, 0.92827, 0.94682, 1.03226, 0.97299}},
        {{0.75000, 0.85793, 0.98380, 1.05654, 0.98024},
            {0.75000, 0.98701, 1.01373, 1.13378, 1.03825},
            {0.80000, 0.94738, 1.01238, 1.09127, 0.99984},
            {0.80000, 0.91455, 0.90857, 0.99919, 0.91523},
            {0.77854, 0.80059, 0.79907, 0.90218, 0.85156},
            {0.68019, 0.31741, 0.50768, 0.38891, 0.64671},
            {0.79492, 0.91278, 0.96083, 1.05711, 0.94795}},
        {{0.75000, 0.83389, 0.86753, 1.05989, 0.93284},
            {0.97970, 0.97147, 0.99551, 1.06849, 1.03015},
            {0.85885, 0.98792, 1.04322, 1.10870, 1.04490},
            {0.80240, 0.95511, 0.91166, 1.04507, 0.94447},
            {0.88489, 0.76621, 0.88539, 0.85907, 0.81819},
            {0.61568, 0.70000, 0.85000, 0.62462, 0.66930},
            {0.83557, 0.94615, 0.97709, 1.04935, 0.97997}},
        {{0.68922, 0.80960, 0.90000, 0.78950, 0.85399},
            {0.85466, 0.85284, 0.93820, 0.92311, 0.95501},
            {0.93860, 0.93298, 1.01039, 1.04395, 1.04164},
            {0.84362, 0.98130, 0.95159, 0.94610, 0.96633},
            {0.69474, 0.81469, 0.57265, 0.40000, 0.72683},
            {0.21137, 0.67178, 0.41634, 0.29729, 0.49805},
            {0.84354, 0.88233, 0.91176, 0.89842, 0.96021}}},
    {{{1.05488, 1.07521, 1.06846, 1.15337, 1.06922},
        {1.00000, 1.06222, 1.01347, 1.08817, 1.04620},
        {0.88509, 0.99353, 0.94259, 1.05499, 1.01274},
        {0.92000, 0.95000, 0.97872, 1.02028, 0.98444},
        {0.85000, 0.90850, 0.83994, 0.98557, 0.96218},
        {0.80000, 0.80000, 0.81008, 0.95000, 0.96155},
        {1.03859, 1.06320, 1.03444, 1.11278, 1.03780}},
        {{1.01761, 1.02836, 1.05896, 1.13318, 1.04562},
            {0.92000, 0.99897, 1.03359, 1.08903, 1.02206},
            {0.91237, 0.94993, 0.97977, 1.02042, 0.98177},
            {0.84716, 0.93530, 0.93054, 0.95505, 0.94656},
            {0.88026, 0.86711, 0.87413, 0.97265, 0.88342},
            {0.62715, 0.62715, 0.70000, 0.77407, 0.84513},
            {0.97370, 1.00624, 1.02619, 1.07196, 1.01724}},
        {{1.02871, 1.01757, 1.02590, 1.08179, 1.02424},
            {0.92498, 0.98550, 1.01410, 1.09221, 0.99961},
            {0.82857, 0.93492, 0.99495, 1.02459, 0.94971},
            {0.90081, 0.90133, 0.92883, 0.97957, 0.91310},
            {0.76103, 0.84515, 0.80536, 0.93679, 0.85346},
            {0.62640, 0.54675, 0.73050, 0.85000, 0.68905},
            {0.95763, 0.98548, 0.99179, 1.05022, 0.98790}},
        {{0.99273, 0.99388, 1.01715, 1.05912, 1.01745},
            {0.97561, 0.98716, 1.02682, 1.07544, 1.00725},
            {0.87109, 0.93319, 0.97469, 0.97984, 0.95273},
            {0.82875, 0.86809, 0.83492, 0.90551, 0.87153},
            {0.78154, 0.78247, 0.76791, 0.76414, 0.79589},
            {0.74346, 0.69339, 0.51487, 0.63015, 0.71566},
            {0.93476, 0.95787, 0.95964, 0.97251, 0.98164}},
        {{0.96584, 0.94124, 0.98710, 1.02254, 1.01116},
            {0.98863, 0.99477, 0.97659, 0.95000, 1.03484},
            {0.95820, 1.01808, 0.97448, 0.92000, 0.98987},
            {0.81172, 0.86909, 0.81202, 0.85000, 0.82105},
            {0.68203, 0.67948, 0.63245, 0.74658, 0.73855},
            {0.66829, 0.44586, 0.50000, 0.67892, 0.69651},
            {0.92694, 0.95335, 0.95905, 0.87621, 0.99149}},
        {{0.94894, 0.99776, 0.85000, 0.82652, 0.99847},
            {1.01786, 0.97000, 0.85000, 0.70000, 0.98856},
            {1.00000, 0.95000, 0.85000, 0.60624, 0.94726},
            {1.00000, 0.74614, 0.75174, 0.59839, 0.72523},
            {0.92221, 0.50000, 0.37680, 0.51711, 0.54863},
            {0.50000, 0.45000, 0.42997, 0.40449, 0.53994},
            {0.96043, 0.88163, 0.77564, 0.59635, 0.93768}}},
    {{{1.03000, 1.04000, 1.00000, 1.00000, 1.04951},
        {1.05000, 0.99000, 0.99000, 0.95000, 0.99653},
        {1.05000, 0.99000, 0.99000, 0.82000, 0.97194},
        {1.05000, 0.79000, 0.88000, 0.82000, 0.95184},
        {1.00000, 0.53000, 0.44000, 0.71000, 0.92873},
        {0.54000, 0.47000, 0.50000, 0.55000, 0.77395},
        {1.03827, 0.92018, 0.91093, 0.82114, 1.03456}},
        {{1.04102, 0.99752, 0.96160, 1.00000, 1.03578},
            {0.94803, 0.98000, 0.90000, 0.95036, 0.97746},
            {0.95000, 0.97725, 0.86927, 0.80000, 0.95168},
            {0.95187, 0.85000, 0.74877, 0.70000, 0.88385},
            {0.90000, 0.82319, 0.72745, 0.60000, 0.83987},
            {0.85000, 0.80502, 0.69231, 0.50000, 0.78841},
            {1.01009, 0.89527, 0.77303, 0.81628, 1.01168}},
        {{1.02245, 1.00460, 0.98365, 1.00000, 1.03294},
            {0.94396, 0.99924, 0.98392, 0.90599, 0.97815},
            {0.93624, 0.94648, 0.85000, 0.85000, 0.93032},
            {0.81642, 0.88500, 0.64495, 0.81765, 0.86531},
            {0.74296, 0.76569, 0.56152, 0.70000, 0.82714},
            {0.64387, 0.59671, 0.47446, 0.60000, 0.65120},
            {0.97174, 0.94056, 0.71488, 0.86438, 1.00165}},
        {{0.99526, 0.97701, 1.00000, 1.00000, 1.03525},
            {0.93981, 0.97525, 0.93998, 0.95000, 0.98255},
            {0.87687, 0.87944, 0.85000, 0.90000, 0.91781},
            {0.87348, 0.87345, 0.75147, 0.85000, 0.86304},
            {0.76147, 0.70236, 0.63877, 0.75000, 0.78312},
            {0.73408, 0.65000, 0.60000, 0.65000, 0.71566},
            {0.94216, 0.91910, 0.77034, 0.73117, 0.99518}},
        {{0.95256, 0.91678, 0.92000, 0.90000, 1.00588},
            {0.92862, 0.99442, 0.90000, 0.90000, 0.98372},
            {0.91307, 0.85000, 0.85000, 0.80000, 0.92428},
            {0.86809, 0.80717, 0.82355, 0.60000, 0.84452},
            {0.76957, 0.71987, 0.65000, 0.55000, 0.73350},
            {0.58025, 0.65000, 0.60000, 0.50000, 0.62885},
            {0.90477, 0.85265, 0.70837, 0.49373, 0.94903}},
        {{0.91197, 0.80000, 0.80000, 0.80000, 0.95632},
            {0.91262, 0.68261, 0.75000, 0.70000, 0.95011},
            {0.65345, 0.65933, 0.70000, 0.60000, 0.85611},
            {0.64844, 0.60000, 0.64112, 0.50000, 0.69578},
            {0.57000, 0.55000, 0.59880, 0.40000, 0.56015},
            {0.47523, 0.50000, 0.51864, 0.33997, 0.52023},
            {0.74344, 0.59219, 0.60306, 0.31693, 0.79439}}}};

	/*logStream << "alt : " << setprecision(3) << setw(5) << *alt
	<< ", doy : " << setprecision(3) << setw(3) << *doy
	<< ", td : " << setprecision(3) << setw(6) <<  *td
	<< ", z : " << setprecision(3) << setw(5) << z[0] << " " << setprecision(3) << setw(5) << z[1] << " " << setprecision(3) << setw(5) << z[2]
	<< ", g : " << setprecision(3) << setw(5) << g[0] << " " << setprecision(3) << setw(5) << g[1] << " " << setprecision(3) << setw(5) << g[2] << endl << flush;
	*/

    /* Initialized data */
    static const float	ktbin[6] = {0.24, 0.4, 0.56, 0.7, 0.8, 99999.99},
    zbin[6] = {25.0, 40.0, 55.0, 70.0, 80.0, 99999.99},
    dktbin[6] = {0.015, 0.035, 0.07, 0.15, 0.3, 99999.99},
    wbin[4] = {1.0, 2.0, 3.0, 99999.99},
    zbin2[6] = {19.0, 32.5, 47.5, 62.5, 75, 82.5},
    wbin2[4] = {0.75, 1.5, 2.5, 3.5};


    /* System generated locals */
    int			i__1;
    float		ret_val, d__1, d__2, d__3;

    /* Local variables */
    static float	bmax, a, b, c, w, ktpam[3], am[3],
    am2, am3, am4, kt2, kt3, coef,
    io, cz[3], kt[3], zenith[3], kt1[3], knc, dkt1;
    static int		i, j, k, l;


    /* Parameter adjustments */
    --z;
    --g;

    /* Function Body */
    if (g[2] < 1.0)
        return 0.0;

    // extraterrestrial irradiation is taken here as a fixed value: 1367 W/m^2
    io = (cos((*doy)*0.0172142)*0.033 + 1.0)*1367.0;
    j = 1;
    k = 3;
    if (g[1] == -999.0 || z[1] == -999.0) {
        j = 2;
        kt1[0] = -999.0;
    }
    if (g[3] == -999.0 || z[3] == -999.0) {
        k = 2;
        kt1[2] = -999.0;
    }

    i__1 = k;
    for (i=j; i <= i__1; ++i) {
        cz[i-1] = cos(z[i]);

        if (cz[i-1] < 0.0) {
            kt1[i-1] = -999.0;
        } else {
            zenith[i-1] = z[i] * (180./M_PI);

            /* Computing MAX */
            d__1 = 0.065, d__2 = cz[i-1];
            kt[i-1] = g[i] / (io * max(d__1, d__2));

            /* Computing MIN */
            d__3 = 93.9 - zenith[i-1];
            d__1 = 15.25, d__2 = 1.0 / (cz[i-1] + pow(d__3, (float) -1.253) * 0.15);
            am[i-1] = min(d__1, d__2);

            ktpam[i-1] = am[i-1] * exp(*alt * -1.184e-4);
            kt1[i-1] = kt[i-1] /
            (exp(-1.4 / (9.4 / ktpam[i-1] + 0.9)) * 1.031 + 0.1);
        }
    }

    if (cz[1] < 0.0)
        return 0.0;

	//logStream << "fraction kt : " << kt[0] << " " << kt[1] << " " << kt[2] << endl << flush;
    kt2 = kt[1]*kt[1];
    kt3 = kt2*kt[1];

    if (kt[1] <= 0.6) {
        a = 0.512 - kt[1]*1.56 + kt2*2.286 - kt3*2.22;
        b = kt[1]*0.962 + 0.37;
        c = kt[1]*0.932 - 0.28 - kt2*2.048;
    } else {
        a = -5.743 + kt[1]*21.77 - kt2*27.49 + kt3*11.56;
        b = 41.4 - kt[1]*118.5 + kt2*66.05 + kt3*31.9;
        c = kt[1]*184.2 - 47.01 - kt2*222.0 + kt3*73.81;
    }

    am2 = am[1]*am[1];
    am3 = am2*am[1];
    am4 = am3*am[1];
    knc = 0.866 - am[1]*0.122 + am2*0.0121 - am3*6.53e-4 + am4*1.4e-5;

	//logStream << "Calcul de bmax : io = " << io << ", knc = " << knc << ", a = " << a << ", b = " << b << ", c = " << c << endl << flush;
    bmax = io * (knc - (a + b * exp(c * am[1]))); // ***************** problem before here ***********

	//logStream << "bmax : " << bmax << endl << flush;

    if (kt1[0] == -999.0 && kt1[2] == -999.0) {
        k = 6;
    } else {
        if (kt1[0] == -999.0 || zenith[0] >= 85.0) {
            dkt1 = fabs(kt1[2] - kt1[1]);
        } else if (kt1[2] == -999.|| zenith[2] >= 85.0) {
            dkt1 = fabs(kt1[1] - kt1[0]);
        } else {
            dkt1 = 0.5*(fabs(kt1[1] - kt1[0]) + fabs(kt1[2] - kt1[1]));
        }

        for (k=0; dkt1 > dktbin[k]; k++) {
        }
    }


    for (i=0; kt1[1] > ktbin[i]; i++) {
    }

    for (j=0; zenith[1] > zbin[j]; j++) {
    }

    if (*td == -999.0) {
        l = 4;
    } else {
        w = exp(*td * 0.07 - 0.075);

        for (l=0; w > wbin[l]; l++) {
        }
    }



    /*
     ** If kt is very high, do some bin interpolating
     */
    if (i >= 4) {
        float		aak, bbk, rz, aaz, bbz, rw, aaw, bbw,
        c11, c12, c13, c14, c1, c2;

        if (kt1[1] >= 0.75) {
            aak = min(1.0, (kt1[1] - 0.75)/0.07);
            bbk = 1.0 - aak;
        } else {
            aak = 0.0;
            bbk = 1.0;
        }
        rz = zenith[1] - zbin2[j];
        if (rz >= 0.0) {
            if (j == 5) {
                aaz = 1.0;
                bbz = 0.0;
            } else {
                aaz = (zbin2[j+1] - zenith[1]) / (zbin2[j+1] - zbin2[j]);
                bbz = 1.0 - aaz;
            }
        } else {
            if (j == 0) {
                aaz = 1.0;
                bbz = 0.0;
            } else {
                j--;
                aaz = (zbin2[j+1] - zenith[1]) / (zbin2[j+1] - zbin2[j]);
                bbz = 1.0 - aaz;
            }
        }

        rw = w -wbin2[l];
        if (l == 4) {
            aaw = 1.0;
            bbw = 0.0;
        } else {
            if (rw >= 0.0) {
                if (l == 3) {
                    aaw = 1.0;
                    bbw = 0.0;
                } else {
                    aaw = (wbin2[l+1] - w) / (wbin2[l+1] - wbin2[l]);
                    bbw = 1.0 - aaw;
                }
            } else {
                if (l == 0) {
                    aaw = 1.0;
                    bbw = 0.0;
                } else {
                    l--;
                    aaw = (wbin2[l+1] - w) / (wbin2[l+1] - wbin2[l]);
                    bbw = 1.0 - aaw;
                }
            }
        }

        c11 = aak*cm[5][min(5, j+1)][k][l] + bbk*cm[4][min(5, j+1)][k][l];
        c12 = aak*cm[5][j][k][l] + bbk*cm[4][j][k][l];
        c13 = aak*cm[5][min(5, j+1)][k][min(4, l+1)] +
        bbk*cm[4][min(5, j+1)][k][min(4, l+1)];
        c14 = aak*cm[5][j][k][min(4, l+1)] + bbk*cm[4][j][k][min(4, l+1)];
        c1 = aaz*c12 + bbz*c11;
        c2 = aaz*c14 + bbz*c13;
        coef = aaw*c1 + bbw*c2;
    } else
        coef = cm[i][j][k][l];

    if ((ret_val = bmax*coef) < 0.0){
		//logStream << " ******** ret_val < 0, bmax : " << bmax << endl << flush;
        ret_val = 0.0;
	}

    return ret_val;
}				/* dirint */

double Climate::dirint_(float *g, float *z__, float *td, float *alt, float *i0)
{
    /*     FUNCTION NAME:    DIRINT */

    /*     PROGRAMMED BY:    HOWARD M. BISNER AFTER PEREZ ET AL. */
    /*     REVISED BY   :    ANTOINE J.F. ZELENKA */
    /*     ARGUMENTS: G(3) - GLOBAL IRRADIANCE (WATTS / SQ. METER) */
    /*                Z(3) - SOLAR ZENITH ANGLE (RADIANS) */
    /*                  TD - DEW POINT TEMPERATURE ( DEGREES C) */
    /*                  I0 - TOA NORMAL BEAM OF DATE (WATT / SQ. METER) */
    /*                 ALT - ALTITUDE OF SITE (METERS) */
    /*     RETURNS: DIRMAX - BEAM IRRADIANCE (WATTS / SQ. METER) */
    /*     NOTES:            THE FUNCTION DIRINT USES THE DISC BEAM MODEL TO */
    /*                       CALCULATE THE BEAM IRRADIANCE RETURNED. THE */
    /*                       ARGUMENT G IS AN ARRAY OF 3 VALUES. THE VALUES */
    /*                       ARE THE GLOBAL IRRADIANCE OF THE PREVIOUS READING, */
    /*                       THE CURRENT READING ,AND THE NEXT READING IN THAT */
    /*                       ORDER. THE ARGUMENT Z USES THE SAME FORMAT, EXCEPT */
    /*                       THE VALUES ARE THE RESPECTIVE SOLAR ZENITH ANGLES. */
    /*                       IF ANY OF THE G OR Z VALUES ARE NOT AVAILABLE OR THE */
    /*                       PREVIOUS OR NEXT READINGS DID NOT OCCUR WITHIN 1.5 */
    /*                       HOURS OF THE CURRENT READING THEN THE APPROPRIATE */
    /*                       VALUE OR VALUES SHOULD BE REPLACED WITH A -999.0. */
    /*                       IF THE ARGUMENT TD IS MISSING THEN THE VALUE */
    /*                       -999.0 SHOULD BE USED IN PLACE OF THE MISSING ARGUMENT. */
    /*                       THE CURRENT GLOBAL IRRADIANCE (G(2)) MUST HAVE A */
    /*                       VALUE. IF THE DEW POINT TEMPERATURE (TD) IS MISSING */
    /*                       THEN TD IS NOT USED TO FIND AN INDEX INTO THE */
    /*                       CORRECTION MATRIX (CM), INSTEAD A SPECIAL COLUMN IN */
    /*                       THE MATRIX IS USED. IF THE PREVIOUS GLOBAL */
    /*                       IRRADIANCE (G(1)) OR SOLAR ZENITH ANGLE (Z(1)) AND */
    /*                       THE NEXT GLOBAL IRRADIANCE (G(3)) OR SOLAR ZENITH ANGLE */
    /*                       (Z(3)) ARE MISSING THEN DELTA KT' (DKT1) IS NOT USED TO */
    /*                       FIND AN INDEX INTO THE CORRECTION MATRIX (CM), INSTEAD */
    /*                       A SPECIAL COLUMN IN THE MATRIX IS USED. */

    /*                       ADDED A 0.82 CAP ON KTPRIME */

    /* Table of constant values */
    static double c_b3 = -1.253;

    /* Initialized data */
    static float ktbin[6] = { .24f,.4f,.56f,.7f,.8f,1.f };
    static float zbin[6] = { 25.f,40.f,55.f,70.f,80.f,90.f };
    static float dktbin[6] = { .015f,.035f,.07f,.15f,.3f,1.f };
    static float wbin[4] = { 1.f,2.f,3.f,1e33f };
    static float zbin2[6] = { 19.f,32.5f,47.5f,62.5f,75.f,82.5f };
    static float wbin2[4] = { .75f,1.5f,2.5f,3.5f };
    static float rtod = 57.295779513082316f;
    static float cm[1260]	/* was [6][6][7][5] */ = { .38523f,.33744f,
	    1.24221f,1.f,1.05488f,1.03f,.13128f,.3f,.7f,.95f,1.01761f,
	    1.04102f,.58751f,.37887f,.6f,.85f,1.02871f,1.02245f,.12697f,
	    .31059f,.78211f,.75f,.99273f,.99526f,.12697f,.597f,.79175f,.75f,
	    .96584f,.95256f,.12697f,.597f,.5f,.68922f,.94894f,.91197f,.33839f,
	    .33744f,.05698f,.96058f,1.f,1.05f,.22371f,.21987f,.8863f,.80412f,
	    .92f,.94803f,.30621f,.10521f,.60445f,.81853f,.92498f,.94396f,
	    .81082f,.97519f,.89448f,.75f,.97561f,.93981f,.81082f,.31479f,
	    .56614f,.9797f,.98863f,.92862f,.81082f,.31479f,.5f,.85466f,
	    1.01786f,.91262f,.23568f,.33744f,.08909f,.87147f,.88509f,1.05f,
	    .22997f,.38665f,.90218f,.73754f,.91237f,.95f,.22402f,.3129f,
	    .45585f,.76538f,.82857f,.93624f,3.24168f,.17558f,.48746f,.8f,
	    .87109f,.87687f,3.24168f,.16651f,.64871f,.85885f,.9582f,.91307f,
	    3.24168f,.16651f,.5f,.9386f,1.f,.65345f,.83013f,.58404f,1.05385f,
	    1.20159f,.92f,1.05f,.0901f,.74673f,1.0953f,1.03298f,.84716f,
	    .95187f,.42154f,.11907f,.52658f,.77561f,.90081f,.81642f,4.f,
	    .71928f,.70931f,.8f,.82875f,.87348f,4.f,.40102f,.63763f,.8024f,
	    .81172f,.86809f,4.f,.40102f,.55471f,.84362f,1.f,.64844f,.54801f,
	    .33744f,1.15154f,1.06501f,.85f,1.f,.13153f,.57542f,1.20166f,.9f,
	    .88026f,.9f,.70668f,.78161f,1.03611f,1.00099f,.76103f,.74296f,
	    12.49417f,.42624f,.86392f,.77854f,.78154f,.76147f,12.49417f,
	    .40036f,.73638f,.88489f,.68203f,.76957f,12.49417f,.40036f,.71251f,
	    .69474f,.92221f,.57f,.54801f,.33744f,1.49498f,1.06501f,.8f,.54f,
	    1.11625f,1.31967f,1.52585f,.6f,.62715f,.85f,4.8644f,1.15229f,
	    1.04844f,.90237f,.6264f,.64387f,21.74424f,2.50112f,1.28035f,
	    .68019f,.74346f,.73408f,21.74424f,3.3153f,1.18097f,.61568f,
	    .66829f,.58025f,21.74424f,3.3153f,1.89926f,.21137f,.5f,.47523f,
	    .58269f,.33744f,1.01845f,1.07157f,1.03859f,1.03827f,.0901f,
	    .66519f,1.28822f,.9268f,.9737f,1.01009f,.39208f,.42466f,.81741f,
	    .85658f,.95763f,.97174f,3.24168f,.49164f,.72542f,.79492f,.93476f,
	    .94216f,3.24168f,.53079f,.70056f,.83557f,.92694f,.90477f,3.24168f,
	    .20434f,.65388f,.84354f,.96043f,.74344f,.38523f,.33744f,1.24221f,
	    1.f,1.07521f,1.04f,.13128f,.3f,.7f,.97339f,1.02836f,.99752f,.13f,
	    .97406f,1.02991f,1.02971f,1.01757f,1.0046f,.12697f,.71441f,
	    .56428f,.85793f,.99388f,.97701f,.12697f,.5f,.65404f,.83389f,
	    .94124f,.91678f,.12697f,.5f,.5f,.8096f,.99776f,.8f,.33839f,
	    .33744f,.05698f,.96058f,1.06222f,.99f,.22371f,.21987f,.8863f,
	    .91387f,.99897f,.98f,.12983f,.26347f,1.02991f,.96001f,.9855f,
	    .99924f,.81082f,.36342f,.68073f,.98701f,.98716f,.97525f,.81082f,
	    .33631f,.94899f,.97147f,.99477f,.99442f,.81082f,.33631f,.5f,
	    .85284f,.97f,.68261f,.23568f,.33744f,.08909f,.87147f,.99353f,.99f,
	    .22997f,.38665f,.90218f,.93597f,.94993f,.97725f,.26062f,.34524f,
	    .75058f,.9535f,.93492f,.94648f,2.5f,.19625f,.81895f,.94738f,
	    .93319f,.87944f,2.5f,.46044f,.63773f,.98792f,1.01808f,.85f,2.5f,
	    .46044f,.5f,.93298f,.95f,.65933f,.83013f,.58404f,1.05385f,
	    1.20159f,.95f,.79f,.18458f,.39983f,1.07506f,1.03454f,.9353f,.85f,
	    .75397f,.36512f,.93231f,.90961f,.90133f,.885f,3.f,.69862f,.87278f,
	    .91455f,.86809f,.87345f,3.f,.55911f,.76761f,.95511f,.86909f,
	    .80717f,3.f,.55911f,.73473f,.9813f,.74614f,.6f,.54801f,.33744f,
	    1.15154f,1.06501f,.9085f,.53f,.13153f,.9367f,1.20166f,.97721f,
	    .86711f,.82319f,.37353f,.83739f,1.10069f,.88188f,.84515f,.76569f,
	    12.49417f,1.46484f,.94777f,.80059f,.78247f,.70236f,12.49417f,
	    .75083f,.94606f,.76621f,.67948f,.71987f,12.49417f,.75083f,
	    1.20599f,.81469f,.5f,.55f,.54801f,.33744f,1.49498f,1.06501f,.8f,
	    .47f,1.11625f,4.00257f,1.52585f,.75f,.62715f,.80502f,.11739f,
	    1.15229f,1.65272f,.87596f,.54675f,.59671f,21.74424f,1.78913f,
	    .86672f,.31741f,.69339f,.65f,21.74424f,1.51038f,.85f,.7f,.44586f,
	    .65f,21.74424f,1.51038f,1.55971f,.67178f,.45f,.5f,.58269f,.33744f,
	    1.01845f,1.07157f,1.0632f,.92018f,.237f,.67891f,1.08281f,.96503f,
	    1.00624f,.89527f,.49329f,.52955f,.97616f,.92827f,.98548f,.94056f,
	    12.49417f,.67761f,.86997f,.91278f,.95787f,.9191f,12.49417f,
	    .74585f,.80144f,.94615f,.95335f,.85265f,12.49417f,1.15774f,
	    .79312f,.88233f,.88163f,.59219f,.38523f,.96911f,1.24221f,1.05f,
	    1.06846f,1.f,.38546f,.7f,1.02346f,.85252f,1.05896f,.9616f,.4f,.5f,
	    .85989f,.9611f,1.0259f,.98365f,.12697f,.25245f,.6f,.9838f,
	    1.01715f,1.f,.12697f,.3f,.48317f,.86753f,.9871f,.92f,.12697f,.3f,
	    .58677f,.9f,.85f,.8f,.22127f,.96911f,.65699f,1.05953f,1.01347f,
	    .99f,.19356f,.52653f,1.33362f,.98099f,1.03359f,.9f,.20446f,
	    .40704f,.85989f,.99645f,1.0141f,.98392f,.81082f,.5f,.54199f,
	    1.01373f,1.02682f,.93998f,.81082f,.4f,.97182f,.99551f,.97659f,.9f,
	    .81082f,.4f,1.05622f,.9382f,.85f,.75f,.24128f,1.53059f,1.04043f,
	    .99586f,.94259f,.99f,.27502f,.11932f,.95433f,.99994f,.97977f,
	    .86927f,.33408f,1.14418f,.80493f,.94826f,.99495f,.85f,2.29144f,
	    .47636f,.84183f,1.01238f,.97469f,.85f,2.29144f,.55257f,.87051f,
	    1.04322f,.97448f,.85f,2.29144f,.55257f,.63183f,1.01039f,.85f,.7f,
	    .17197f,.84725f,1.39969f,.99361f,.97872f,.88f,.2605f,.47097f,
	    1.17649f,.96846f,.93054f,.74877f,.75066f,.56052f,.90862f,.9278f,
	    .92883f,.64495f,2.f,.65777f,.90848f,.90857f,.83492f,.75147f,2.f,
	    .40363f,.92567f,.91166f,.81202f,.82355f,2.f,.40363f,.98582f,
	    .95159f,.75174f,.64112f,.478f,.31024f,1.11829f,.82866f,.83994f,
	    .44f,.37019f,1.6492f,1.4382f,.94596f,.87413f,.72745f,1.24567f,
	    1.27042f,.84838f,.87595f,.80536f,.56152f,8.f,.67855f,.87622f,
	    .79907f,.76791f,.63877f,8.f,.84264f,1.11759f,.88539f,.63245f,.65f,
	    8.f,.84264f,.90951f,.57265f,.3768f,.5988f,1.f,1.01501f,1.7f,
	    .62369f,.81008f,.5f,.92803f,1.27639f,1.86916f,.75f,.7f,.69231f,
	    .26518f,1.49208f,.9f,.80799f,.7305f,.47446f,21.74424f,1.38709f,
	    .76979f,.50768f,.51487f,.6f,21.74424f,2.44365f,1.05f,.85f,.5f,.6f,
	    21.74424f,2.44365f,1.f,.41634f,.42997f,.51864f,.22972f,.96911f,
	    1.1536f,.95807f,1.03444f,.91093f,.30004f,1.01236f,1.28637f,
	    .96852f,1.02619f,.77303f,.65156f,.96691f,.8613f,.94682f,.99179f,
	    .71488f,1.62076f,.68561f,.86881f,.96083f,.95964f,.77034f,1.62076f,
	    .69305f,.96197f,.97709f,.95905f,.70837f,1.62076f,2.00308f,.90332f,
	    .91176f,.77564f,.60306f,.46288f,1.09719f,1.24221f,1.17038f,
	    1.15337f,1.f,.51107f,1.1f,.7f,1.0922f,1.13318f,1.f,.53721f,
	    .49188f,.55f,1.05567f,1.08179f,1.f,.12697f,.5f,.6f,1.05654f,
	    1.05912f,1.f,.12697f,.31005f,.409f,1.05989f,1.02254f,.9f,.12697f,
	    .31005f,.47055f,.7895f,.82652f,.8f,.31673f,1.11603f,.65699f,
	    1.17903f,1.08817f,.95f,.30456f,.80961f,.8f,1.09458f,1.08903f,
	    .95036f,.5f,.55346f,.6567f,1.08197f,1.09221f,.90599f,.81082f,.4f,
	    .8f,1.13378f,1.07544f,.95f,.81082f,.4f,.65357f,1.06849f,.95f,.9f,
	    .81082f,.4f,1.26014f,.92311f,.7f,.7f,.15783f,1.02442f,1.23248f,
	    1.14191f,1.05499f,.82f,.31273f,.57612f,1.12669f,1.05649f,1.02042f,
	    .8f,.50104f,.85479f,.823f,1.05211f,1.02459f,.85f,2.29144f,
	    1.07247f,.87254f,1.09127f,.97984f,.9f,2.29144f,1.f,.8606f,1.1087f,
	    .92f,.8f,2.29144f,1.f,.84262f,1.04395f,.60624f,.6f,.84107f,
	    .91494f,1.08464f,1.10938f,1.02028f,.82f,.68748f,.98653f,1.13947f,
	    1.03208f,.95505f,.7f,3.70684f,.79372f,.98352f,.9878f,.97957f,
	    .81765f,.97543f,1.19084f,.95329f,.99919f,.90551f,.85f,.97543f,
	    1.01671f,.99031f,1.04507f,.85f,.6f,.97543f,1.01671f,.91564f,
	    .9461f,.59839f,.5f,.96688f,1.43502f,1.53164f,.93997f,.98557f,.71f,
	    1.38035f,1.49584f,1.25628f,1.00884f,.97265f,.6f,.86486f,1.53798f,
	    1.03527f,.9491f,.93679f,.7f,5.08352f,1.15773f,1.07875f,.90218f,
	    .76414f,.75f,5.08352f,1.8026f,1.02934f,.85907f,.74658f,.55f,
	    5.08352f,1.8026f,1.07826f,.4f,.51711f,.4f,3.01237f,1.09719f,
	    1.80081f,1.11962f,.95f,.55f,3.52549f,2.64455f,1.98541f,.84471f,
	    .77407f,.5f,.35918f,1.24537f,2.35041f,.94241f,.85f,.6f,21.74424f,
	    2.39418f,1.07875f,.38891f,.63015f,.65f,21.74424f,1.63882f,.95f,
	    .62462f,.67892f,.5f,21.74424f,1.63882f,1.15f,.29729f,.40449f,
	    .33997f,.89271f,1.14573f,1.32189f,1.11413f,1.11278f,.82114f,
	    .81247f,1.19994f,1.16617f,1.04491f,1.07196f,.81628f,1.93278f,
	    1.03346f,.97478f,1.03226f,1.05022f,.86438f,1.37525f,1.0824f,
	    .95119f,1.05711f,.97251f,.73117f,1.37525f,1.45804f,.90614f,
	    1.04935f,.87621f,.49373f,1.37525f,2.62208f,.94407f,.89842f,
	    .59635f,.31693f,.31744f,1.11608f,1.24221f,1.17809f,1.06922f,
	    1.04951f,.12794f,.79694f,.94583f,1.09659f,1.04562f,1.03578f,
	    .83249f,.66529f,.8136f,1.0097f,1.02424f,1.03294f,.12697f,.6076f,
	    .66574f,.98024f,1.01745f,1.03525f,.12697f,.41351f,.59718f,.93284f,
	    1.01116f,1.00588f,.12697f,.80092f,.62979f,.85399f,.99847f,.95632f,
	    .50365f,.6239f,.92516f,1.13169f,1.0462f,.99653f,.19394f,.6493f,
	    1.06662f,1.04242f,1.02206f,.97746f,.68164f,.58259f,.92884f,
	    1.03647f,.99961f,.97815f,.81082f,.5028f,.66914f,1.03825f,1.00725f,
	    .98255f,.81082f,.44246f,.71855f,1.03015f,1.03484f,.98372f,.81082f,
	    .23704f,.65814f,.95501f,.98856f,.95011f,.26944f,.90848f,1.2053f,
	    1.1146f,1.01274f,.97194f,.24461f,.68546f,1.09731f,1.05006f,
	    .98177f,.95168f,.35047f,.61228f,.911f,1.00014f,.94971f,.93032f,
	    2.29144f,.49051f,.70904f,.99984f,.95273f,.91781f,2.29144f,.46161f,
	    .6943f,1.0449f,.98987f,.92428f,2.29144f,.58199f,.58278f,1.04164f,
	    .94726f,.85611f,.45737f,1.2893f,1.23334f,1.12632f,.98444f,.95184f,
	    .57944f,.78537f,1.09611f,1.01578f,.94656f,.88385f,.98379f,.8026f,
	    .98809f,.9521f,.9131f,.86531f,1.96557f,.68111f,.84435f,.91523f,
	    .87153f,.86304f,1.96557f,.67149f,.84767f,.94447f,.82105f,.84452f,
	    1.96557f,.89857f,.89826f,.96633f,.72523f,.69578f,1.03637f,
	    1.85283f,1.41184f,1.01793f,.96218f,.92873f,1.05227f,1.33559f,
	    1.19806f,.96996f,.88342f,.83987f,1.99263f,1.29295f,1.04238f,
	    .89369f,.85346f,.82714f,8.79239f,.97843f,.93691f,.85156f,.79589f,
	    .78312f,8.79239f,1.02383f,.94702f,.81819f,.73855f,.7335f,8.79239f,
	    3.40039f,.88561f,.72683f,.54863f,.56015f,1.97654f,2.11723f,
	    1.6716f,1.13226f,.96155f,.77395f,2.31692f,2.51867f,1.91159f,
	    .8991f,.84513f,.78841f,3.31082f,2.1771f,1.08295f,.91792f,.68905f,
	    .6512f,21.74424f,2.39418f,.97513f,.64671f,.71566f,.71566f,
	    21.74424f,2.13399f,.88858f,.6693f,.69651f,.62885f,21.74424f,
	    2.50878f,1.12039f,.49805f,.53994f,.52023f,.56995f,1.4764f,
	    1.29467f,1.12711f,1.0378f,1.03456f,.66497f,.98658f,1.11933f,
	    1.03231f,1.01724f,1.01168f,.89873f,.95873f,1.00458f,.97299f,
	    .9879f,1.00165f,2.33162f,.73541f,.82922f,.94795f,.98164f,.99518f,
	    2.33162f,.8045f,.82388f,.97997f,.99149f,.94903f,2.33162f,1.40938f,
	    .79613f,.96021f,.93768f,.79439f };

    /* System generated locals */
    int i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    float ret_val, r__1, r__2;
    double d__1;

    /* Local variables */
    static float a, b, c__;
    static int i__, j, k, l;
    static float w, c1, c2, c11, c12, c13, c14, am[3], cz[3], kt[3], rw, rz,
	    kt1[3], aak, bbk, aim, aaw, bbw, knc, aaz, bbz, xkt, dkt1, coef,
	    calt, bmax, zenith[3];

    /* Parameter adjustments */
    --z__;
    --g;

    /* Function Body */
    if (g[2] < 1.f || z__[2] > 90.f / rtod) {
	ret_val = 0.f;
	return ret_val;
    }
    calt = exp(*alt * -1.184e-4f);
    j = 1;
    k = 3;
    if (g[1] == -999.f || z__[1] == -999.f) {
	j = 2;
	kt1[0] = -999.f;
    }
    if (g[3] == -999.f || z__[3] == -999.f) {
	k = 2;
	kt1[2] = -999.f;
    }
    i__1 = k;
    for (i__ = j; i__ <= i__1; ++i__) {
	cz[i__ - 1] = cos(z__[i__]);
	if (cz[i__ - 1] < 0.f) {
	    kt1[i__ - 1] = -999.f;
	} else {
	    zenith[i__ - 1] = z__[i__] * rtod;
/* Computing MAX */
	    r__1 = .0521791f, r__2 = cz[i__ - 1];
	    kt[i__ - 1] = g[i__] / (*i0 * max(r__1,r__2));
/* Computing MIN */
	    d__1 = (double) (93.885f - zenith[i__ - 1]);
	    r__1 = 15.25f, r__2 = 1.f / (cz[i__ - 1] + pow(d__1, c_b3) *
		    .15f);
	    am[i__ - 1] = min(r__1,r__2);
	    kt1[i__ - 1] = kt[i__ - 1] / (exp(-1.4f / (9.4f / (am[i__ - 1] *
		    calt) + .9f)) * 1.031f + .1f);
	}
	if (kt1[i__ - 1] >= .82f) {
	    kt1[i__ - 1] = .82f;
	}
/* L100: */
    }
    aim = am[1];
    xkt = kt[1];
    knc = aim * (aim * (aim * (aim * 1.4e-5f - 6.53e-4f) + .0121f) - .122f) +
	    .866f;
    if (xkt <= .6f) {
	a = xkt * (xkt * (2.286f - xkt * 2.222f) - 1.56f) + .512f;
	b = xkt * .962f + .37f;
	c__ = xkt * (.932f - xkt * 2.048f) - .28f;
    } else {
	a = xkt * (xkt * (xkt * 11.56f - 27.49f) + 21.77f) - 5.743f;
	b = xkt * (xkt * (xkt * 31.9f + 66.05f) - 118.5f) + 41.4f;
	c__ = xkt * (xkt * (xkt * 73.81f - 222.f) + 184.2f) - 47.01f;
    }
    bmax = knc - (a + b * exp(c__ * aim));
    if (bmax < 0.f) {
	ret_val = 0.f;
	return ret_val;
    }
    if (kt1[0] == -999.f && kt1[2] == -999.f) {
	k = 7;
    } else {
	if (kt1[0] == -999.f || zenith[0] >= 85.f) {
	    dkt1 = (r__1 = kt1[2] - kt1[1], abs(r__1));
	} else {
	    if (kt1[2] == -999.f || zenith[2] >= 85.f) {
		dkt1 = (r__1 = kt1[1] - kt1[0], abs(r__1));
	    } else {
		dkt1 = ((r__1 = kt1[1] - kt1[0], abs(r__1)) + (r__2 = kt1[2]
			- kt1[1], abs(r__2))) / 2.f;
	    }
	}
	for (k = 1; k <= 6; ++k) {
	    if (k == 6 || dkt1 < dktbin[k - 1]) {
		goto L111;
	    }
/* L110: */
	}
    }
L111:
    for (i__ = 1; i__ <= 6; ++i__) {
	if (i__ == 6 || kt1[1] < ktbin[i__ - 1]) {
	    goto L121;
	}
/* L120: */
    }
L121:
    for (j = 1; j <= 6; ++j) {
	if (j == 6 || zenith[1] < zbin[j - 1]) {
	    goto L131;
	}
/* L130: */
    }
L131:
    if (*td == -999.f) {
	l = 5;
    } else {
	w = exp(*td * .07f - .075f);
	for (l = 1; l <= 4; ++l) {
	    if (l == 4 || w < wbin[l - 1]) {
		goto L141;
	    }
/* L140: */
	}
    }
L141:
    if (i__ >= 5) {
	if (kt1[1] >= .75f) {
/* Computing MIN */
	    r__1 = 1.f, r__2 = (kt1[1] - .75f) / .07f;
	    aak = min(r__1,r__2);
	    bbk = 1.f - aak;
	} else {
	    aak = 0.f;
	    bbk = 1.f - aak;
	}
	rz = zenith[1] - zbin2[j - 1];
	if (rz >= 0.f) {
	    if (j == 6) {
		aaz = 1.f;
		bbz = 1.f - aaz;
	    } else {
		aaz = (zbin2[j] - zenith[1]) / (zbin2[j] - zbin2[j - 1]);
		bbz = 1.f - aaz;
	    }
	} else {
	    if (j == 1) {
		aaz = 1.f;
		bbz = 1.f - aaz;
	    } else {
		--j;
		aaz = (zbin2[j] - zenith[1]) / (zbin2[j] - zbin2[j - 1]);
		bbz = 1.f - aaz;
	    }
	}
	if (l == 5) {
	    aaw = 1.f;
	    bbw = 1.f - aaw;
	} else {
	    rw = w - wbin2[l - 1];
	    if (rw >= 0.f) {
		if (l == 4) {
		    aaw = 1.f;
		    bbw = 1.f - aaw;
		} else {
		    aaw = (wbin2[l] - w) / (wbin2[l] - wbin2[l - 1]);
		    bbw = 1.f - aaw;
		}
	    } else {
		if (l == 1) {
		    aaw = 1.f;
		    bbw = 1.f - aaw;
		} else {
		    --l;
		    aaw = (wbin2[l] - w) / (wbin2[l] - wbin2[l - 1]);
		    bbw = 1.f - aaw;
		}
	    }
	}
/* Computing MIN */
	i__1 = 6, i__2 = j + 1;
/* Computing MIN */
	i__3 = 6, i__4 = j + 1;
	c11 = aak * cm[(min(i__1,i__2) + (k + l * 7) * 6) * 6 - 289] + bbk *
		cm[(min(i__3,i__4) + (k + l * 7) * 6) * 6 - 290];
	c12 = aak * cm[(j + (k + l * 7) * 6) * 6 - 289] + bbk * cm[(j + (k +
		l * 7) * 6) * 6 - 290];
/* Computing MIN */
	i__1 = 6, i__2 = j + 1;
/* Computing MIN */
	i__3 = 5, i__4 = l + 1;
/* Computing MIN */
	i__5 = 6, i__6 = j + 1;
/* Computing MIN */
	i__7 = 5, i__8 = l + 1;
	c13 = aak * cm[(min(i__1,i__2) + (k + min(i__3,i__4) * 7) * 6) * 6 -
		289] + bbk * cm[(min(i__5,i__6) + (k + min(i__7,i__8) * 7) *
		6) * 6 - 290];
/* Computing MIN */
	i__1 = 5, i__2 = l + 1;
/* Computing MIN */
	i__3 = 5, i__4 = l + 1;
	c14 = aak * cm[(j + (k + min(i__1,i__2) * 7) * 6) * 6 - 289] + bbk *
		cm[(j + (k + min(i__3,i__4) * 7) * 6) * 6 - 290];
	c1 = aaz * c12 + bbz * c11;
	c2 = aaz * c14 + bbz * c13;
	coef = aaw * c1 + bbw * c2;
    } else {
	coef = cm[i__ + (j + (k + l * 7) * 6) * 6 - 295];
    }
    ret_val = *i0 * bmax * coef;
    return ret_val;
} /* dirint_ */

float Climate::getTgroundCelsius(unsigned int day, unsigned int hour, float depth, float alpha, float depthMax) {

    // case with a fixed depth, if depthMax is not given, it is put to 0
    if (depthMax <= depth) {

        // if no surface ground temperature is given, or if depth is given
        if (Tground.empty() || depth > 0.f) {

            // returns a calculated ground temperature from Darren Robinson's equation
            return meanAnnualTemperature - (meanDailyTemperature[hotDay]-meanDailyTemperature[coolDay])/2.f //Modified by Max (added the divided by 2 for the amplitude)
                                           *exp(-depth*sqrt(M_PI/(365.f*alpha)))
                                           *cos(2.f*M_PI*(static_cast<float>(day)-static_cast<float>(coolDay+1)-depth*0.5f*sqrt(365.f/(M_PI*alpha)))/365.f);

        }
        else return Tground[(day-1)*24 + hour -1];

    }
    else {
    // case with the average temperature

        double integralTemperature = depthMax*meanAnnualTemperature - depth*meanAnnualTemperature
                                     - (sqrt(365.f/M_PI)*(meanDailyTemperature[hotDay]-meanDailyTemperature[coolDay])
                                       *(( cos((-2.f*day*M_PI + 2.f*static_cast<float>(coolDay+1)*M_PI + depth*sqrt(365.f*M_PI)*sqrt(1.f/alpha))/365.f)
                                          +sin((M_PI*(2.f*day - 2.f*static_cast<float>(coolDay+1) - depth*sqrt(365.f/M_PI)*sqrt(1.f/alpha)))/365.f))
                                         /exp(depth*sqrt(M_PI/365.f)*sqrt(1.f/alpha)) +
                                         (-cos((-2.f*day*M_PI + 2.f*static_cast<float>(coolDay+1)*M_PI + depthMax*sqrt(365.f*M_PI)*sqrt(1.f/alpha))/365.f)
                                          +sin((-2.f*day*M_PI + 2.f*static_cast<float>(coolDay+1)*M_PI + depthMax*sqrt(365.f*M_PI)*sqrt(1.f/alpha))/365.f))
                                         /exp(depthMax*sqrt(M_PI/365.f)*sqrt(1.f/alpha))))/(2.f*sqrt(1.f/alpha));

        return integralTemperature/(depthMax-depth);

    }

}

void Climate::exportCliFile(string filename){
    ofstream file(filename.c_str());
    if (!file.is_open()) throw(string("Error creating model .xml file: "+filename));

    string loc=location;
    loc.erase(loc.find_last_not_of(" \n\r\t")+1);
    file << loc << "\r\n" ;
    file << latitudeN << "," << longitudeE << "," << altitude << "," << meridianE << "\r\n" << "\r\n";
    if(Tground.empty())
        file << "dm	m	h	G_Dh	G_Bn	Ta	FF	DD	RH	RR	N" << "\r\n";
    else
        file << "dm	m	h	G_Dh	G_Bn	Ta  Tg	FF	DD	RH	RR	N" << "\r\n";

    unsigned int month=1, day=1, hour=1;
    for (int i=0; i<8760; ++i){
        if(Tground.empty())
            file << day << "\t" << month << "\t" << hour << "\t" << Idh[i] << "\t" << Ibn[i] << "\t" << Tout[i] << "\t" << windSpeed[i] << "\t" << windDirection[i] << "\t" << relativeHumidity[i] << "\t" << Prec[i] << "\t" << cloudiness[i] << "\r\n";
        else
            file << day << "\t" << month << "\t" << hour << "\t" << Idh[i] << "\t" << Ibn[i] << "\t" << Tout[i] << "\t" << Tground[i] << "\t" << windSpeed[i] << "\t" << windDirection[i] << "\t" << relativeHumidity[i] << "\t" << Prec[i] << "\t" << cloudiness[i] << "\r\n";
    }
    file << "\r\n";
    file.close();
}

void Climate::loadCli2(string filename) {

    // Climate2 filename opening
    fstream input(filename.c_str(), ios::in | ios::binary);
    if (!input.is_open()) {
        logStream << "No .cli2 weather file." << endl << flush;
        return;
    }
    else logStream << "Loading: " << filename << endl << flush;

    // variables used to store the elements
    string buffer;
    vector<pair<unsigned int, string> > headerLine;

    getline(input,buffer); // gets a line in the buffer

    size_t pos1=0, pos2=0;
    do {
        pos1 = buffer.find(":", pos2)+1;
        pos2 = buffer.find(":", pos1);
        //cout << "pos1: " << pos1 << "\tpos2: " << pos2 << endl;
        //cout << buffer.substr(pos1,pos2-pos1) << endl;
        unsigned int surfaceId = to<unsigned int>(buffer.substr(pos1,pos2-pos1));
        pos1 = ++pos2;
        pos2 = buffer.find("\t", pos1);
        //cout << buffer.substr(pos1, pos2-pos1) << endl;
        string header = buffer.substr(pos1, pos2-pos1);
        if (header == "KeCoeff1")
            KeCoeff1.insert(pair<unsigned int, vector<float> > (surfaceId, vector<float>()));
        else if (header == "KeCoeff2")
            KeCoeff2.insert(pair<unsigned int, vector<float> > (surfaceId, vector<float>()));
        else if (header == "KeCoeff3")
            KeCoeff3.insert(pair<unsigned int, vector<float> > (surfaceId, vector<float>()));
        else if (header.substr(0,3) == "Tla")
            Tla.insert(pair<unsigned int, vector<float> > (surfaceId, vector<float>()));
        else throw(string("Error in the .cli2 header."));
        headerLine.push_back(pair<unsigned int, string>(surfaceId,header));
    } while ( pos2 != std::string::npos);

    //cout << "headerline contains: " << headerLine.size() << endl;

    while (!input.eof()) {
        // read the data line by line
        for (size_t index=0; index < headerLine.size();) {
            input >> buffer; if (buffer=="") break;
            KeCoeff1[headerLine[index++].first].push_back(to<float>(buffer));
            input >> buffer; if (buffer=="") break;
            KeCoeff2[headerLine[index++].first].push_back(to<float>(buffer));
            input >> buffer; if (buffer=="") break;
            KeCoeff3[headerLine[index++].first].push_back(to<float>(buffer));
            input >> buffer; if (buffer=="") break;
            Tla[headerLine[index++].first].push_back(to<float>(buffer));
        }
    }

    input.close();

    return;

}

void Climate::loadCli3(string filename) {

     // Climate3 filename opening
    fstream input(filename.c_str(), ios::in | ios::binary);
    if (!input.is_open()) {
        logStream << "No .cli3 weather file." << endl << flush;
        return;
    }
    else logStream << "Loading: " << filename << endl << flush;

    // variables used to store the elements
    string buffer;
    getline(input,buffer);

    stringstream line(buffer);
    vector<string> tokens{istream_iterator<string>{line}, istream_iterator<string>{}};

    // read the values
    bool firstLine = true;
    while (!input.eof()) {
        for (size_t i=0;i<tokens.size();++i) {
            // prepare the maps to store the data
            string header = tokens.at(i).substr(0,tokens.at(i).find_first_of("("));
            float height = to<float>(tokens.at(i).substr(tokens.at(i).find_first_of("(")+1,tokens.at(i).find_first_of(")")-tokens.at(i).find_first_of("(")-1));

            if (header == "FF") {
                if (firstLine) FF.insert(pair<float, vector<float> > (height, vector<float>()));
                input >> buffer; if (buffer=="") break;
                FF[height].push_back(to<float>(buffer));
            }
            else if (header == "DD") {
                if (firstLine) DD.insert(pair<float, vector<float> > (height, vector<float>()));
                input >> buffer; if (buffer=="") break;
                DD[height].push_back(to<float>(buffer));
            }
            else if (header == "Ta") {
                if (firstLine) Ta.insert(pair<float, vector<float> > (height, vector<float>()));
                input >> buffer; if (buffer=="") break;
                Ta[height].push_back(to<float>(buffer));
            }
            else throw(string("Error in the .cli3 header."));
        }
        // header was read
        firstLine = false;
    }

    return;

}
