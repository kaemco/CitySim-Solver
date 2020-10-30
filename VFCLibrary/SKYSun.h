#ifndef _INC_SKYSUN_INCLUDED
#define _INC_SKYSUN_INCLUDED

/**
 *  Class to describe movement of sun
 *
 */

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "SKYSiteLocation.h"
#include "GENPoint.h"

class SKYSun
{
public:
	SKYSun (const SKYSiteLocation& location=SKYSiteLocation()); // latitude in radians, day is integer 1-365

    double getAperture() { return 0.533; } // return the diameter of the sun seen from earth in degrees
	double getSolidAngle() { return 2.*M_PI*(1.-cos((getAperture()/2.)*M_PI/180.)); } // return the solid angle in sr

	// Difference between solar and clock time (add result to clock time to get solar)
	// return value in hours
	float TimeDiffHours() const;

	bool SetDay(int day); // day from [1,365]
	int GetDay() const;

	bool SetClockTime(float hour); // hour from [0,24[
	bool SetClockTime1(float hour); // hour from [1,24]

	bool SetSolarTime(float hour);
	float GetSolarTime() const;

	const GENPoint& GetPosition() const; // solar position (in x,y,z form)

	const GENAngle& GetSunriseHourAngle() const;

	// True if sun is above horizon
	bool SunUp() const;

	bool SetHourAngle(const GENAngle &hourangle);

private:

	void Update();
	GENAngle CalculateAltitude();
	GENAngle CalculateAzimuth(const GENAngle &altitude);
	void CalculateSunrise();
	void CalculatePosition();

	GENAngle m_latitudeN;
	GENAngle m_longitudeE;
	GENAngle m_meridianE;
	GENAngle m_northAngleAnticlockwiseFromY;

	// solar time
	GENAngle m_dayAngle;
	GENAngle m_hourAngle;

	GENAngle m_declination;
	GENAngle m_sunriseHourAngle;
	GENPoint m_solarPosition;
};

#endif //_INC_SKYSUN_INCLUDED
