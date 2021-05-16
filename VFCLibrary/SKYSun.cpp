#include "SKYSun.h"
#include "SKYSiteLocation.h"

SKYSun::SKYSun (const SKYSiteLocation& location) :
		m_latitudeN(location.latitudeN),
		m_longitudeE(location.longitudeE),
		m_meridianE(location.meridianE),
		m_northAngleAnticlockwiseFromY(location.northAngleW),
		m_solarPosition(GENPoint::Origin())
{
	this->SetDay(1);
	this->SetHourAngle(GENAngle::Radians((float)M_PI/2));

	Update();
}

float SKYSun::TimeDiffHours() const
{
	const float Et = 229.2f * (0.000075f +
					0.001868f*cos(m_dayAngle) - 0.032077f*sin(m_dayAngle) -
					0.014615f*cos(m_dayAngle*2.f) - 0.04089f*sin(m_dayAngle*2.f));

	return (4.f/60.f)*(m_longitudeE-m_meridianE).degrees() + (Et/60.f);
}

bool SKYSun::SetDay(int day)
{
    // take the day within the boundaries, JK - 21.06.2015
    day = (day-1)%365+1;

	if (day >=1 && day <=365)
	{
		m_dayAngle=GENAngle::Radians(2.0f*(float)M_PI*(day-1.0f)/365.0f);
		m_declination = GENAngle::Radians(
			0.006918f - 0.399912f*cos(m_dayAngle) + 0.070257f*sin(m_dayAngle)
				- 0.006758f*cos(m_dayAngle*2.f) + 0.000907f*sin(m_dayAngle*2.f)
				- 0.002697f*cos(m_dayAngle*3.f) + 0.00148f*sin(m_dayAngle*3.f));

		Update();
		return true;
	}
	else
		return false;
}

int SKYSun::GetDay() const
{
	return (int((m_dayAngle.radians()*365/(2*M_PI)) +365.5) % 365) + 1;
}

bool SKYSun::SetHourAngle(const GENAngle &hourangle)
{
	m_hourAngle=hourangle;

	Update();

	if (hourangle.radians() >= m_sunriseHourAngle.radians() || hourangle.radians() < - m_sunriseHourAngle.radians())
	{
		return true;
	}
	else
		// return false if the user tries to set the time outside of the solar
        // day (but still let them do it)
		return false;
}

bool SKYSun::SetSolarTime(float hour)
{
	return SetHourAngle(GENAngle::Radians(hour*(float)M_PI/12));
}

float SKYSun::GetSolarTime()  const
{
	const float solarTime=m_hourAngle.radians()*(12.f/(float)M_PI);
	return solarTime<0 ? solarTime+24 : solarTime;
}

bool SKYSun::SetClockTime(float hour)
{
	return SetSolarTime(hour+TimeDiffHours());
}

bool SKYSun::SetClockTime1(float hour)
{
	return SetSolarTime(hour-0.5f+TimeDiffHours()); // MeteoNorm uses an average value on the hour -> middle position for the sun
}

void SKYSun::Update()
{
    CalculateSunrise();
	CalculatePosition();
}

GENAngle SKYSun::CalculateAltitude()
{
	return GENAngle::Radians(asin(
			sin(m_latitudeN)*sin(m_declination) -
			cos(m_latitudeN)*cos(m_declination)*cos(m_hourAngle)
		));
}

GENAngle SKYSun::CalculateAzimuth(const GENAngle &altitude)
{
  const float temp = (-sin(m_latitudeN)*sin(altitude) + sin(m_declination))
					/(cos(m_latitudeN)*cos(altitude));

  float az=0;
  if (temp > 1)
	 az= 0;
  else if (temp < -1)
	 az= (float)M_PI;
  else if (m_hourAngle.radians()>=0)
     az= acos(temp);
  else
     az= 2*(float)M_PI - acos(temp);

  return GENAngle::Radians(az)-m_northAngleAnticlockwiseFromY;
}

const GENPoint& SKYSun::GetPosition() const
{
	return m_solarPosition;
}

const GENAngle& SKYSun::GetSunriseHourAngle() const
{
	return m_sunriseHourAngle;
}

bool SKYSun::SunUp() const
{
	return m_solarPosition[GENPointCoords::Z]>0.f;
}

void SKYSun::CalculateSunrise()
{
    //const float t1=tan(m_latitudeN.radians())*tan(m_declination.radians()); // approximation equation
	// JK - better equation for sunrise and sunset
	const float t=-(sin(-0.83f/180.f*M_PI)-sin(m_latitudeN.radians())*sin(m_declination.radians()))/(cos(m_latitudeN.radians())*cos(m_declination.radians()));

	if (t>=1)
		m_sunriseHourAngle=GENAngle();
	else if (t<=-1)
		m_sunriseHourAngle=GENAngle::Radians((float)M_PI);
	else
		m_sunriseHourAngle=GENAngle::Radians(acos(t));
}

void SKYSun::CalculatePosition()
{
	const GENAngle altitude=CalculateAltitude();
	const GENAngle azimuth=CalculateAzimuth(altitude);
	m_solarPosition=GENPoint::Polar(altitude, azimuth);
}
