#include "./DATAInsolationFactors.h"

#include "GENAssert.h"

#include "SKYSun.h"
#include "GENPoint.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

using std::max;
using std::min;

DATAInsolationFactors::DATAInsolationFactors()
{
	SetInterpolationScheme(12,24);
}

void DATAInsolationFactors::SetInterpolationScheme(int daysPerYear, int stepsPerDay)
{
	m_interpolationDaysPerYear=daysPerYear;
	m_storedInterpolationDaysPerYear=1+m_interpolationDaysPerYear/2;
	m_interpolationPeriodsPerDay=stepsPerDay;
	m_insolationFactors.resize(m_storedInterpolationDaysPerYear*m_interpolationPeriodsPerDay,0.f);
}

void DATAInsolationFactors::GetSunPositionsToRender(const SKYSiteLocation &location,
												std::vector<GENPoint> &sunPositions) const
{
	SKYSun Sun(location);

	sunPositions.resize(m_storedInterpolationDaysPerYear*m_interpolationPeriodsPerDay);
	for (int i=0; i<m_storedInterpolationDaysPerYear; ++i)
	{
		// start with 356, and go through to 173
		const int day=dayOfYearFromIndex(((i-1)+m_interpolationDaysPerYear) %m_interpolationDaysPerYear);
		Sun.SetDay(day);
		for (int j=0; j<m_interpolationPeriodsPerDay; ++j)
		{
			Sun.SetSolarTime(hourFromIndex(j));
			sunPositions[i*m_interpolationPeriodsPerDay+j]=Sun.GetPosition();
		}
	}
}

float DATAInsolationFactors::GetIndexedViewFactor(unsigned int dayIndex, unsigned int hourIndex) const
{
	GENAssert(dayIndex>=0 && dayIndex<m_interpolationDaysPerYear);
	GENAssert(hourIndex>=0 && hourIndex<m_interpolationPeriodsPerDay);

	unsigned int arrayIndex=(dayIndex+1)%m_interpolationDaysPerYear;
	if (arrayIndex>=m_storedInterpolationDaysPerYear) arrayIndex=m_interpolationDaysPerYear-arrayIndex;

	return m_insolationFactors[arrayIndex*m_interpolationPeriodsPerDay+hourIndex];
}

float DATAInsolationFactors::GetIndexedViewFactor(unsigned int sunIndex) const
{
	return m_insolationFactors[sunIndex];
}

void DATAInsolationFactors::SetInsolationFactor(int sunPositionIndex, float vf)
{
	GENAssert(sunPositionIndex>=0 && sunPositionIndex<m_insolationFactors.size());
	m_insolationFactors[sunPositionIndex]=vf;
}

float DATAInsolationFactors::GetInsolationFactor(const SKYSun &Sun) const
{
	int hourIndex[3];
	hourIndex[1]=indexForTimePointBefore(Sun.GetSolarTime());

	if (hourIndex[1]<1) hourIndex[1]=1;
	if (hourIndex[1]>=m_interpolationPeriodsPerDay-1) hourIndex[1]=m_interpolationPeriodsPerDay-2;

	hourIndex[0]=hourIndex[1]-1;
	hourIndex[2]=hourIndex[1]+1;

	int dayIndex[3];
	dayIndex[1]=indexForDayPointBefore(Sun.GetDay());
	dayIndex[0]=dayIndex[1]-1;
	dayIndex[2]=dayIndex[1]+1;

	if (dayIndex[0] < 0) dayIndex[0]+=m_interpolationDaysPerYear;
	if (dayIndex[2]>m_interpolationDaysPerYear-1) dayIndex[2]-=m_interpolationDaysPerYear;

	int interpolationDays[3];

	for (unsigned int i=0; i<3; ++i)
		interpolationDays[i]=dayOfYearFromIndex(dayIndex[i]);

	if (0==dayIndex[1])
		interpolationDays[0]-=365;
	else if (m_interpolationDaysPerYear-1==dayIndex[1])
		interpolationDays[2]+=365;

	if (Sun.GetDay()<dayOfYearFromIndex(0))
	{
		for (unsigned int i=0; i<3; ++i)
			interpolationDays[i]-=365;
	}

	GENAssert(interpolationDays[1]>interpolationDays[0]);
	GENAssert(interpolationDays[2]>interpolationDays[1]);
	GENAssert(interpolationDays[1]-interpolationDays[0]<1.2*365/m_interpolationDaysPerYear);
	GENAssert(interpolationDays[2]-interpolationDays[1]<1.2*365/m_interpolationDaysPerYear);

	GENAssert(Sun.GetDay()>interpolationDays[0]);
	GENAssert(Sun.GetDay()<interpolationDays[2]);

	// interpolate across the hours
	float interpolatedFactors[3];
	for (unsigned int i=0; i<3; i++)
	{
		interpolatedFactors[i]=LagrangianInterpolation(Sun.GetSolarTime(),
						hourFromIndex(hourIndex[0]),GetIndexedViewFactor(dayIndex[i],hourIndex[0]),
						hourFromIndex(hourIndex[1]),GetIndexedViewFactor(dayIndex[i],hourIndex[1]),
						hourFromIndex(hourIndex[2]),GetIndexedViewFactor(dayIndex[i],hourIndex[2]));

		// Check that result is not out of range
		interpolatedFactors[i]=min(1.0f,interpolatedFactors[i]);
		interpolatedFactors[i]=max(0.0f,interpolatedFactors[i]);
	}

	float insolationFactor=LagrangianInterpolation((float)Sun.GetDay(),
										(float)interpolationDays[0],interpolatedFactors[0],
										(float)interpolationDays[1],interpolatedFactors[1],
										(float)interpolationDays[2],interpolatedFactors[2]);

	insolationFactor=min(1.0f,insolationFactor);
	insolationFactor=max(0.0f,insolationFactor);

	return insolationFactor;
}

float DATAInsolationFactors::GetInsolationFactor_linear_interp(const SKYSun &Sun) const
{
	int hourIndex[2];
	hourIndex[0]=indexForTimePointBefore(Sun.GetSolarTime());
	hourIndex[1]=hourIndex[0]+1;
	if (hourIndex[1]>=m_interpolationPeriodsPerDay) hourIndex[1]-=m_interpolationPeriodsPerDay;

	int interpolationDays[2];
	int dayIndex[2];

	dayIndex[0]=indexForDayPointBefore(Sun.GetDay());
	dayIndex[1]=(dayIndex[0]+1)%m_interpolationDaysPerYear;

	interpolationDays[0]=dayOfYearFromIndex(dayIndex[0]);
	interpolationDays[1]=(dayIndex[0]==m_interpolationDaysPerYear-1) ? dayOfYearFromIndex(0)+365 : dayOfYearFromIndex(dayIndex[1]);
	if (Sun.GetDay()<dayOfYearFromIndex(0))
	{
		interpolationDays[0]-=365;
		interpolationDays[1]-=365;
	}

	GENAssert(dayIndex[0]>=0 && dayIndex[0]<m_interpolationDaysPerYear);
	GENAssert(dayIndex[1]>=0 && dayIndex[1]<m_interpolationDaysPerYear);

	GENAssert(Sun.GetDay()>=interpolationDays[0]);
	GENAssert(Sun.GetDay()<=interpolationDays[1]);
	GENAssert(interpolationDays[1]-interpolationDays[0]<1.2*365/m_interpolationDaysPerYear);

	// interpolate across the hours
	float interpolatedViewFactors[2];
	for (unsigned int i=0; i<2; i++)
	{
		const float f=(Sun.GetSolarTime()-hourFromIndex(hourIndex[0]))/(24.f/m_interpolationPeriodsPerDay);
		GENAssert(f>=0 && f<=1);
		interpolatedViewFactors[i]=
			f*GetIndexedViewFactor(dayIndex[i],hourIndex[1])+(1-f)*GetIndexedViewFactor(dayIndex[i],hourIndex[0]);

		// Check that result is not out of range
		interpolatedViewFactors[i]=min(1.0f,interpolatedViewFactors[i]);
		interpolatedViewFactors[i]=max(0.0f,interpolatedViewFactors[i]);
	}

	const float f=((float)Sun.GetDay()-interpolationDays[0])/(interpolationDays[1]-interpolationDays[0]);
	GENAssert(f>=0 && f<=1);
	float insolationFactor=(1-f)*interpolatedViewFactors[0]+f*interpolatedViewFactors[1];

	insolationFactor=min(1.0f,insolationFactor);
	insolationFactor=max(0.0f,insolationFactor);

	return insolationFactor;
}

float DATAInsolationFactors::LagrangianInterpolation(float x, float x1, float y1, float x2, float y2, float x3, float y3) const
{
	// perform Lagrangian interpolation

	// some intermediate values
	const float xmx1=x-x1;
	const float xmx2=x-x2;
	const float xmx3=x-x3;

	// work out the multiplication factor for each y value
	const float scale1=xmx2*xmx3/((x1-x2)*(x1-x3));
	const float scale2=xmx1*xmx3/((x2-x1)*(x2-x3));
	const float scale3=xmx1*xmx2/((x3-x1)*(x3-x2));

	// return the result
	return scale1*y1+scale2*y2+scale3*y3;
}


float DATAInsolationFactors::hourFromIndex(int hourIndex) const
{
	GENAssert(hourIndex>=0 && hourIndex<DATAInsolationFactors::m_interpolationPeriodsPerDay);
	return 24.f*hourIndex/DATAInsolationFactors::m_interpolationPeriodsPerDay;
}
int DATAInsolationFactors::indexForTimePointBefore(float solarTime) const
{
	return static_cast<int>(solarTime*DATAInsolationFactors::m_interpolationPeriodsPerDay/24.);
}

int DATAInsolationFactors::dayOfYearFromIndex(int dayIndex) const
{
	GENAssert(dayIndex>=0 && dayIndex<DATAInsolationFactors::m_interpolationDaysPerYear);
	int day=356+static_cast<int>((dayIndex+1)*182./(DATAInsolationFactors::m_storedInterpolationDaysPerYear-1));
	if (day>365) day-=365;

	return day;
}

int DATAInsolationFactors::indexForDayPointBefore(int dayOfYear) const
{
	GENAssert(dayOfYear>=1 && dayOfYear<=365);

	if (dayOfYear>=dayOfYearFromIndex(DATAInsolationFactors::m_interpolationDaysPerYear-1))
		return DATAInsolationFactors::m_interpolationDaysPerYear-1;

    int i;
    // corrige le int i sorti de la boucle
	for (i=-1; dayOfYearFromIndex(i+1)<=dayOfYear; ++i)
	{
	}
	return (i+DATAInsolationFactors::m_interpolationDaysPerYear)%DATAInsolationFactors::m_interpolationDaysPerYear;
}
