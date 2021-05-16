#ifndef _INC_DATAINSOLATIONFACTORS_INCLUDED
#define _INC_DATAINSOLATIONFACTORS_INCLUDED

#include <algorithm>
#include <vector>

#include "GENPoint.h"

class SKYSiteLocation;
class SKYSun;

/**
 * This class stores view factors to the sun for one surface.
 *
 */
class DATAInsolationFactors
{
public:
	DATAInsolationFactors();

	void GetSunPositionsToRender(const SKYSiteLocation &location,
										std::vector<GENPoint> &sunPositions) const;

	void SetInsolationFactor(int sunPositionIndex, float vf);
	float GetInsolationFactor(const SKYSun &Sun) const;
	float GetInsolationFactor_linear_interp(const SKYSun &Sun) const;

	// Public for CIBSE tests - not ideal
    float GetIndexedViewFactor(int dayIndex, int hourIndex) const;
	float GetIndexedViewFactor(unsigned int sunIndex) const;

	// use to change number of days/time points used in view factor interpolation
	void SetInterpolationScheme(int daysPerYear, int stepsPerDay);

private:

	float LagrangianInterpolation(float x, float x1, float y1, float x2, float y2, float x3, float y3) const;
	float hourFromIndex(int hourIndex) const;
	int indexForTimePointBefore(float solarTime) const;
	int dayOfYearFromIndex(int dayIndex) const;
	int indexForDayPointBefore(int dayOfYear) const;

	std::vector<float> m_insolationFactors;

	int m_interpolationDaysPerYear;
	int m_storedInterpolationDaysPerYear;
	int m_interpolationPeriodsPerDay;
};


#endif //_INC_DATAINSOLATIONFACTORS_INCLUDED
