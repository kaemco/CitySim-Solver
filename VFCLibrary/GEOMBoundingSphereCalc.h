#ifndef _INC_GEOMBOUNDINGSPHERECALC_INCLUDED
#define _INC_GEOMBOUNDINGSPHERECALC_INCLUDED

#include <list>
#include <vector>

#include "GENPoint.h"
#include "GEOMSphere.h"

/*
 * Calculate a bounding sphere using algorithm from "Smallest enclosing disks
 * (balls and ellipsoids)", by Emo Welzl
 */
class  GEOMBoundingSphereCalc
{
public:
	typedef std::list<GENPoint> List_t;
	typedef std::list<GENPoint>::iterator PointIt_t;
	typedef std::vector<const GENPoint*> Vector_t; 

	template <class ITERATOR>
	GEOMSphere Calculate(const ITERATOR pointBegin, const ITERATOR pointEnd);

private:
	GEOMSphere Calculate(List_t::iterator endPoint);
	void MoveToFrontOfList(List_t::iterator p);
	void Recurse(List_t::iterator endPoint);

	List_t m_points;
	GEOMSphericalBasis m_basis;	
	GEOMSphere m_currentSphere;
};

// ..........................................................................

template <class ITERATOR>
GEOMSphere GEOMBoundingSphereCalc::Calculate(const ITERATOR pointBegin, const ITERATOR pointEnd)
{
	for(ITERATOR it = pointBegin; it !=pointEnd ; ++it)
	{
		m_points.push_back(*it);
	}

	GEOMSphere MB = Calculate(m_points.end());

	return MB;
}


#endif // _INC_GEOMBOUNDINGSPHERECALC_INCLUDED
