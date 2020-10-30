#include "GEOMBoundingSphereCalc.h"

#include "GENAssert.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

// ..........................................................................

GEOMBoundingSphereCalc::List_t::iterator previous(const GEOMBoundingSphereCalc::List_t::iterator &i)
{
	GEOMBoundingSphereCalc::List_t::iterator result=i;
	--result;
	return result;
}

GEOMSphere GEOMBoundingSphereCalc::Calculate(List_t::iterator endPoint)
{
	Recurse(endPoint);
	return m_currentSphere;
}

// ..........................................................................

void GEOMBoundingSphereCalc::Recurse(List_t::iterator endPoint)
{
	m_currentSphere=GEOMSphere(m_basis);

	GENAssert(m_basis.size()<=4);
	if (endPoint==m_points.begin() || 4==m_basis.size())
		return;

	for (List_t::iterator it=m_points.begin(); it!=endPoint;)
	{
		// it may be invalidated by MoveToFrontOfList, so increment here
		List_t::iterator thisIt=it;
		++it;

		if (!m_currentSphere.PointIsInside(*thisIt))
		{
			m_basis.push_back(*thisIt);
			Recurse(thisIt);
			m_basis.pop_back();
			MoveToFrontOfList(thisIt);
		}
	}
}

// ..........................................................................

void GEOMBoundingSphereCalc::MoveToFrontOfList(List_t::iterator p)
{
	 m_points.splice (m_points.begin(), m_points, p);
}

