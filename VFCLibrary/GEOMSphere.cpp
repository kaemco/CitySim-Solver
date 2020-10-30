#include "GEOMSphere.h"


#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

// ..........................................................................

GEOMSphere::GEOMSphere(const GEOMSphericalBasis &basis) :
	m_sphericalBasis(basis)
{}

GEOMSphere::GEOMSphere()
{}


// ..........................................................................

GENPoint GEOMSphere::Centre() const 
{ 
	return m_sphericalBasis.Centre(); 
}

// ..........................................................................

float GEOMSphere::Radius() const 
{ 
	return m_sphericalBasis.Radius(); 
}

// ..........................................................................

bool GEOMSphere::PointIsInside(const GENPoint &p, double distanceSquaredTolerance) const
{
	const GENPoint dp=p-m_sphericalBasis.Centre();
	const float rpsq=dot_product(dp,dp);
	const float radius=m_sphericalBasis.Radius();
	return rpsq<=radius*radius+distanceSquaredTolerance;
}

