#ifndef _INC_GEOMSPHERE_INCLUDED
#define _INC_GEOMSPHERE_INCLUDED

#include <vector>

#include "GENPoint.h"
#include "GEOMSphericalBasis.h"

class GEOMSphere
{
public:
	GEOMSphere();
	GEOMSphere(const GEOMSphericalBasis &basis);

	GENPoint Centre() const;
	float Radius() const;

	bool PointIsInside(const GENPoint &p, double distanceSquaredTolerance=1e-3) const;
private:
	GEOMSphericalBasis m_sphericalBasis;
};

#endif // _INC_GEOMSPHERE_INCLUDED
