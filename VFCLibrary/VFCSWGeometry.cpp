#include "./VFCSWGeometry.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER


void VFCSWGeometry::Initialise(const GENPoint &externalWallCentroid, 
									const GENPoint &externalWallNormal)
{
	m_externalWallNormal=externalWallNormal;
	m_externalWallCentroid=externalWallCentroid;
}

// ..........................................................................

void VFCSWGeometry::GetViewPoint(unsigned int /*viewPointNo*/, 
									  GENPoint &viewPointDir, 
									  GENPoint &viewPointPos) const
{
	viewPointDir=m_externalWallNormal;
	viewPointPos=m_externalWallCentroid;
}

// ..........................................................................

const GENPoint& VFCSWGeometry::GetViewNormal() const
{
	return m_externalWallNormal;
}

// ..........................................................................

void VFCSWGeometry::GetRenderDirections(std::vector<GENPoint>& viewDirections) const
{
	viewDirections.resize(RenderCount());

	GENPoint up=GENPoint::Cartesian(0,0,1);

	// special case of a horizontal surface
	// TODO NOTE ALSO NEED CASE OF SURFACE FACING GROUND
	if (ApproxEqual(m_externalWallNormal,up) )
	{
		up=GENPoint::Cartesian(0,1,0);
	}

	// set a to be 'left or right'
	GENPoint a = cross_product(up, m_externalWallNormal);

	// now make up perpendicular to a and surf norm
	up=cross_product(m_externalWallNormal,a);

	// for each of up and down
	int i=0;
	for (int updown=-1; updown<2; updown+=2)
	{
		// for each of left and right
		for (int leftright=-1; leftright<2; leftright+=2)
		{

			viewDirections[i]=GENPoint::UnitVector((m_externalWallNormal +a*(float)leftright));
			viewDirections[i]+=up*(float)updown;
			++i;
		}
	}
	return;
}

// ..........................................................................

const GENPoint& VFCSWGeometry::GetRenderPoint() const
{
	return m_externalWallCentroid;
}
