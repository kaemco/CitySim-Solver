#include "./VFCVISGeometry.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

VFCVISGeometry::VFCVISGeometry(void) :
	m_roomMidPoint(GENPoint::Origin()),
	m_externalWallNormal(GENPoint::Origin()),
	m_externalWallCentroid(GENPoint::Origin()),
	m_s(GENPoint::Origin()),
	m_t(GENPoint::Origin())
{
}

VFCVISGeometry::~VFCVISGeometry(void)
{
}

void VFCVISGeometry::Initialise(const GENPoint &externalWallCentroid, const GENPoint &externalWallNormal)
{
	// store the external surface normal
	m_externalWallCentroid=externalWallCentroid;
	m_externalWallNormal=GENPoint::UnitVector(externalWallNormal);

	GENPoint tmp=m_externalWallNormal;
	tmp[GENPointCoords::Z]=0;
	if (dot_product(tmp,tmp)==0)
	{
		// TODO view direction is vertically up, just do something vaguely "sensible" for now
		tmp[GENPointCoords::Y]=-1;
	}
	GENPoint horizAwayFromWall=GENPoint::UnitVector(tmp);

	// now find the s and t vectors
	GENPoint up=GENPoint::Cartesian(0,0,1);
	GENPoint right = cross_product(up, m_externalWallNormal);

	m_t = cross_product(m_externalWallNormal,right);
	m_s=-right;

	// calculate where the view points are
	// first clear any existing points
	m_viewPoints.clear();
	m_lamdanumerators.clear();

	// find the midpoint of the room
	m_roomMidPoint=m_externalWallCentroid - horizAwayFromWall*3.f;
	m_roomMidPoint[GENPointCoords::Z]-=.7f;

	for (float viewPointDistance=1.5; viewPointDistance<5; viewPointDistance+=3)
	{
		m_viewPoints.push_back(m_externalWallCentroid - horizAwayFromWall*viewPointDistance);

		// centroid of room's external wall is at surface centroid, => viewpoint is .7 m below this (i.e. .8 m from floor)
		m_viewPoints[m_viewPoints.size()-1][GENPointCoords::Z]-=.7f;

		// finally pre-calculate the numerator for the calculation of plane intersection point
		float d = dot_product(m_externalWallCentroid, m_externalWallNormal); // plane eqn ax+by+cz=d
		m_lamdanumerators.push_back(d-(float)dot_product(m_externalWallNormal,m_viewPoints[m_viewPoints.size()-1]));
	}

//	std::cerr << "ext wall centroid (x,y,z): (" << m_externalWallCentroid[0] << "," << m_externalWallCentroid[1] << "," << m_externalWallCentroid[2] << ")" << std::endl;
//	std::cerr << "room mid point (x,y,z): (" << m_roomMidPoint[0] << "," << m_roomMidPoint[1] << "," << m_roomMidPoint[2] << ")" << std::endl;
//	std::cerr << "1st view point (x,y,z): (" << m_viewPoints[0][0] << "," << m_viewPoints[0][1] << "," << m_viewPoints[0][2] << ")" << std::endl;
//	std::cerr << "2nd view point (x,y,z): (" << m_viewPoints[1][0] << "," << m_viewPoints[1][1] << "," << m_viewPoints[1][2] << ")" << std::endl;
//
//    std::string essai;
//    std::cin >> essai;

}

bool VFCVISGeometry::IsValid(const GENPoint &rayDir, unsigned int viewPointNo) const
{
	// first find point where ray intersects with plane of the window
	float lamdadenominator=(float)dot_product(m_externalWallNormal,rayDir);

	/* Test not really required, since we can't call this routine unless a patch is being tested,
	 * and we don't test patches that are not visible to the external surface
	// denominator = zero implies no intersection point
	// denominator < zero implies patch/cell is in opposite direction to window wall
	if (lamdadenominator<=0)
		return false;*/

	float lamda = m_lamdanumerators[viewPointNo]/lamdadenominator;
	GENPoint inter = m_viewPoints[viewPointNo] + rayDir * lamda;

	// now find alpha and beta - the components of s and t required to move from window centroid to intersection point
	GENPoint imc=inter-m_externalWallCentroid;

	const float alpha = dot_product(imc,m_s);
	const float beta = dot_product(imc,m_t);

	// return true if ray hits within bounds of window
	return ( (-0.7f <= beta) && (beta <= 1.5f) && (-4.5f <= alpha) && ( alpha <= 4.5f) );
}

void VFCVISGeometry::GetRenderDirections(std::vector<GENPoint>& viewDirections) const
{
	viewDirections.resize(RenderCount());

	GENPoint up=m_externalWallCentroid - m_viewPoints[1];
	GENPoint viewNormal=GENPoint::Cartesian(0,0,1);

	// set a to be 'left or right'
	GENPoint a = cross_product(up, viewNormal);

	// now make up perpendicular to a and surf norm
	up=cross_product(viewNormal,a);

	int i=0;
	// for each of left and right
	for (int leftright=-1; leftright<2; leftright+=2)
	{
		/// \todo - this fragment is common with SWGeometry
		viewDirections[i]=GENPoint::UnitVector((viewNormal +a*(float)leftright));
		viewDirections[i]+=up;
		++i;
	}
}

void VFCVISGeometry::GetViewPoint(unsigned int viewPointNo,
									GENPoint &viewPointDir,
									GENPoint &viewPointPos) const
{
	viewPointDir = GENPoint::Cartesian(0,0,1);
	viewPointPos = m_viewPoints[viewPointNo];
}

GENPoint VFCVISGeometry::GetViewNormal() const
{
	return GENPoint::Cartesian(0,0,1);
}

const GENPoint& VFCVISGeometry::GetRenderPoint() const
{
	return m_viewPoints[1];
}
