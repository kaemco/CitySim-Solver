#ifndef _INC_VFCVISGEOMETRY_INCLUDED
#define _INC_VFCVISGEOMETRY_INCLUDED

#include <vector>
#include "GENPoint.h"

/*
 * This class performs all of the geometrical tasks necessary for the visible diffuse
 * view factor calculation.  This class and VFCSWGeometry are used to
 * separate out the shortwave and daylight specific aspects of diffuse view
 * factor calculations from the main algorithm.
 *
 */
class VFCVISGeometry
{
public:
	VFCVISGeometry(void);
	~VFCVISGeometry(void);

	/// initialise the geometry of the window
	void Initialise(const GENPoint &externalWallCentroid, const GENPoint &externalWallNormal);

	/// test if the ray in the given direction passes through the window (starting from the viewpoint)
	bool IsValid(const GENPoint &rayDir, unsigned int viewPointNo) const;

	/// Calculate the view directions which must be rendered
	void GetRenderDirections(std::vector<GENPoint>& viewDirections) const;

	/// Return a given view point within the room
	void GetViewPoint(unsigned int viewPointNo, GENPoint &viewPointDir, GENPoint &viewPointPos) const;

	/// Get the view point normal
	GENPoint GetViewNormal() const;

	/// Return the point from which views should be rendered
	const GENPoint& GetRenderPoint() const;

	/// Return the number of view points
	static unsigned int ViewPointCount() { return 2; }

	/// Return the number of render directions needed for this calc
	static unsigned int RenderCount() { return 2; }

private:
	// the view points (in world coordinates)
	std::vector<GENPoint> m_viewPoints;
	GENPoint m_roomMidPoint;

	// external surface normal and centroid
	GENPoint m_externalWallNormal;
	GENPoint m_externalWallCentroid;

	// s and t - t is the external wall's 'up' vector, s is 'right' (when viewed from inside room)
	GENPoint m_s,m_t;

	// used to find intersection point of ray with window plane
	// (precalculated for efficiency)
	std::vector<float> m_lamdanumerators;
};

#endif //_INC_VFCVISGEOMETRY_INCLUDED
