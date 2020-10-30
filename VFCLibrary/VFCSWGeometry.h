#ifndef _INC_VFCSWGEOMETRY_INCLUDED
#define _INC_VFCSWGEOMETRY_INCLUDED

#include <vector>
#include "GENPoint.h"

/*
 * This class performs all of the geometrical tasks necessary for the shortwave diffuse
 * view factor calculation.  This class and VFCVISGeometry are used to
 * separate out the shortwave and daylight specific aspects of diffuse view
 * factor calculations from the main algorithm.

 */
class VFCSWGeometry
{
public:
	VFCSWGeometry(void) :
	  m_externalWallNormal(GENPoint::Origin()),
	  m_externalWallCentroid(GENPoint::Origin())
	{}

	/// initialise the geometry
	void Initialise(const GENPoint &externalWallCentroid,
					const GENPoint &externalWallNormal);

	/// All rays are valid in the SW calculation
	bool IsValid(const GENPoint /*&rayDir*/, unsigned int /*viewPointNo*/) const { return true; }

	/// Calculate the view directions which must be rendered
	void GetRenderDirections(std::vector<GENPoint>& viewDirections) const;

	/// Return a given view point within the room
	void GetViewPoint(unsigned int viewPointNo,
						GENPoint &viewPointDir,
						GENPoint &viewPointPos) const;

	/// Return the point from which views should be rendered
	const GENPoint& GetRenderPoint() const;

	/// Get the view point normal
	const GENPoint& GetViewNormal() const;

	/// Return the number of view points
	static unsigned int ViewPointCount() { return 1; }

	/// Return the number of render directions needed for this calc
	static unsigned int RenderCount() { return 4; }

private:
	// external surface normal and centroid
	GENPoint m_externalWallNormal;
	GENPoint m_externalWallCentroid;
};

#endif //_INC_VFCSWGEOMETRY_INCLUDED
