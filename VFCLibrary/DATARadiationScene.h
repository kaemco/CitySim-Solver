#ifndef _INC_DATARadiationScene_INCLUDED
#define _INC_DATARadiationScene_INCLUDED

// SYSTEM INCLUDES
//
#include <vector>
#include <memory>
#include <list>
#include <map>
#include <limits>

// PROJECT INCLUDES
//


// LOCAL INCLUDES
//
#include "SKYSiteLocation.h"
#include "GENHandle.h"
#include "DATASurfaceDelegateABC.h"
#include "GENPoint.h"

// FORWARD REFERENCES
//
class RENIndexedFaceSet;
class DATASurface;
class RENBoundingSphereABC;
class DATASurfaceIterator;
class DATASurfaceBuildingIterator;

/*
 * Stores a set of DATASurface objects for the radiation model and also the
 * mesh used in the view factor calculations.  Methods are provided to access
 * either all of the DATASurfaces or just the DATASurfaceBuilding
 * objects.
 *
 */
class DATARadiationScene
{
public:
	DATARadiationScene();
	virtual ~DATARadiationScene();

	void AddBuildingSurface(GENHandle<DATASurfaceDelegateABC> surfaceSegment);
	void AddGroundSurface(GENHandle<DATASurfaceDelegateABC> surfaceSegment);

	/// Add a point for which incident diffuse radiation will be calculated, but which does not cause shading, etc.
	void AddDiffuseSamplingPoint(GENHandle<DATASurfaceDelegateABC> surfaceSegment);
	// JK - 18.04.2015 - sets m_surfaceIsPartOfRadiationMesh to false

	const DATASurface& GetSurface(unsigned int index) const;
	int GetMeshIndex(const DATASurface &surface);

    // gets smallest surface, added by JK - 24/03/09
    const DATASurface& GetSmallestSurface() const;

	DATASurfaceBuildingIterator GetBuildingSurfaces() const;
	DATASurfaceIterator GetAllSurfaces() const;

	unsigned int SurfaceCount() const;
	unsigned int BuildingSurfaceCount() const;

    std::unique_ptr<RENBoundingSphereABC> GetBoundingSphere() const;
    std::shared_ptr<RENIndexedFaceSet> IndexedFaceSet() const;

	void SetLocation(const SKYSiteLocation& location);
	const SKYSiteLocation& GetLocation() const;

private:
	DATARadiationScene(const DATARadiationScene &s);
	DATARadiationScene& operator=(const DATARadiationScene &s);

	void DeleteStoredSurfaces();

	unsigned int m_buildingSurfaceCount;

	std::vector<DATASurface*> m_surfaces;
	std::vector<bool> m_surfaceIsPartOfRadiationMesh;

	SKYSiteLocation m_siteLocation;
};

#endif //_INC_DATARadiationScene_INCLUDED
