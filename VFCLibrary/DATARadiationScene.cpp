#include <vector>

#include "./DATARadiationScene.h"

#include <list>
#include "DATASurfaceDelegateABC.h"
#include "DATASurface.h"
#include "DATASurfaceIterator.h"

#include "GEOMPolygonInfo.h"

#include "RENIndexedFaceSetBuilder.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

using namespace std;

// ..........................................................................

DATARadiationScene::DATARadiationScene() :
	m_buildingSurfaceCount(0)
{
}

// ..........................................................................

DATARadiationScene::~DATARadiationScene()
{
	DeleteStoredSurfaces();
}

// ..........................................................................

DATARadiationScene::DATARadiationScene(const DATARadiationScene &)
{
}

// ..........................................................................

DATARadiationScene& DATARadiationScene::operator=(const DATARadiationScene &/*s*/)
{
	throw std::string("not implemented");
}


// ..........................................................................

struct VertexSet : public std::vector<GENPoint>, public DATASurfaceDelegateABC::VertexVisitor
{
	void operator()(const GENPoint &p)
	{
		push_back(p);
	}
};

GEOMPolygonInfo CalculatePolygonInfo(GENHandle<DATASurfaceDelegateABC> surfaceDelegate)
{
	VertexSet vertices;
	vertices.reserve(surfaceDelegate->vertexCount());
	surfaceDelegate->sendVertices(vertices);

	return GEOMPolygonInfo(vertices.begin(), vertices.end(),surfaceDelegate->normal());
}

// ..........................................................................

void DATARadiationScene::AddBuildingSurface(GENHandle<DATASurfaceDelegateABC> surfaceSegment)
{
	const GEOMPolygonInfo& info=CalculatePolygonInfo(surfaceSegment);

	m_surfaces.push_back(new DATASurface(surfaceSegment,
										info.normal(),
										info.centroid(),
										true));
	m_surfaceIsPartOfRadiationMesh.push_back(true);
	++m_buildingSurfaceCount;
}

// ..........................................................................

void DATARadiationScene::AddGroundSurface(GENHandle<DATASurfaceDelegateABC> surfaceSegment)
{
	const GEOMPolygonInfo& info=CalculatePolygonInfo(surfaceSegment);

	m_surfaces.push_back(new DATASurface(surfaceSegment,
										info.normal(),
										info.centroid(),
										false));
	m_surfaceIsPartOfRadiationMesh.push_back(true);
}

// ..........................................................................

void DATARadiationScene::AddDiffuseSamplingPoint(GENHandle<DATASurfaceDelegateABC> surfaceDelegate)
{
	VertexSet vertices;
	vertices.reserve(surfaceDelegate->vertexCount());
	surfaceDelegate->sendVertices(vertices);

	GEOMPolygonInfo info(vertices.begin(), vertices.end());

	m_surfaces.push_back(new DATASurface(surfaceDelegate,
										surfaceDelegate->normal(),
										info.centroid(),
										false));
	m_surfaceIsPartOfRadiationMesh.push_back(false);
}

// ..........................................................................

std::shared_ptr<RENIndexedFaceSet> DATARadiationScene::IndexedFaceSet() const
{
    std::shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet);
    RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);

	unsigned int colour=1;
	std::vector<bool>::const_iterator addToMeshIt=m_surfaceIsPartOfRadiationMesh.begin();
	for (DATASurfaceIterator surfaces=GetAllSurfaces();
		!surfaces.isAtEnd();
		++surfaces, ++colour, ++addToMeshIt)
	{
		if (*addToMeshIt && (*surfaces).Area()>1e-4)
			meshBuilder.AddSurface(surfaces->SurfaceDelegate(),colour);
	}
	return indexedFaceSet;
}

// ..........................................................................

std::unique_ptr<RENBoundingSphereABC> DATARadiationScene::GetBoundingSphere() const
{
    std::shared_ptr<RENIndexedFaceSet> faceSet(IndexedFaceSet());
	return faceSet->GetBoundingSphere();
}

// ..........................................................................

DATASurfaceBuildingIterator DATARadiationScene::GetBuildingSurfaces() const
{
	return DATASurfaceBuildingIterator(m_surfaces.begin(), m_surfaces.end());
}

// ..........................................................................

DATASurfaceIterator DATARadiationScene::GetAllSurfaces() const
{
	return DATASurfaceIterator(m_surfaces.begin(), m_surfaces.end());
}

// ..........................................................................

unsigned int DATARadiationScene::SurfaceCount() const
{
	return static_cast<unsigned int>(m_surfaces.size());
}

// ..........................................................................

unsigned int DATARadiationScene::BuildingSurfaceCount() const
{
	return m_buildingSurfaceCount;
}

// ..........................................................................

void DATARadiationScene::SetLocation(const SKYSiteLocation& location)
{
	m_siteLocation=SKYSiteLocation(location);
}

// ..........................................................................

const SKYSiteLocation& DATARadiationScene::GetLocation() const
{
	return m_siteLocation;
}

// ..........................................................................

void DATARadiationScene::DeleteStoredSurfaces()
{
	for (unsigned int i=0; i<m_surfaces.size(); ++i)
	{
		delete m_surfaces[i];
	}
}

// ..........................................................................

int DATARadiationScene::GetMeshIndex(const DATASurface &surface)
{
	for (unsigned int i=0; i<m_surfaces.size(); ++i)
	{
		if (&surface==m_surfaces[i])
			return i;
	}
	return -1;
}

// ..........................................................................

const DATASurface& DATARadiationScene::GetSurface(unsigned int index) const
{
	return *(m_surfaces[index]);
}

// ..........................................................................

const DATASurface& DATARadiationScene::GetSmallestSurface() const
{
    unsigned int index=numeric_limits<unsigned int>::signaling_NaN();
    double smallestSurface = std::numeric_limits<double>::max();
	for (unsigned int i=0; i<m_surfaces.size(); ++i)
	{
        if (m_surfaces[i]->SurfaceDelegate()->getRadius() < smallestSurface) { smallestSurface = m_surfaces[i]->SurfaceDelegate()->getRadius(); index = i; }
	}
	return *(m_surfaces[index]);
}
