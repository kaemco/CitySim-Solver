

#include "DATASurface.h"
#include "GENPoint.h"

DATASurface::DATASurface(GENHandle<DATASurfaceDelegateABC> surface,
								   const GENPoint& normal,
								   const GENPoint& centroid,
								   bool isBuildingSurface) :
	m_surface(surface),
	m_normal(normal),
	m_centroid(centroid),
	m_isBuildingSurface(isBuildingSurface)
{
}

DATASurface::DATASurface(const DATASurface& )
{
	GENAssert(true==false);
}

DATASurface::~DATASurface()
{
}

const GENPoint& DATASurface::Normal() const
{
	return m_normal;
}

const GENPoint& DATASurface::Centroid() const
{
	return m_centroid;
}

double DATASurface::Area() const
{
	return m_surface->getArea();
}

DATASurfaceDelegateABC* DATASurface::SurfaceDelegate()
{
	return &*m_surface;
}

const DATASurfaceDelegateABC* DATASurface::SurfaceDelegate() const 
{
	return &*m_surface;
}

