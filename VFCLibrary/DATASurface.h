#ifndef _INC_DATASURFACE_INCLUDED
#define _INC_DATASURFACE_INCLUDED

#include "GENHandle.h"
#include "GENPoint.h"
#include "DATASurfaceDelegateABC.h"
#include "DATAViewFactorSetSparse.h"
#include "DATAInsolationFactors.h"

/* The view factor calc uses this class to represent a surface - it stores view
 * factors in this class and access info (vertices, area, etc) through this
 * class
 */
class DATASurface
{
public:
	DATASurface(GENHandle<DATASurfaceDelegateABC> surface, const GENPoint& normal, const GENPoint& centroid, bool isBuildingSurface);

	virtual ~DATASurface();

	const GENPoint& Normal() const;
	const GENPoint& Centroid() const;
	double Area() const;

	bool IsBuildingSurface() const { return m_isBuildingSurface; }

	virtual DATASurfaceDelegateABC* SurfaceDelegate();
	virtual const DATASurfaceDelegateABC* SurfaceDelegate() const;

	DATAViewFactorSetSparse& DaylightViewFactors(unsigned int viewPoint) { return m_daylightVFs[viewPoint]; }
	const DATAViewFactorSetSparse& DaylightViewFactors(unsigned int viewPoint) const  { return m_daylightVFs[viewPoint]; }

	DATAViewFactorSetSparse& SWViewFactors(unsigned int /*viewPoint*/) { return m_SWVFs; }
	const DATAViewFactorSetSparse& SWViewFactors() const { return m_SWVFs; }

	DATAInsolationFactors& InsolationFactors() { return m_insolationFactors; }
	const DATAInsolationFactors& InsolationFactors() const { return m_insolationFactors; }

private:
	DATASurface(const DATASurface& s);

	DATAViewFactorSetSparse m_daylightVFs[2];
	DATAViewFactorSetSparse m_SWVFs;
	DATAInsolationFactors m_insolationFactors;

	GENHandle<DATASurfaceDelegateABC> m_surface;

	GENPoint m_normal;
	GENPoint m_centroid;

	bool m_isBuildingSurface;
};

#endif //_INC_DATASURFACE_INCLUDED
