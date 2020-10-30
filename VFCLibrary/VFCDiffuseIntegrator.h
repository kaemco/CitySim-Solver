#ifndef _INC_VFCDIFFUSEINTEGRATOR_INCLUDED
#define _INC_VFCDIFFUSEINTEGRATOR_INCLUDED

#include <set> // has a set to keep track of which elements of a large array need re-zeroing
#include <vector>

#include "SKYCell.h"
#include "GENWrappedArray.h"
#include "DATAViewFactorSetSparse.h"

/*
 * Provides the logic required for integrating cos(theta)*dw across sky patches.
 * Makes use of a GEOMETRY template class to handle the differences between the
 * visible and short-wave calculations (see e.g. VFCVISGeometry and
 * VFCSWGeometry)
 */
template <class GEOMETRY>
class VFCDiffuseIntegrator
{
public:
	VFCDiffuseIntegrator(size_t surfaceCount, const GEOMETRY *geomClass);

	/// re-zeros arrays (for new calculation)
	inline void StartCalculation();

	/// call for each cell of a patch that is obstructed
	inline void HitObstruction(const SKYCell &cell, int ObstructingSurface);

	/// call for each cell of a patch that is unobstructed
	inline void HitSky(const SKYCell &cell);

	/// save the results into patch patchNum of the RESULTSOBJECT
	void GetResults(unsigned int viewPoint, DATAViewFactorSetSparse &viewFactors, unsigned int patchNum);

private:
	 /* Since there will be FAR fewer surfaces providing obstructions to a patch than there are surfaces
	  * in the world, it is more efficient to keep track of which elements of
	  * m_cosZetadwObstructedPerSurface have been changed (and therefore need to be re-zeroed) than re-zeroing the whole lot
	  */
	std::vector< std::set<int> > m_changeTracker;

 	GENWrappedArray<float> m_obstructedVF;   // the summed obstructed view factor for this patch
	GENWrappedArray<float> m_unobstructedVF; // the summed unobstructed view factor for this patch
	GENWrappedArray<int> m_mainObstructingSurface;   // the main obstructing surface
	GENWrappedArray<GENWrappedArray<float> > m_obstructedVFPerSurface; // an array storing the obstructed v.f. contribution for each surface
									 // (in order to determine main obstructing surface)
	const GEOMETRY *m_geometry;
};

template <class GEOMETRY>
VFCDiffuseIntegrator<GEOMETRY>::VFCDiffuseIntegrator(size_t surfaceCount, const GEOMETRY *geomClass) :
	m_obstructedVF(geomClass->ViewPointCount()),
	m_unobstructedVF(geomClass->ViewPointCount()),
	m_mainObstructingSurface(geomClass->ViewPointCount()),
	m_obstructedVFPerSurface(geomClass->ViewPointCount()),
	m_geometry(geomClass)
{
	const int viewPointCount = geomClass->ViewPointCount();

	m_changeTracker.resize(viewPointCount);

	for (int i=0; i<viewPointCount; ++i)
	{
		m_obstructedVF[i]=0;
		m_unobstructedVF[i]=0;
		m_mainObstructingSurface[i]=0;
	}

	for (int i=0; i<viewPointCount; ++i)
	{
		m_obstructedVFPerSurface[i].resize(surfaceCount);
		for (unsigned int j=0; j<surfaceCount; ++j)
		{
			m_obstructedVFPerSurface[i][j]=0;
		}
	}
}

// ..........................................................................

template <class GEOMETRY>
void VFCDiffuseIntegrator<GEOMETRY>::GetResults(unsigned int viewPoint, DATAViewFactorSetSparse &viewFactors, unsigned int patchNum)
{
	viewFactors.SetViewFactors(patchNum,m_obstructedVF[viewPoint],m_unobstructedVF[viewPoint],m_mainObstructingSurface[viewPoint]);
}

// ..........................................................................

template <class GEOMETRY>
void VFCDiffuseIntegrator<GEOMETRY>::StartCalculation()
{
	// zero the running totals of the view factors
	for (unsigned int i=0; i<m_geometry->ViewPointCount(); ++i)
	{
		m_obstructedVF[i]=0;
		m_unobstructedVF[i]=0;
		m_mainObstructingSurface[i]=0;

		for (std::set<int>::iterator it=m_changeTracker[i].begin(); it!=m_changeTracker[i].end(); ++it)
		{
			m_obstructedVFPerSurface[i][*it]=0;
		}
		m_changeTracker[i].clear();
	}
}

// ..........................................................................

template <class GEOMETRY>
void VFCDiffuseIntegrator<GEOMETRY>::HitObstruction(const SKYCell &cell, int ObstructingSurface)
{
	// called when a cell is obstructed

	const float coszetadw = cell.formFactor(m_geometry->GetViewNormal());
	for (unsigned int i=0; i<m_geometry->ViewPointCount(); ++i)
	{
		if (m_geometry->IsValid(cell.direction(), i))
		{
			// add the elemental change in view factor to the running total and keep track of which surface
			// caused it
			m_obstructedVF[i]+=coszetadw;
			m_obstructedVFPerSurface[i][ObstructingSurface]+=coszetadw;

			// make a note which element of m_obstructedVFPerSurface has changed (for subsequent zeroing)
			m_changeTracker[i].insert(ObstructingSurface);

			// check whether the main obstructing surface has been changed by this hit
			if (m_obstructedVFPerSurface[i][ObstructingSurface] > m_obstructedVFPerSurface[i][m_mainObstructingSurface[i]])
				m_mainObstructingSurface[i]=ObstructingSurface;
		}
	}
}

// ..........................................................................

template <class GEOMETRY>
void VFCDiffuseIntegrator<GEOMETRY>::HitSky(const SKYCell &cell)
{
	// called when a surface has an unobstructed view of a cell
	const float coszetadw = cell.formFactor(m_geometry->GetViewNormal());
	for (unsigned int i=0; i<m_geometry->ViewPointCount(); ++i)
	{
		if (m_geometry->IsValid(cell.direction(), i))
		{
			m_unobstructedVF[i]+=coszetadw;
		}
	}
}

#endif //_INC_VFCDIFFUSEINTEGRATOR_INCLUDED
