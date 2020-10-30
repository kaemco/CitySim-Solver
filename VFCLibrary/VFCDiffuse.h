#ifndef _INC_VFCDIFFUSE_INCLUDED
#define _INC_VFCDIFFUSE_INCLUDED

#include <vector>

#include "VFCFishEyeRenderer.h"
#include "SKYVaultABC.h"		// gives the sky vault geometry
#include "VFCDiffuseIntegrator.h"

#ifdef OPENGL
#include "RENOpenGLRenderer.h"
#define RENRenderer RENOpenGLRenderer
#else
#include "RENSoftwareRenderer.h"
#define RENRenderer RENSoftwareRenderer
#endif

/*
 * Uses VFCDiffuseRenderer and VFCDiffuseIntegrator to calculate view
 * factors.  The GEOMETRY_T template argument deals with the differences between
 * the visible and shortwave calculates (the differences are all geometrical in
 * terms of are we claculating from inside a room, etc., etc.).  See for example
 * VFCSWGeometry and VFCVISGeometry.
 *
 */
template<class GEOMETRY_T, class RENDERER_T=RENRenderer>
class VFCDiffuse
{
public:
    VFCDiffuse(std::shared_ptr<SKYVault> skyGeom, size_t surfaceCount, std::shared_ptr<RENDERER_T> renderer);

	/// Calculate view factors for given centroid and normal
	void  Calculate(const GENPoint &surfaceCentroid,
				    const GENPoint &surfaceNormal,
					std::vector<DATAViewFactorSetSparse*> &retVals);

//private:
	/* Attention! Construction order here! In VFCDiffuse constructor
     * m_geometry gets passed to m_integrator's constructor.  Currently this is
     * not a problem as integrator only calls static methods on geometry, but if
     * this changes then geometry must be constructed first
     */
	GEOMETRY_T m_geometry;

	// gives the layout of sky patches and cells, not owned by this class
    std::shared_ptr<SKYVault> m_theSky;

	VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T> m_fisheyeRenderer;
	VFCDiffuseIntegrator<GEOMETRY_T> m_integrator;
};

template<class GEOMETRY_T, class RENDERER_T>
 VFCDiffuse<GEOMETRY_T,RENDERER_T>::VFCDiffuse(std::shared_ptr<SKYVault> skyGeom,
														size_t surfaceCount,
                                                         std::shared_ptr<RENDERER_T> renderer) :
    m_theSky(skyGeom),
	m_fisheyeRenderer(renderer),
    m_integrator(surfaceCount, &m_geometry)
{
}

// ..........................................................................
template<class GEOMETRY_T, class RENDERER_T>
void  VFCDiffuse<GEOMETRY_T,RENDERER_T>::Calculate(
				const GENPoint &surfaceCentroid,
				const GENPoint &surfaceNormal,
				std::vector<DATAViewFactorSetSparse*> &retVals)
{
	/* ALGORITHM:
	 *	Initialise ViewFactorSet (tell it how many patches there are, clear any existing values)
	 *	Render views
	 *	For each sky patch
	 *		Initialise Integrator
	 *		Do the integrating
	 *		Store the result
	 *	Loop
	 */
	m_geometry.Initialise(surfaceCentroid, surfaceNormal);
	m_fisheyeRenderer.Render(surfaceCentroid,surfaceNormal);

	const GENPoint viewNormal=m_geometry.GetViewNormal();

	// for each patch
	const unsigned int patchCount=m_theSky->PatchCount();
	for (unsigned int i=0; i<patchCount; i++)
	{
		// re-initialise the integrator
		m_integrator.StartCalculation();

		const SKYPatch& thisPatch=*m_theSky->GetPatch(i);

		// for each Cell in the patch
		const unsigned int cellCount=thisPatch.cellCount();
		for (unsigned int j=0; j<cellCount; j++)
		{
			const SKYCell& thisCell=thisPatch.getCell(j);

			// is the cell potentially visible?
			if ( dot_product(thisCell.direction(),
								viewNormal) > 0)
			{
				const GENPoint& cellDir = thisCell.direction();

				const unsigned int pixelColour=m_fisheyeRenderer.GetPixelInDirection(cellDir);

				if (pixelColour>0)
				{
					// subtract 1 from r for consistency with other (surfaces
                    // are numbered from 0, but here 0 represents "no
                    // obstruction")
					m_integrator.HitObstruction(thisCell, (int)pixelColour-1);
				}
				else
				{
					m_integrator.HitSky(thisCell);
				}
			}
		}
		for (unsigned int vp=0; vp<m_geometry.ViewPointCount(); ++vp)
			m_integrator.GetResults(vp,*retVals[vp],i);
	}
}

#endif //_INC_VFCDIFFUSE_INCLUDED
