#ifndef _INC_VFCFISHEYERENDERER_INCLUDED
#define _INC_VFCFISHEYERENDERER_INCLUDED

#include <vector>

#include "GENPoint.h"
#include "RENRenderTarget.h"  // DP modified when refactoring VFCViewFactorCalculation
#include "RENVertex.h"   // DP modified when refactoring VFCViewFactorCalculation

//class RENRenderTarget; // DP modified when refactoring VFCViewFactorCalculation

/*
 * This class uses a standard perspective renderer to implement a fish-eye view
 * of the world.  The GEOMETRY_T template argument describes what direction to
 * render views in (if it is known that a full 180 degree view is not required
 * then only views covering the required region need to be rendered).
 *
 * The reason for using a perspective renderer is that fisheye renderers are
 * very difficult to implement - straight lines need to come out as curves.
 * Instead, render enough perspective views to cover the region of interest then
 * look up in the correct image to get requested data points.
 */
template<class GEOMETRY_T, class RENDERER_T>
class VFCFishEyeRenderer
{
public:
    VFCFishEyeRenderer(std::shared_ptr<RENDERER_T> renderer);

	void  Render(const GENPoint &surfaceCentroid,
				    const GENPoint &surfaceNormal);

	unsigned int GetPixelInDirection(const GENPoint& direction);

private:
	enum { PERSPECTIVEMODE_XSIZE=256, PERSPECTIVEMODE_YSIZE=256 };


	void PrepareStorage();
	void PrepareRenderer();

	void CalculateViewDirections(const GENPoint &surfaceCentroid,
								const GENPoint &surfaceNormal,
								std::vector<GENPoint> &viewDirections);

	void RenderViews(const std::vector<GENPoint> &viewDirections);

	bool IsValidDepth(double z)
	{
		return (z>=0) && (z<=1);
	}

	GEOMETRY_T m_geometry;

    std::shared_ptr<RENDERER_T> m_renderer;

	std::vector<RENRenderTarget> m_renderTargets;
	std::vector<typename RENDERER_T::Matrix_t> m_viewMatrices;

	unsigned int m_currentView;
};

// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::VFCFishEyeRenderer(std::shared_ptr<RENDERER_T> renderer) :
    m_renderer(renderer),
	m_currentView(0)
{
	PrepareStorage();
}

// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
void  VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::PrepareStorage()
{
	m_renderTargets.resize(GEOMETRY_T::RenderCount());
	for (unsigned int i=0; i<GEOMETRY_T::RenderCount(); ++i)
	{
		m_renderTargets[i].SetSize(PERSPECTIVEMODE_XSIZE, PERSPECTIVEMODE_YSIZE);
	}

	m_viewMatrices.resize(GEOMETRY_T::RenderCount());

}
// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
void  VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::Render(
				const GENPoint &surfaceCentroid,
				const GENPoint &surfaceNormal)
{
	std::vector<GENPoint> viewDirections;
	CalculateViewDirections(surfaceCentroid,surfaceNormal,viewDirections);

	PrepareRenderer();
	RenderViews(viewDirections);
}

// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
void  VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::PrepareRenderer()
{
	m_renderer->SetPerspective(120.f,1.f,1.f,10000.f);
	m_renderer->SetEyePosition(m_geometry.GetRenderPoint());
	m_renderer->EnableBackFaceCulling(true);
}

// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
void  VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::CalculateViewDirections(
								const GENPoint &surfaceCentroid,
								const GENPoint &surfaceNormal,
								std::vector<GENPoint> &viewDirections)
{
	// Note this object's GEOMETRY_T object is different to the VFCDiffuse one
	m_geometry.Initialise(surfaceCentroid, surfaceNormal);
	m_geometry.GetRenderDirections(viewDirections);
}

// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
void  VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::RenderViews(const std::vector<GENPoint> &viewDirections)
{
	const size_t viewDirectionCount=viewDirections.size();
	for (unsigned int vd=0; vd<viewDirectionCount; ++vd)
	{
		m_renderer->SetViewTarget(viewDirections[vd]+m_geometry.GetRenderPoint());

		m_renderer->Render(m_renderTargets[vd]);
		m_viewMatrices[vd]=m_renderer->GetTransformMatrix();
	}
}

// ..........................................................................

template<class GEOMETRY_T, class RENDERER_T>
unsigned int VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::GetPixelInDirection(const GENPoint& direction)
{
	const unsigned int viewDirectionCount=static_cast<unsigned int>(m_renderTargets.size());

	const GENPoint& viewPoint=m_geometry.GetRenderPoint();
	RENVertex samplePoint=RENVertex ::CreateCartesian(direction[GENPointCoords::X]*100+viewPoint[GENPointCoords::X],
																direction[GENPointCoords::Y]*100+viewPoint[GENPointCoords::Y],
																direction[GENPointCoords::Z]*100+viewPoint[GENPointCoords::Z]);

	unsigned int loopcounter=0;

	bool newview=false;
	while(loopcounter<viewDirectionCount)
	{
		// loop counter makes sure we don't get stuck in this loop!
		++loopcounter;

		// select the next saved view making use of the fact that if previously
		// specified direction was in a given view then it is likely that
        // current direction is also in that view => i.e. assume that caller is
		// traversing through directions such that successive directions are
        // next to each other (not randomly all over the place)
		if (newview)
		{
			m_currentView=(m_currentView+1)%viewDirectionCount;
			newview = false;
		}

		// convert the ray direction into screen coordinates
		RENVertex winpos;
		m_renderer->ProjectTransform(samplePoint,m_viewMatrices[m_currentView],PERSPECTIVEMODE_XSIZE, PERSPECTIVEMODE_YSIZE,winpos);

		// round position to nearest pixel
		const int winx=static_cast<int>(winpos.m_position[0]+0.5f);
		const int winy=static_cast<int>(winpos.m_position[1]+0.5f);

		// get the pixel at the correct screen place if we hit inside the viewport
		if (m_renderTargets[0].GetViewPort().IsInside(winx,winy) &&
			IsValidDepth(winpos.m_position[2]) )
		{
			return m_renderTargets[m_currentView].GetPixel(winx,winy);
		}
		else
		{
			newview=true;
		}
	}
	// means we couldn't find the specified direction in any of the rendered views if we get here
	// this is okay as the cell might be blocked by a window
	//throw std::runtime_error("Strange error in VFCFishEyeRenderer<GEOMETRY_T,RENDERER_T>::GetPixelInDirection");
	return 0;
}


#endif //_INC_VFCFISHEYERENDERER_INCLUDED
