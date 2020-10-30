#include <algorithm>
#include <string.h>

#include "./VFCDirect.h"
#include "../surface.h"

#ifdef OPENGL
#include "RENOpenGLRenderer.h"
#define RENRenderer RENOpenGLRenderer
#else
#include "RENSoftwareRenderer.h"
#define RENRenderer RENSoftwareRenderer
#endif

#include "DATARadiationScene.h"
#include "DATASurface.h"
#include "DATASurfaceIterator.h"

#include "RENBoundingSphereABC.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

#include <cstdio>
#include <cstring>

using std::vector;
using std::list;
using std::min;

// ..........................................................................

void VFCDirect::Initialise(RenderClass_t *renderer,
					 DATARadiationScene &scene)
{
	m_scene=&scene;
	m_renderer=renderer;

	// array to hold the number of pixels taken up by each surface
	m_pixelSum.resize(m_scene->SurfaceCount());

	// set up the view parameters
    std::unique_ptr<RENBoundingSphereABC> boundingSphere=m_scene->GetBoundingSphere();
	const float modelRadius=boundingSphere->Radius();
	const GENPoint renderCentre=boundingSphere->Centre();
	std::cerr << "Scene bounding sphere centre: (" << renderCentre << ")" << " radius: " << modelRadius << std::endl;

    // get the smallest surface given by vertices, JK - 24.03.09 in order to compute pixels per metre
    if (m_scene->GetSmallestSurface().SurfaceDelegate()->getRadius()<=0.f) {
        throw (string("In surface id=") + toString(((Surface*)m_scene->GetSmallestSurface().SurfaceDelegate())->getId()) + ": smallest surface radius is zero.");
	}
	#ifndef WS
    std::cerr << "Smallest surface radius: " << m_scene->GetSmallestSurface().SurfaceDelegate()->getRadius() << std::endl;
    #endif

    // previous value was fixed to 12 pixelsPerMetre, not adaptive to the smallestSurface
	const static float pixelsPerMetre=12.f/m_scene->GetSmallestSurface().SurfaceDelegate()->getRadius(); //12; // aim for 12 pixels per smallest surface
	#ifndef WS
    std::cerr << "Pixels per metre: " << pixelsPerMetre << std::endl;
    #endif

	m_renderWindowDim=std::min(static_cast<unsigned int>(pixelsPerMetre*modelRadius),512u); // TODO: change to 2048u
    std::cerr << "Direct renderer window dimension: " << m_renderWindowDim << std::endl;

	m_renderTarget.SetSize(m_renderWindowDim,m_renderWindowDim);
	m_renderer->SetOrthogonal(modelRadius*2,modelRadius*2,1+2*modelRadius,1);
	m_renderer->EnableBackFaceCulling(true);

	m_viewingSystem.reset(new VFCDirectViewSystem(renderCentre, modelRadius));
	m_renderer->SetViewTarget(m_viewingSystem->ViewTarget());

	// area of one pixel
	m_pixelArea=(modelRadius*modelRadius*4)/(m_renderWindowDim*m_renderWindowDim);

}

// ..........................................................................

void VFCDirect::CalculateInsolationFactors(const GENPoint &sunPosition,
				const unsigned int sunPositionIndex)
{
	if (sunPosition[GENPointCoords::Z]<=0)
	{
		ZeroInsolationFactors(sunPositionIndex);
		return;
	}

	UpdatePixelSum(sunPosition);

	for (DATASurfaceIterator surfaces=m_scene->GetAllSurfaces();
		!surfaces.isAtEnd();
		++surfaces)
	{
		DATASurface &surface=*surfaces;
		const int surfaceIndex=m_scene->GetMeshIndex(surface);

		if (surfaceIndex>=0)
		{
			const float projectedArea=(float)surface.Area() * dot_product(sunPosition,surface.Normal());
			if (projectedArea>0)
			{
				const float insolationFraction=min(m_pixelArea*m_pixelSum[surfaceIndex]/projectedArea,1.f);
				surface.InsolationFactors().SetInsolationFactor(sunPositionIndex,insolationFraction);

				//std::cerr << "surfaceIndex: " << surfaceIndex << "\tprojectedArea: " << projectedArea << "\tinsolationFraction: " << insolationFraction << "\tm_pixelSum[surfaceIndex]: " << m_pixelSum[surfaceIndex] << std::endl;


			}
			else
			{
				surface.InsolationFactors().SetInsolationFactor(sunPositionIndex,0);
			}
		}
	}
}

// ..........................................................................

void VFCDirect::ZeroInsolationFactors(const unsigned int sunPositionIndex)
{
	for (DATASurfaceIterator surfaceIt=(*m_scene).GetAllSurfaces();
		!surfaceIt.isAtEnd();
		++surfaceIt)
	{
		surfaceIt->InsolationFactors().SetInsolationFactor(
							sunPositionIndex,0);
	}
}

// ..........................................................................

void VFCDirect::UpdatePixelSum(const GENPoint &sunPosition)
{
	Render(sunPosition);
	ProcessRenderedImage(m_renderTarget.GetColourBuffer());
}

// ..........................................................................

void VFCDirect::Render(const GENPoint &sunPosition)
{
	m_viewingSystem->SetViewDirection(-sunPosition);
	m_renderer->SetEyePosition(m_viewingSystem->EyePosition());
	m_renderer->Render(m_renderTarget);
}


// ..........................................................................

void VFCDirect::ProcessRenderedImage(const unsigned int* pixels)
{
	if (0!=m_pixelSum.size())
	{
	    //std::cerr << "m_pixelSum.size(): " << m_pixelSum.size() << std::endl;

		memset(&m_pixelSum[0],0,sizeof(int)*m_pixelSum.size());

		for (unsigned int j=0; j<m_renderWindowDim*m_renderWindowDim; ++j)
		{
		    //std::cerr << "pixel(" << j << "): surface id: " << pixels[j] << std::endl;
			if (pixels[j]>0)
			{
				// note return value is in range 1->n, we want 0->n-1 => subtract 1
				++m_pixelSum[pixels[j]-1];
			}
		}


		//for (unsigned int i=0; i<m_pixelSum.size(); i++) std::cerr << "m_pixelSum[" << i << "]: " << m_pixelSum[i] << std::endl;
	}
}

