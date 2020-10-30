#include "VFCViewFactorCalculation.h"

#include "SKYTregenza.h"

#ifdef OPENGL
#include "RENOpenGLRenderer.h"
#define RENRenderer RENOpenGLRenderer
#else
#include "RENSoftwareRenderer.h"
#define RENRenderer RENSoftwareRenderer
#endif


#include "VFCDirect.h"

#include "VFCSWGeometry.h"
#include "VFCVISGeometry.h"
#include "DATARadiationScene.h"

#include "DATASurface.h"
#include "DATASurfaceIterator.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER


using std::list;
using std::vector;

// ..........................................................................

VFCViewFactorCalculation::VFCViewFactorCalculation(void) :
    m_theSky(new SKYTregenza),
	m_renderer(new RENRenderer)
{
}

// ..........................................................................

VFCViewFactorCalculation::~VFCViewFactorCalculation(void)
{
}

// ..........................................................................

void VFCViewFactorCalculation::CalculateViewFactors(DATARadiationScene &scene)
{
    //std::cout << "VFCViewFactorCalculation" << std::endl;
	m_renderer->SetModel(scene.IndexedFaceSet());

	CalculateInsolationFactors(scene);
    //std::cout << "CalculateDiffuseFactors<VFCSWGeometry>" << std::endl;
    CalculateDiffuseFactors<VFCSWGeometry>(scene.GetAllSurfaces(),scene.SurfaceCount(),&DATASurface::SWViewFactors);
    //std::cout << "CalculateDiffuseFactors<VFCVISGeometry>" << std::endl;
    CalculateDiffuseFactors<VFCVISGeometry,DATASurfaceBuildingIterator,DATAViewFactorSetSparse>(scene.GetBuildingSurfaces(),scene.SurfaceCount(),&DATASurface::DaylightViewFactors);
}

// ..........................................................................

void VFCViewFactorCalculation::CalculateInsolationFactors(DATARadiationScene &scene)
{
	if (scene.GetAllSurfaces().end()==scene.GetAllSurfaces())
		return; // can't get the relevant sun positions without a
                // DATAInsolationFactors object, can't get a
                // DATAInsolationFactors object without a surface...

	VFCDirect solarCalc;
	solarCalc.Initialise(m_renderer.get(), scene);

	vector<GENPoint> sunPositions;
//	DATAInsolationFactors::GetSunPositionsToRender(scene.GetLocation(),sunPositions);
	(*scene.GetAllSurfaces()).InsolationFactors().GetSunPositionsToRender(scene.GetLocation(),sunPositions);

	// for each sun position
	for (unsigned int i=0; i<sunPositions.size(); ++i)
	{
		solarCalc.CalculateInsolationFactors(sunPositions[i],i);
	}
}

/*
template <class GEOMETRY, class SURFACE_ITERATOR, class RESULT_OBJECT>
void VFCViewFactorCalculation::CalculateDiffuseFactors(SURFACE_ITERATOR surfaces,
                             size_t surfaceCount,
                             RESULT_OBJECT&(DATASurface::*resultObjectAccessor)(unsigned int))
{
    VFCDiffuse<GEOMETRY> calc(*m_theSky,surfaceCount,m_renderer.get());
    unsigned int tempCounter=0;
    for (;
        !surfaces.isAtEnd();
        ++surfaces, ++tempCounter)
    {
        std::vector<RESULT_OBJECT*> resultObjects;
        for (unsigned int i=0; i<calc.m_geometry.ViewPointCount(); ++i)
            resultObjects.push_back(&((*surfaces).*resultObjectAccessor)(i));

        calc.Calculate((*surfaces).Centroid(),
                            (*surfaces).Normal(),
                            resultObjects);
    }
}

*/
