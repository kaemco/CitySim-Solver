#ifndef _INC_VFCVIEWFACTORCALCULATION_INCLUDED
#define _INC_VFCVIEWFACTORCALCULATION_INCLUDED

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <memory>

#include "DATARadiationScene.h"
#include "VFCDiffuse.h"
#include "RENRenderTarget.h"
class SKYVault;
class RENRenderer;
class DATARadiationScene;

/*
 * This class is the main class for the view factor calculation, but most of the
 * work is done elsewhere
 */
class VFCViewFactorCalculation
{
public:
	VFCViewFactorCalculation(void);
	~VFCViewFactorCalculation(void);

	void CalculateViewFactors(DATARadiationScene &surfaces);

    // DP: added to decompose pre-process simulation steps in Pro
    std::shared_ptr<RENRenderer> getM_renderer(){ return m_renderer; }

    // DP: these used to be private, but are needed in Pro
    void CalculateInsolationFactors(DATARadiationScene &surfaces);

    template <class GEOMETRY, class SURFACE_ITERATOR, class RESULT_OBJECT>
    void CalculateDiffuseFactors(SURFACE_ITERATOR surfaces,
                                 size_t surfaceCount,
                                 RESULT_OBJECT&(DATASurface::*resultObjectAccessor)(unsigned int))
    {
        VFCDiffuse<GEOMETRY> calc(m_theSky,surfaceCount,m_renderer);
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

private:

    std::shared_ptr<SKYVault> m_theSky;
    std::shared_ptr<RENRenderer> m_renderer;
};


#endif //_INC_VFCVIEWFACTORCALCULATION_INCLUDED
