#ifndef _INC_VFCDIRECT_INCLUDED
#define _INC_VFCDIRECT_INCLUDED

#include <vector>
#include <list>

#include "GENPoint.h"
#include "RENRenderTarget.h"
#include "GEOMBoundingSphereCalc.h"

#ifdef OPENGL
#define RENRenderer RENOpenGLRenderer
#else
#define RENRenderer RENSoftwareRenderer
#endif

class RENRenderer;
class DATARadiationScene;
class VFCDirectViewSystem;

/*
 * This class calculates the view factors used by the solar radiation model.
 * The basic algorithm is as follows:
 * - For each sun position
 *    - Render the model from the point of view of the sun in orthographic view,
 *      with each surface uniquely coloured
 *    - Count how many pixels each surface covers
 *    - Convert number of pixels into an area
 *    - Divide by surface area to get fraction of surface that is insolated.
 */
class VFCDirect
{
	typedef RENRenderer RenderClass_t;
public:
    VFCDirect() : m_viewingSystem(nullptr) {}

	void Initialise(RenderClass_t *renderer,
					 DATARadiationScene &surfaces);


	void CalculateInsolationFactors(const GENPoint &sunPosition,
				const unsigned int sunPositionIndex);

private:
	void ZeroInsolationFactors(const unsigned int sunPositionIndex);

	void UpdatePixelSum(const GENPoint &sunPosition);
	void ProcessRenderedImage(const unsigned int* pixels);
	void Render(const GENPoint &sunPosition);

	RENRenderTarget m_renderTarget;
	RenderClass_t *m_renderer;
    std::unique_ptr<VFCDirectViewSystem> m_viewingSystem;
	unsigned int m_renderWindowDim;
	float m_pixelArea;
	std::vector<unsigned int> m_pixelSum;

	DATARadiationScene *m_scene;


};

// ..........................................................................

class VFCDirectViewSystem
{
public:
	VFCDirectViewSystem(const GENPoint &modelCentre, float modelRadius) :
		m_modelCentre(modelCentre),
        m_viewDirection(GENPoint::Origin()),
		m_modelRadius(modelRadius)
	{}

	void SetViewDirection(const GENPoint &direction)
	{
		m_viewDirection=direction;
	}

	GENPoint EyePosition() const
	{
		return m_modelCentre-m_viewDirection*(m_modelRadius+1);
	}

	GENPoint ViewTarget() const
	{
		return m_modelCentre;
	}

private:
	GENPoint m_modelCentre;
	GENPoint m_viewDirection;
	float m_modelRadius;
};

#endif //_INC_VFCDIRECT_INCLUDED
