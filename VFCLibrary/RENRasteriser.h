#ifndef _INC_RENRASTERISER_INCLUDED
#define _INC_RENRASTERISER_INCLUDED

#include "RENVertex.h"

class RENRenderTarget;

/*
 * Generic interface for a rasteriser for use by RENSoftwareRenderer.  The
 * rasteriser could e.g. flat shade polygons, wireframe render them, etc.
 */
class RENRasteriser
{
public:
	/// Rasterise the triangle described by v1,v2,v3 (in screen co-ordinates)
	virtual void RenderTriangle(const RENVertex *v1, const RENVertex *v2, const RENVertex *v3,
						const unsigned int colour, RENRenderTarget &target) = 0;

	//void DrawLine(const RENVertex *v1, const RENVertex *v2,
	//				const unsigned int colour, RENRenderTarget &target);

protected:
	void DrawHorizontalLine(unsigned int y,
							unsigned int x1,
							unsigned int x2,
							RENRenderTarget &target, unsigned int colour);

	void DrawVerticalLine(unsigned int x,
							unsigned int y1,
							unsigned int y2,
							RENRenderTarget &target, unsigned int colour);

	void DrawNonAxisAlignedLine(const RENVertex *v1, const RENVertex *v2,
				   RENRenderTarget &target, unsigned int colour);

	void Bresenham(int x1, int y1, int x2, int y2,
					RENRenderTarget &target,
					unsigned int colour);
};

#endif //_INC_RENRASTERISER_INCLUDED
