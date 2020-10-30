#ifndef _INC_RENRASTERISERFLATSHADING_INCLUDED
#define _INC_RENRASTERISERFLATSHADING_INCLUDED

#include "RENRasteriser.h"

/*
 * Rasterises a triangle using flat shading
 */
class RENRasteriserFlatShading : public RENRasteriser
{
public:
	RENRasteriserFlatShading(void);
	virtual ~RENRasteriserFlatShading(void);

	/// Rasterise the triangle described by v1,v2,v3 (in screen co-ordinates)
	virtual void RenderTriangle(const RENVertex *v1, const RENVertex *v2, const RENVertex *v3,
						const unsigned int colour, RENRenderTarget &target);

private:
	void ProcessTrapezoid(const float ystart,
						const float ystop,
						const float xl,
						const float zl,
						const float xr,
						const float zr,
						const float dxldy,
						const float dxrdy,
						const float dzldy,
						const float dzrdy,
						const unsigned int renderColour,
						RENRenderTarget &target);


	void RasteriseHalfTriangle(const unsigned int ystart,
					const unsigned int yend, // one past last row to be rendered
				   const float dxldy,
				   const float dxrdy,
				   const float dzldy,
				   const float dzrdy,
				   const unsigned int colour,
				   float xl,
				   float xr,
				   float zl,
				   float zr,
				   RENRenderTarget &target);
};

#endif //_INC_RENRASTERISERFLATSHADING_INCLUDED
