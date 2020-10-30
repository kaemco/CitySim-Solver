#include "./RENRasteriserFlatShading.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "GENAssert.h"
#include "RENRenderTarget.h"

RENRasteriserFlatShading::RENRasteriserFlatShading(void)
{
}

// ..........................................................................

RENRasteriserFlatShading::~RENRasteriserFlatShading(void)
{
}
	
// ..........................................................................

void ScanLine(const int y,
			  const unsigned int stride,
			  const unsigned int colour,
			  const int xmin, 
			  const int xmax, 
			  const float xl, 
			  const float xr, 
			  const float zl, 
			  const float zr,RENRenderTarget &target,
			  float* depthBuffer) 
{
	const float dzdx=(zr-zl)/(xr-xl);
	const float deltax=(static_cast<float>(xmin)-xl);

	float z=zl+deltax*dzdx;

	for (int x=xmin; x<=xmax; ++x)
	{
		if (z<depthBuffer[x+y*stride])
		{
			depthBuffer[x+y*stride] = z;
			target.SetPixel(x,y,colour);
		}
		z+=dzdx;
	}
}

// ..........................................................................

/*
 * Previously this method worked directly on two arrays (one for depth buffer,
 * one for colour buffer) as I had assumed that this would be quicker than
 * calling a .setPixel method on the render target.  As it turns out, that is
 * not the case, and the .setPixel method is just as quick.  The .setPixel
 * approach is preferrable as
 * - the render target can implement something like "dirty rectangles" to keep
 *   track of which bits have been drawn on.
 * - the render target does not have to use contiguous arrays for the depth and
 *   colour buffers if it doesn't want to.
 */
void RENRasteriserFlatShading::RasteriseHalfTriangle(const unsigned int ystart,
					const unsigned int yend, 
				   const float dxldy,
				   const float dxrdy,
				   const float dzldy,
				   const float dzrdy,
				   const unsigned int colour,
				   float xl,
				   float xr,
				   float zl,
				   float zr,
				   RENRenderTarget &target)
{
	float* depthBuffer = target.GetDepthBuffer();
	const int stride=target.GetStride();

	for(unsigned int y = ystart; y<yend; ++y)
    {
		const int xmin=static_cast<int>(xl+1);
		const int xmax=static_cast<int>(xr);

		GENAssert(xmin>=0);
		if (xmax>=xmin)
		{
			ScanLine(y,stride,colour,xmin,xmax,xl,xr,zl,zr,target,depthBuffer);
		}

		xl+=dxldy;
		xr+=dxrdy;

		zl+=dzldy;
		zr+=dzrdy;
	}
}

float gradient(const float x1, const float y1, const float x2, const float y2)
{
	return (x2-x1)/(y2-y1);
}

void RENRasteriserFlatShading::ProcessTrapezoid(const float ystart,
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
													RENRenderTarget &target)
{
	const int yfirst=static_cast<int>(ceil(ystart));
	const int ylast=static_cast<int>(ceil(ystop));

	// calculate xl and xr - the x-coordinates of the left & right edges at the
	// first integer value of y after y1
	const float deltay=(static_cast<float>(yfirst)-ystart);
	const float xlfirst=xl+dxldy*deltay;
	const float xrfirst=xr+dxrdy*deltay;

	const float zlfirst=zl+dzldy*deltay;
	const float zrfirst=zr+dzrdy*deltay;

	RasteriseHalfTriangle(yfirst,
							ylast,
							dxldy,
							dxrdy,
							dzldy,
							dzrdy,
							renderColour,
							xlfirst,
							xrfirst,
							zlfirst,
							zrfirst,
							target);
}

// ..........................................................................

// Expects coordinates to be in pixel coordinate system (i.e. whole numbers
// correspond to a centre of a pixel)
void RENRasteriserFlatShading::RenderTriangle(const RENVertex *v1, 
										const RENVertex *v2, 
										const RENVertex *v3,
										const unsigned int renderColour, 
										RENRenderTarget &target)
{
	// sort vertices bottom to top
	if(v1->m_position[1] > v3->m_position[1]) 
		std::swap(v1, v3);
	if(v2->m_position[1] > v3->m_position[1]) 
		std::swap(v2, v3);
	if(v1->m_position[1] > v2->m_position[1]) 
		std::swap(v1, v2);

	// coordinates of points
    const float y1 = v1->m_position[1];
    const float y2 = v2->m_position[1];
    const float y3 = v3->m_position[1];

    const float x1 = v1->m_position[0];
    const float x2 = v2->m_position[0];
    const float x3 = v3->m_position[0];

    const float z1 = v1->m_position[2];
    const float z2 = v2->m_position[2];
    const float z3 = v3->m_position[2];

	if (y3==y1) // triangle is a horizontal line
		return;

	// gradients on the long edge
	const float dx13dy=gradient(x1,y1,x3,y3);
	const float dz13dy=gradient(z1,y1,z3,y3);

	// gradients on the first short edge
	const float dx12dy=gradient(x1,y1,x2,y2);
	const float dz12dy=gradient(z1,y1,z2,y2);

	// gradients on the second short edge
	const float dx23dy=gradient(x2,y2,x3,y3);
	const float dz23dy=gradient(z2,y2,z3,y3);

	const float x2b=x1+dx13dy*(y2-y1);
	const float z2b=z1+dz13dy*(y2-y1);

	// figure out which way around the edges are (is 1-2 or 1-3 furthest left?)
	if (dx12dy>dx13dy)
	{
		// 1-3 is left edge
		ProcessTrapezoid(y1,y2,x1,z1,x1,z1,dx13dy,dx12dy,dz13dy,dz12dy,renderColour,target);
		ProcessTrapezoid(y2,y3,x2b,z2b,x2,z2,dx13dy,dx23dy,dz13dy,dz23dy,renderColour,target);
	}
	else
	{
		// 1-3 is right edge
		ProcessTrapezoid(y1,y2,x1,z1,x1,z1,dx12dy,dx13dy,dz12dy,dz13dy,renderColour,target);
		ProcessTrapezoid(y2,y3,x2,z2,x2b,z2b,dx23dy,dx13dy,dz23dy,dz13dy,renderColour,target);
	}
}



