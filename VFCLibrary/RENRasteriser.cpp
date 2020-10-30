#include "./RENRasteriser.h"

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


#include "RENRenderTarget.h"

// ..........................................................................

/*void RENRasteriser::DrawLine(const RENVertex *v1, const RENVertex *v2, unsigned int colour, RENRenderTarget &target)
{
	const unsigned int y1 = static_cast<unsigned int>(v1->m_position[1]);
    const unsigned int y2 = static_cast<unsigned int>(v2->m_position[1]);

	const unsigned int x1 = static_cast<unsigned int>(v1->m_position[0]);
    const unsigned int x2 = static_cast<unsigned int>(v2->m_position[0]);

	if (y1==y2)
		DrawHorizontalLine(y1,x1,x2, target,colour);
	else if (x1==x2)
		DrawVerticalLine(x1,y1,y2, target,colour);
	else
		DrawNonAxisAlignedLine(v1,v2,target,colour);
}*/

// ..........................................................................

void RENRasteriser::DrawHorizontalLine(unsigned int y, 
											unsigned int x1, 
											unsigned int x2,
											RENRenderTarget &target,
											unsigned int colour)
{
	if (x1>x2)
		std::swap(x1,x2);

	for (unsigned int x=x1; x<x2; ++x)
	{
		target.SetPixel(x,y,colour);
	}
}

// ..........................................................................

void RENRasteriser::DrawVerticalLine(unsigned int x, 
										unsigned int y1, 
										unsigned int y2,
										RENRenderTarget &target,
										unsigned int colour)
{
	if (y1>y2)
		std::swap(y1,y2);

	for (unsigned int y=y1; y<y2; ++y)
	{
		target.SetPixel(x,y,colour);
	}
}

// ..........................................................................

void RENRasteriser::DrawNonAxisAlignedLine(const RENVertex *v1, 
												const RENVertex *v2,
												RENRenderTarget &target,
												unsigned int colour)
{
	const unsigned int y1 = static_cast<unsigned int>(v1->m_position[1]);
    const unsigned int y2 = static_cast<unsigned int>(v2->m_position[1]);

	const unsigned int x1 = static_cast<unsigned int>(v1->m_position[0]);
    const unsigned int x2 = static_cast<unsigned int>(v2->m_position[0]);
	Bresenham(x1,y1,x2,y2,target,colour);
}

// ..........................................................................

/* Implementation of Bresenham line drawing algorithm from Graphic Gems
 */
void RENRasteriser::Bresenham(int x1, int y1, int x2, int y2, 
									RENRenderTarget &target,
								   unsigned int colour)
{
	const int dx=x2-x1;
	const int ax=abs(dx)*2;
	const int sx= (dx>=0) ? 1 : -1;

	const int dy=y2-y1;
	const int ay=abs(dy)*2;
	const int sy= (dy>=0) ? 1 : -1;

	int x=x1;
	int y=y1;

	if (ax>ay)
	{
		int d=ay-ax/2;

        for (;;)
		{
			target.SetPixel(x,y,colour);
			if (x==x2) return;
			if (d>=0) {
				y+=sy;
				d-=ax;
			}
			x+=sx;
			d+=ay;
		}
	}
	else
	{
		int d=ax-ay/2;

		for (;;)
		{
			target.SetPixel(x,y,colour);
			if (y==y2) return;
			if (d>=0) {
				x+=sx;
				d-=ay;
			}
			y+=sy;
			d+=ax;
		}
	}
}


