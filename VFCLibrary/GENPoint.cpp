#include "GENPoint.h"


#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

using namespace GENPointCoords;

// ..........................................................................
bool operator<(const GENPoint &a, const GENPoint &b)  
{
	const static float tolerance=1e-5f;

	// first test x
	if (a[X]<b[X]-tolerance)
		return true;
	if (a[X]>b[X]+tolerance)
		return false;

	// then test y
	if (a[Y]<b[Y]-tolerance)
		return true;
	if (a[Y]>b[Y]+tolerance)
		return false;

	// finally test z
	if (a[Z]<b[Z]-tolerance)
		return true;
	if (a[Z]>b[Z]+tolerance)
		return false;

	return false;
}

// ..........................................................................


/**
 * Given the line a,b there is a point x0 that is the closest point on the line
 * to [point].  This function determines where that point is, returning lamda
 * such that:
 *
 * x0 = a + (b-a) * lamda
 *
 * see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
 */
float closest_point_on_line(const GENPoint &point, const GENPoint &a, const GENPoint &b)
{
	return -dot_product(a-point,b-a)/length_squared(b-a);
}

// ..........................................................................

// Tests for collinearity of a, b and c
bool collinear(const GENPoint &a,const GENPoint &b,const GENPoint &c)
{
	const static float tolerance=1e-6f;

	if (ApproxEqual(a,b) || ApproxEqual(b,c) || ApproxEqual(a,c))
		return true;

	const float l1=length_squared(a-b);
	const float l2=length_squared(a-c);
	const float l3=length_squared(b-c);

	// choose the two points that are furthest away from each other to define a line
	float distanceSqFromLine;
	float lineLengthSq;
	if (l1>l2 && l1>l3)
	{
		lineLengthSq=l1;
		distanceSqFromLine=distance_from_line_squared(c,a,b-a);
	}
	else if (l2>l3)
	{
		lineLengthSq=l2;
		distanceSqFromLine=distance_from_line_squared(b,a,c-a);
	}
	else
	{
		lineLengthSq=l3;
		distanceSqFromLine=distance_from_line_squared(a,b,c-b);
	}

	if (distanceSqFromLine<tolerance*lineLengthSq)
	{
		return true;
	}
	else 
	{
		return false;
	}
}

