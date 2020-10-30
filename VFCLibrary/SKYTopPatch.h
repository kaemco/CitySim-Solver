#ifndef _INC_SKYTOPPATCH_INCLUDED
#define _INC_SKYTOPPATCH_INCLUDED

/*
 *  Implements a round sky patch (for the top of a dome) - a circular region of
 *  sky bounded by a constant azimuth band
 *
 */

#include "SKYPatch.h"
#include "GENFixedVector.h"

class SKYTopPatch :
	public SKYPatch
{
public:
	SKYTopPatch(float   AltRange=0, bool IsTop=true);
	virtual ~SKYTopPatch(void);

	virtual inline const float& Alt() const;

protected:
	// geometry variables
	float  m_Alt;

	// spacing of cells
	static const int m_altBands=10;
	static const int m_azBands=60;
	static const int m_numCells=m_altBands*m_azBands;

	// IsTop = true for top patch, false for bottom patch
	bool m_IsTop;
};

#endif // _INC_SKYTOPPATCH_INCLUDED
