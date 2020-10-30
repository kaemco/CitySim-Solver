#ifndef _INC_SKYREGULARPATCH_INCLUDED
#define _INC_SKYREGULARPATCH_INCLUDED


/*
 * Implements a sky patch with constant altitude and constant azimuth sides - a
 * patch of sky delimited by bands of constant altitude and azimuth.
 *
 */

#include "SKYPatch.h" // this classes base class
#include "GENFixedVector.h"

class SKYRegularPatch :
	public SKYPatch
{
public:
	SKYRegularPatch (float Alt=0, float Az=0, float AltRange=0, float AzRange=0);
	virtual ~SKYRegularPatch (void);

protected:
	// spacing of cells
	static const int m_altBands=10;
	static const int m_azBands=10;
	static const int m_numCells=m_altBands*m_azBands;
};

#endif // _INC_SKYREGULARPATCH_INCLUDED
