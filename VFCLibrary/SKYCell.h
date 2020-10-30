#ifndef _INC_SKYCELL_INCLUDED
#define _INC_SKYCELL_INCLUDED

/**
 * Implements a 'cell' - an elemental region of a sky patch.  A cell has a
 * direction and a solid angle. A collection of these can be used to describe a
 * sky patch.
 *
 */

#include "GENPoint.h"

class SKYCell
{
public:
	SKYCell() {}
	SKYCell(const GENPoint &direction, float solidAngle);

	const GENPoint& direction() const { return m_direction; }
	float solidAngle() const { return m_solidAngle; }

	float formFactor(const GENPoint& viewDir) const;

private:
	GENPoint m_direction;
	float m_solidAngle;
};

#endif //_INC_SKYCELL_INCLUDED
