#include <algorithm>

#include "./SKYCell.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

SKYCell::SKYCell(const GENPoint &direction, float solidAngle) :
	m_direction(direction),
	m_solidAngle(solidAngle)
{
}

// ..........................................................................

float SKYCell::formFactor(const GENPoint& viewDir) const
{
	return (std::max)(0.f,dot_product(viewDir, m_direction)*m_solidAngle);
}
