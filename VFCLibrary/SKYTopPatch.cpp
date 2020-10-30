#include "SKYTopPatch.h"
#include "GENAssert.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER


SKYTopPatch::SKYTopPatch(float AltRange, bool IsTop) :
	SKYPatch(GENPoint::Cartesian(0,0,IsTop ? 1.f : -1.f)),
    m_Alt(IsTop ? static_cast<float>(M_PI/2.) : -static_cast<float>(M_PI/2.)),
    m_IsTop(IsTop)
{
	// create the sub cells
	float  ThisAlt, ThisAz;

	int up;
	if (IsTop)
		up=1;
	else
		up=-1;

	const float AltInc=up*AltRange/static_cast<float>(m_altBands);
	const float AzInc=2.f*static_cast<float>(M_PI)/static_cast<float>(m_azBands);

	int cellCounter=0;
    for (ThisAlt=up*(static_cast<float>(M_PI/2.)-AltRange)+AltInc/2;
         up*ThisAlt < M_PI/2.;
		 ThisAlt+=AltInc)
	{
		for (ThisAz=0+AzInc/2; ThisAz < 2*M_PI; ThisAz+=AzInc)
		{
			float dw=fabs(AzInc * (sin(ThisAlt+AltInc/2.f)-sin(ThisAlt-AltInc/2.f)));

			// add a new cell to this patch's cells list
			addCell(SKYCell(GENPoint::Polar(GENAngle::Radians(ThisAlt),GENAngle::Radians(ThisAz)),dw));
			++cellCounter;
		}
	}
	GENAssert(cellCounter==600);

}

SKYTopPatch::~SKYTopPatch(void)
{
}

inline const float& SKYTopPatch::Alt() const
{
	return m_Alt;
}

