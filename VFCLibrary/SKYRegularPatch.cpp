#include "./SKYRegularPatch.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

SKYRegularPatch ::SKYRegularPatch (float  alt, float az, float altRange, float azRange) :
	SKYPatch(GENPoint::Polar(GENAngle::Radians(alt),GENAngle::Radians(az)))
{
	const float altInc=altRange/(float)m_altBands; 
	const float azInc=azRange/(float)m_azBands;

	// Create the cells that make up this patch
	int cellID=0;
	for (float thisAlt=alt-altRange/2.f+altInc/2.f; thisAlt < alt+altRange/2.f; thisAlt+=altInc)
	{
		for (float thisAz=az-azRange/2.f+azInc/2.f; thisAz < az+azRange/2.f; thisAz+=azInc, ++cellID)
		{
			float dw=azInc * (sin(thisAlt+altInc/2.f)-sin(thisAlt-altInc/2.f));

			// add a new cell to this patch's cells list
			addCell(SKYCell(GENPoint::Polar(GENAngle::Radians(thisAlt),GENAngle::Radians(thisAz)),dw));
		}
	}
}

// ..........................................................................

SKYRegularPatch ::~SKYRegularPatch (void)
{
}
