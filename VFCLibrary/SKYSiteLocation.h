#ifndef _INC_SKYSITELOCATION_INCLUDED
#define _INC_SKYSITELOCATION_INCLUDED

#include "GENAngle.h"

class SKYSiteLocation
{
public:
	SKYSiteLocation() :
	  latitudeN(GENAngle::Degrees(52))
		{}

	SKYSiteLocation(const GENAngle &latN, const GENAngle &longE, const GENAngle &merE, const GENAngle &northW) :
		latitudeN(latN),
		longitudeE(longE),
		meridianE(merE),
		northAngleW(northW)
		{}

	GENAngle  latitudeN;
	GENAngle   longitudeE;
	GENAngle  meridianE;
	GENAngle  northAngleW; // west of north
};

#endif //_INC_SKYSITELOCATION_INCLUDED
