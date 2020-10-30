#ifndef _INC_SKYTREGENZA_INCLUDED
#define _INC_SKYTREGENZA_INCLUDED

/*
 *   defines a sky discretised according to Tregenza's method (with an inverted sky dome for the ground)
 */

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <algorithm>

#include "./SKYVaultABC.h"
#include "GENFixedVector.h"

class SKYTregenza :
	public SKYVault
{
public:

	enum { PATCHCOUNT=290,BANDS=8,ALT=12 };

	SKYTregenza(void);
	virtual ~SKYTregenza(void);

	virtual void CreatePatches();
	virtual unsigned int PatchCount() const { return PATCHCOUNT; }

	// return a pointer to a particular sky patch
	virtual const SKYPatch* GetPatch(unsigned int i) const { return m_skyPatches[i]; }

    // return information about the present configuration
    unsigned int getBands() const { return BANDS; }
    double getDeltaAltitude() const { return ALT*M_PI/180.; }
    double getDeltaAzimuth(unsigned int band);
    unsigned int getPatchesPerBand(unsigned int band) const { return patchesPerBand[band%BANDS]; }
    unsigned int getPatchIndex(unsigned int band, int zone) const { unsigned int index=0; for (unsigned int i=0; i<band; ++i) index += patchesPerBand[i%BANDS]; index += zone; return index; }
    unsigned int getPatchBand(unsigned int patch) const { unsigned int band=0; while ( patch >= getPatchIndex(band+1,0) ) ++band; return band; }
    double getPatchCenterAltitude(unsigned int patch);
    double getPatchCenterAzimuth(unsigned int patch);


    // computes the spherical area of a regular patch and a triangle
    static double zoneSphericalArea(double phi1, double theta1, double phi2, double theta2);
    static double triangleSphericalArea(double phi1, double theta1, double phi2, double theta2);

private:
	// can only store pointers to sky patches since we want to store a mixture of "regular" and "top" patches
	GENFixedVector<SKYPatch*, PATCHCOUNT> m_skyPatches;
    static const unsigned int patchesPerBand[];

};

#endif //_INC_SKYTREGENZA_INCLUDED
