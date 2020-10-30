#include "SKYTregenza.h"

#include "SKYRegularPatch.h"
#include "SKYTopPatch.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER


using std::vector;

const unsigned int SKYTregenza::patchesPerBand[8]={30,30,24,24,18,12,6,1};

SKYTregenza::SKYTregenza(void)
{
	CreatePatches();
}

SKYTregenza::~SKYTregenza(void)
{
	for (int i=0; i<(int)m_skyPatches.size(); ++i)
	{
		delete m_skyPatches[i];
	}
}

void SKYTregenza::CreatePatches()
{
	int i,j;

	const int AzIncDegrees[7]={12,12,15,15,20,30,60};

	// generate the sky patches
	int patchCounter=0;

	const float AltRange=12.f*(float)M_PI/180.f;

	float Alt=6.f*(float)M_PI/180.f;
	for (i=0; i<7; i++)
	{
		int RowCount=360/AzIncDegrees[i];

		const float AzRange=2.f*(float)M_PI/RowCount;
		float Az=0;
		for (j=0; j<RowCount; j++)
		{
			// add a new sky patch to the vector
			m_skyPatches[patchCounter] = new SKYRegularPatch(Alt,Az,AltRange, AzRange);
			++patchCounter;

			Az+=AzRange;
		}
		Alt+=12.f*(float)M_PI/180.f;
	}

	// add the top patch
	m_skyPatches[patchCounter] = new SKYTopPatch(6.f*(float)M_PI/180.f,true);
	++patchCounter;

	// now the inverted sky
	Alt=6.f*(float)M_PI/180.f;
	for (i=0; i<7; i++)
	{
		int RowCount=360/AzIncDegrees[i];

		const float AzRange=2.f*(float)M_PI/RowCount;
		float Az=0;
		for (j=0; j<RowCount; j++)
		{
			m_skyPatches[patchCounter] = new SKYRegularPatch(-Alt,Az,AltRange, AzRange);
			++patchCounter;
			Az+=AzRange;
		}
		Alt+=12.f*(float)M_PI/180.f;
	}

	// add the bottom patch
	m_skyPatches[patchCounter] = new SKYTopPatch(6.f*(float)M_PI/180.f,false);
	++patchCounter;

}

double SKYTregenza::getDeltaAzimuth(unsigned int band) { return 2.*M_PI/patchesPerBand[band]; }

double SKYTregenza::getPatchCenterAltitude(unsigned int patch) {

    if ( patch == 144 ) return M_PI/2.;
    else return asin( (sin( ALT*M_PI/180.*getPatchBand(patch) ) + sin( ALT*M_PI/180.*(getPatchBand(patch)+1) ) ) / 2. );

}

double SKYTregenza::getPatchCenterAzimuth(unsigned int patch) {

    return fmod(m_skyPatches[patch]->centroid().Azimuth().radians()+2.*M_PI,2.*M_PI);

}

double SKYTregenza::zoneSphericalArea(double phi1, double theta1, double phi2, double theta2) {

    return (sin(theta2)-sin(theta1))*(phi2-phi1);

}

double SKYTregenza::triangleSphericalArea(double phi1, double theta1, double phi2, double theta2) {

    if (theta2 < theta1) {
        return -((phi1-phi2)*(-cos(theta1)+cos(theta2)+(theta2-theta1)*sin(theta2)))/(theta1-theta2);
    }
    else if (theta2 > theta1) {
        return ((phi1-phi2)*(-cos(theta1)+cos(theta2)+(theta2-theta1)*sin(theta1)))/(theta2-theta1);
    }
    else return 0.;

}
