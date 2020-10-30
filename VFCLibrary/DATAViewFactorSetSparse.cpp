#include "./DATAViewFactorSetSparse.h"
#include "GENAssert.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath> // modified by JK - math.c
#include <string> // added by JK - string usage in error handling

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

using std::find;

// ..........................................................................

bool operator<(const SDATAViewFactorTuple &a, const SDATAViewFactorTuple &b)
{
	return (a.patchNo<b.patchNo);
}

bool operator==(const SDATAViewFactorTuple &a, const SDATAViewFactorTuple &b)
{
	return (a.patchNo==b.patchNo);
}

bool operator==(const SDATAViewFactorTuple *a, const SDATAViewFactorTuple &b)
{
	return (a->patchNo==b.patchNo);
}

/**
 *
 * \param patchNum - the number of the patch for which view factors are being stored
 * \param obstructedVF - the view factor to the obstructed region of sky
 * \param unobstructedVF - the view factor to the unobstructed region of sky
 * \param mainObstructingSurface - the main obstructing surface of the patch
 */
void DATAViewFactorSetSparse::SetViewFactors(
				unsigned int patchNum,
				float obstructedVF,
				float unobstructedVF,
				unsigned int mainObstructingSurface)
{
	if (unobstructedVF>0 || obstructedVF>0)
	{
		GENEnforce(m_storedValues.end()==find(m_storedValues.begin(),m_storedValues.end(),SDATAViewFactorTuple(patchNum)),"View factors have already been set for this patch");

		m_storedValues.push_back(
			SDATAViewFactorTuple(patchNum,unobstructedVF, obstructedVF, mainObstructingSurface));

		std::sort(m_storedValues.begin(), m_storedValues.end());
	}
}


// ..........................................................................
