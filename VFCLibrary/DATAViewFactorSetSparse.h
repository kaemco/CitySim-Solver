#ifndef _INC_DATAVIEWFACTORSETDAYLIGHT_INCLUDED
#define _INC_DATAVIEWFACTORSETDAYLIGHT_INCLUDED

#include <vector>
#include <algorithm>

struct SDATAViewFactorTuple
{
	SDATAViewFactorTuple(int patchNo_, float unobstructed_=0, float obstructed_=0, unsigned int mainobstructing_=0) :
		patchNo(patchNo_),
		unobstructed(unobstructed_),
		obstructed(obstructed_),
		mainobstructing(mainobstructing_) {}
	int patchNo;
	float unobstructed;
	float obstructed;
	unsigned int mainobstructing;
};

/*
 * Stores a set of view factors for the diffuse calc in a sparse format.
 *
 * NOTE: storage in a sparse format is not really worth it for diffuse solar sky patch v.f.s
 * since approximately 40% of matrix is non-zero, but it is worth it for daylight externally
 * refelected component v.f.s
 *
 */
class DATAViewFactorSetSparse
{
public:
	typedef  std::vector<SDATAViewFactorTuple>::const_iterator const_iterator;

	void SetViewFactors(
				unsigned int patchNum,
				float obstructedVF,
				float unobstructedVF,
				unsigned int mainObstructingSurface);

	inline const_iterator GetVFs() const
	{ return m_storedValues.begin(); }

	inline const_iterator GetLastVF() const
	{ return m_storedValues.end(); }

private:
	std::vector<SDATAViewFactorTuple> m_storedValues;
};

#endif //_INC_DATAVIEWFACTORSETDAYLIGHT_INCLUDED

