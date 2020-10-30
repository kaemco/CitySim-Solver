#ifndef _INC_GEOMSPHERICALBASIS_INCLUDED
#define _INC_GEOMSPHERICALBASIS_INCLUDED

#include <vector>
#include <set>
#include <list>

#include "GENPoint.h"

/* Part of the bounding sphere calculation
 */
class GEOMSphericalBasis
{
public:
	typedef std::vector<GENPoint> Vector_t; 
	typedef std::list<GENPoint> Basis_t; 

	GEOMSphericalBasis()  { Update();}

	size_t size() const;

	bool push_back(const GENPoint& p);
	void pop_back();

	GENPoint Centre() const;
	float Radius() const;

private:
	void AddToBasis(const GENPoint &p);
	void RemoveFromBasis(const GENPoint &p);

	bool Update();

	bool Update0(const Basis_t &basis);
	bool Update1(const Basis_t &basis);
	bool Update2(const Basis_t &basis);
	bool Update3(const Basis_t &basis);
	bool Update4(const Basis_t &basis);
	
	GENPoint m_centre;
	float m_radius;

	Vector_t m_points;	
	Basis_t m_basis;
};

#endif // _INC_GEOMSPHERICALBASIS_INCLUDED
