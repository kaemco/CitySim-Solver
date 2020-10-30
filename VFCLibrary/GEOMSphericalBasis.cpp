#include "GEOMSphericalBasis.h"
#include "GENAssert.h"
#include "../util.h"
#include <algorithm>
#include <cassert>

// ..........................................................................

size_t GEOMSphericalBasis::size() const {
	return m_basis.size();
}

// ..........................................................................

bool GEOMSphericalBasis::push_back(const GENPoint& p) {
	m_points.push_back(p);
	AddToBasis(p);

	if (!Update())
	{
		this->pop_back();
		return false;
	}
	return true;
}

// ..........................................................................

void GEOMSphericalBasis::pop_back() {
	RemoveFromBasis(m_points.back());
	m_points.pop_back();

	Update();
}

// ..........................................................................

void GEOMSphericalBasis::AddToBasis(const GENPoint &p)
{
	if (m_basis.end()==std::find(m_basis.begin(),m_basis.end(),p))
		m_basis.push_back(p);
}

// ..........................................................................

void GEOMSphericalBasis::RemoveFromBasis(const GENPoint &p)
{
	if (!(m_basis.end()==std::find(m_basis.begin(),m_basis.end(),p)))
	{
		m_basis.remove(p);
	}
}

// ..........................................................................

GENPoint GEOMSphericalBasis::Centre() const
{
	return m_centre;
}

// ..........................................................................

float GEOMSphericalBasis::Radius() const
{
	return m_radius;
}

// ..........................................................................

bool GEOMSphericalBasis::Update()
{
	bool result=true;
	switch(m_basis.size())
	{
	case 0:
		result=Update0(m_basis);
		break;
	case 1:
		result=Update1(m_basis);
		break;
	case 2:
		result=Update2(m_basis);
		break;
	case 3:
		result=Update3(m_basis);
		break;
	case 4:
		result=Update4(m_basis);
		break;
	default:
		throw std::logic_error("GEOMSphericalBasis:: Invalid size for basis function!");
		break;
	}
	return result;
}

// ..........................................................................

bool GEOMSphericalBasis::Update0(const Basis_t &)
{
	m_centre=GENPoint::Origin();
	m_radius=0;
	return true;
}

// ..........................................................................

bool GEOMSphericalBasis::Update1(const Basis_t &basis)
{
	m_centre=basis.front();
	m_radius=0;
	return true;
}

// ..........................................................................

bool GEOMSphericalBasis::Update2(const Basis_t &basis)
{
	GENPoint diff=basis.back()-basis.front();

	m_centre=basis.front()+diff*0.5f;
	m_radius=diff.Radius()/2.f;
	return true;
}

// ..........................................................................

bool GEOMSphericalBasis::Update3(const Basis_t &basis)
{
	Basis_t::const_iterator it=basis.begin();
    const GEN::Point<double> basis0=*it;
	++it;
    const GEN::Point<double> basis1=*it;
	++it;
    const GEN::Point<double> basis2=*it;

    GEN::Point<double>  a = basis1-basis0;
    GEN::Point<double>  b = basis2-basis0;

    GEN::Point<double> crossProdab=cross_product(a,b);
    double denominator = 2. * (dot_product(crossProdab,crossProdab));

	// denominator==0 implies that all three points are collinear, which should
    // never happen during the bounding sphere calculation.
    if (denominator==0.) throw std::string("Three collinear points in bounding sphere calculation: (" + toString(basis0) + ") and (" + toString(basis1) + ") and (" + toString(basis2) + ")");

    GEN::Point<double> o = (cross_product(crossProdab,a)*dot_product(b,b) +
                            cross_product(b,crossProdab)*dot_product(a,a)) / denominator;

	m_radius=o.Radius();
	m_centre=basis0+o;
	return true;
}

// ..........................................................................

// matrix determinant
float Determinant3x3(float m11, float m12, float m13,
		       float m21, float m22, float m23,
		       float m31, float m32, float m33)
{
	return m11 * (m22 * m33 - m32 * m23) -
	       m21 * (m12 * m33 - m32 * m13) +
	       m31 * (m12 * m23 - m22 * m13);
}

// ..........................................................................

bool GEOMSphericalBasis::Update4(const Basis_t &basis)
{
	Basis_t::const_iterator it=basis.begin();
	const GENPoint basis0=*it;
	++it;
	const GENPoint basis1=*it;
	++it;
	const GENPoint basis2=*it;
	++it;
	const GENPoint basis3=*it;

	GENPoint a = basis1-basis0;
	GENPoint b = basis2-basis0;
	GENPoint c = basis3-basis0;

	float denominator = 2.0f * Determinant3x3(a[GENPointCoords::X], a[GENPointCoords::Y], a[GENPointCoords::Z],
											b[GENPointCoords::X], b[GENPointCoords::Y], b[GENPointCoords::Z],
		                                    c[GENPointCoords::X], c[GENPointCoords::Y], c[GENPointCoords::Z]);

	if (fabs(denominator)<1e-4)
	{
		// denominator=0 implies all four points are on the sphere.  We know
        // that each point is unique, so can use just first 3 to establish the
        // sphere.
		m_basis.pop_back();
		return Update3(basis);
	}

	GENPoint tmp = (cross_product(a,b)*dot_product(c,c) +
						cross_product(c,a)*dot_product(b,b) +
		                cross_product(b,c)*dot_product(a,a));
	GENPoint o=tmp/ denominator;

	m_radius = o.Radius();
	m_centre = basis0 + o;
	return true;
}

