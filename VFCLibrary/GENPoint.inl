#include "GENPoint.h"


#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

// ..........................................................................

template<class NUMBER>
Point<NUMBER> Point<NUMBER>::Cartesian(NUMBER x,NUMBER y,NUMBER z)
{
	return Point<NUMBER>(x,y,z);
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER> Point<NUMBER>::Polar(const Angle<NUMBER> &altitude, const Angle<NUMBER> &azimuth, NUMBER radius)
{
	return Point<NUMBER>(altitude,azimuth,radius);
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER> Point<NUMBER>::UnitVector(const Point<NUMBER> &A)
{
	if (A.Radius()!=0)
		return Cartesian(A.m_coords[GENPointCoords::X]/A.Radius(),A.m_coords[GENPointCoords::Y]/A.Radius(),A.m_coords[GENPointCoords::Z]/A.Radius());
	else
		return Origin();
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER> Point<NUMBER>::CartesianNormalised(NUMBER x, NUMBER y, NUMBER z)
{
	return UnitVector(Cartesian(x,y,z));
}

// ..........................................................................

template<class NUMBER>
const Point<NUMBER>& Point<NUMBER>::Origin()
{
	// a little optimisation here - origin point is often required as a
    // temporary in calcs, so we create a static const version to return as a
    // const reference
	static const Point<NUMBER> origin(0,0,0);
	return origin;
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>::Point() :
	m_sphericalCalculated(false)
{
	std::fill(m_coords,m_coords+3,0.f);
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>::Point(const Angle<NUMBER> &altitude, const Angle<NUMBER> &azimuth, NUMBER radius) :
	m_altitude(altitude),
	m_azimuth(azimuth),
	m_radius(radius),
	m_sphericalCalculated(true)
{
	m_coords[GENPointCoords::X]=sin(azimuth)*cos(altitude)*m_radius;
	m_coords[GENPointCoords::Y]=cos(azimuth)*cos(altitude)*m_radius;
	m_coords[GENPointCoords::Z]=sin(altitude)*m_radius;
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>::Point(NUMBER Xcoord, NUMBER Ycoord, NUMBER Zcoord) :
	m_sphericalCalculated(false)
{
	m_coords[GENPointCoords::X]=Xcoord;
	m_coords[GENPointCoords::Y]=Ycoord;
	m_coords[GENPointCoords::Z]=Zcoord;
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>::Point(const Point<NUMBER> &A) :
	m_altitude(A.m_altitude),
	m_azimuth(A.m_azimuth),
	m_radius(A.m_radius),
	m_sphericalCalculated(A.m_sphericalCalculated)
{
	m_coords[0]=A.m_coords[0];
	m_coords[1]=A.m_coords[1];
	m_coords[2]=A.m_coords[2];
}

// ..........................................................................

template<class NUMBER>
template <class N>
Point<NUMBER>::operator Point<N>() const
{
	return Point<N>::Cartesian(m_coords[0],m_coords[1],m_coords[2]);
}


// ..........................................................................

template<class NUMBER>
Point<NUMBER>& Point<NUMBER>::operator=(const Point<NUMBER> &a)
{
	if (this!=&a)
	{
		this->m_coords[0]=a.m_coords[0];
		this->m_coords[1]=a.m_coords[1];
		this->m_coords[2]=a.m_coords[2];

		this->m_altitude=a.m_altitude;
		this->m_azimuth=a.m_azimuth;
		this->m_radius=a.m_radius;
		this->m_sphericalCalculated=a.m_sphericalCalculated;
	}

	return *this;
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>& Point<NUMBER>::operator-=(const Point<NUMBER> &a)
{
	this->m_coords[0]-=a.m_coords[0];
	this->m_coords[1]-=a.m_coords[1];
	this->m_coords[2]-=a.m_coords[2];

	this->m_sphericalCalculated=false;

	return *this;
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>& Point<NUMBER>::operator+=(const Point<NUMBER> &a)
{
	this->m_coords[0]+=a.m_coords[0];
	this->m_coords[1]+=a.m_coords[1];
	this->m_coords[2]+=a.m_coords[2];

	this->m_sphericalCalculated=false;

	return *this;
}

// ..........................................................................

template<class NUMBER>
Point<NUMBER>& Point<NUMBER>::operator*=(NUMBER k)
{
	this->m_coords[0]*=k;
	this->m_coords[1]*=k;
	this->m_coords[2]*=k;

	this->m_sphericalCalculated=false;

	return *this;
}

template<class NUMBER>
Point<NUMBER>& Point<NUMBER>::operator/=(NUMBER k)
{
	this->m_coords[0]/=k;
	this->m_coords[1]/=k;
	this->m_coords[2]/=k;

	this->m_sphericalCalculated=false;

	return *this;
}
// ..........................................................................

template<class NUMBER>
NUMBER& Point<NUMBER>::operator[](unsigned int i)
{
	this->m_sphericalCalculated=false;
	return this->m_coords[i];
}

// ..........................................................................

template<class NUMBER>
const NUMBER& Point<NUMBER>::operator[](unsigned int i) const
{
	return this->m_coords[i];
}

// ..........................................................................

template<class NUMBER>
const Angle<NUMBER>& Point<NUMBER>::Altitude() const
{
	if (!m_sphericalCalculated)
		CalculateSphericalCoordinates();

	return m_altitude;
}

// ..........................................................................

template<class NUMBER>
const Angle<NUMBER>& Point<NUMBER>::Azimuth() const
{
	if (!m_sphericalCalculated)
		CalculateSphericalCoordinates();

	return m_azimuth;
}

// ..........................................................................

template<class NUMBER>
NUMBER Point<NUMBER>::Radius() const
{
	if (!m_sphericalCalculated)
		CalculateSphericalCoordinates();

	return m_radius;
}

// ..........................................................................

template<class NUMBER>
void Point<NUMBER>::CalculateSphericalCoordinates() const
{
	m_radius=sqrt(dot_product(*this,*this));

	CalculateAltitude();
	CalculateAzimuth();
	m_sphericalCalculated=true;
}

// ..........................................................................

template<class NUMBER>
void Point<NUMBER>::CalculateAltitude() const
{
	if (0==m_radius)
	{
		m_altitude=Angle<NUMBER>();
	}
	else
	{
		m_altitude=Angle<NUMBER>::Radians(asin(m_coords[GENPointCoords::Z]/m_radius));
	}
}

// ..........................................................................

template<class NUMBER>
void Point<NUMBER>::CalculateAzimuth() const
{
	// shouldn't do atan2(0,0) (occurs when altitude = +/- 90 degrees)
	if (m_coords[GENPointCoords::X]==0 && m_coords[GENPointCoords::Y]==0)
	{
		m_azimuth=Angle<NUMBER>();
	}
	else
	{
		m_azimuth=Angle<NUMBER>::Radians(atan2(m_coords[GENPointCoords::X],m_coords[GENPointCoords::Y]));
	}
}


// ..........................................................................

template<class NUMBER>
const Point<NUMBER> operator+(const Point<NUMBER> &a, const Point<NUMBER> &b)
{
	Point<NUMBER> c=a;
	c+=b;
	return c;
}

// ..........................................................................

template<class NUMBER>
const Point<NUMBER> operator-(const Point<NUMBER> &a, const Point<NUMBER> &b)
{
	Point<NUMBER> c=a;
	c-=b;
	return c;
}

// ..........................................................................

template<class NUMBER, class SCALAR>
const Point<NUMBER> operator/(const Point<NUMBER> &a, SCALAR k)
{
	Point<NUMBER> tmp(a);
	tmp/=k;
	return tmp;
}

// ..........................................................................

template<class NUMBER, class SCALAR>
const Point<NUMBER> operator*(const Point<NUMBER> &a, SCALAR k)
{
	Point<NUMBER> tmp(a);
	tmp*=k;
	return tmp;
}

// ..........................................................................

template<class NUMBER>
const Point<NUMBER> operator-(const Point<NUMBER> &a)
{
	return Point<NUMBER>::Origin()-a;
}

// ..........................................................................

template <class NUMBER>
bool ApproxEqual(const Point<NUMBER> &a, const Point<NUMBER> &b)
{
	const static NUMBER EPSILON=static_cast<NUMBER>(1e-4);
	const NUMBER distanceABSquared=length_squared(a-b);
	return distanceABSquared<EPSILON;
}

// ..........................................................................

template <class NUMBER>
NUMBER dot_product(const Point<NUMBER> &a, const Point<NUMBER> &b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ..........................................................................

template <class NUMBER>
Point<NUMBER> cross_product(const Point<NUMBER> &a,const  Point<NUMBER> &b)
{
	return Point<NUMBER>::Cartesian(a[GENPointCoords::Y]*b[GENPointCoords::Z]-b[GENPointCoords::Y]*a[GENPointCoords::Z],
											-(a[GENPointCoords::X]*b[GENPointCoords::Z]-b[GENPointCoords::X]*a[GENPointCoords::Z]),
											a[GENPointCoords::X]*b[GENPointCoords::Y]-b[GENPointCoords::X]*a[GENPointCoords::Y]);

}

// ..........................................................................

template <class NUMBER>
NUMBER length_squared(const Point<NUMBER> &a)
{
	return dot_product(a,a);
}

// ..........................................................................

template <class NUMBER>
NUMBER distance_squared(const Point<NUMBER> &a, const Point<NUMBER> &b)
{
	return length_squared(a-b);
}

// ..........................................................................
template <class NUMBER>
NUMBER distance_from_line_squared(const Point<NUMBER> &point,
							  const Point<NUMBER> &line_origin,
							  const Point<NUMBER> &line_direction)
{
	return length_squared(cross_product(line_direction,line_origin-point))/length_squared(line_direction);
}


// ..........................................................................

template <class NUMBER>
std::ostream& operator<<(std::ostream &o, const Point<NUMBER>& A)
{
	o << A[0] << "," << A[1] << "," << A[2];
	return o;
}
