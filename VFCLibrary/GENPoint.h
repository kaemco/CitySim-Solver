#ifndef _INC_GENPOINT_INCLUDED
#define _INC_GENPOINT_INCLUDED

// SYSTEM INCLUDES
//
#include <iostream>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include "GENAngle.h"

namespace GENPointCoords
{
	// helper enum for accessing x, y and z coordinates in GENPoint
	enum { X=0, Y=1, Z=2};
}

namespace GEN
{
	/* Represents a three dimensional point, handles conversions between spherical and cartesian
	 * coordinates as necessary
	 */
	template <class NUMBER>
	class  Point
	{
	public:
		typedef NUMBER number_type;

		Point(); ///< default constructor creates point at origin
		Point& operator=(const Point &a);
		Point(const Point &A);

		template <class N>
		operator Point<N>() const;

		static Point UnitVector(const Point &A);
		static Point CartesianNormalised(NUMBER Xcoord, NUMBER Ycoord, NUMBER Zcoord);
		static Point Cartesian(NUMBER Xcoord, NUMBER Ycoord, NUMBER Zcoord);
		static Point Polar(const Angle<NUMBER> &altitude, const Angle<NUMBER> &azimuthFromNorth, NUMBER radius=1);
		static const Point& Origin();

		Point& operator+=(const Point &a);
		Point& operator-=(const Point &a);
		Point& operator*=(NUMBER k);
		Point& operator/=(NUMBER k);

		NUMBER& operator[](unsigned int i);			///< access element i
		const NUMBER& operator[](unsigned int i) const;///< access element i

		const Angle<NUMBER>& Altitude() const;
		const Angle<NUMBER>& Azimuth() const;
		NUMBER Radius() const;

	private:
		// construct polar point
		Point(const Angle<NUMBER> &altitudeRads, const Angle<NUMBER> &azimuthRads, NUMBER radius);

		// construct cartesian point
		Point(NUMBER Xcoord, NUMBER Ycoord, NUMBER Zcoord);

		void CalculateSphericalCoordinates() const;
		void CalculateAltitude() const;
		void CalculateAzimuth() const;

		/// store the actual co-ordinates (3 dimensional points)
		NUMBER m_coords[3];

		mutable Angle<NUMBER> m_altitude;
		mutable Angle<NUMBER> m_azimuth;
		mutable NUMBER m_radius;
		mutable bool m_sphericalCalculated; // note - lazy evaluation of alt and az is used
	};

	template <class NUMBER>
	const Point<NUMBER> operator+(const Point<NUMBER> &a, const Point<NUMBER> &b);
	template <class NUMBER>
	const Point<NUMBER> operator-(const Point<NUMBER> &a, const Point<NUMBER> &b);
	template <class NUMBER, class SCALAR>
	const Point<NUMBER> operator/(const Point<NUMBER> &a, SCALAR b);
	template <class NUMBER, class SCALAR>
	const Point<NUMBER> operator*(const Point<NUMBER> &a, SCALAR b);
	template <class NUMBER>
	const Point<NUMBER> operator-(const Point<NUMBER> &a);

	template <class NUMBER>
	Point<NUMBER> cross_product(const Point<NUMBER> &a, const Point<NUMBER> &b);
	template <class NUMBER>
	NUMBER dot_product(const Point<NUMBER> &a, const Point<NUMBER> &b);

	template <class NUMBER>
	NUMBER length_squared(const Point<NUMBER> &a);
	template <class NUMBER>
	NUMBER distance_squared(const Point<NUMBER> &a, const Point<NUMBER> &b);
	template <class NUMBER>
	NUMBER distance_from_line_squared(const Point<NUMBER> &point,
								const Point<NUMBER> &line_origin,
								const Point<NUMBER> &line_direction);

	template <class NUMBER>
	std::ostream& operator<<(std::ostream &o, const Point<NUMBER>& a);

	template <class NUMBER>
	bool ApproxEqual(const Point<NUMBER> &a, const Point<NUMBER> &b);

	template <class NUMBER>
	bool operator==(const Point<NUMBER> &a, const Point<NUMBER> &b)
	{
		return ApproxEqual(a,b);
	}

	template <class NUMBER>
	bool operator!=(const Point<NUMBER> &a, const Point<NUMBER> &b)
	{
		return !(a==b);
	}


	#include "GENPoint.inl"
}

// we get a significant increase in speed by using float instead of double points
typedef GEN::Point<float> GENPoint;


bool operator<(const GENPoint &a, const GENPoint &b);

float closest_point_on_line(const GENPoint &point, const GENPoint &a, const GENPoint &b);

// returns true if points are (roughly) collinear
bool collinear(const GENPoint &a,const GENPoint &b,const GENPoint &c);


#endif // _INC_GENPOINT_INCLUDED
