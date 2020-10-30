#ifndef _INC_GENANGLE_INCLUDED
#define _INC_GENANGLE_INCLUDED

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
//#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif // M_PI

namespace GEN
{
	/* This class is used to represent the concept of an angle - basically it makes it difficult to
	 * accidentally treat a quantity in radians as a quantity in degrees!
	 * To create:
	 * Angle<double>::Degrees(90)
	 * Angle<float>::Radians(M_PI)
	 *
	 * To use:
	 * Angle<double> a1=Angle<double>::Degrees(90)
	 * double inRadians=a1.radians()
	 *
	 * Arithmetic operations also exist:
	 * Angle<double> a2=Angle<double>::Radians(M_PI)
	 * Angle<double> a3=a1+a2   <-- result is 270 degrees
	 */
	template <class NUMBER>
	class Angle
	{
	public:
		static Angle Radians(NUMBER radians) {
			return Angle(radians);
		}
		static Angle Degrees(NUMBER degrees) {
            return Angle(static_cast<NUMBER>(degrees*M_PI/180.));
		}
		Angle() : m_radians(0) {}

		NUMBER radians() const { return m_radians; }
        NUMBER degrees() const { return static_cast<NUMBER>(m_radians*180/M_PI); }

		template <class OTHER_NUMBER>
		operator Angle<OTHER_NUMBER>() const {
			return Angle<OTHER_NUMBER>::Radians(m_radians);
		}

		Angle<NUMBER>& operator+=(const Angle<NUMBER>& a) {
			m_radians+=a.m_radians;
			normalise();
			return *this;
		}
		Angle<NUMBER>& operator-=(const Angle<NUMBER>& a) {
			m_radians-=a.m_radians;
			normalise();
			return *this;
		}
		Angle<NUMBER>& operator*=(NUMBER c) {
			m_radians*=c;
			normalise();
			return *this;
		}
		Angle<NUMBER>& operator/=(NUMBER c) {
			m_radians/=c;
			normalise();
			return *this;
		}

	private:
		void normalise() {
			while (m_radians>(NUMBER)M_PI)
				m_radians-=2*(NUMBER)M_PI;
			while (m_radians<=-(NUMBER)M_PI)
				m_radians+=2*(NUMBER)M_PI;
		}
		Angle(NUMBER radians) : m_radians(radians) {
			 normalise();
		}

		NUMBER m_radians;
	};

	template <class NUMBER>
	Angle<NUMBER> operator+(const Angle<NUMBER>& a, const Angle<NUMBER>& b)
	{
		Angle<NUMBER> r=a;
		r+=b;
		return r;
	}
	template <class NUMBER>
	Angle<NUMBER> operator-(const Angle<NUMBER>& a, const Angle<NUMBER>& b)
	{
		Angle<NUMBER> r=a;
		r-=b;
		return r;
	}
	template <class NUMBER>
	Angle<NUMBER> operator*(const Angle<NUMBER>& a, NUMBER c)
	{
		Angle<NUMBER> r=a;
		r*=c;
		return r;
	}
	template <class NUMBER>
	Angle<NUMBER> operator/(const Angle<NUMBER>& a, NUMBER c)
	{
		Angle<NUMBER> r=a;
		r/=c;
		return r;
	}

	template <class NUMBER>
	NUMBER sin(const GEN::Angle<NUMBER> &a) {
		return ::sin(a.radians());
	}
	template <class NUMBER>
	NUMBER cos(const GEN::Angle<NUMBER> &a) {
		return ::cos(a.radians());
	}
}

// need to bring these into global scope so that sin(Angle a) doesn't hide the global sin function
using GEN::sin;
using GEN::cos;

// we get a significant increase in speed by using float instead of double points
typedef GEN::Angle<float> GENAngle;


#endif // _INC_GENANGLE_INCLUDED
