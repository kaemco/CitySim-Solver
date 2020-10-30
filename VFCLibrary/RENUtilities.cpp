#include "RENUtilities.h"

#include <math.h>

namespace REN
{
	RENMatrix BuildPerspectiveProjectionMatrix(float fov, float aspect, float znear, float zfar)
	{
		RENMatrix matrix;
		const float h = 1/tan(fov/2);
		const float w = h*aspect;
		const float dz=znear-zfar;

		matrix(0,0)=w;
		matrix(0,1)=0;
		matrix(0,2)=0;
		matrix(0,3)=0;

		matrix(1,0)=0;
		matrix(1,1)=h;
		matrix(1,2)=0;
		matrix(1,3)=0;

		matrix(2,0)=0;
		matrix(2,1)=0;
		matrix(2,2)=zfar/dz;
		matrix(2,3)=znear*zfar/dz;

		matrix(3,0)=0;
		matrix(3,1)=0;
		matrix(3,2)=-1;
		matrix(3,3)=0;

		return matrix;
	}

	// ..........................................................................

	/**
	* Calculates the projection matrix (same as D3DXMatrixOrthoRH function (needs transpose?)):
	* 2/w  0    0           0
	* 0    2/h  0           0
	* 0    0    1/(zn-zf)   0
	* 0    0    zn/(zn-zf)  1
	*
	* \param w - width of view
	* \param h - height of view
	* \param znear - location of near clipping plane
	* \param zfar - location of far clipping plane
	*/
	RENMatrix BuildOrthogonalProjectionMatrix(float w, float h, float znear, float zfar)
	{
		RENMatrix matrix;

		matrix(0,0)=2/w;
		matrix(0,1)=0;
		matrix(0,2)=0;
		matrix(0,3)=0;

		matrix(1,0)=0;
		matrix(1,1)=2/h;
		matrix(1,2)=0;
		matrix(1,3)=0;

		matrix(2,0)=0;
		matrix(2,1)=0;
		matrix(2,2)=1/(znear-zfar);
		matrix(2,3)=znear/(znear-zfar);

		matrix(3,0)=0;
		matrix(3,1)=0;
		matrix(3,2)=0;
		matrix(3,3)=1;

		return matrix;
	}

		RENMatrix BuildLookAtMatrix(
							const GENPoint &eye,
							const GENPoint &at,
							const GENPoint &up)
	{
		GENPoint z=GENPoint::UnitVector(eye-at);

		GENPoint x=GENPoint::UnitVector(cross_product(up,z));

		// y guaranteed to be length 1 since z and x are
		GENPoint y=cross_product(z, x);

		RENMatrix matrix;
		matrix(0,0) = x[GENPointCoords::X];
		matrix(0,1) = x[GENPointCoords::Y];
		matrix(0,2) = x[GENPointCoords::Z];
		matrix(0,3) = -dot_product(x,eye);

		matrix(1,0) = y[GENPointCoords::X];
		matrix(1,1) = y[GENPointCoords::Y];
		matrix(1,2) = y[GENPointCoords::Z];
		matrix(1,3) = -dot_product(y,eye);

		matrix(2,0) = z[GENPointCoords::X];
		matrix(2,1) = z[GENPointCoords::Y];
		matrix(2,2) = z[GENPointCoords::Z];
		matrix(2,3) = -dot_product(z,eye);

		matrix(3,0) = 0;
		matrix(3,1) = 0;
		matrix(3,2) = 0;
		matrix(3,3) = 1;

		return matrix;
	}
}
