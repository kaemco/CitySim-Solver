#ifndef _INC_RENUtilities_INCLUDED
#define _INC_RENUtilities_INCLUDED

#include "RENMatrix.h"

namespace REN
{
	RENMatrix BuildPerspectiveProjectionMatrix(float fov, float aspect, float znear, float zfar);

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
	RENMatrix BuildOrthogonalProjectionMatrix(float w, float h, float znear, float zfar);

	/**
	* Calculate a "look at" matrix (equivalent to D3DXMatrixLookAtRH (needs transpose?)):
	* zaxis = normal(Eye - At)
	* xaxis = normal(cross(Up, zaxis))
	* yaxis = cross(zaxis, xaxis)
	*
	* xaxis.x           yaxis.x           zaxis.x          0
	* xaxis.y           yaxis.y           zaxis.y          0
	* xaxis.z           yaxis.z           zaxis.z          0
	* -dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1
	*
	* \param eye - location of view point
	* \param at - point to look at
	* \param up
	*/
	RENMatrix BuildLookAtMatrix(
							const GENPoint &eye,
							const GENPoint &at,
							const GENPoint &up);

}
#endif //_INC_RENUtilities_INCLUDED
