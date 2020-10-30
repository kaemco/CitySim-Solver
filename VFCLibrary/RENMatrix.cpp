#include "RENMatrix.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER



RENMatrix::RENMatrix(void)
{
}

// ..........................................................................

RENMatrix::~RENMatrix(void)
{
}

// ..........................................................................

/// Inverse
RENMatrix RENMatrix::operator!() const
{
	const RENMatrix &M = *this;
	RENMatrix I;

	float M2233 = M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3);
	float M1233 = M(1, 2) * M(3, 3) - M(3, 2) * M(1, 3);
	float M1223 = M(1, 2) * M(2, 3) - M(2, 2) * M(1, 3);
	float M2133 = M(2, 1) * M(3, 3) - M(3, 1) * M(2, 3);
	float M1133 = M(1, 1) * M(3, 3) - M(3, 1) * M(1, 3);
	float M1123 = M(1, 1) * M(2, 3) - M(2, 1) * M(1, 3);
	float M2132 = M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2);
	float M1132 = M(1, 1) * M(3, 2) - M(3, 1) * M(1, 2);
	float M1122 = M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2);
	float M0233 = M(0, 2) * M(3, 3) - M(3, 2) * M(0, 3);
	float M0223 = M(0, 2) * M(2, 3) - M(2, 2) * M(0, 3);
	float M0133 = M(0, 1) * M(3, 3) - M(3, 1) * M(0, 3);
	float M0123 = M(0, 1) * M(2, 3) - M(2, 1) * M(0, 3);
	float M0132 = M(0, 1) * M(3, 2) - M(3, 1) * M(0, 2);
	float M0122 = M(0, 1) * M(2, 2) - M(2, 1) * M(0, 2);
	float M0213 = M(0, 2) * M(1, 3) - M(1, 2) * M(0, 3);
	float M0113 = M(0, 1) * M(1, 3) - M(1, 1) * M(0, 3);
	float M0112 = M(0, 1) * M(1, 2) - M(1, 1) * M(0, 2);

	// Adjoint Matrix
	I(0, 0) =  M(1, 1) * M2233 - M(2, 1) * M1233 + M(3, 1) * M1223;
	I(1, 0) = -M(1, 0) * M2233 + M(2, 0) * M1233 - M(3, 0) * M1223;
	I(2, 0) =  M(1, 0) * M2133 - M(2, 0) * M1133 + M(3, 0) * M1123;
	I(3, 0) = -M(1, 0) * M2132 + M(2, 0) * M1132 - M(3, 0) * M1122;

	I(0, 1) = -M(0, 1) * M2233 + M(2, 1) * M0233 - M(3, 1) * M0223;
	I(1, 1) =  M(0, 0) * M2233 - M(2, 0) * M0233 + M(3, 0) * M0223;
	I(2, 1) = -M(0, 0) * M2133 + M(2, 0) * M0133 - M(3, 0) * M0123;
	I(3, 1) =  M(0, 0) * M2132 - M(2, 0) * M0132 + M(3, 0) * M0122;

	I(0, 2) =  M(0, 1) * M1233 - M(1, 1) * M0233 + M(3, 1) * M0213;
	I(1, 2) = -M(0, 0) * M1233 + M(1, 0) * M0233 - M(3, 0) * M0213;
	I(2, 2) =  M(0, 0) * M1133 - M(1, 0) * M0133 + M(3, 0) * M0113;
	I(3, 2) = -M(0, 0) * M1132 + M(1, 0) * M0132 - M(3, 0) * M0112;

	I(0, 3) = -M(0, 1) * M1223 + M(1, 1) * M0223 - M(2, 1) * M0213;
	I(1, 3) =  M(0, 0) * M1223 - M(1, 0) * M0223 + M(2, 0) * M0213;
	I(2, 3) = -M(0, 0) * M1123 + M(1, 0) * M0123 - M(2, 0) * M0113;
	I(3, 3) =  M(0, 0) * M1122 - M(1, 0) * M0122 + M(2, 0) * M0112;

	// Division by determinant
	I /= M(0, 0) * I(0, 0) +
		    M(1, 0) * I(0, 1) +
		    M(2, 0) * I(0, 2) +
		    M(3, 0) * I(0, 3);

	return I;
}

// ..........................................................................

RENMatrix& RENMatrix::operator/=(float a)
{
	for (int i=0; i<4; ++i)
	{
		for (int j=0; j<4; ++j)
		{
			m[i][j]/=a;
		}
	}
	return *this;
}

// ..........................................................................

RENMatrix& RENMatrix::operator*=(const RENMatrix &N)
{
	RENMatrix M = *this;

	for (int i=0; i<4; ++i)
	{
		for (int j=0; j<4; ++j)
		{
			(*this)(i,j)=0;
			for (int k=0; k<4; ++k)
			{
				(*this)(i,j)+=M(i,k)*N(k,j);
			}
		}
	}

	return *this;
}

// ..........................................................................

void RENMatrixUtils::Identity(RENMatrix &M)
{
	M(0,0)=1;
	M(0,1)=0;
	M(0,2)=0;
	M(0,3)=0;

	M(1,0)=0;
	M(1,1)=1;
	M(1,2)=0;
	M(1,3)=0;

	M(2,0)=0;
	M(2,1)=0;
	M(2,2)=1;
	M(2,3)=0;

	M(3,0)=0;
	M(3,1)=0;
	M(3,2)=0;
	M(3,3)=1;
}

// ..........................................................................

RENMatrix RENMatrixUtils::operator*(const RENMatrix &M, const RENMatrix &N)
{
	RENMatrix result;

	result.m[0][0]=M.m[0][0]*N.m[0][0] + M.m[0][1]*N.m[1][0] + M.m[0][2]*N.m[2][0] +M.m[0][3]*N.m[3][0];
	result.m[0][1]=M.m[0][0]*N.m[0][1] + M.m[0][1]*N.m[1][1] + M.m[0][2]*N.m[2][1] +M.m[0][3]*N.m[3][1];
	result.m[0][2]=M.m[0][0]*N.m[0][2] + M.m[0][1]*N.m[1][2] + M.m[0][2]*N.m[2][2] +M.m[0][3]*N.m[3][2];
	result.m[0][3]=M.m[0][0]*N.m[0][3] + M.m[0][1]*N.m[1][3] + M.m[0][2]*N.m[2][3] +M.m[0][3]*N.m[3][3];

	result.m[1][0]=M.m[1][0]*N.m[0][0] + M.m[1][1]*N.m[1][0] + M.m[1][2]*N.m[2][0] +M.m[1][3]*N.m[3][0];
	result.m[1][1]=M.m[1][0]*N.m[0][1] + M.m[1][1]*N.m[1][1] + M.m[1][2]*N.m[2][1] +M.m[1][3]*N.m[3][1];
	result.m[1][2]=M.m[1][0]*N.m[0][2] + M.m[1][1]*N.m[1][2] + M.m[1][2]*N.m[2][2] +M.m[1][3]*N.m[3][2];
	result.m[1][3]=M.m[1][0]*N.m[0][3] + M.m[1][1]*N.m[1][3] + M.m[1][2]*N.m[2][3] +M.m[1][3]*N.m[3][3];

	result.m[2][0]=M.m[2][0]*N.m[0][0] + M.m[2][1]*N.m[1][0] + M.m[2][2]*N.m[2][0] +M.m[2][3]*N.m[3][0];
	result.m[2][1]=M.m[2][0]*N.m[0][1] + M.m[2][1]*N.m[1][1] + M.m[2][2]*N.m[2][1] +M.m[2][3]*N.m[3][1];
	result.m[2][2]=M.m[2][0]*N.m[0][2] + M.m[2][1]*N.m[1][2] + M.m[2][2]*N.m[2][2] +M.m[2][3]*N.m[3][2];
	result.m[2][3]=M.m[2][0]*N.m[0][3] + M.m[2][1]*N.m[1][3] + M.m[2][2]*N.m[2][3] +M.m[2][3]*N.m[3][3];

	result.m[3][0]=M.m[3][0]*N.m[0][0] + M.m[3][1]*N.m[1][0] + M.m[3][2]*N.m[2][0] +M.m[3][3]*N.m[3][0];
	result.m[3][1]=M.m[3][0]*N.m[0][1] + M.m[3][1]*N.m[1][1] + M.m[3][2]*N.m[2][1] +M.m[3][3]*N.m[3][1];
	result.m[3][2]=M.m[3][0]*N.m[0][2] + M.m[3][1]*N.m[1][2] + M.m[3][2]*N.m[2][2] +M.m[3][3]*N.m[3][2];
	result.m[3][3]=M.m[3][0]*N.m[0][3] + M.m[3][1]*N.m[1][3] + M.m[3][2]*N.m[2][3] +M.m[3][3]*N.m[3][3];

	return result;
}

// ..........................................................................

RENVertex RENMatrixUtils::operator*(const RENMatrix &M, const RENVertex &V)
{
	RENVertex result;
	result.m_colour=V.m_colour;

	result.m_position[0]=	M.m[0][0]*V.m_position[0] +
							M.m[0][1]*V.m_position[1] +
							M.m[0][2]*V.m_position[2] +
							M.m[0][3]*V.m_position[3];

	result.m_position[1]=	M.m[1][0]*V.m_position[0] +
							M.m[1][1]*V.m_position[1] +
							M.m[1][2]*V.m_position[2] +
							M.m[1][3]*V.m_position[3];

	result.m_position[2]=	M.m[2][0]*V.m_position[0] +
							M.m[2][1]*V.m_position[1] +
							M.m[2][2]*V.m_position[2] +
							M.m[2][3]*V.m_position[3];

	result.m_position[3]=	M.m[3][0]*V.m_position[0] +
							M.m[3][1]*V.m_position[1] +
							M.m[3][2]*V.m_position[2] +
							M.m[3][3]*V.m_position[3];

	return result;
}


