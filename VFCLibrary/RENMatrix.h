#ifndef _INC_RENMATRIX_INCLUDED
#define _INC_RENMATRIX_INCLUDED

#include "RENVertex.h"

class RENMatrix;
namespace RENMatrixUtils
{
	/// Set M to be the identity matrix
	void Identity(RENMatrix &M);

	RENMatrix operator*(const RENMatrix &M, const RENMatrix &N);
	RENVertex operator*(const RENMatrix &M, const RENVertex &V);
}

/*
 * Implementation of an optimised 4x4 matrix for use in vertex transformations, etc.
 */
class RENMatrix
{
public:
	RENMatrix(void);
	~RENMatrix(void);

	/// Access element (row, col), starting with (0,0)
	float& operator()(int i, int j) { return m[i][j]; }
	const float& operator()(int i, int j) const { return m[i][j]; }

	/// multiplication by another matrix (this=this*A)
	RENMatrix& operator*=(const RENMatrix& A);

	/// Inverse
	RENMatrix operator!() const;

	/// Division by a constant
	RENMatrix & operator/=(float a);

	friend RENMatrix RENMatrixUtils::operator*(const RENMatrix &M, const RENMatrix &N);
	friend RENVertex RENMatrixUtils::operator*(const RENMatrix &M, const RENVertex &V);
private:
	// Row major order
	float m[4][4];

};



#endif //_INC_RENMATRIX_INCLUDED
