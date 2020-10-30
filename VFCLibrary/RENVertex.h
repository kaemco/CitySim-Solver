#ifndef _INC_RENVERTEX_INCLUDED
#define _INC_RENVERTEX_INCLUDED

#include <ostream>

#include "GENPoint.h"

/*
 * A very simple class to represent a vertex within the software renderer (a vertex
 * has position (x,y,z,w) and colour)
 */
class RENVertex
{
public:

	RENVertex() {}

	static RENVertex CreateCartesian(float x, float y, float z, unsigned int colour=0)
	{
		return RENVertex(x,y,z,1,colour);
	}

	static RENVertex CreateHomogenous(float x, float y, float z, float w, unsigned int colour=0)
	{
		return RENVertex(x,y,z,w,colour);
	}

	operator GENPoint() const {

		return GENPoint::Cartesian(m_position[0]/m_position[3],
												m_position[1]/m_position[3],
												m_position[2]/m_position[3]);
	}

	/// The colour of the vertex
	unsigned int m_colour;

	/// The position of the vertex (x,y,z,w)
	float m_position[4];

private:
	RENVertex(float x, float y, float z, float w, unsigned int colour=0) :
		m_colour(colour)
	{
		m_position[0]=x;
		m_position[1]=y;
		m_position[2]=z;
		m_position[3]=w;
	}
};

bool operator==(const RENVertex &a, const RENVertex &b);
std::ostream& operator<<(std::ostream &o, const RENVertex &a);

#endif //_INC_RENVERTEX_INCLUDED
