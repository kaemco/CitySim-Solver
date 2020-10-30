#include <math.h>

#include "./RENVertex.h"

bool operator==(const RENVertex &a, const RENVertex &b)
{
	if (a.m_colour!=b.m_colour) 
		return false;
	for (unsigned int i=0; i<4; ++i)
	{
		if (fabs(a.m_position[i]-b.m_position[i])>1e-5)
			return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream &o, const RENVertex &a)
{
	o << "{" << a.m_position[0] << "," << a.m_position[1] <<"," << a.m_position[2] << "," << a.m_position[3] << " col=" << a.m_colour;
	return o;
}
