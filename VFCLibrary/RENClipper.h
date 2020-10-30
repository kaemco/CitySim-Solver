#ifndef _INC_RENCLIPPER_INCLUDED
#define _INC_RENCLIPPER_INCLUDED

#include "RENVertex.h"
#include <vector>

/*
 * Clip a triangle in homogenous perspective space - i.e. view volume is bounded
 * by a cuboid from -1<x<1, -1<y<1, 0<z<1
 */
class RENClipper
{
public:
	RENClipper() {}
	~RENClipper();

	/// Clip a triangle
	RENVertex** ClipTriangle(RENVertex &V1,
									RENVertex &V2,
									RENVertex &V3);

	/// Clip a line
	RENVertex** ClipLine(RENVertex &V1,
								 RENVertex &V2);

	/// Return the resulting number of vertices in the last shape that was clipped
	int GetNumVertices() const;

private:
	RENClipper(const RENClipper& ) {
		throw std::string(std::string("logic_error: " + std::string("not implemented")));
	}
	RENClipper& operator=(const RENClipper& ) {
		throw std::string(std::string("logic_error: " + std::string("not implemented")));
	}

	void clipAll();

	typedef float(RENClipper::*distance_func)(const RENVertex *V) const;
	float dNear(const RENVertex *V) const;
	float dFar(const RENVertex *V) const;
	float dLeft(const RENVertex *V) const;
	float dRight(const RENVertex *V) const;
	float dTop(const RENVertex *V) const;
	float dBottom(const RENVertex *V) const;

	void clip(distance_func distanceFromClipEdge);

	void clearStore();
	void addToStore(RENVertex* vertex);
	std::vector<RENVertex*> m_vertexStore;

	std::vector<RENVertex*> m_P[7];   // Pointers to polygon's vertices, indexed by stage

	int m_numVertices;
	int m_stage;
};

#endif //_INC_RENCLIPPER_INCLUDED
