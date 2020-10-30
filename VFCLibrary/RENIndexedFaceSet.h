#ifndef _INC_RENINDEXEDFACESET_INCLUDED
#define _INC_RENINDEXEDFACESET_INCLUDED

#include <vector>
#include <memory>
#include <utility>

#include "RENVertex.h"
#include "GENPoint.h"
#include "RENBoundingSphereABC.h"

/* Store geometry in indexed face set format
 */
class RENIndexedFaceSet
{
public:
	RENIndexedFaceSet(void);
	virtual ~RENIndexedFaceSet(void);

	const RENVertex& GetVertex(int n) const;

	const std::vector<unsigned int>& GetIndices() const;
	unsigned int VertexCount() const;

	const GENPoint& GetNormal(unsigned int polygon);
	float EvaluatePlaneEquation(unsigned int polygon, const GENPoint& point);

	unsigned int AddVertex(const RENVertex &vertex);
	void AddPlaneEquation(const GENPoint &normal, float offset);
	void AddIndex(unsigned int i);
	unsigned int TriangleCount();

    std::unique_ptr<RENBoundingSphereABC> GetBoundingSphere() const;
	std::pair<float,float> getBoundingZ() const;

private:
	std::vector<RENVertex> m_vertices;
	std::vector<unsigned int> m_indices;

	// describe the plane equation (dot_product(normal,point)=offset)
	std::vector<GENPoint> m_normals;
	std::vector<float> m_offsets;
};


#endif //_INC_RENINDEXEDFACESET_INCLUDED
