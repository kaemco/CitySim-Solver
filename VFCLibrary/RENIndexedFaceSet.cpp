#include "./RENIndexedFaceSet.h"

#include <memory>
#include <iostream>
#include <cstring>
#include <sstream>
#include <limits>
#include <cmath>

#include "GEOMBoundingSphereCalc.h"
class WEMSphere : public RENBoundingSphereABC
{
public:
	virtual GENPoint Centre() const {return centre;}
	virtual float Radius() const {return radius; }

	GENPoint& Centre() {return centre;}
	float& Radius() {return radius; }

private:
	GENPoint centre;
	float radius;
};



RENIndexedFaceSet::RENIndexedFaceSet(void)
{
}

// ..........................................................................

RENIndexedFaceSet::~RENIndexedFaceSet(void)
{
}

// ..........................................................................

const RENVertex& RENIndexedFaceSet::GetVertex(int n) const
{
	return m_vertices[n];
}

// ..........................................................................

unsigned int RENIndexedFaceSet::VertexCount() const
{
	return static_cast<unsigned int>(m_vertices.size());
}

// ..........................................................................

const std::vector<unsigned int>& RENIndexedFaceSet::GetIndices() const
{
	return m_indices;
}

// ..........................................................................
const GENPoint& RENIndexedFaceSet::GetNormal(unsigned int polygon)
{
	return m_normals[polygon];
}

// ..........................................................................

float RENIndexedFaceSet::EvaluatePlaneEquation(unsigned int polygon, const GENPoint& point)
{
	return dot_product(point,m_normals[polygon])-m_offsets[polygon];
}

// ..........................................................................

unsigned int RENIndexedFaceSet::AddVertex(const RENVertex &vertex) {
	m_vertices.push_back(vertex);
	return static_cast<unsigned int>(m_vertices.size()-1);
}

// ..........................................................................

void RENIndexedFaceSet::AddPlaneEquation(const GENPoint &normal, float offset) {
	m_normals.push_back(normal);
	m_offsets.push_back(offset);
}

// ..........................................................................

void RENIndexedFaceSet::AddIndex(unsigned int i) {
	m_indices.push_back(i);
}

// ..........................................................................

unsigned int RENIndexedFaceSet::TriangleCount() {
	return static_cast<unsigned int>(m_indices.size()/3);
}

// ..........................................................................

std::unique_ptr<RENBoundingSphereABC> RENIndexedFaceSet::GetBoundingSphere() const
{
	GEOMBoundingSphereCalc boundingSphereCalculator;
    GEOMSphere s;

    try {
        s=boundingSphereCalculator.Calculate(m_vertices.begin(),m_vertices.end());
    }
    catch (std::string msg) {
        std::cerr << msg << " " << std::endl;
    }

    std::unique_ptr<WEMSphere> boundingSphere(new WEMSphere);
	boundingSphere->Radius()=s.Radius();
	boundingSphere->Centre()=s.Centre();
    return boundingSphere;
}
// ..........................................................................

std::pair<float,float> RENIndexedFaceSet::getBoundingZ() const
{
    // define the min and max in Z of the scene
    float Zmin=std::numeric_limits<float>::max(), Zmax=std::numeric_limits<float>::min();

    for (size_t i=0;i<m_vertices.size();++i) {
        Zmin=fmin(GENPoint(m_vertices.at(i))[GENPointCoords::Z],Zmin);
        Zmax=fmax(GENPoint(m_vertices.at(i))[GENPointCoords::Z],Zmax);
    }

    return std::pair<float,float>(Zmin,Zmax);
}
