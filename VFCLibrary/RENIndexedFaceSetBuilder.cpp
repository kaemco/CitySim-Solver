#include "./RENIndexedFaceSetBuilder.h"
#include "DATASurfaceDelegateABC.h"
#include "GEOMPolygonInfo.h"

// ..........................................................................

RENIndexedFaceSetBuilder::RENIndexedFaceSetBuilder(std::shared_ptr<RENIndexedFaceSet> model) :
	m_model(model)
{
	m_gluTriangulator.RegisterCallbacks(this);
}
// ..........................................................................

RENIndexedFaceSetBuilder::~RENIndexedFaceSetBuilder() {
	// free up any vertices still stored in the array (only required if this
	// object is destructed part way through a tesselation, i.e. an error
	// occured)
	for (TemporaryVertexStore_t::iterator it=m_temporaryVertexStore.begin();
		it!=m_temporaryVertexStore.end();
		++it)
	{
		delete *it;
	}
}

// ..........................................................................

// coordinates sorted first on x, then y, then z
bool operator<(RENVertex a, RENVertex b)
{
	const float TOLERANCE=1e-7f;

	for (unsigned int i=0; i<=2; ++i)
	{
		if (a.m_position[i]<b.m_position[i]-TOLERANCE)
			return true;
		if (a.m_position[i]>b.m_position[i]+TOLERANCE)
			return false;
	}

	return false;
}

// ..........................................................................

void RENIndexedFaceSetBuilder::BeginPolygon()
{
	m_vertNo=0;
}

// ..........................................................................

void RENIndexedFaceSetBuilder::CombineData(double coords[3],
                    void *vertex_data[4],
                    float /*weight*/[4], void **dataOut)
{

	RENVertex *vertex = new RENVertex(RENVertex::CreateCartesian((float)coords[0],
														(float)coords[1],
														(float)coords[2],
														static_cast<RENVertex*>(vertex_data[0])->m_colour));
	m_temporaryVertexStore.insert(vertex);

	*dataOut=vertex;
}

// ..........................................................................

void RENIndexedFaceSetBuilder::AddNextVertex(void *newVert)
{
	RENVertex *vertex=static_cast<RENVertex*>(newVert);
	m_temporaryVertexStore.insert(vertex);

	bool firstVert=(m_vertNo==0) ? true : false;

	if (firstVert)
		m_model->AddPlaneEquation(m_currentNormal,m_currentOffset);

	m_vertNo=(m_vertNo+1)%3;

	unsigned int vertIndx=AddVertex(*vertex, firstVert);
	m_model->AddIndex(vertIndx);
}

// ..........................................................................

unsigned int RENIndexedFaceSetBuilder::AddVertex(const RENVertex& newVert, bool IsFirst)
{
	// basically we want a method to avoid duplicate vertices, so if the vertex we are trying to add
	// already exists, return its index, otherwise add a new one to the list.
	// NOTE:  because we are using flat shading, if this is the first vertex of a triangle the colour must
	//        be correct, otherwise colour doesn't matter

	static std::multimap<RENVertex, unsigned int, std::less<RENVertex> >::iterator it;

	it=m_vertexMap.find(newVert);
	if (it!=m_vertexMap.end())
	{
		if (!IsFirst || ((*it).first.m_colour == newVert.m_colour) )
		{
			return (*it).second;
		}
		++it;
		while(!(newVert<(*it).first) && it!=m_vertexMap.end())
		{
			if ((*it).first.m_colour == newVert.m_colour)
			{
				return (*it).second;
			}
			++it;
		}
	}

	unsigned int vertIndx=m_model->AddVertex(newVert);
	m_vertexMap.insert(std::make_pair(newVert,vertIndx));

	return vertIndx;
}

// ..........................................................................

void RENIndexedFaceSetBuilder::SetCurrentPlaneEquation(const GENPoint &normal, float offset)
{
	m_currentNormal=normal;
	m_currentOffset=offset;
}

// ..........................................................................

void RENIndexedFaceSetBuilder::AddSurface(DATASurfaceDelegateABC *surface, unsigned int colour)
{
	struct Vertices : public DATASurfaceDelegateABC::VertexVisitor
	{
		virtual void operator()(const GENPoint& v)
		{
			vertices.push_back(v);
		}

		std::vector<GENPoint> vertices;
	} vs;
	surface->sendVertices(vs);

	GEOMPolygonInfo info(vs.vertices.begin(),
		vs.vertices.end(),surface->normal());

	SetCurrentPlaneEquation(info.normal(),dot_product(info.normal(), info.centroid()));

	if (info.verticesAreOrderedAnticlockwise())
	{
		AddSurface(vs.vertices.begin(),vs.vertices.end(),colour);
	}
	else
	{
		AddSurface(vs.vertices.rbegin(),vs.vertices.rend(),colour);
	}
}
