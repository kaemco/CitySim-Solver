#ifndef _INC_RENINDEXEDFACESETBUILDER_INCLUDED
#define _INC_RENINDEXEDFACESETBUILDER_INCLUDED

#include <map>
#include <set>

#include "GENGluTriangulator.h"
#include "RENIndexedFaceSet.h"

/* Builds an indexed face set representation of geometry - simply call
 * AddSurface repeatedly until you are done
 */
class DATASurfaceDelegateABC;

class RENIndexedFaceSetBuilder : public GENGluTriangulatorCallbacks
{
public:
    RENIndexedFaceSetBuilder(std::shared_ptr<RENIndexedFaceSet> model);
	~RENIndexedFaceSetBuilder();

//  commented by JK - not used in this code 02/02/09
//	template <class SURFACE_T>
//	void AddSurface(SURFACE_T *surface, unsigned int colour);

	void AddSurface(DATASurfaceDelegateABC *surface, unsigned int colour);

    template <class VERTEX_ITERATOR>
    void AddSurface(VERTEX_ITERATOR begin, const VERTEX_ITERATOR end, unsigned int colour)
    {
        m_gluTriangulator.BeginPolygon();
            m_gluTriangulator.BeginContour();

                for (unsigned int vertpointer=0;
                        begin!=end;
                        ++vertpointer,++begin)
                {
                    const GENPoint& vertex=*begin;
                    m_gluTriangulator.AddVertex(vertex[0],vertex[1],vertex[2],
                                                new RENVertex(RENVertex::CreateCartesian(vertex[0],
                                                              vertex[1],
                                                              vertex[2],
                                                              colour)));
                }

            m_gluTriangulator.EndContour();
        m_gluTriangulator.EndPolygon();
    }

private:

	void SetCurrentPlaneEquation(const GENPoint &normal, float offset);

	void BeginPolygon();
	void AddNextVertex(void *newVert);
	void CombineData(double coords[3],
                     void *vertex_data[4],
                     float weight[4], void **dataOut);

	unsigned int AddVertex(const RENVertex& newVert, bool IsFirst);

	// used to provide fast searching of model vertex vector during triangulation
	std::multimap<RENVertex, unsigned int, std::less<RENVertex> > m_vertexMap;

	typedef std::set<RENVertex*> TemporaryVertexStore_t;
	TemporaryVertexStore_t m_temporaryVertexStore; // used to store the vertices produced during the triangulation routine
										// and which should be deleted afterwards
	int m_vertNo;

    std::shared_ptr<RENIndexedFaceSet> m_model;

	GENPoint m_currentNormal;
	float m_currentOffset;

	GENGluTriangulator m_gluTriangulator;

};

// ..........................................................................

// commented by JK - not used in this code 02/02/09
//template <class SURFACE_T>
//void RENIndexedFaceSetBuilder::AddSurface(SURFACE_T *surface, unsigned int colour)
//{
//	GEOMPolygonInfoImpl<VertexConverter> info(surface->FirstVertex(),
//								surface->LastVertex());
//
//	GENAssert(info.verticesAreOrderedAnticlockwise());
//	SetCurrentPlaneEquation(info.normal(),dot_product(info.normal(), info.centroid()));
//
//	m_gluTriangulator.BeginPolygon();
//		m_gluTriangulator.BeginContour();
//
//            SURFACE_T::VertexIterator_t vertices=surface->FirstVertex();
//			for (unsigned int vertpointer=0;
//					vertpointer < surface->VertexCount();
//					++vertpointer,++vertices)
//			{
//				const GENPoint& vertex=vertices->Data();
//				m_gluTriangulator.AddVertex(vertex[0],vertex[1],vertex[2],
//					new RENVertex(RENVertex::CreateCartesian(vertex[0],
//																		vertex[1],
//																		vertex[2],
//																		colour)));
//			}
//
//		m_gluTriangulator.EndContour();
//	m_gluTriangulator.EndPolygon();
//}
//template <class SURFACE_T>
//void RENIndexedFaceSetBuilder::AddSurface(SURFACE_T *surface, unsigned int colour)
//{
//	GEOMPolygonInfoImpl<VertexConverter> info(bind1st(mem_fun(SURFACE_T::getVertex),surface),surface->vertexCount());
//
//	GENAssert(info.verticesAreOrderedAnticlockwise());
//	SetCurrentPlaneEquation(info.normal(),dot_product(info.normal(), info.centroid()));
//
//	m_gluTriangulator.BeginPolygon();
//		m_gluTriangulator.BeginContour();
//
//			for (unsigned int vertpointer=0;
//					vertpointer<surface->vertexCount();
//					++vertpointer)
//			{
//				const GENPoint& vertex=surface->getVertex(vertpointer);
//				m_gluTriangulator.AddVertex(vertex[0],vertex[1],vertex[2],
//					new RENVertex(RENVertex::CreateCartesian(vertex[0],
//																		vertex[1],
//																		vertex[2],
//																		colour)));
//			}
//
//		m_gluTriangulator.EndContour();
//	m_gluTriangulator.EndPolygon();
//}
#endif //_INC_RENINDEXEDFACESETBUILDER_INCLUDED
