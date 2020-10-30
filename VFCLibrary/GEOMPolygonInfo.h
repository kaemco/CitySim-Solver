#ifndef _INC_GEOMPOLYGONINFO_INCLUDED
#define _INC_GEOMPOLYGONINFO_INCLUDED

#include "GENAssert.h"
#include "GENPoint.h"
#include "RENIndexedFaceSet.h"
#include "RENIndexedFaceSetBuilder.h"

/* This file defines a class used for calculating plane equations of polygons.
 * It is not as easy as you might think because you have to handle things like
 * convex corners.
 *
 */

// a,b,c are vertices anti-clockwise around triangle
template<class POINT, class VECTOR>
float TriangleArea(POINT &a, POINT &b, POINT &c, const VECTOR &normal)
{
	return dot_product(cross_product(b-a,c-a),normal)/2.f;
}

// ..........................................................................

/* Some mucking about with templates in this class because sometimes the vertex
 * iterator needs to reference vertices, some times pointers to vertices
 *
 * VERTEX_CONVERTER is a class that converts from whatever the vertex format is
 * to const GENPoint&, e.g.:
 *
 * const GENPoint &p=VERTEX_CONVERTER(*v)
 *
 */
template <class VERTEX_CONVERTER>
class GEOMPolygonInfoImpl
{
public:
	GEOMPolygonInfoImpl() {}

	/* VERTEX_ITERATOR must be able to be dereferenced to either "const
     * GENPoint &" or "const GENPoint*"
     */
	template <class VERTEX_ITERATOR>
	GEOMPolygonInfoImpl(VERTEX_ITERATOR firstVertex,
						VERTEX_ITERATOR lastVertex)
	{
		m_normal=SignedNormal(VertexIteratorPair<VERTEX_ITERATOR>(firstVertex,lastVertex));

		calculateData(VertexIteratorPair<VERTEX_ITERATOR>(firstVertex,lastVertex),m_normal);
		if (!m_verticesAreAnticlockwise)
		{
			m_normal=-m_normal;
			m_verticesAreAnticlockwise=true;
		}
	}

	template <class VERTEX_ITERATOR>
	GEOMPolygonInfoImpl(VERTEX_ITERATOR firstVertex,
						VERTEX_ITERATOR lastVertex,
						const GENPoint &normal) :
		m_normal(normal)
	{
		m_normal=SignedNormal(VertexIteratorPair<VERTEX_ITERATOR>(firstVertex,lastVertex));
		if (dot_product(normal,m_normal)<0)
			m_normal=-m_normal;

		calculateData(VertexIteratorPair<VERTEX_ITERATOR>(firstVertex,lastVertex),m_normal);
		if (SignedNormal(VertexIteratorPair<VERTEX_ITERATOR>(firstVertex,lastVertex)).Radius()==0)
		{
			// this means that polygon either has zero area or is very close to
            // being a thin strip
			m_area=0;
		}
	}

	float area() const {
		return m_area;
	}
	float volume() const {
		return m_volume;
	}
	bool verticesAreOrderedAnticlockwise() const {
		return m_verticesAreAnticlockwise;
	}

	const GENPoint& normal() const {
		return m_normal;
	}

	const GENPoint& centroid() const {
		return m_centroid;
	}

    // returns the normal from three points (not normalised)
	template<class POINT>
	static GENPoint NormalFromThreePoints(const POINT &a, const POINT &b, const POINT &c)
	{
		const GENPoint v1=(b-a);
		const GENPoint v2=(c-b);
        return cross_product(v1,v2); // was previously normalised
	}

private:
	template <class VERTEX_ITERATOR>
	struct VertexIteratorPair
	{

        VERTEX_ITERATOR first;
		VERTEX_ITERATOR last;

	    VertexIteratorPair(const VERTEX_ITERATOR &first, const VERTEX_ITERATOR &last) : first(first), last(last) {}

		const GENPoint& get() {

	return VERTEX_CONVERTER(dereferencer<typename std::iterator_traits<VERTEX_ITERATOR>::value_type>().deref(first));

            // modified by JK, made it shorter 02/02/09
			//return dereference(first);
			//const GENPoint& dereference(VERTEX_ITERATOR it)
            //{
            //    return VERTEX_CONVERTER(dereferencer<VERTEX_ITERATOR,std::iterator_traits<VERTEX_ITERATOR>::value_type>().deref(it));
            //}

		}
		void next() {
			++first;
		}
		bool isAtEnd() const {
			return first==last;
		}

		/* There are two versions of the struct "dereferencer", one works with a
		* vertex iterator that returns pointers to points, one with an iterator
		* that returns references to points.  The dereference function will
		* automatically use the correct one based on the vertex iterator's
		* ::value_type
		*/
		template<class POINT> // modified by JK, changed the class template name to avoid shadowing 1/02/09
		struct dereferencer
		{
			const POINT& deref(VERTEX_ITERATOR it) const
			{
				return *it;
			}
		};

		template<class POINT> // modified by JK, changed the class template name to avoid shadowing 1/02/09
		struct dereferencer<POINT*>
		{
			const POINT& deref(VERTEX_ITERATOR it) const
			{
				return **it;
			}
		};

	};

	template<class VERTEX_SET>
	void calculateDataDegeneratePolygon(VERTEX_SET vertices)
	{
		m_area=0;
		m_normal=GENPoint::Origin();

		if (vertices.isAtEnd())
		{
			m_centroid=GENPoint::Origin();
			return;
		}
		GENPoint centroidtmp=GENPoint::Origin();

		float count=0;
		for (; !vertices.isAtEnd(); vertices.next(), ++count)
		{
			centroidtmp+=vertices.get();
		}

		m_centroid=centroidtmp/count;
	}

	template<class VERTEX_SET>
	void calculateData(const VERTEX_SET &verts, const GENPoint &normal)
	{
		if (verts.isAtEnd())
		{
			calculateDataDegeneratePolygon(verts);
			return;
		}

		VERTEX_SET vertices=verts;
		const GENPoint &vert1 = vertices.get();

		float summedArea=0;
		GENPoint centroidtmp=GENPoint::Origin();
		unsigned int triangleCount=0;

		if (!vertices.isAtEnd())
			vertices.next();

		GENPoint previous;
		if (!vertices.isAtEnd())
		{
			previous=vertices.get();
			vertices.next();
		}

		for (;
			!vertices.isAtEnd();
			vertices.next(), ++triangleCount)
		{
			GENPoint current=vertices.get();

			// form a triangle from 1st vertex, vertex i and vertex im1
			const GENPoint &vertim1=previous;
			const GENPoint &verti=current;

			const float triangleArea = TriangleArea(vert1,
												vertim1,
												verti,
												normal);

			summedArea+=triangleArea;
			centroidtmp+=(vert1+vertim1+verti)*triangleArea;

			previous=current;
		}

		if (0!=summedArea)
		{
			m_verticesAreAnticlockwise=summedArea>=0;
			m_area=fabs(summedArea);
			m_centroid=centroidtmp/(3*summedArea);
		}
		else
			calculateDataDegeneratePolygon(verts);
	}

	// Calculate unsigned normal - i.e. no guarantees whether it points in or out of polygon
	template <class VERTEX_SET>
	GENPoint UnsignedNormal(VERTEX_SET vertices)
	{

		if (vertices.isAtEnd())
			return GENPoint::Origin();

		// need to find three vertices that are not collinear
		const GENPoint* vertexarr[3];

		vertexarr[0] = &vertices.get();
		vertices.next();

		while (!vertices.isAtEnd() && ApproxEqual(vertices.get(),*(vertexarr[0])))
		{
			vertices.next();
		}

		if (vertices.isAtEnd())
			return GENPoint::Origin();

		vertexarr[1] = &vertices.get();
		vertices.next();

		vertexarr[2]=NULL;

		for (; !vertices.isAtEnd(); vertices.next())
		{
			const GENPoint* vert2=&vertices.get();

			if (!collinear(*(vertexarr[0]),*(vertexarr[1]),*vert2))
			{
				vertexarr[2]=vert2;
				break;
			}

			// if vertexarr[0] and vertexarr[1] are too close then it is v.
			// difficult to find another vertex which doesn't count as collinear
            // (within the tolerance), so keep track of the two which are further apart
			if (length_squared(*vert2-*vertexarr[0]) > length_squared(*vertexarr[1]-*vertexarr[0]))
				vertexarr[1]=&(*vert2);
		}

		if (NULL==vertexarr[2])
		{
			// couldn't find 3 non-collinear points => zero area polygon
			return GENPoint::Origin();
		}

		const GENPoint result=GENPoint::UnitVector(NormalFromThreePoints(*vertexarr[0],*vertexarr[1],*vertexarr[2]));
		GENAssert(result.Radius()>0.999 && result.Radius()<1.001);

		return result;
	}

	// Calculate signed normal - i.e. the mean normal correctly oriented
	// This needs a triangulation of the polygon
	template <class VERTEX_SET>
	GENPoint SignedNormal(VERTEX_SET vertices)
	{

		if (vertices.isAtEnd())
			return GENPoint::Origin();

        std::shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
        RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
        meshBuilder.AddSurface(vertices.first,vertices.last,0);
        // sets the normal to zero
        GENPoint normal = GENPoint::Origin();
        const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
        for (unsigned int i=0; i<indexedFaceSet->TriangleCount()*3; i+=3) {
            GENPoint a = GENPoint(indexedFaceSet->GetVertex(pointIndices[i]));
            GENPoint b = GENPoint(indexedFaceSet->GetVertex(pointIndices[i+1]));
            GENPoint c = GENPoint(indexedFaceSet->GetVertex(pointIndices[i+2]));
            normal += NormalFromThreePoints(a,b,c);
            // JK - 18.04.2015 - added the volume calculation
            m_volume += (a[GENPointCoords::Z] + b[GENPointCoords::Z] + c[GENPointCoords::Z]) / 3. * TriangleArea(a, b, c, GENPoint::Cartesian(0,0,1));
        }
        // returns the value normalised
        return GENPoint::UnitVector(normal);
	}

	GENPoint m_normal;
	GENPoint m_centroid;
	float m_area;
	float m_volume = 0.f; // JK - 18.04.2015 - added volume
	bool m_verticesAreAnticlockwise;
};

/* A vertex converter class suitable for any vertex classes that define their
 * own conversions to GENPoint
 */
struct DummyConverter
{
	template <class VERTEX>
	DummyConverter(const VERTEX &v) : m_p(v) {}

	operator const GENPoint&() const {
		return m_p;
	}

private:
	DummyConverter& operator=(const DummyConverter& ) { return *this;}
	const GENPoint &m_p;
};

typedef GEOMPolygonInfoImpl<DummyConverter> GEOMPolygonInfo;


/* A vertex converter class suitable for classes where GENPoint is access via a .Data method
 */
struct VertexConverter
{
	template <class VERTEX>
	VertexConverter(const VERTEX &v) : m_p(v.Data()) {}

	operator const GENPoint&() const {
		return m_p;
	}

private:
	VertexConverter& operator=(VertexConverter &) { return *this; }
	const GENPoint &m_p;
};

#endif //_INC_GEOMPOLYGONINFO_INCLUDED
