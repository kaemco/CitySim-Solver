#ifndef _INC_DATASURFACEDELEGATEABC_INCLUDED
#define _INC_DATASURFACEDELEGATEABC_INCLUDED

#include <vector>

#include "GENPoint.h"

/* This class provides an interface through which the view factor calculation
 * can interface with surface objects from the main application.  At some point
 * during the calculation each of these methods are called by the view factor
 * calculation procedures, these methods are implemented by the surface class of
 * the main application.
 */
class DATASurfaceDelegateABC
{
public:
	virtual ~DATASurfaceDelegateABC() {}

	// return the number of vertices
	virtual unsigned int vertexCount() const = 0;

	struct VertexVisitor
	{
		virtual void operator()(const GENPoint&) = 0;
	};

	/* for each vertex in the surface apply the VertexVistor v, i.e. the implementor needs to do something like:
	 * for each vertex vi in polygon:
	 *		v(vi)
	 */
	virtual void sendVertices(VertexVisitor &v) const = 0;

	/* report the surface area
	 */
	virtual float getArea() const = 0;

	/* report the surface normal
	 */
	virtual GENPoint normal() const = 0;

    /* report the surface radius
     */
    virtual float getRadius() const = 0;
};

#endif //_INC_DATASURFACEDELEGATEABC_INCLUDED
