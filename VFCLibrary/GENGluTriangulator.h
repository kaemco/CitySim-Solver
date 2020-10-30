#ifndef _INC_GENGLUTRIANGULATOR_INCLUDED
#define _INC_GENGLUTRIANGULATOR_INCLUDED

/**
 * Class that wraps the GLU triangulator API into a slightly more user friendly form
 *
 */

class GLUtesselator;

class GENGluTriangulatorCallbacks
{
public:

    GENGluTriangulatorCallbacks() {};
    virtual ~GENGluTriangulatorCallbacks() {};

	virtual void BeginPolygon() = 0;

	// AddNextVertex is called for each new vertex with the data passed into
    // GENGluTriangulator::AddVertex
	virtual void AddNextVertex(void *newVert) = 0;

	// CombineData us called whenever a new vertex is created, vertex_data
    // contains the data passed to AddVertex for each of the new interpolated
    // vertices.  This function is supposed to generate vertex data for the new
    // vertex that will then be passed to AddNextVertex
	virtual void CombineData(double coords[3],
		                     void *vertex_data[4],
		                     float weight[4], void **dataOut) = 0;
};

class GENGluTriangulator
{
public:
	GENGluTriangulator();
	~GENGluTriangulator();

	void RegisterCallbacks(GENGluTriangulatorCallbacks *callbacks);

	void BeginPolygon();
	void BeginContour();
	void AddVertex(double x, double y, double z, void *data);

	void EndContour();
	void EndPolygon();

private:
	GLUtesselator *m_triangulator;
	void *m_callbacks;
};


#endif //_INC_GENGLUTRIANGULATOR_INCLUDED
