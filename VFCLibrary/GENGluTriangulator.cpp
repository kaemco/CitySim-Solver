#include "./GENGluTriangulator.h"

// VC_EXTRALEAN removes in windows.h the extra features such as sound, etc.
#define VC_EXTRALEAN

#ifdef _MSC_VER
 #define NOMINMAX
 #include <windows.h>
#endif

#ifndef CALLBACK
# if defined(__MINGW32__)
#  define CALLBACK __attribute__((__stdcall__))
# elif defined(__CYGWIN__)
#  define CALLBACK APIENTRY
# elif defined(__APPLE_CC__)
#  define CALLBACK
# else
#  define CALLBACK APIENTRY
# endif
#endif

#ifdef __APPLE_CC__
 #include <OpenGL/glu.h>
#else
 #include <GL/glu.h>
#endif

#include <string>
#include <stdexcept>

namespace GENGluTriangulatorUtils
{

	void CALLBACK errorCallBack(GLenum /*errorCode*/)
	{
		throw std::runtime_error("Tesselation error!");
	}

	void CALLBACK combineDataCallback(double coords[3],
						void *vertex_data[4],
						float weight[4], void **dataOut, GENGluTriangulatorCallbacks * callbacks )
	{
		callbacks->CombineData(coords,vertex_data,weight,dataOut);
	}

	void CALLBACK beginCallback(GLenum /*type*/, GENGluTriangulatorCallbacks * callbacks)
	{
		callbacks->BeginPolygon();
	}

	void CALLBACK vertexCallback(void * thisVertex, GENGluTriangulatorCallbacks * callbacks)
	{
		callbacks->AddNextVertex(thisVertex);
	}

	void CALLBACK edgeFlagCallback (GLboolean /*flag*/)
	{
	}
};

using namespace GENGluTriangulatorUtils;

GENGluTriangulator::GENGluTriangulator() :
	m_triangulator(gluNewTess())
{
	// set up winding rule
	gluTessProperty(m_triangulator, GLU_TESS_WINDING_RULE,
					GLU_TESS_WINDING_POSITIVE);

}

GENGluTriangulator::~GENGluTriangulator()
{
	gluDeleteTess(m_triangulator);
}

void GENGluTriangulator::BeginPolygon() {
	gluTessBeginPolygon(m_triangulator, m_callbacks);
}
void GENGluTriangulator::BeginContour() {
	gluTessBeginContour(m_triangulator);
}

void GENGluTriangulator::AddVertex(double x, double y, double z, void *data) {
	double coords[3] = { x,y,z};
	gluTessVertex(m_triangulator,coords,data);
}

void GENGluTriangulator::EndContour() {
	gluTessEndContour(m_triangulator);
}
void GENGluTriangulator::EndPolygon() {
	gluTessEndPolygon(m_triangulator);
}
void GENGluTriangulator::RegisterCallbacks(GENGluTriangulatorCallbacks *callbacks)
{
	m_callbacks=callbacks;
	gluTessCallback(m_triangulator, GLU_TESS_VERTEX_DATA,
                (void (CALLBACK*) ()) &vertexCallback);
	gluTessCallback(m_triangulator, GLU_TESS_BEGIN_DATA,
                (void (CALLBACK*) ()) &beginCallback);
	gluTessCallback(m_triangulator, GLU_TESS_ERROR,
                (void (CALLBACK*) ()) &errorCallBack);
	gluTessCallback(m_triangulator, GLU_TESS_COMBINE_DATA,
                (void (CALLBACK*) ()) &combineDataCallback);

	// registering this callback forces the triangulator to generate only triangles (not triangle strips or fans)
	gluTessCallback(m_triangulator, GLU_TESS_EDGE_FLAG,
                (void (CALLBACK*) ()) &edgeFlagCallback);

}
