#ifndef _INC_RENHARDWARERENDERER_INCLUDED
#define _INC_RENHARDWARERENDERER_INCLUDED

// SYSTEM INCLUDES
//
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <memory>

#include "RENIndexedFaceSet.h"

#include "RENVertexProcessor.h"

#include "GENPoint.h"

#ifdef _MSC_VER
 #define NOMINMAX
 #include <windows.h>
 #include <GL/gl.h>
 #include <GL/glu.h>
#elif __APPLE_CC_
 #include <OpenGL/gl.h>
 #include <OpenGL/glu.h>
 #include <OpenGL/glext.h>
#else
 #include <GL/gl.h>
 #include <GL/glu.h>
 #define GL_GLEXT_PROTOTYPES    1
 #include <GL/glext.h>
#endif

class RENRenderTarget;


/*
 * An OpenGL rendering engine
 */
class RENOpenGLRenderer
{
public:
	typedef RENMatrix Matrix_t;

// LIFECYCLE
	RENOpenGLRenderer();
	virtual ~RENOpenGLRenderer(void);

// OPERATORS

// OPERATIONS

	/// Set the mesh to be renderered
	void SetModel(std::auto_ptr<RENIndexedFaceSet> W);

	/// Set the position of the camera
	void SetEyePosition(const GENPoint &eye) { m_eyePosition=eye; }

	/// Set the point to look at
	void SetViewTarget(const GENPoint &target) { m_viewTarget=target; }

	/// Render view and return pointer to colour buffer
	void Render(RENRenderTarget &target);

	/// Set orthogonal mode
	void SetOrthogonal(float xsize, float ysize, float zmax, float zmin);

	/// Set perspective mode
	void SetPerspective(const float fieldOfViewDegrees,
						const float aspectRatio,
						const float znear,
						const float zfar);

	void EnableBackFaceCulling(bool state) { m_backFaceCullingEnabled=state; }

	/// Project the given world coordinates into screen space using the provided view matrix
	void Project(const RENVertex& worldCoords,
				const Matrix_t &viewMatrix,
				unsigned int windowWidth,
				unsigned int windowHeight,
				RENVertex& screenCoord);

	/// Project the given world coordinates into screen space using the provided transformation matrix
	void ProjectTransform(const RENVertex& worldCoords,
						const Matrix_t &transformMatrix,
						unsigned int windowWidth,
						unsigned int windowHeight,
						RENVertex& screenCoord);

// ACCESS
	/// Get the transform matrix
	const Matrix_t& GetTransformMatrix();

// INQUIRY

private:
	void ConvertPerspectiveToWindow(RENVertex& vertex, unsigned int windowWidth, unsigned int windowHeight);

	/// Calculate the view matrix
	void UpdateViewMatrix();

    // variables for the view
    float m_xsize, m_ysize, m_zmax, m_zmin;

	bool m_OrthogonalMode;
	bool m_backFaceCullingEnabled;

	GENPoint m_eyePosition;
	GENPoint m_viewTarget;

	std::auto_ptr<RENIndexedFaceSet> m_model;

	RENVertexProcessor m_vertexProcessor;

    GLuint listObjects;
    unsigned int maxSurfaceId;

    GLuint frameBuffer;

};

#endif //_INC_RENSOFTWARERENDERER_INCLUDED
