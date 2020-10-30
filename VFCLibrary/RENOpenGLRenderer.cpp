// compile this file only if explicitely with OPENGL define
#ifdef OPENGL

#include "./RENOpenGLRenderer.h"
#include "RENRenderTarget.h"
#include "RENUtilities.h"
#include "RENIndexedFaceSetBuilder.h"

#include "GENAssert.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

RENOpenGLRenderer::RENOpenGLRenderer(void) :
	m_backFaceCullingEnabled(true),
	m_eyePosition(GENPoint::Origin()),
	m_viewTarget(GENPoint::Origin()),
	maxSurfaceId(0)
{
    // for information purposes
    std::cerr << "GL_VENDOR: " << glGetString(GL_VENDOR) << "\nGL_RENDERER: " << glGetString(GL_RENDERER) << "\nGL_VERSION: " << glGetString(GL_VERSION) << std::endl;

    // disable lighting
    glDisable( GL_LIGHTING ) ;

    // the clear color and depth
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);

    // create the texture for storing the image
    //glGenFramebuffersEXT(1, &frameBuffer);

}

// ..........................................................................

RENOpenGLRenderer::~RENOpenGLRenderer(void)
{
    glDeleteLists(listObjects, 1);
}

// ..........................................................................

void RENOpenGLRenderer::SetModel(std::auto_ptr<RENIndexedFaceSet> model)
{
	m_model=model;
	m_vertexProcessor.SetVertexBuffer(m_model.get());

    listObjects = glGenLists(1);
    glNewList(listObjects, GL_COMPILE);

    // loads the model in triangles (after triangulation)
	const std::vector<unsigned int>& pointIndices=m_model->GetIndices();

    // loop on triangles to find the maxSurfaceId
    for (unsigned int i=0; i<pointIndices.size()/3; i++) {

        if (maxSurfaceId < m_model->GetVertex(pointIndices[i*3]).m_colour) maxSurfaceId = m_model->GetVertex(pointIndices[i*3]).m_colour;
    }

    // loop on triangles
    for (unsigned int i=0; i<pointIndices.size()/3; i++) {

        unsigned int firstVertexColour = m_model->GetVertex(pointIndices[i*3]).m_colour;
        std::cerr << "Triangle(" << i << ") colour: " << firstVertexColour << "\tfloat: " << (float)firstVertexColour << std::endl;
        glColor3f(static_cast<float>(firstVertexColour)/static_cast<float>(maxSurfaceId), 0.f, 0.f);
        GENPoint vertex;
        glBegin(GL_TRIANGLES);
        vertex = m_model->GetVertex(pointIndices[i*3]);
        glVertex3f(vertex[0], vertex[1], vertex[2]);
        vertex = m_model->GetVertex(pointIndices[i*3+1]);
        glVertex3f(vertex[0], vertex[1], vertex[2]);
        vertex = m_model->GetVertex(pointIndices[i*3+2]);
        glVertex3f(vertex[0], vertex[1], vertex[2]);
        glEnd();

    }

    glEndList();

}

/**
 * Renders the currently set up view, returns a pointer to the resulting colour buffer
 *
 * \return - pointer to the colour buffer
 */
void RENOpenGLRenderer::Render(RENRenderTarget &target)
{
    // clear the target storing the colour and depth
	target.ClearColourBuffer();
	target.ClearDepthBuffer();
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // makes a test to know which face comes first from a view point, but does not store the depth
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); // stores the pixel if the depth is less from the viewpoint
    glDepthMask(GL_TRUE); // the depth buffer is enabled

    // flat normal mode
    glShadeModel(GL_FLAT);

    // back face culling
    if (m_backFaceCullingEnabled)   { glEnable(GL_CULL_FACE); glCullFace(GL_BACK); }
    else                            glDisable(GL_CULL_FACE);

	// Select rasterization mode
	glPolygonMode(GL_FRONT, GL_FILL);

	// Scale viewport according to image resolution
	glViewport( 0, 0, target.GetWidth(), target.GetHeight() );
	std::cerr << "Viewport width: " << target.GetWidth() << " height: " << target.GetHeight() << std::endl;

    // creates the view matrix
	UpdateViewMatrix();

    // calls the list of objects
    glCallList(listObjects);

    // flushes the view and waits for its finish
    glFinish();

    // read the output
	float *pixels = new float[target.GetWidth()*target.GetHeight()];
	glReadPixels( 0, 0, target.GetWidth(), target.GetHeight(), GL_RED, GL_FLOAT, pixels);
    target.SetPixels(pixels, maxSurfaceId);
}

// ..........................................................................

void RENOpenGLRenderer::SetOrthogonal(float xsize, float ysize, float zmax, float zmin)
{
	m_OrthogonalMode=true;

    m_xsize = xsize;
    m_ysize = ysize;
    m_zmax  = zmax;
    m_zmin  = zmin;
}

// ..........................................................................

void RENOpenGLRenderer::SetPerspective(const float fieldOfViewDegrees,
												const float aspectRatio,
												const float znear,
												const float zfar)
{
	m_vertexProcessor.SetProjectionMatrix(
				REN::BuildPerspectiveProjectionMatrix(fieldOfViewDegrees*(float)M_PI/180.f,aspectRatio,znear,zfar)
			);

	m_OrthogonalMode=false;
}

// ..........................................................................

/**
 * This method projects a point in world coordinates into screen coordinates
 * using the provided view transform matrix + window size and the current projection matrix
 *
 * \param worldCoords - the point to be projected
 * \param &viewMatrix - the view matrix
 * \param screenCoord - the results
 */
void RENOpenGLRenderer::Project(const RENVertex& worldCoords,
									   const Matrix_t &viewMatrix,
										unsigned int windowWidth,
										unsigned int windowHeight,
									   RENVertex& screenCoord)
{
	// project the point to homogenous clip co-ords
	m_vertexProcessor.Project(worldCoords, viewMatrix, screenCoord);
	ConvertPerspectiveToWindow(screenCoord, windowWidth,windowHeight);
}

// ..........................................................................

/**
 * This method projects a point in world coordinates into screen coordinates
 * using the provided transform matrix (which includes the projection matrix) and window size
 *
 * \param worldCoords - the point to be projected
 * \param &transformMatrix - the view matrix
 * \param screenCoord - the results
 */
void RENOpenGLRenderer::ProjectTransform(const RENVertex& worldCoords,
												const Matrix_t &transformMatrix,
												unsigned int windowWidth,
												unsigned int windowHeight,
												RENVertex& screenCoord)
{
	// project the point to homogenous clip co-ords
	m_vertexProcessor.ProjectTransform(worldCoords, transformMatrix, screenCoord);
	ConvertPerspectiveToWindow(screenCoord,windowWidth,windowHeight);
}

// ..........................................................................

/* Take a point in homogenous canonical perspective space and convert it to pixel (window) coordinates */
void RENOpenGLRenderer::ConvertPerspectiveToWindow(RENVertex& vertex, unsigned int windowWidth, unsigned int windowHeight)
{
	// after perspective divide, points within the viewport map to -1<x<1, -1<y<1
	const unsigned int W = windowWidth/2;
	const unsigned int H = windowHeight/2;
	const float RHW = 1.0f / vertex.m_position[3];

	// offset pixel positions by 0.5 pixels - we want to convert
    // from a system where the bottom left of the screen is position
    // 0,0 and top right is 128,128 to where the centre of the bottom
    // left pixel is 0,0 and centre of the top left pixel is 127,127
	vertex.m_position[0] = W * (vertex.m_position[0]*RHW+1) - 0.5f;
	vertex.m_position[1] = H * (1-vertex.m_position[1]*RHW) - 0.5f;
	vertex.m_position[2] = vertex.m_position[2] * RHW;
	vertex.m_position[3]= 1;
}

// ..........................................................................

/**
 *
 * \return
 */
void RENOpenGLRenderer::UpdateViewMatrix()
{
    // here we save the transformation matrix (for VFCFishEye)
	GENPoint up=GENPoint::Cartesian(0,0,1);
	//if (m_OrthogonalMode)
	{
		GENPoint z=GENPoint::UnitVector(m_eyePosition-m_viewTarget);
		if (dot_product(z,up)>0.999)
		{
			up=GENPoint::Cartesian(0,1,0);
		}
	}

	m_vertexProcessor.SetViewMatrix(REN::BuildLookAtMatrix(m_eyePosition,m_viewTarget,up));

    // using OpenGL for projection, this sets up the view

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();

	// Applying orthographic projection is the same as dividing by bounding sphere radius
 	glOrtho(-m_xsize/2, m_xsize/2, -m_ysize/2, m_ysize/2, -(m_zmax-m_zmin), (m_zmax-m_zmin)); // TODO: verifiy the matrix from gluLookAt and save matrix

    // Create the lookat matrix for projection
    gluLookAt(m_eyePosition[0], m_eyePosition[1], m_eyePosition[2], m_viewTarget[0], m_viewTarget[1], m_viewTarget[2], up[0], up[1], up[2]);

}

// ..........................................................................

const RENOpenGLRenderer::Matrix_t& RENOpenGLRenderer::GetTransformMatrix()
{
	return m_vertexProcessor.GetTransformMatrix();
}

#endif
