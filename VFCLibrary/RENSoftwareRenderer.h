#ifndef _INC_RENSOFTWARERENDERER_INCLUDED
#define _INC_RENSOFTWARERENDERER_INCLUDED

// SYSTEM INCLUDES
//
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <memory>

#include "RENIndexedFaceSet.h"

#include "RENVertexProcessor.h"
#include "RENClipper.h"

#include "GENPoint.h"

class RENRenderTarget;
class RENRasteriser;


/*
 * A software-implemented rendering engine.
 */
class RENSoftwareRenderer
{
public:
	typedef RENMatrix Matrix_t;

// LIFECYCLE
	RENSoftwareRenderer();
	virtual ~RENSoftwareRenderer(void);

// OPERATORS

// OPERATIONS
    void SetRasteriser(std::shared_ptr<RENRasteriser> rasteriser);

	/// Set the mesh to be renderered
    void SetModel(std::shared_ptr<RENIndexedFaceSet> W);

	/// Set the position of the camera
	void SetEyePosition(const GENPoint &eye);

	/// Set the point to look at
	void SetViewTarget(const GENPoint &target);

	/// Render view and return pointer to colour buffer
	void Render(RENRenderTarget &target);

	/// Draw a line on the view (points specified in world coordinates, no depth testing)
	//void Draw3DLine(const GENPoint& v1WC, const GENPoint& v2WC,
	//								const unsigned int colour, RENRenderTarget &target);

	/// Draw a line on the view (points specified in screen coordinates, no depth testing)
	//void Draw2DLine(const GENPoint& v1SC, const GENPoint& v2SC,
	//								const unsigned int colour, RENRenderTarget &target);

	/// Set orthogonal mode
	void SetOrthogonal(float xsize, float ysize, float zmax, float zmin);

	/// Set perspective mode
	void SetPerspective(const float fieldOfViewDegrees,
						const float aspectRatio,
						const float znear,
						const float zfar);

	void EnableBackFaceCulling (bool state);

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

	/// Test if a face should be culled
	bool CullFace(const unsigned int faceIndex);

	/// Indicates current matrix mode
	bool m_OrthogonalMode;

	bool m_backFaceCullingEnabled;

	GENPoint m_eyePosition;
	GENPoint m_viewTarget;

    std::shared_ptr<RENIndexedFaceSet> m_model;

	enum {BATCHSIZE=16 };
	RENVertex m_triangleBatch[BATCHSIZE*3];

	RENVertexProcessor m_vertexProcessor;
	RENClipper m_clipper;

    std::shared_ptr<RENRasteriser> m_rasteriser;
};

#endif //_INC_RENSOFTWARERENDERER_INCLUDED
