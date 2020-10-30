#include "./RENSoftwareRenderer.h"
#include "RENRenderTarget.h"
#include "RENRasteriserFlatShading.h"
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

RENSoftwareRenderer::RENSoftwareRenderer(void) :
	m_backFaceCullingEnabled(true),
	m_eyePosition(GENPoint::Origin()),
	m_viewTarget(GENPoint::Origin()),
	m_rasteriser(new RENRasteriserFlatShading)
{
}

// ..........................................................................

RENSoftwareRenderer::~RENSoftwareRenderer(void)
{
}

// ..........................................................................

void RENSoftwareRenderer::SetRasteriser(std::shared_ptr<RENRasteriser> rasteriser)
{
	m_rasteriser=rasteriser;
}

// ..........................................................................
void RENSoftwareRenderer::EnableBackFaceCulling(bool state)
{
	m_backFaceCullingEnabled=state;
}

// ..........................................................................

/**
 * Renders the currently set up view, returns a pointer to the resulting colour buffer
 *
 * \return - pointer to the colour buffer
 */
void RENSoftwareRenderer::Render(RENRenderTarget &target)
{
	UpdateViewMatrix();

	target.ClearColourBuffer();
	target.ClearDepthBuffer();

	const std::vector<unsigned int>& faces=m_model->GetIndices();
	const unsigned int numFaces=static_cast<unsigned int>(faces.size()/3);

	for(unsigned int i = 0; i < numFaces; i += BATCHSIZE)
	{
		for(unsigned int j = 0; j < BATCHSIZE; ++j)
		{
			const unsigned int k = i + j;
			if(k >= numFaces) break;

			// check if face can be culled
			if (!CullFace(k))
			{
				m_triangleBatch[3 * j + 0] = m_vertexProcessor.ProcessVertex(faces[k*3]);
				m_triangleBatch[3 * j + 1] = m_vertexProcessor.ProcessVertex(faces[k*3+1]);
				m_triangleBatch[3 * j + 2] = m_vertexProcessor.ProcessVertex(faces[k*3+2]);

				RENVertex **V = m_clipper.ClipTriangle(m_triangleBatch[3 * j + 0],
													m_triangleBatch[3 * j + 1],
													m_triangleBatch[3 * j + 2]);
				unsigned int n = m_clipper.GetNumVertices();
				if(n < 3)
					continue;

				for(unsigned int k = 0; k < n; ++k)
				{
					ConvertPerspectiveToWindow(*V[k], target.GetWidth(), target.GetHeight());
				}

				for (unsigned int k=2; k<n; ++k)
				{
					// NOTE: We want colour of rendered triangle to be colour of
					// first vertex, of original (unclipped) triangle, i.e. m_triangleBatch[3 * j + 0].m_colour
					m_rasteriser->RenderTriangle(V[0],V[k-1],V[k],m_triangleBatch[3 * j + 0].m_colour, target);
				}
			}
		}
	}
}

// ..........................................................................

/// Draw a line on the view
/*void RENSoftwareRenderer::Draw3DLine(const GENPoint& v1WC,
														 const GENPoint& v2WC,
														 unsigned int colour,
														 RENRenderTarget &target)
{
	UpdateViewMatrix();

	RENVertex persp1=m_vertexProcessor.ProcessVertex(v1WC);
	RENVertex persp2=m_vertexProcessor.ProcessVertex(v2WC);

	RENVertex **V = m_clipper.ClipLine(persp1,persp2);

	unsigned int n = m_clipper.GetNumVertices();

	if(n >= 2)
	{
		GENAssert(n==2);

		for(unsigned int k = 0; k < n; ++k)
		{
			ConvertPerspectiveToWindow(*V[k], target.GetWidth(), target.GetHeight());
		}

		m_rasteriser->DrawLine(V[0],V[1],colour, target);
	}
}*/

// ..........................................................................

/// Draw a line on the view
/*void RENSoftwareRenderer::Draw2DLine(const GENPoint& v1SC,
														 const GENPoint& v2SC,
														 unsigned int colour,
														 RENRenderTarget &target)
{
	RENVertex v1=RENVertex::CreateCartesian(v1SC[GENPointCoords::X],v1SC[GENPointCoords::Y],0,1);
	RENVertex v2=RENVertex::CreateCartesian(v2SC[GENPointCoords::X],v2SC[GENPointCoords::Y],0,1);

	m_rasteriser->DrawLine(&v1,&v2,colour, target);
}
*/
// ..........................................................................

void RENSoftwareRenderer::SetModel(std::shared_ptr<RENIndexedFaceSet> model)
{
    m_model=model;
    m_vertexProcessor.SetVertexBuffer(m_model);
}

// ..........................................................................

void RENSoftwareRenderer::SetEyePosition(const GENPoint &eye)
{
	m_eyePosition=eye;
}

// ..........................................................................

void RENSoftwareRenderer::SetViewTarget(const GENPoint &target)
{
	m_viewTarget=target;
}

// ..........................................................................

void RENSoftwareRenderer::SetOrthogonal(float xsize, float ysize, float zmax, float zmin)
{
	m_vertexProcessor.SetProjectionMatrix(
				REN::BuildOrthogonalProjectionMatrix(xsize,ysize,zmin,zmax)
			);

	m_OrthogonalMode=true;
}

// ..........................................................................

void RENSoftwareRenderer::SetPerspective(const float fieldOfViewDegrees,
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
void RENSoftwareRenderer::Project(const RENVertex& worldCoords,
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
void RENSoftwareRenderer::ProjectTransform(const RENVertex& worldCoords,
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
void RENSoftwareRenderer::ConvertPerspectiveToWindow(RENVertex& vertex, unsigned int windowWidth, unsigned int windowHeight)
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
void RENSoftwareRenderer::UpdateViewMatrix()
{
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
}

// ..........................................................................

bool RENSoftwareRenderer::CullFace(const unsigned int faceIndex)
{
	if (m_backFaceCullingEnabled)
	{
		if (m_OrthogonalMode)
		{
			// cull if surface normal does not point towards view direction
			const GENPoint viewDirection=m_viewTarget-m_eyePosition;
			if (dot_product(viewDirection,m_model->GetNormal(faceIndex))>=0) return true;
		}
		else
		{
			// cull if view point is in negative half space of plane that polygon is in
			if (m_model->EvaluatePlaneEquation(faceIndex,m_eyePosition) <= 1e-5) return true;
		}
	}
	return false;
}

// ..........................................................................

const RENSoftwareRenderer::Matrix_t& RENSoftwareRenderer::GetTransformMatrix()
{
	return m_vertexProcessor.GetTransformMatrix();
}

