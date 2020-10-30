#ifndef _INC_RENVERTEXPROCESSOR_INCLUDED
#define _INC_RENVERTEXPROCESSOR_INCLUDED

#include <memory>
#include "RENMatrix.h"
#include "RENVertex.h"

class RENIndexedFaceSet;

/*
 * This class performs the transformation of vertices from world coordinates to
 * screen coordinates.
 */
class RENVertexProcessor
{
public:
	RENVertexProcessor();
	~RENVertexProcessor(void);

        void SetVertexBuffer(std::shared_ptr<RENIndexedFaceSet> vb);

	/// Process vertex #vertexId from the vertex buffer
	RENVertex ProcessVertex(int vertexId);

	/// Process the provided vertex
	RENVertex ProcessVertex(const GENPoint& worldCoordinates);

	/// Process the specified vertex using the specified view matrix and current projection matrix
	void Project(const RENVertex& worldCoords, const RENMatrix &viewMatrix,
				RENVertex& perspCoord);

	/// Process the specified vertex using the specified transform matrix
	void ProjectTransform(const RENVertex& worldCoords, const RENMatrix &transformMatrix,
				RENVertex& perspCoord);

	/// Set the projection matrix
	void SetProjectionMatrix(const RENMatrix &projectionMatrix);

	/// Set the view matrix
	void SetViewMatrix(const RENMatrix &viewMatrix);

	/// Get the current transform matrix (view and perspective combined)
	const RENMatrix& GetTransformMatrix() { return m_transform; }
private:
	/// Transform vertex using m_transform
	RENVertex Transform(const RENVertex& inputVertex);

	/// Make sure matrices are current, etc.
	void Update();

	RENMatrix m_projectionMatrix;
	RENMatrix m_viewMatrix;
	RENMatrix m_transform;

	/// Pointer to a vertex buffer (not owned by this class)
    std::shared_ptr<RENIndexedFaceSet> m_vertexBuffer;

	bool m_needsUpdate;
};

#endif //_INC_RENVERTEXPROCESSOR_INCLUDED
