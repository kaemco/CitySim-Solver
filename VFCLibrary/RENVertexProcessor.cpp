
#include "./RENVertexProcessor.h"



#include "RENIndexedFaceSet.h"
#include "GENPoint.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

using namespace RENMatrixUtils;

RENVertexProcessor::RENVertexProcessor() :
	m_vertexBuffer(NULL)
{
	Identity(m_projectionMatrix);
	Identity(m_viewMatrix);

	m_needsUpdate=false;
}

// ..........................................................................

RENVertexProcessor::~RENVertexProcessor(void)
{
}

// ..........................................................................

void RENVertexProcessor::SetVertexBuffer(std::shared_ptr<RENIndexedFaceSet> vb)
{
	m_vertexBuffer=vb;
}


// ..........................................................................

RENVertex RENVertexProcessor::ProcessVertex(int vertexId)
{
	Update();
	return Transform(m_vertexBuffer->GetVertex(vertexId));
}

// ..........................................................................

RENVertex RENVertexProcessor::ProcessVertex(const GENPoint& worldCoordinates)
{
	Update();
	RENVertex input=RENVertex::CreateCartesian(worldCoordinates[GENPointCoords::X],
														worldCoordinates[GENPointCoords::Y],
														worldCoordinates[GENPointCoords::Z]);

	input.m_colour=0;

	return Transform(input);
}

// ..........................................................................

void RENVertexProcessor::Project(const RENVertex& worldCoords, const RENMatrix &viewMatrix, 
				RENVertex& perspCoord)
{
	perspCoord=m_projectionMatrix * viewMatrix * worldCoords;
}

// ..........................................................................

void RENVertexProcessor::ProjectTransform(const RENVertex& worldCoords, const RENMatrix &transformMatrix, 
				RENVertex& perspCoord)
{
	perspCoord=transformMatrix * worldCoords;
}

// ..........................................................................

void RENVertexProcessor::SetProjectionMatrix(const RENMatrix &projectionMatrix)
{
	m_projectionMatrix=projectionMatrix;
	m_needsUpdate=true;
}

// ..........................................................................

void RENVertexProcessor::SetViewMatrix(const RENMatrix &viewMatrix)
{
	m_viewMatrix=viewMatrix;
	m_needsUpdate=true;
}

// ..........................................................................

RENVertex RENVertexProcessor::Transform(const RENVertex& inputVertex)
{
	return m_transform*inputVertex;
}

// ..........................................................................

void RENVertexProcessor::Update()
{
	if (m_needsUpdate)
	{
		m_transform=m_projectionMatrix;
		m_transform*=m_viewMatrix;
		m_needsUpdate=false;
	}
}



