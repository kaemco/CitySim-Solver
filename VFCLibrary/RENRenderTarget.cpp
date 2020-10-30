#include "./RENRenderTarget.h"
#include <memory.h>
#include "GENAssert.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

// ..........................................................................

RENRenderTarget::RENRenderTarget(void)
{
	m_colourBuffer = 0;
	m_depthBuffer=0;
}

// ..........................................................................

RENRenderTarget::~RENRenderTarget(void)
{
	delete[] m_depthBuffer;
	delete[] m_colourBuffer;
}

// ..........................................................................

void RENRenderTarget::ClearColourBuffer(unsigned int colour)
{
	GENAssert(m_colourBuffer);

	const unsigned int stop=GetWidth() * GetHeight();

	unsigned int* ptColourBuffer=m_colourBuffer;
	for (unsigned int i=0; i<stop; ++i, ++ptColourBuffer)
	{
		*ptColourBuffer=colour;
	}
}

// ..........................................................................

void RENRenderTarget::ClearColourBuffer()
{
	GENAssert(m_colourBuffer);

	const unsigned int stop=GetWidth() * GetHeight();

	memset(m_colourBuffer,0,stop*sizeof(unsigned int));
}

// ..........................................................................

void RENRenderTarget::ClearDepthBuffer(float depth)
{
	GENAssert(m_depthBuffer);

	const unsigned int stop=GetWidth() * GetHeight();
	float* ptDepthBuffer=m_depthBuffer;
	for(unsigned int i = 0; i < stop; ++i, ++ptDepthBuffer)
	{
		*ptDepthBuffer = depth;
	}
}

// ..........................................................................

unsigned int* RENRenderTarget::GetColourBuffer()
{
	return m_colourBuffer;
}

// ..........................................................................

float* RENRenderTarget::GetDepthBuffer()
{
	return m_depthBuffer;
}

// ..........................................................................

int RENRenderTarget::GetStride() const
{
	return GetWidth();
}

// ..........................................................................

unsigned int RENRenderTarget::GetWidth() const
{
	return m_viewport.width;
}

// ..........................................................................

unsigned int RENRenderTarget::GetHeight() const
{
	return m_viewport.height;
}

// ..........................................................................

void RENRenderTarget::SetSize(unsigned int width, unsigned int height)
{
	delete[] m_depthBuffer;
	delete[] m_colourBuffer;

	m_depthBuffer=new float[width*height];
	m_colourBuffer=new unsigned int[width*height];

	m_viewport.width=width;
	m_viewport.height=height;
}

