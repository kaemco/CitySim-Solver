#include "./DATASurfaceIterator.h"
#include "DATASurface.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

using std::vector;

DATASurfaceBuildingIterator::DATASurfaceBuildingIterator(
							const std::vector< DATASurface* >::const_iterator &startpoint,
							const std::vector< DATASurface* >::const_iterator &endpoint) :
		DATASurfaceIterator(startpoint,endpoint)
{
	while (!(m_ptr==m_end) && !(*m_ptr)->IsBuildingSurface()  )
		++m_ptr;
}

// ..........................................................................

DATASurfaceBuildingIterator& DATASurfaceBuildingIterator::operator++() {
	if (m_ptr<m_end)
	{
		++m_ptr;
		while (!(m_ptr==m_end) && !(*m_ptr)->IsBuildingSurface()  )
			++m_ptr;
	}

	return *this;
}

// ..........................................................................
// ..........................................................................

DATASurfaceIterator::DATASurfaceIterator(
								const std::vector< DATASurface* >::const_iterator &startpoint,
								const std::vector< DATASurface* >::const_iterator &endpoint) :
		m_ptr(startpoint),
		m_end(endpoint)
{}

// ..........................................................................

DATASurface& DATASurfaceIterator::operator*() {
	return static_cast<DATASurface&>(**m_ptr);
}

// ..........................................................................

DATASurface* DATASurfaceIterator::operator->() {
	return static_cast<DATASurface*>(*m_ptr);
}

// ..........................................................................

bool DATASurfaceIterator::isAtEnd() const {
	return m_ptr==m_end;
}


// ..........................................................................

DATASurfaceIterator& DATASurfaceIterator::operator++() {
	if (m_ptr!=m_end)
		++m_ptr;

	return *this;
}

// ..........................................................................

bool DATASurfaceIterator::operator==(const DATASurfaceIterator& _Right) const
{	// test for iterator equality
	return (m_ptr == _Right.m_ptr);
}

// ..........................................................................

bool DATASurfaceIterator::operator!=(const DATASurfaceIterator& _Right) const
{	// test for iterator inequality
	return (!(*this == _Right));
}

// ..........................................................................

DATASurfaceIterator DATASurfaceIterator::end() const
{
	return DATASurfaceIterator(m_end,m_end);
}
