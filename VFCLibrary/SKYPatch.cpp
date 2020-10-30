#include "./SKYPatch.h"


#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

SKYPatch::SKYPatch(const GENPoint &patchDir) :
	m_patchDir(patchDir)
{
}

SKYPatch::~SKYPatch(void)
{
}

/**
 * Note that this function calculates the centroid of the potentially visible
 * cells, if all cells are visible the centroid is *NOT* the same as the patch
 * direction specified at construction, this is due to a couple of reasons:
 * - solid angle of cells decrease as altitude increases
 * - even if dw was the same for each cell, the centroid works out a bit lower
 *   than the average angle of altitude of the cells
 *
 */
GENPoint SKYPatch::averageVisibleDirection(const GENPoint& viewDir) const
{
	GENPoint sum=GENPoint::Origin();

	for (std::vector<SKYCell>::const_iterator cell=m_Cells.begin();
		cell!=m_Cells.end();
		++cell)
	{
		const GENPoint& thisCellDir=cell->direction();
		if (dot_product(viewDir, thisCellDir) > 0)
		{
			sum+=thisCellDir*cell->solidAngle();
		}
	}

	return GENPoint::UnitVector(sum);
}

// ..........................................................................

float SKYPatch::formFactor(const GENPoint& viewDir) const
{
	float viewfactor=0;
	for (std::vector<SKYCell>::const_iterator cell=m_Cells.begin();
		cell!=m_Cells.end();
		++cell)
	{
		viewfactor+=cell->formFactor(viewDir);
	}

	return viewfactor;
}

// ..........................................................................

void SKYPatch::addCell(const SKYCell &cell)
{
	m_Cells.push_back(cell);
}

// ..........................................................................

const SKYCell& SKYPatch::getCell(int i) const
{
	return m_Cells[i];
}

// ..........................................................................

unsigned int SKYPatch::cellCount() const
{
	return (unsigned int)m_Cells.size();
}

// ..........................................................................

const GENPoint& SKYPatch::centroid() const
{
	return m_patchDir;
}

// ..........................................................................

float SKYPatch::solidAngle() const
{
	float solidAngle=0.f;
	for (std::vector<SKYCell>::const_iterator cell=m_Cells.begin();
		cell!=m_Cells.end();
		++cell)
	{
		solidAngle+=cell->solidAngle();
	}

	return solidAngle;
}
