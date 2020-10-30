#ifndef _INC_SKYPATCH_INCLUDED
#define _INC_SKYPATCH_INCLUDED

#include <vector>

#include "./SKYCell.h" // patch is made up of elemental cells
#include "GENPoint.h"

class SKYPatch
{
public:
	SKYPatch(const GENPoint &patchDir);
	virtual ~SKYPatch(void);

	const SKYCell& getCell(int i) const;
	unsigned int cellCount() const;

	/// Return direction of centroid of visible part (if view direction is as specified)
	GENPoint averageVisibleDirection(const GENPoint& viewDir) const;

	float formFactor(const GENPoint& viewDir) const;
	const GENPoint& centroid() const;
    float solidAngle() const;

protected:
	void addCell(const SKYCell &cell);

private:
	float  m_solidAngle;
	GENPoint m_patchDir;

	std::vector<SKYCell> m_Cells;
};

#endif // _INC_SKYPATCH_INCLUDED

