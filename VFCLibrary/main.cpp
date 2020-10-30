#include "./VFCViewFactorCalculation.h"
#include "./DATARadiationScene.h"
#include "./DATASurfaceDelegateABC.h"
#include "./DATASurfaceIterator.h"
#include "./DATAViewFactorSetSparse.h"
#include "./DATASurface.h"
#include "./DATAInsolationFactors.h"
#include "./SKYSun.h"

using namespace std;

class MySurfaceDelegate : public DATASurfaceDelegateABC
{
public:
	virtual ~MySurfaceDelegate() {}

	virtual unsigned int vertexCount() const
	{
		return vertices.size();
	}

	virtual void sendVertices(VertexVisitor &v) const
	{
		for (std::vector<GENPoint>::const_iterator it=vertices.begin();
			it!=vertices.end();
			++it)
			v(*it);
	}

	virtual float getArea() const
	{
		return area;
	}

	virtual GENPoint normal() const
	{
		return norm;
	}

	double area;
	GENPoint norm;
	std::vector<GENPoint> vertices;
};



#ifdef TEST_VFCLIBRARY

int main(void)
{
	GENHandle<MySurfaceDelegate> s1(new MySurfaceDelegate);
	s1->vertices.push_back(GENPoint::Cartesian(0,0,0));
	s1->vertices.push_back(GENPoint::Cartesian(1,0,0));
	s1->vertices.push_back(GENPoint::Cartesian(1,1,0));
	s1->vertices.push_back(GENPoint::Cartesian(0,1,0));
	s1->norm=GENPoint::Cartesian(0,0,1);
	s1->area=1;

	GENHandle<MySurfaceDelegate> s2(new MySurfaceDelegate);
	s2->vertices.push_back(GENPoint::Cartesian(0,1,0));
	s2->vertices.push_back(GENPoint::Cartesian(0,1,1));
	s2->vertices.push_back(GENPoint::Cartesian(1,1,1));
	s2->vertices.push_back(GENPoint::Cartesian(1,1,0));
	s2->norm=GENPoint::Cartesian(0,-1,0);
	s2->area=1;

	DATARadiationScene scene;
	scene.AddBuildingSurface(s1);
	scene.AddBuildingSurface(s2);

	VFCViewFactorCalculation v;
	v.CalculateViewFactors(scene);

	for (DATASurfaceIterator it=scene.GetAllSurfaces();
		!it.isAtEnd();
		++it)
	{
		for (DATAViewFactorSetSparse::const_iterator factors=it->SWViewFactors().GetVFs();
			factors!=it->SWViewFactors().GetLastVF();
			++factors)
		{
			std::cout << "Unobstructed factor to patch " << factors->patchNo << ": " << factors->unobstructed << "\n";
		}
		SKYSun sun;
		sun.SetDay(30);
		sun.SetClockTime(10.5);
		std::cout << "Insolation factor: " << it->InsolationFactors().GetInsolationFactor(sun) << "\n";
	}

	return 0;
}
#endif
