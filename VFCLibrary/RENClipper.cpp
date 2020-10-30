#include "./RENClipper.h"

#ifdef _MSC_VER
#if defined(DEBUG_MEM)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__, __LINE__)
#endif
#endif // _MSC_VER

RENVertex interpolate(const RENVertex &Vi, const RENVertex &Vj, float d)
{
	RENVertex W;

	W.m_position[0] = Vi.m_position[0] + d * (Vj.m_position[0] - Vi.m_position[0]);
	W.m_position[1] = Vi.m_position[1] + d * (Vj.m_position[1] - Vi.m_position[1]);
	W.m_position[2] = Vi.m_position[2] + d * (Vj.m_position[2] - Vi.m_position[2]);
	W.m_position[3] = Vi.m_position[3] + d * (Vj.m_position[3] - Vi.m_position[3]);

	W.m_colour=Vi.m_colour;

	return W;
}

// ..........................................................................

RENClipper::~RENClipper()
{
	clearStore();
}

// ..........................................................................

/**
 *
 * \param &V1 
 * \param &V2 
 * \param &V3 
 * \return 
 * \todo - change return type to a std::vector
 */
RENVertex** RENClipper::ClipTriangle(RENVertex &V1, 
											   RENVertex &V2, 
											   RENVertex &V3)
{
	clearStore();
	m_P[0].clear();

	m_P[0].push_back(&V1);
	m_P[0].push_back(&V2);
	m_P[0].push_back(&V3);

	m_stage = 0;
	m_numVertices = 3;

	clipAll();
	if (0!=GetNumVertices())
		return &m_P[m_stage][0];
	else
		return NULL;
}

// ..........................................................................

RENVertex** RENClipper::ClipLine(RENVertex &V1, 
										 RENVertex &V2)
{
	clearStore();
	m_P[0].clear();

	m_P[0].push_back(&V1);
	m_P[0].push_back(&V2);

	m_stage = 0;
	m_numVertices = 2;

	clipAll();
	return &m_P[m_stage][0];
}

// ..........................................................................

void RENClipper::clipAll()
{
	clip(&RENClipper::dNear);
	clip(&RENClipper::dFar);
	clip(&RENClipper::dLeft);
	clip(&RENClipper::dRight);
	clip(&RENClipper::dTop);
	clip(&RENClipper::dBottom);
}

// ..........................................................................

int RENClipper::GetNumVertices() const
{
	return m_numVertices;
}

// ..........................................................................

void RENClipper::clip(distance_func distanceFromClipEdge)
{
	std::vector<RENVertex*> *preClipVertices = &m_P[m_stage];
	std::vector<RENVertex*> *postClipVertices = &m_P[m_stage + 1];
	postClipVertices->clear();

	int t = 0;

	for(int current = 0; current < m_numVertices; current++)
	{
		const int next = (current+1)%m_numVertices;

		const float di = (this->*distanceFromClipEdge)((*preClipVertices)[current]);
		const float dj = (this->*distanceFromClipEdge)((*preClipVertices)[next]);

		if(di >= 0) 
		{
			postClipVertices->push_back((*preClipVertices)[current]);
			t++;
			if(dj < 0)
			{
				RENVertex* newVertex=new RENVertex(interpolate(*(*preClipVertices)[current], *(*preClipVertices)[next], di / (di - dj)));
				addToStore(newVertex);
				postClipVertices->push_back(newVertex);
				t++;
			}
		}
		else
		{
			if(dj > 0)
			{
				RENVertex* newVertex=new RENVertex(interpolate(*(*preClipVertices)[next], *(*preClipVertices)[current], dj / (dj - di)));
				addToStore(newVertex);
				postClipVertices->push_back(newVertex);
				t++;
			}
		}
	}

	m_numVertices = t;
	m_stage += 1;
}

// ..........................................................................

float RENClipper::dNear(const RENVertex *V) const
{
	return V->m_position[2];
}

// ..........................................................................

float RENClipper::dFar(const RENVertex *V) const
{
	return V->m_position[3] - V->m_position[2];
}

// ..........................................................................

float RENClipper::dLeft(const RENVertex *V) const
{
	return V->m_position[3] + V->m_position[0];
}

// ..........................................................................

float RENClipper::dRight(const RENVertex *V) const
{
	return V->m_position[3] - V->m_position[0];
}

// ..........................................................................

float RENClipper::dTop(const RENVertex *V) const
{
	return V->m_position[3] - V->m_position[1];
}

// ..........................................................................

float RENClipper::dBottom(const RENVertex *V) const
{
	return V->m_position[3] + V->m_position[1];
}

// ..........................................................................

void RENClipper::clearStore()
{
	for (std::vector<RENVertex*>::iterator it=m_vertexStore.begin();
		it!=m_vertexStore.end();
		++it)
	{
		delete *it;
	}
	m_vertexStore.clear();
}

// ..........................................................................

void RENClipper::addToStore(RENVertex* vertex)
{
	m_vertexStore.push_back(vertex);
}

