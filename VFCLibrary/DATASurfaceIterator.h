#ifndef _INC_DATASURFACEBUILDINGITERATOR_INCLUDED
#define _INC_DATASURFACEBUILDINGITERATOR_INCLUDED


#include <vector>

class DATASurface;

/* An iterator that goes through each surface in a radiation scene */
class DATASurfaceIterator
{
public:
	typedef std::forward_iterator_tag iterator_category;
	typedef DATASurface value_type;
	typedef value_type* pointer;
	typedef value_type& reference;
	typedef unsigned int difference_type;

	DATASurfaceIterator(const std::vector< DATASurface* >::const_iterator &startpoint, const std::vector< DATASurface* >::const_iterator &endpoint);

	DATASurface& operator*();
	DATASurface* operator->();
	DATASurfaceIterator& operator++();

	bool operator==(const DATASurfaceIterator& _Right) const;

	bool operator!=(const DATASurfaceIterator& _Right) const;

	DATASurfaceIterator end() const;

	bool isAtEnd() const;
protected:
	DATASurfaceIterator& operator=(DATASurfaceIterator &A);

	std::vector< DATASurface* >::const_iterator m_ptr;
	const std::vector< DATASurface* >::const_iterator m_end;
};

/* An iterator that goes through each building surface in a radiation scene */
class DATASurfaceBuildingIterator : public DATASurfaceIterator
{
public:
	DATASurfaceBuildingIterator(const std::vector< DATASurface* >::const_iterator &startpoint,
									const std::vector< DATASurface* >::const_iterator &endpoint);

	DATASurfaceBuildingIterator& operator++();


private:
	DATASurfaceBuildingIterator& operator=(DATASurfaceBuildingIterator &A);
};

#endif //_INC_DATASURFACEBUILDINGITERATOR_INCLUDED
