#ifndef _INC_GENWRAPPEDARRAY_INCLUDED
#define _INC_GENWRAPPEDARRAY_INCLUDED

#include <cstdlib>
#include <algorithm>
#include <exception>

/* Just wraps a simple dynamically allocated array, but using this means you don't need to
 * worry about deleting the array - it will be automatically deleted when this class goes out of scope
 */

template <class STORED_TYPE>
class GENWrappedArray
{
public:
	typedef STORED_TYPE& reference;
	typedef const STORED_TYPE& const_reference;

	typedef STORED_TYPE* iterator;
	typedef const STORED_TYPE* const_iterator;

	GENWrappedArray(void) :
        m_arraySize(0),
		m_array(NULL)
	{}

	GENWrappedArray(const GENWrappedArray &A) {
		SafeAlloc(&m_array,A.size());
		std::copy(A.m_array,A.m_array+A.m_arraySize,this->m_array);
	}

	GENWrappedArray<STORED_TYPE>& operator=(const GENWrappedArray<STORED_TYPE> &A) {
		if (A.size()!=this->size())
		{
			this->resize(A.size());
		}

		std::copy(A.m_array,A.m_array+A.m_arraySize,this->m_array);

		return *this;
	}

	GENWrappedArray(size_t size)
	{
		SafeAlloc(&m_array,size);
	}
	~GENWrappedArray(void) {
		delete[] m_array;
	}

	STORED_TYPE* get() {
		return m_array;
	}
	const STORED_TYPE* get() const {
		return m_array;
	}

	void resize(size_t newSize) {
		STORED_TYPE *tmp=m_array;
		SafeAlloc(&m_array,newSize);
		delete[] tmp;
	}

	const_reference operator[](size_t i) const {
		return m_array[i];
	}

	reference operator[](size_t i) {
		return m_array[i];
	}

	size_t size() const {
		return m_arraySize;
	}

private:
	void SafeAlloc(STORED_TYPE** target, size_t size) {
		STORED_TYPE* tmp=NULL;
		try {
			tmp=new STORED_TYPE[size];
			m_arraySize=size;
        } catch (std::exception& e) {
			delete[] tmp;
			throw e;
		}
		*target=tmp;
	}

	size_t m_arraySize;
	STORED_TYPE *m_array;
};

#endif //_INC_GENWRAPPEDARRAY_INCLUDED
