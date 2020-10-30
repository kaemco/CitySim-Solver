#ifndef _INC_GENFIXEDVECTOR_INCLUDED
#define _INC_GENFIXEDVECTOR_INCLUDED

#include <algorithm>


/*
 * This class implements a fixed size vector (whose size is set at
 * compile-time). The interface is very similar to that of the STL vector class,
 * the only difference is that storage is a fixed size and because the size is
 * known at compile time it can be allocated on the stack rather than on the
 * heap (through new or whatever), so that access is quicker.
 *
 */
template <typename Type, size_t dimension>
class GENFixedVector
{
public:
	/// Typedefs for STL compatibility
	typedef Type value_type;
	typedef Type* iterator;
	typedef const Type * const_iterator;
	typedef Type& reference;
	typedef const Type& const_reference;

	/// Default constructor (does nothing)
	GENFixedVector() {};

	/// Constructor that initialises all values to initVal
	GENFixedVector(Type initVal);

	/// Default destructor (does nothing)
	~GENFixedVector() {};

	/// Copy constructor
	GENFixedVector(const GENFixedVector<Type,dimension>& V);

	/// Assignment operator
	GENFixedVector& operator=(const GENFixedVector<Type,dimension>& V);

	/// Returns size of the vector
	inline size_t size() const {return dimension; }

	/// Provides access to individual elements of the vector
	inline reference operator[](size_t indx) { return _myvals[indx]; }

	/// Provides access to individual elements of the vector (const version)
	inline const_reference operator[](size_t indx) const { return _myvals[indx];}

	/// Iterators
	iterator begin() {return _myvals;}
	const_iterator begin() const {return _myvals;}

	iterator end() {return _myvals+dimension;}
	const_iterator end() const {return _myvals+dimension;}

private:
	/// Stores the data
	Type _myvals[dimension];
};

// ..........................................................................

template <typename Type, size_t dimension>
GENFixedVector<Type,dimension>::GENFixedVector(Type initVal)
{
	for (unsigned int i=0; i<dimension; ++i)
	{
		_myvals[i]=initVal;
	}
}

// ..........................................................................

template <typename Type, size_t dimension>
GENFixedVector<Type,dimension>::GENFixedVector(const GENFixedVector<Type,dimension>& V)
{
	std::copy(V._myvals, V._myvals+dimension, _myvals);
}

// ..........................................................................

template <typename Type, size_t dimension>
GENFixedVector<Type,dimension>& GENFixedVector<Type,dimension>::operator=(const GENFixedVector<Type,dimension>& V)
{
	std::copy(V._myvals, V._myvals+dimension, _myvals);
	return *this;
}

#endif //_INC_GENFIXEDVECTOR_INCLUDED
