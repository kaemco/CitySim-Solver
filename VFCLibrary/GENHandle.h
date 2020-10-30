#ifndef _INC_GENHANDLE_INCLUDED
#define _INC_GENHANDLE_INCLUDED


#include "GENAssert.h"

#ifndef NULL
#define NULL 0
#endif

/* GENHandle is basically a class to make memory management a little easier -
 * you give the GENHandle a pointer to an object, you can copy GENHandles to
 * share that pointer in various places of the program, then when the last
 * handle is deleted the object that it points to is also deleted
 *
 */
namespace GENHandleUtils {
	class SimpleCounter
	{
	public:
		SimpleCounter(unsigned int count) :
			m_count(new unsigned int(count))
		{}

		virtual ~SimpleCounter()
		{
			delete m_count;
		}

		void Decrement() {
            GENAssert(*m_count>0); --*m_count;
		}
		void Increment() {
			++*m_count;
		}

		unsigned int Count() { return *m_count; }
	private:
		SimpleCounter(const SimpleCounter &X):
			m_count(X.m_count)
		{
			throw "Shouldn't be called";
		}

		unsigned int *m_count;
	};


	struct DeleteViaDestructor
	{
		template <class T>
		void DeleteObject(T* obj)
		{
			delete obj;
		}

	protected:
		~DeleteViaDestructor() {};
	};

	// ..........................................................................

	//template <class T>
	//struct DeleteViaMemberFunction
	//{
	//public:
	//	void DeleteObject(T* obj)
	//	{
	//		obj->Destruct();
	//	}

	//protected:
	//	~DeleteViaMemberFunction() {};

	//};

	// ..........................................................................

	template <class T, void(*DELETIONFUNCTION)(T*)>
	struct DeleteViaUnaryFunction
	{
	public:
		void DeleteObject(T* obj)
		{
			DELETIONFUNCTION(obj);
		}

	protected:
		~DeleteViaUnaryFunction() {};

	};

	// ..........................................................................

	//template <class T>
	//struct DontDelete
	//{
	//public:
	//	void DeleteObject(T* obj)
	//	{
	//	}

	//protected:
	//	~DontDelete() {};

	//};
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY=GENHandleUtils::DeleteViaDestructor>
class GENHandle : protected DELETIONPOLICY
{
    // commented by JK 30/01/09
	//friend class GENHandle;
public:
	GENHandle() : m_target(NULL), m_referenceCounter(NULL) {}
	explicit GENHandle(TARGET* obj);
	GENHandle(const GENHandle<TARGET,DELETIONPOLICY> &handle);

	virtual ~GENHandle();

	GENHandle<TARGET,DELETIONPOLICY>& operator=(const GENHandle<TARGET,DELETIONPOLICY>& handle);

	template <class PT> operator GENHandle<PT>();

	TARGET* get();
	TARGET* operator->();
	TARGET& operator*();
	const TARGET* operator->() const;
	const TARGET& operator*() const ;

	bool IsNull() const;

	explicit GENHandle(TARGET* obj, GENHandleUtils::SimpleCounter *referenceCounter);

private:

	void DecrementRefCount();

	TARGET *m_target;
	GENHandleUtils::SimpleCounter *m_referenceCounter;
};


// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
GENHandle<TARGET,DELETIONPOLICY>::GENHandle(TARGET* obj) :
	m_target(obj)
{
	if (!IsNull())
	{
		m_referenceCounter=new GENHandleUtils::SimpleCounter(1);
	}
	else
	{
		m_referenceCounter=NULL;
	}
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
GENHandle<TARGET,DELETIONPOLICY>::GENHandle(const GENHandle &handle) :
		DELETIONPOLICY(handle),
		m_target(handle.m_target),
		m_referenceCounter(handle.m_referenceCounter)

{
	if (!IsNull())
		m_referenceCounter->Increment();
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
GENHandle<TARGET,DELETIONPOLICY>& GENHandle<TARGET,DELETIONPOLICY>::operator=(const GENHandle& handle)
{
	if (m_target==handle.m_target) return *this;

	DecrementRefCount();

	m_target=handle.m_target;
	m_referenceCounter=handle.m_referenceCounter;

	if (!IsNull())
		m_referenceCounter->Increment();
	return *this;

}

// ..........................................................................

/// \bug - broken! - no checking that DELETIONPOLICYs are compatible
template <class TARGET, class DELETIONPOLICY>
template <class PT>
GENHandle<TARGET,DELETIONPOLICY>::operator GENHandle<PT>()
{
	return GENHandle<PT>(m_target, m_referenceCounter);
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
GENHandle<TARGET,DELETIONPOLICY>::~GENHandle()
{
	DecrementRefCount();
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
TARGET* GENHandle<TARGET,DELETIONPOLICY>::operator->()
{
	return m_target;
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
TARGET& GENHandle<TARGET,DELETIONPOLICY>::operator*()
{
	return *m_target;
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
const TARGET* GENHandle<TARGET,DELETIONPOLICY>::operator->() const
{
	return m_target;
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
TARGET* GENHandle<TARGET,DELETIONPOLICY>::get()
{
	return m_target;
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
const TARGET& GENHandle<TARGET,DELETIONPOLICY>::operator*() const
{
	return *m_target;
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
bool GENHandle<TARGET,DELETIONPOLICY>::IsNull() const
{
	return NULL==m_target;
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
void GENHandle<TARGET,DELETIONPOLICY>::DecrementRefCount()
{
	if (!IsNull())
	{
		m_referenceCounter->Decrement();
		if (0==m_referenceCounter->Count())
		{
			delete m_referenceCounter;
			m_referenceCounter=NULL;

			DELETIONPOLICY::DeleteObject(m_target);
			m_target=NULL;
		}
	}
}

// ..........................................................................

template <class TARGET, class DELETIONPOLICY>
GENHandle<TARGET,DELETIONPOLICY>::GENHandle(TARGET* obj, GENHandleUtils::SimpleCounter *referenceCounter) :
	m_target(obj),
	m_referenceCounter(referenceCounter)
{
	m_referenceCounter->Increment();
}


#endif //_INC_GENHANDLE_INCLUDED
