#ifndef _INC_GENASSERT_INCLUDED
#define _INC_GENASSERT_INCLUDED

// Define this constant if optional debug test blocks should be executed in the code
#define GENAssert_PERFORM_OPTIONAL_DEBUG_TESTS

// Define to disable assertions (e.g. for release builds)
//#define GEN_DISABLE_ASSERTIONS

typedef void(AssertFailureHandler)(const char* condition, const char* file, unsigned int lineNumber);
void RegisterAssertFailureHandler(AssertFailureHandler f);

#ifdef GEN_DISABLE_ASSERTIONS

#define GENAssert(exp) ((void)0)

#else

void GENHandleFailedAssert(const char* condition, const char* file, unsigned int lineNumber);
#define GENAssert(exp) (void)((exp) || (GENHandleFailedAssert(#exp, __FILE__, __LINE__),0) )

#endif

// commented by JK - 30.01.09
#include <stdexcept>
//class _enforce_failure : public std::exception
//{
//public:
//	explicit _enforce_failure(const std::string& msg)
//		: m_msg(msg)
//	{}
//
//	virtual ~_enforce_failure()
//	{}
//
//	virtual const char *what() const
//	{
//		return (m_msg.c_str());
//	}
//
//private:
//	std::string m_msg;
//};


// commented by JK, original code
/*#define GENEnforce(condition,message)		\
	(void)((condition) || (throw _enforce_failure(message),0))*/

#define GENEnforce(condition,message)		\
	(void)((condition) || (throw std::string(std::string("Enforce failure: ") + std::string(message)),0))

#endif //_INC_GENASSERT_INCLUDED
