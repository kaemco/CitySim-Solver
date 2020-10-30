#include "./GENAssert.h"

#include <sstream>
#include <fstream>
#include <cassert>

// ..........................................................................

/* If GENAssert_SHOULD_THROW is defined then whenever assert(XYZ) fails an
 * exception will be thrown (in both debug and release builds).  This allows for
 * a kind of release build with debugging checks still in place
 */

#define GENAssert_SHOULD_THROW

// Define this constant if the assert should write the error message to a file
#define GENAssert_SHOULD_DUMP_MESSAGE_TO_FILE
#define GENAssert_MESSAGE_FILE "c:/radiation_error.txt"

void GENHandleFailedAssertImpl(const char* condition, const char* file, unsigned int lineNumber)
{
	std::ostringstream message;
	message << "Assertion failed:\n";
	message << std::string(condition) << " in " << std::string(file) << " line number " << lineNumber;

#ifdef GENAssert_SHOULD_DUMP_MESSAGE_TO_FILE
	std::ofstream errorFile(GENAssert_MESSAGE_FILE);
	errorFile << message.str();
	errorFile.close();
#endif

#ifdef GENAssert_SHOULD_THROW
	throw std::runtime_error(message.str());
#else
	_assert(condition,file,lineNumber);
#endif

}

// ..........................................................................

namespace GEN
{
	static AssertFailureHandler *failedAssertHandler=&GENHandleFailedAssertImpl;

};

void RegisterAssertFailureHandler(AssertFailureHandler f)
{

	GEN::failedAssertHandler=NULL!=f ? f : &GENHandleFailedAssertImpl;
}

void GENHandleFailedAssert(const char* condition, const char* file, unsigned int lineNumber)
{
	(*GEN::failedAssertHandler)(condition,file,lineNumber);
}
