#ifndef _INC_SKYVAULT_INCLUDED
#define _INC_SKYVAULT_INCLUDED

/*
 *  Defines a generic sky vault made up of patches.  The theory is that by using
 *  this base class the view factor calculation is independent of the sky
 *  discretisation scheme, so you could easily implement a different
 *  discretisation.
 */

#include "./SKYPatch.h" // sky is made up of sky patches

// JK, unused class commented 24.02.2009
class SKYVault
{
public:
	virtual ~SKYVault(void) {}

	// return all of the patches that make up the sky vault
	// NOTE: used to have a method to return only the patches that are visible to surface, but the overhead
	//		 in calculating these, plus difficulties in outputing results meant that it wasn't worth it
	virtual const SKYPatch* GetPatch(unsigned int i) const = 0;
 	virtual unsigned int PatchCount() const = 0;

protected:

};

#endif // _INC_SKYVAULT_INCLUDED
