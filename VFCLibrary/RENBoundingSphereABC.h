#ifndef _INC_GEOMSPHEREABC_INCLUDED
#define _INC_GEOMSPHEREABC_INCLUDED

// SYSTEM INCLUDES
// 

// PROJECT INCLUDES
//


// LOCAL INCLUDES
//
//class GENPoint;

// FORWARD REFERENCES
//

class RENBoundingSphereABC
{
public:
	virtual ~RENBoundingSphereABC() {}

	virtual GENPoint Centre() const = 0;
	virtual float Radius() const = 0;

private:

};

#endif // _INC_GEOMSPHEREABC_INCLUDED
