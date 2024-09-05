#ifndef FLUIDELEMENT
#define FLUIDELEMENT

#include "DiffDoubClass.h"

//dup1
class DiffDoub0FlPrereq {
public:
	DiffDoub0* globNds;
	DiffDoub0* globVel;
	DiffDoub0* flDen;
	DiffDoub0* flVel;
	DiffDoub0* flTemp;
	DiffDoub0 refVisc;
	DiffDoub0 denVisCoef;
	DiffDoub0 tempVisCoef;
	DiffDoub0 gradVisCoef;
	DiffDoub0 refEnth;
	DiffDoub0 denEnthCoef;
	DiffDoub0 denEnthExp;
	DiffDoub0 denPresCoef;
	DiffDoub0 denPresExp;
	DiffDoub0 refDen;
	DiffDoub0 refTemp;
	DiffDoub0 thermCond;
	DiffDoub0 specHeat;
	DiffDoub0 iGConst;
	DiffDoub0 denMax;
	
	DiffDoub0FlPrereq();

	~DiffDoub0FlPrereq();
};
//end dup

class FluidElement {
private:

public:
};

#endif
