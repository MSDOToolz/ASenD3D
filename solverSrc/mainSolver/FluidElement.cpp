#include "FluidElement.h"

using namespace std;

//dup1
DiffDoub0FlPrereq::DiffDoub0FlPrereq() {
	globNds = new DiffDoub0[30];
	globVel = new DiffDoub0[30];
	flDen = new DiffDoub0[10];
	flVel = new DiffDoub0[30];
	flTemp = new DiffDoub0[10];

	return;
}

DiffDoub0FlPrereq::~DiffDoub0FlPrereq() {
	delete[] globNds;
	delete[] globVel;
	delete[] flDen;
	delete[] flVel;
	delete[] flTemp;
	return;
}

//end dup