#include "DesignVariableClass.h"
#include <string>
#include <vector>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "SetClass.h"
using namespace std;

const int max_int = 2000000000;

DesignVariable::DesignVariable() {
	category = "";
	component = 1;
	layer = 0;
	elSetName = "";
	elSetPtr = max_int;
	ndSetName = "";
	ndSetPtr = max_int;
	coefs.clear();
	activeTime[0] = 0.0;
	activeTime[1] = 1.0e+100;
	value.setVal(0.0);
	diffVal.setVal(0.0);
	compElList.clear();
}

void DesignVariable::setActiveTime(double newAt[]) {
	activeTime[0] = newAt[0];
	activeTime[1] = newAt[1];
}

void DesignVariable::getValue(DiffDoub0& inp) {
	inp.setVal(value);
	return;
}

void DesignVariable::getValue(DiffDoub1& inp) {
	inp.setVal(diffVal);
	return;
}

void DesignVariable::addCompEl(int eli) {
	for (auto& el : compElList) {
		if (el == eli) {
			return;
		}
	}
	compElList.push_back(eli);
	return;
}
