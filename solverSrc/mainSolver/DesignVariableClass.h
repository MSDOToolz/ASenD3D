#ifndef DVARCLASS
#define DVARCLASS
#include <string>
#include <vector>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "SetClass.h"

class DesignVariable {
	public:
	    std::string category;
		int component;
		int layer;
		double activeTime[2];
		std::string elSetName;
		int elSetPtr;
		std::string ndSetName;
		int ndSetPtr;
		std::list<double> coefs;
		DiffDoub0 value;
		DiffDoub1 diffVal;
		std::list<int> compElList;
		
	    DesignVariable();
		
		void setActiveTime(double newAt[]);
		
		void getValue(DiffDoub0& inp);
		
		void getValue(DiffDoub1& inp);

		void addCompEl(int eli);

};

#endif