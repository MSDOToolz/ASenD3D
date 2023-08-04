#ifndef DVARCLASS
#define DVARCLASS
#include <string>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "SetClass.h"

class DesignVariable {
	private:
	    std::string category;
		int component;
		int layer;
		double activeTime[2];
		std::string elSetName;
		std::string ndSetName;
		Set* ndSetPtr;
		DoubList coefs;
		Doub value;
		DiffDoub diffVal;
		IntList compElList;
		DesignVariable *nextDV;
		
    public:
	    DesignVariable(std::string newCat);
		
		void setActiveTime(double newAt[]);
		
		void setLayer(int newLay);
		
		void setElset(std::string newElset);
		
		void setNdset(std::string newNdset);

		void setNdsetPtr(Set* newSet);
		
		void setValue(double newVal);
		
		void setDiffVal(double newVal, double dNewVal);
		
		void setComponent(int newComp);
		
		void setNext(DesignVariable* newNext);
		
		void addCoefficient(double newCoef);

        void addCompEl(int newEl);		
		
		std::string getCategory();
		
		std::string getElSet();
		
		std::string getNdSet();
		
		DoubList* getCoefs();
		
		void getValue(Doub& inp);
		
		void getValue(DiffDoub& inp);
		
		int getComponent();
		
		int getLayer();

		IntListEnt* getFirstEl();

		IntListEnt* getFirstNd();
		
		DesignVariable* getNext();
		
		void destroy();
};

class DVPt {
	public:
	    DesignVariable *ptr;
		
		DVPt();
};

class DesVarList {
	private:
	    DesignVariable *firstDVar;
		DesignVariable *lastDVar;
		int length;
		
	public:
	    DesVarList();
		
		void addDVar(DesignVariable* newDVar);
		
		int getLength();
		
		DesignVariable* getFirst();

		void destroy();
};

#endif