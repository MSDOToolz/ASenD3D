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
		Set* elSetPtr;
		std::string ndSetName;
		Set* ndSetPtr;
		DoubList coefs;
		DiffDoub0 value;
		DiffDoub1 diffVal;
		IntList compElList;
		DesignVariable *nextDV;
		
    public:
	    DesignVariable(std::string newCat);
		
		void setActiveTime(double newAt[]);
		
		void setLayer(int newLay);
		
		void setElset(std::string newElset);

		void setElsetPtr(Set* newSet);
		
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

		Set* getElsetPtr();
		
		std::string getNdSet();

		Set* getNdsetPtr();
		
		DoubList* getCoefs();
		
		void getValue(DiffDoub0& inp);
		
		void getValue(DiffDoub1& inp);
		
		int getComponent();
		
		int getLayer();

		IntListEnt* getFirstEl();

		IntListEnt* getFirstNd();
		
		DesignVariable* getNext();
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

		~DesVarList();
};

#endif