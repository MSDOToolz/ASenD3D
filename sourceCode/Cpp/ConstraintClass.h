#ifndef CONSTRAINTS
#define CONSTRAINTS
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>

class ConstraintTerm {
	private:
	    std::string nodeSet;
		int dof;
		double coef;
		ConstraintTerm *nextTerm;
		
	public:
	    ConstraintTerm(std::string newNSet, int newDof, double newCoef);
		
		std::string getSetName();
		
		int getDof();
		
		double getCoef();
		
		void setNext(ConstraintTerm *newNext);
};

class Constraint {
	private:
	    ConstraintTerm *firstTerm;
		ConstraintTerm *lastTerm;
		double rhs;
		SparseMat mat;
		Constraint *nextConst;
		
	public:
	    Constraint();
		
		void addTerm(string nodeSet, int dof, double coef);
		
		Constraint* getNext();
		
		void buildMat(Set *firstSet, NdPt *ndAr);
		
		int getMatDim();
		
		MatrixEnt* getMatFirst(int row);
};

class ConstraintList {
	private:
	    Constraint *firstConst;
		Constraint *lastConst;
	
	public:
	    ConstraintList();
		
		Constraint* getFirst();
};

#endif