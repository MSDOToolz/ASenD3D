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
	    ConstraintTerm(std::string newNSet);
		
		void setDof(int newDof);
		
		void setCoef(double newCoef);
		
		std::string getSetName();
		
		int getDof();
		
		double getCoef();
		
		ConstraintTerm* getNext();
		
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
		
		void addTerm(ConstraintTerm *newTerm);
		
		void setRhs(double newRhs);
		
		void setNext(Constraint *newNext);
		
		Constraint* getNext();
		
		void buildMat(Set *firstSet, NdPt ndAr[]);
		
		int getMatDim();
		
		MatrixEnt* getMatFirst(int row);
		
		void destroy();
};

class ConstraintList {
	private:
	    Constraint *firstConst;
		Constraint *lastConst;
	
	public:
	    ConstraintList();
		
		void addConstraint(Constraint *newConst);
		
		Constraint* getFirst();
		
		void destroy();
};

#endif