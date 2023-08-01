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
		std::string type;
		double rhs;
		SparseMat mat;
		double scaleFact;
		Constraint *nextConst;
		
	public:
	    Constraint();

		void setType(std::string newType);
		
		void addTerm(ConstraintTerm *newTerm);
		
		void setRhs(double newRhs);

		void setScaleFact(double newSFact);
		
		void setNext(Constraint *newNext);
		
		Constraint* getNext();
		
		void buildMat(Set *firstSet, NdPt ndAr[]);
		
		int getMatDim();

		double getScaleFact();
		
		MatrixEnt* getMatFirst(int row);

		void getLoad(double cLd[], double uVec[], double qVec[],int resDim);
		
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

		void scaleElastic(double newSF);

		void getTotalLoad(double cLd[], double uVec[], double qVec[], int resDim);
		
		void destroy();
};

#endif