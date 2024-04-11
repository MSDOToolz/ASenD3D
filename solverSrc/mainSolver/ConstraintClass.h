#ifndef CONSTRAINTS
#define CONSTRAINTS
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>
#include <fstream>

class ConstraintTerm {
	private:
	    std::string nodeSet;
		Set* nsPtr;
		int dof;
		double coef;
		ConstraintTerm *nextTerm;
		
	public:
	    ConstraintTerm(std::string newNSet);

		void setNsPtr(Set* newSet);
		
		void setDof(int newDof);
		
		void setCoef(double newCoef);
		
		std::string getSetName();

		Set* getSetPtr();
		
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
		Constraint* nextConst;
		
	public:
	    Constraint();

		void setType(std::string newType);
		
		void addTerm(ConstraintTerm *newTerm);
		
		void setRhs(double newRhs);

		void setScaleFact(double newSFact);
		
		void setNext(Constraint *newNext);
		
		Constraint* getNext();
		
		void buildMat(Node* ndAr[]);

		ConstraintTerm* getFirst();
		
		int getMatDim();

		double getScaleFact();
		
		MatrixEnt* getMatFirst(int row);

		void fullVecMultiply(double prod[], double vec[], double tmpV[]);

		void getLoad(double cLd[], double uVec[], double qVec[],int resDim);

		void writeToFile(std::ofstream& outFile);
		
		~Constraint();
};

class ConstraintList {
	private:
	    Constraint *firstConst;
		Constraint *lastConst;
	
	public:
	    ConstraintList();
		
		void addConstraint(Constraint *newConst);
		
		Constraint* getFirst();

		void setScaleFact(double newSF);

		void getTotalVecMult(double prod[], double vec[], double tmpV[]);

		void getTotalLoad(double cLd[], double uVec[], double qVec[], int resDim);

		void writeAllToFile(std::string fileName);
		
		~ConstraintList();
};

#endif