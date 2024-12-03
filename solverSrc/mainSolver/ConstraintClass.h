#ifndef CONSTRAINTS
#define CONSTRAINTS
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include <string>
#include <fstream>

class ConstraintTerm {
	public:
	    std::string nodeSet;
		int nsPtr;
		int dof;
		double coef;
		
	    ConstraintTerm();
		
};

class Constraint {
	public:
		std::string type;
		std::list<ConstraintTerm> terms;
		double rhs;
		SparseMat mat;
		double scaleFact;
		
	    Constraint();
		
		void buildMat(std::vector<Node>& ndAr, std::vector<Set>& setAr);

		void fullVecMultiply(std::vector<double>& prod, std::vector<double>& vec, std::vector<double>& tmpV);

		void getLoad(std::vector<double>& cLd, std::vector<double>& uVec, std::vector<double>& qVec,int resDim);

		void writeToFile(std::ofstream& outFile);
};

class ConstraintList {
	public:
		std::vector<Constraint> constVec;
	
	    ConstraintList();

		void setScaleFact(double newSF);

		void getTotalVecMult(std::vector<double>& prod, std::vector<double>& vec, std::vector<double>& tmpV);

		void getTotalLoad(std::vector<double>& cLd, std::vector<double>& uVec, std::vector<double>& qVec, int resDim);

		void writeAllToFile(std::string fileName);
};

#endif