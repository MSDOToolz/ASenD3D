#ifndef OBJECTIVE
#define OBJECTIVE
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"

class ObjectiveTerm {
	public:
	    std::string category;
		std::string optr;
		double activeTime[2];
		int component;
		int layer;
		double coef;
		double expnt;
		std::string elSetName;
		int elSetPtr;
		std::string ndSetName;
		int ndSetPtr;
		std::string tgtTag;
		std::list<double> tgtVals;
		double value;

		std::vector<double> qVec;
		std::vector<double> elVolVec;
		std::vector<double> tgtVec;
		std::vector<double> errNormVec;

		int qLen;

		SparseMat dQdU;
		SparseMat dQdV;
		SparseMat dQdA;
		SparseMat dQdT;
		SparseMat dQdTdot;
		SparseMat dQdD;
		SparseMat dVdD;
		
	    ObjectiveTerm();
		
		void setActiveTime(double newAt[]);

		void allocateObj(std::vector<Set>& ndSets, std::vector<Set>& elSets);

		void allocateObjGrad();

		double getPowerNorm();

		void dPowerNormdU(std::vector<double>& dLdU, std::vector<double>& dLdV, std::vector<double>& dLdA, std::vector<double>& dLdT, std::vector<double>& dLdTdot);

		void dPowerNormdD(std::vector<double>& dLdD);

		double getVolIntegral();

		void dVolIntegraldU(std::vector<double>& dLdU, std::vector<double>& dLdV, std::vector<double>& dLdA, std::vector<double>& dLdT, std::vector<double>& dLdTdot);

		void dVolIntegraldD(std::vector<double>& dLdD);

		double getVolAverage();

		void dVolAveragedU(std::vector<double>& dLdU, std::vector<double>& dLdV, std::vector<double>& dLdA, std::vector<double>& dLdT, std::vector<double>& dLdTdot);

		void dVolAveragedD(std::vector<double>& dLdD);

		void getObjVal(double time, bool nLGeom, std::vector<Node>& ndAr, std::vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre);

		void getdLdU(std::vector<double>& dLdU, std::vector<double>& dLdV, std::vector<double>& dLdA, std::vector<double>& dLdT, std::vector<double>& dLdTdot, double time, bool nLGeom, std::vector<Node>& ndAr,  std::vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre);

		void getdLdD(std::vector<double>& dLdD, double time, bool nLGeom, std::vector<Node>& ndAr, std::vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr, DiffDoub1StressPrereq& stPre);
};

class Objective {
	public:
		std::vector<ObjectiveTerm> terms;
		
	    Objective();

		void clearValues();

		void calculateTerms(double time, bool nLGeom, std::vector<Node>& ndAr,  std::vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre);

		void calculatedLdU(std::vector<double>& dLdU, std::vector<double>& dLdV, std::vector<double>& dLdA, std::vector<double>& dLdT, std::vector<double>& dLdTdot, double time, bool nLGeom, std::vector<Node>& ndAr, std::vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr, DiffDoub0StressPrereq& stPre);

		void calculatedLdD(std::vector<double>& dLdD, double time, bool nLGeom, std::vector<Node>& ndAr,  std::vector<Element>& elAr, std::vector<Set>& ndSets, std::vector<Set>& elSets, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr, DiffDoub1StressPrereq& stPre);

};

#endif
