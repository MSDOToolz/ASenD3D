#ifndef OBJECTIVE
#define OBJECTIVE
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"

class ObjectiveTerm {
	private:
	    std::string category;
		std::string optr;
		double activeTime[2];
		int component;
		int layer;
		double coef;
		double expnt;
		std::string elSetName;
		Set* elSetPtr;
		std::string ndSetName;
		Set* ndSetPtr;
		std::string tgtTag;
		DoubList tgtVals;
		ObjectiveTerm* next;

		double* qVec;
		double* elVolVec;
		double* tgtVec;
		double* errNormVec;

		int qLen;

		SparseMat dQdU;
		SparseMat dQdV;
		SparseMat dQdA;
		SparseMat dQdT;
		SparseMat dQdTdot;
		SparseMat dQdD;
		SparseMat dVdD;
		
	public:
	    ObjectiveTerm(std::string newCat);
		
		void setOperator(std::string newOp);
		
		void setActiveTime(double newAt[]);
		
		void setComponent(int newComp);
		
		void setLayer(int newLay);
		
		void setCoef(double newCoef);
		
		void setExponent(double newExp);
		
		void setElset(std::string newElset);

		void setElsetPtr(Set* newPtr);
		
		void setNdset(std::string newNdset);

		void setNdsetPtr(Set* newPtr);
		
		void setTgtTag(std::string newTag);

		void setNext(ObjectiveTerm* newNext);
		
		void addTargetValue(double newTgt);

		std::string getElsetName();

		std::string getNdsetName();

		ObjectiveTerm* getNext();

		void allocateObj();

		void allocateObjGrad();

		double getPowerNorm();

		void dPowerNormdU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[]);

		void dPowerNormdD(double dLdD[]);

		double getVolIntegral();

		void dVolIntegraldU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[]);

		void dVolIntegraldD(double dLdD[]);

		double getVolAverage();

		void dVolAveragedU(double dLdU[], double dLdV[], double dLdA[], double dLdT[], double dLdTdot[]);

		void dVolAveragedD(double dLdD[]);

		double getObjVal(NdPt ndAr[], ElPt elAr[], DVPt dvAr[]);
};

class Objective {
	private:
	    ObjectiveTerm *firstTerm;
		ObjectiveTerm *lastTerm;
		int length;
		
	public:
	    Objective();
		
		void addTerm(ObjectiveTerm *newTerm);
		
		int getLength();
		
		ObjectiveTerm* getFirst();
};

#endif
