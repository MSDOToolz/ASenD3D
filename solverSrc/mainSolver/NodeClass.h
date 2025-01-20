#ifndef NODECLASS
#define NODECLASS
#include <string>
#include <iostream>
#include <vector>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"

class Node {
	public:
	    int label;
		bool fluid;
		int numDof;
		int dofIndex[6];
		int sortedRank;
	    double coord[3];
		double displacement[6];
		double velocity[6];
		double flVel[3];
		double flVelDot[3];
		double acceleration[6];
		double temperature;
		double tempChangeRate;
		double flDen;
		double flDenDot;
		double turbE;
		double turbEDot;
		double prevDisp[6];
		double prevVel[6];
		double prevFlVel[3];
		double pFlVelLF[3];
		double prevFlVelDot[3];
		double prevAcc[6];
		double prevTemp;
		double pTempLF;
		double prevTdot;
		double prevFlDen;
		double pFlDenLF;
		double prevFlDenDot;
		double prevTurbE;
		double prevTurbEDot;
		double pTurbELF;
		double initialDisp[6];
		double initialVel[6];
		double initialFlVel[3];
		double initialFlVelDot[3];
		double initialAcc[6];
		double initialTemp;
		double initialTdot;
		double initialFlDen;
		double initialFlDenDot;
		double initialTurbE;
		double initialTurbEDot;
		std::list<IDCapsule> dVarLst;
		
	    Node();
		
		void setCrd(double newCrd[]);
		
		void setDisplacement(double newDisp[]);

		void setVelocity(double newVel[]);

		void setAcceleration(double newAcc[]);

		void setFlVel(double newVel[]);

		void setFlVelDot(double newFlVDot[]);
		
		void addToDisplacement(double delDisp[]);

		void addToFlVel(double delFlVel[]);
		
		void setInitialDisp(double newDisp[]);
		
		void setInitialVel(double newVel[]);
		
		void setInitialAcc(double newAcc[]);

		void setInitialFlVel(double newFlVel[]);

		void setInitialFlVDot(double newFlVDot[]);

		void setPrevDisp(double newDisp[]);

		void setPrevVel(double newVel[]);

		void setPrevAcc(double newAcc[]);

		void setPrevFlVel(double newFlVel[]);

		void setPFlVelLF(double newVel[]);

		void setPrevFlVDot(double newFlVDot[]);
		
		void initializeDisp();

		void initializeFlVel();

		void initializeTemp();

		void initializeFlDen();

		void initializeTurbE();
		
		void updateVelAcc(double nmBeta, double nmGamma, double delT);

		void updateFlVelDot(double nmGamma, double delT);
		
		void updateTdot(double nmGamma, double delT);

		void updateFlDenDot(double nmGamma, double delT);

		void updateTurbEDot(double nmGamma, double delT);
		
		void advanceDisp();

		void advanceFlVel();

		void advanceTemp();

		void advanceFlDen();

		void advanceTurbE();

		void backstepDisp();

		void backstepFlVel();

		void backstepTemp();

		void backstepFlDen();

		void backstepTurbE();
		
		void addDesignVariable(int dIndex, double coef);
		
//dup1
		void getCrd(DiffDoub0 crdOut[], std::vector<DesignVariable>& dvAr);
		
		void getDisp(DiffDoub0 disp[]);
		
		void getElasticDVLoad(DiffDoub0 ld[], std::vector<DesignVariable>& dvAr);

		void getThermalDVLoad(DiffDoub0& ld, std::vector<DesignVariable>& dvAr);
//end dup	
 
//skip 
 
//DiffDoub1 versions: 
//dup1
		void getCrd(DiffDoub1 crdOut[], std::vector<DesignVariable>& dvAr);
		
		void getDisp(DiffDoub1 disp[]);
		
		void getElasticDVLoad(DiffDoub1 ld[], std::vector<DesignVariable>& dvAr);

		void getThermalDVLoad(DiffDoub1& ld, std::vector<DesignVariable>& dvAr);
//end dup	
 
//end skip 
 
 
};

#endif