#ifndef NODECLASS
#define NODECLASS
#include <string>
#include <iostream>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"

class Node {
	private:
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
		double initialDisp[6];
		double initialVel[6];
		double initialFlVel[3];
		double initialFlVelDot[3];
		double initialAcc[6];
		double initialTemp;
		double initialTdot;
		double initialFlDen;
		double initialFlDenDot;
		IntList dVars;
		DoubList coefs;
		Node *nextNd;
		
    public:
	    Node(int newLab);

		void setFluid(bool newFluid);
		
		void setNumDof(int nDof);
		
		void setCrd(double newCrd[]);
		
		void setDofIndex(int dof, int index);

		void setSortedRank(int newRank);
		
		void setDisplacement(double newDisp[]);

		void setVelocity(double newVel[]);

		void setAcceleration(double newAcc[]);

		void setFlVel(double newVel[]);

		void setFlVelDot(double newFlVDot[]);

		void setTemperature(double newTemp);

		void setTdot(double newTdot);

		void setFlDen(double newFlDen);

		void setFlDenDot(double newFlDDot);
		
		void addToDisplacement(double delDisp[]);

		void addToFlVel(double delFlVel[]);

		void addToTemperature(double delTemp);

		void addToFlDen(double delFlDen);
		
		void setInitialDisp(double newDisp[]);
		
		void setInitialVel(double newVel[]);
		
		void setInitialAcc(double newAcc[]);

		void setInitialFlVel(double newFlVel[]);

		void setInitialFlVDot(double newFlVDot[]);
		
		void setInitialTemp(double newTemp);
		
		void setInitialTdot(double newTdot);

		void setInitialFlDen(double newFlDen);

		void setInitialFlDDot(double newFlDen);

		void setPrevDisp(double newDisp[]);

		void setPrevVel(double newVel[]);

		void setPrevAcc(double newAcc[]);

		void setPrevFlVel(double newFlVel[]);

		void setPFlVelLF(double newVel[]);

		void setPrevFlVDot(double newFlVDot[]);

		void setPrevTemp(double newTemp);

		void setPTempLF(double newTemp);

		void setPrevTdot(double newTdot);

		void setPrevFlDen(double newFlDen);

		void setPFlDenLF(double newDen);

		void setPrevFlDDot(double newFlDDot);
		
		void initializeDisp();

		void initializeFlVel();

		void initializeTemp();

		void initializeFlDen();
		
		void updateVelAcc(double nmBeta, double nmGamma, double delT);

		void updateFlVelDot(double nmGamma, double delT);
		
		void updateTdot(double nmGamma, double delT);

		void updateFlDenDot(double nmGamma, double delT);
		
		void advanceDisp();

		void advanceFlVel();

		void advanceTemp();

		void advanceFlDen();

		void backstepDisp();

		void backstepFlVel();

		void backstepTemp();

		void backstepFlDen();
		
		void addDesignVariable(int dIndex, double coef);
		
		IntList* getDesignVars();
		
		void getCrd(double crdOut[]);
		
		int getLabel();

		bool isFluid();
		
		int getNumDof();
		
		void getDisp(double dispOut[]);
		
		void getVel(double velOut[]);
		
		void getAcc(double accOut[]);

		void getFlVel(double velOut[]);

		void getFlVelDot(double velDotOut[]);
		
		double getTemperature();
		
		double getTdot();

		double getFlDen();

		double getFlDenDot();

		void getPrevDisp(double dispOut[]);

		void getPrevVel(double velOut[]);

		void getPrevAcc(double accOut[]);

		void getPrevFlVel(double flVelOut[]);

		void getPrevFlVelDot(double flVDOut[]);

		double getPrevTemp();

		double getPrevTdot();

		double getPrevFlDen();

		double getPrevFlDenDot();
		
//dup1
		void getCrd(DiffDoub0 crdOut[], DesignVariable* dvAr[]);

		void getDefCrd(DiffDoub0 crdOut[], DesignVariable* dvAr[]);
		
		void getDisp(DiffDoub0 disp[]);
		
		void getElasticDVLoad(DiffDoub0 ld[], DesignVariable* dvAr[]);

		void getThermalDVLoad(DiffDoub0& ld, DesignVariable* dvAr[]);
//end dup		
 
//skip 
 
//DiffDoub1 versions: 
//dup1
		void getCrd(DiffDoub1 crdOut[], DesignVariable* dvAr[]);
		
		void getDisp(DiffDoub1 disp[]);
		
		void getElasticDVLoad(DiffDoub1 ld[], DesignVariable* dvAr[]);

		void getThermalDVLoad(DiffDoub1& ld, DesignVariable* dvAr[]);
//end dup		
 
//end skip 
 
 
		int getDofIndex(int dof);

		int getSortedRank();
		
		Node* getNext();
		
		void setNext(Node *newNext);
};

class NodeList {
	private:
	    Node *firstNode;
		Node *lastNode;
		int length;
		
	public:
	    NodeList();
		
		void addNode(Node *newNd);
		
		int getLength();
		
		Node* getFirst();

		~NodeList();
};

#endif