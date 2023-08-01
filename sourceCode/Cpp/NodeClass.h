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
		int numDof;
		int dofIndex[6];
		int sortedRank;
	    double coord[3];
		double displacement[6];
		double velocity[6];
		double acceleration[6];
		double temperature;
		double tempChangeRate;
		double prevDisp[6];
		double prevVel[6];
		double prevAcc[6];
		double prevTemp;
		double prevTdot;
		double initialDisp[6];
		double initialVel[6];
		double initialAcc[6];
		double initialTemp;
		double initialTdot;
		IntList dVars;
		DoubList coefs;
		Node *nextNd;
		
    public:
	    Node(int newLab);
		
		void setNumDof(int nDof);
		
		void setCrd(double newCrd[]);
		
		void setDofIndex(int dof, int index);

		void setSortedRank(int newRank);
		
		void setDisplacement(double newDisp[]);
		
		void addToDisplacement(double delDisp[]);
		
		void setInitialDisp(double newDisp[]);
		
		void setInitialVel(double newVel[]);
		
		void setInitialAcc(double newAcc[]);
		
		void setInitialTemp(double newTemp);
		
		void setInitialTdot(double newTdot);
		
		void initializeDisp();
		
		void updateVelAcc(double nmBeta, double nmGamma, double delT);
		
		void advanceDisp();
		
		void addDesignVariable(int dIndex, double coef);
		
		IntList* getDesignVars();
		
		void getCrd(double crdOut[]);
		
		int getLabel();
		
		int getNumDof();
		
		void getDisp(double dispOut[]);
		
		void getVel(double velOut[]);
		
		void getAcc(double accOut[]);
		
		double getTemperature();
		
		double getTdot();
		
		//dup1
		void getCrd(Doub crdOut[], DVPt dvAr[]);
		
		void getDisp(Doub disp[]);
		
		void getElasticDVLoad(Doub ld[], DVPt dvAr[]);
		//end dup		
//skip 
 
//DiffDoub versions: 
		void getCrd(DiffDoub crdOut[], DVPt dvAr[]);
		
		void getDisp(DiffDoub disp[]);
		
		void getElasticDVLoad(DiffDoub ld[], DVPt dvAr[]);
 
//end skip 
		
		int getDofIndex(int dof);

		int getSortedRank();
		
		Node* getNext();
		
		void setNext(Node *newNext);
		
		void destroy();
};

class NdPt {
	public:
	    Node *ptr;
		
		NdPt();
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
};

#endif