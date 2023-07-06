#ifndef NODECLASS
#define NODECLASS
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"


class Node {
	private:
	    int label;
		int numDof;
		int dofIndex[6];
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
		
		void setInitialDisp(double newDisp[]);
		
		void setInitialVel(double newVel[]);
		
		void setInitialAcc(double newAcc[]);
		
		void setInitialTemp(double newTemp);
		
		void setInitialTdot(double newTdot);
		
		void addDesignVariable(int dIndex, double coef);
		
		IntList* getDesignVars();
		
		void getCrd(double crdOut[]);
		
		int getLabel();
		
		int getNumDof();
		
		//dup1
		void getCrd(Doub crdOut[], DVPt dvAr[]);
		
		void getDisp(Doub disp[]);
		//end dup		
		
		int getDofIndex(int dof);
		
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