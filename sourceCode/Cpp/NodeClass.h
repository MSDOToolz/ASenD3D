#ifndef NODECLASS
#define NODECLASS
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"


class Node {
	private:
	    int label;
		int dofIndex[6];
	    double coord[3];
		double displacement[6];
		double velocity[6];
		double acceleration[6];
		double temperature;
		double tempChangeRate;
		IntList dVars;
		DoubList coefs;
		Node *nextNd;
		
    public:
	    Node(int newLab, double& newCrd);
		
		void addDesignVariable(int dIndex, double coef);
		
		//dup1
		void getCrd(Doub& crdOut, DVPt& dvAr);
		
		void getDisp(Doub& disp);
		//end dup
		
		//skip
		
		void getCrd(DiffDoub& crdOut, DVPt& dvAr);
		
		//end skip
		
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
		
		void addNode(int newLab, double newCrd[]);
		
		int getLength();
		
		Node* getFirst();
};

#endif