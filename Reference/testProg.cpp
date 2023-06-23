#include <iostream>
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "DiffDoubClass.h"

using namespace std;

int main() {
	DVPt *dvAr = new DVPt[2];
	dvAr[0].ptr = new DesignVariable("nodeCoord");
	dvAr[0].ptr->setValue(1.0);
	dvAr[0].ptr->setComponent(1);
	dvAr[1].ptr = new DesignVariable("nodeCoord");
	dvAr[1].ptr->setValue(2.0);
	dvAr[1].ptr->setComponent(2);
	
	double crd[] = {0.0,0.0,0.0};
	Node myNd = Node(0,crd);
	myNd.addDesignVariable(0,1.0);
	myNd.addDesignVariable(1,2.0);
	
	Doub crdOut[3];
	myNd.getCrd(crdOut,dvAr);
	
	cout << "Coordinates: " << to_string(crdOut[0].val) << "  " << to_string(crdOut[1].val) << "  " << to_string(crdOut[2].val) << endl;
	
	dvAr[0].ptr->setDiffVal(1.0,1.0);
	dvAr[1].ptr->setDiffVal(2.0,1.0);
	DiffDoub dcrdOut[3];
	myNd.getCrd(dcrdOut,dvAr);
	
	cout << "Coordinates: " << to_string(dcrdOut[0].dval) << "  " << to_string(dcrdOut[1].dval) << "  " << to_string(dcrdOut[2].dval) << endl;
	
	cout << "check 1" << endl;
	delete dvAr[0].ptr;
	cout << "check 1.1" << endl;
	delete dvAr[1].ptr;
	cout << "check 1.2" << endl;
	delete[] dvAr;
	
	cout << "check 2" << endl;
	
	myNd.destroy();
	
	cout << "check 3" << endl;
	
	/* double crd[3] = {1.0,2.0,3.0};
	Node myNd = Node(0,crd);
	double newCrd[3] = {0.0,0.0,0.0};
	myNd.getCrd(newCrd);
	cout << "Coordinates: " << to_string(newCrd[0]) << "  " << to_string(newCrd[1]) << "  " << to_string(newCrd[2]) << endl; */
	return 0;
}