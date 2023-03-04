#include <iostream>
#include "NodeClass.cpp"

using namespace std;

int main() {
	DVPt *dvAr = new DVPt[2];
	dvAr[0].ptr = new DesignVariable("nodeCoord",1);
	dvAr[0].ptr->r_setValue(1.0);
	dvAr[1].ptr = new DesignVariable("nodeCoord",2);
	dvAr[1].ptr->r_setValue(2.0);
	
	double crd[] = {0.0,0.0,0.0};
	Node myNd = Node(0,crd);
	myNd.addDesignVariable(0,1.0);
	myNd.addDesignVariable(1,2.0);
	
	myNd.r_getCrd(crd,dvAr);
	
	cout << "Coordinates: " << to_string(crd[0]) << "  " << to_string(crd[1]) << "  " << to_string(crd[2]) << endl;
	
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