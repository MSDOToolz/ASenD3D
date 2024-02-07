#ifndef MESHELEMENT
#define MESHELEMENT
#include "MeshNode.h"

class MeshElement {
private:
	MeshNode* nodes[4];
	MeshElement* next;

public:
	MeshElement();

	void setNdPt(MeshNode* newNds[]);

	void setNext(MeshElement* newNext);

	MeshNode** getNodes();

	void getCentroid(double cent[]);

	double getVolume();

	bool pointIn(double pt[]);

	MeshElement* getNext();
	
};

class MEList {
private:
	MeshElement* first;
	MeshElement* last;
	int length;

public:
	MEList();

	void addEnt(MeshElement* newEl);

	MeshElement* getFirst();

	int getLength();

	~MEList();

};

#endif
