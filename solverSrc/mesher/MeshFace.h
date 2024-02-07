#ifndef MESHFACE
#define MESHFACE
#include "MeshNode.h"
#include "MeshElement.h"

class MeshFace {
private:
	int nodeLabs[3];
	MeshNode* nodes[3];
	MeshElement* elements[2];
	double normDir[3];
	double projDist;
	MeshFace* next;

public:
	MeshFace();

	void setNodeLabs(int newLabs[]);

	void setNodePt(MeshNode* newNds[]);

	void setPtFromLabs(MeshNode* ndAr[]);

	void setElPt(MeshElement* newEls[]);

	void initNormDir();

	void normDirFromElCent(double cent[]);

	void setNext(MeshFace* newNext);

	MeshNode** getNdPt();

	MeshElement** getElPt();

	double* getNormDir();

	void getCentroid(double cent[]);

	double getProjDist();

	bool getIntersection(double outParam[], double pt[], double vec[]);

	int getSharedNodes(MeshNode* ndPts[], bool shared[], MeshFace* fc);

	MeshFace* getNext();
};

class MFList {
private:
	MeshFace* first;
	MeshFace* last;
	int length;

public:
	MFList();

	void addEnt(MeshFace* newFc);

	MeshFace* getFirst();

	int getLength();

	~MFList();
};

#endif