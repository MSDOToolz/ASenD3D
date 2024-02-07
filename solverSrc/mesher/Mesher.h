#ifndef MESHER
#define MESHER
#include "MeshNode.h"
#include "MeshElement.h"
#include "MeshFace.h"
#include "SpatialGrid.h"
#include <string>

class Mesher {
private:
	MNList nodes;
	SpatialGrid nodeGrid;
	MEList elements;
	SpatialGrid elementGrid;
	MFList faces;
	SpatialGrid faceGrid;
	
	MeshEnt** gridOut1;
	MeshEnt** gridOut2;
	int gOLen;

	int numBoundNds;
	double avgProj;
	double globProjWt;

public:
	Mesher();

	void readInput(std::string fileName);

	void initBoundaryNormals();

	void prep();

	bool checkNewEl(MeshElement* newEl, MeshFace* newFaces[]);

	bool addFaceIfAbsent(MeshFace* newFace, MeshElement* newEl);

	bool adoptConnectedNd(MeshFace* thisFc, double tgtPt[], double srchRad);

	bool adoptAnyNd(MeshFace* thisFc, double tgtPt[], double srchRad);

	bool createNewNd(MeshFace* thisFc, double tgtPt[]);

	void generateMesh();

	void distributeNodes();

	void writeOutput(std::string fileName);

	void printCurrentMesh();

	~Mesher();
};

#endif
