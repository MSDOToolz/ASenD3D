#ifndef MESHER
#define MESHER
#include "MeshNode.h"
#include "MeshElement.h"
#include "MeshFace.h"
#include "SpatialGrid.h"
#include <string>
#include <vector>

class Mesher {
private:
	std::vector<MeshNode> nodes;
	int nd_ct;
	int nd_cap;
	SpatialGrid nodeGrid;
	std::vector<MeshElement> elements;
	int el_ct;
	int el_cap;
	SpatialGrid elementGrid;
	std::vector<MeshFace> faces;
	int fc_ct;
	int fc_cap;
	SpatialGrid faceGrid;
	
	std::vector<int> gridOut1;
	std::vector<int> gridOut2;
	int gOLen;

	int numBoundNds;
	double avgProj;
	double maxProj;
	double maxEdgeLen;
	double globProjWt;
	int maxNumEls;

	int newEl;
	std::vector<MeshFace> newElFcs;
	int newNd;

public:
	Mesher();

	void readInput(std::string fileName);

	void initBoundaryNormals();

	void prep();

	bool checkNewEl(MeshElement& newEl, std::vector<MeshFace>& newFaces);

	bool addFaceIfAbsent(int newEl);

	bool adoptConnectedNd(int fc_i, double tgtPt[], double srchRad);

	bool adoptAnyNd(int fc_i, double tgtPt[], double srchRad);

	bool createNewNd(int fc_i, double tgtPt[]);

	bool generateMesh();

	void distributeNodes();

	void writeOutput(std::string fileName);

	void printCurrentMesh();
};

#endif
