#ifndef MESHFACE
#define MESHFACE
#include "MeshNode.h"
#include "MeshElement.h"
#include <vector>

class MeshFace {
public:

	int nodes[3];
	int elements[2];
	double normDir[3];
	double projDist;

	MeshFace();

	void copy_data(MeshFace& f_in);

	void initNormDir(std::vector<MeshNode>& nd_ar);

	void normDirFromElCent(double cent[], std::vector<MeshNode>& nd_ar);

	void getCentroid(double cent[], std::vector<MeshNode>& nd_ar);

	double getLongestEdgeLen(std::vector<MeshNode>& nd_ar);

	bool getIntersection(double outParam[], double pt[], double vec[], std::vector<MeshNode>& nd_ar);

	bool edgesIntersect(int fc, double distTol, std::vector<MeshNode>& nd_ar, std::vector<MeshFace>& fc_ar);

	int getSharedNodes(int ndPts[], bool shared[], int fc, std::vector<MeshNode>& nd_ar, std::vector<MeshFace>& fc_ar);

	void printInfo(std::vector<MeshNode>& nd_ar);
};

#endif