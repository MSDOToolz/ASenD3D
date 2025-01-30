#ifndef MESHELEMENT
#define MESHELEMENT
#include "MeshNode.h"
#include <vector>

class MeshElement {
public:
	int nodes[4];

	MeshElement();

	void getCentroid(double cent[], std::vector<MeshNode>& nd_ar);

	double getVolume(std::vector<MeshNode>& nd_ar);

	bool pointIn(double pt[], std::vector<MeshNode>& nd_ar);
	
};

#endif
