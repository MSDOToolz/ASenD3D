#ifndef meshelement
#define meshelement
#include "MeshNode.h"
#include <vector>

class MeshElement {
public:
	int nodes[4];

	MeshElement();

	void get_centroid(double cent[], std::vector<MeshNode>& nd_ar);

	double get_volume(std::vector<MeshNode>& nd_ar);

	bool point_in(double pt[], std::vector<MeshNode>& nd_ar);
	
};

#endif
