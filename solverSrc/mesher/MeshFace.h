#ifndef meshface
#define meshface
#include "MeshNode.h"
#include "MeshElement.h"
#include <vector>

class MeshFace {
public:

	int nodes[3];
	int elements[2];
	double norm_dir[3];
	double proj_dist;

	MeshFace();

	void copy_data(MeshFace& f_in);

	void init_norm_dir(std::vector<MeshNode>& nd_ar);

	void norm_dir_from_el_cent(double cent[], std::vector<MeshNode>& nd_ar);

	void get_centroid(double cent[], std::vector<MeshNode>& nd_ar);

	double get_longest_edge_len(std::vector<MeshNode>& nd_ar);

	bool get_intersection(double out_param[], double pt[], double vec[], std::vector<MeshNode>& nd_ar);

	bool edges_intersect(int fc, double dist_tol, std::vector<MeshNode>& nd_ar, std::vector<MeshFace>& fc_ar);

	int get_shared_nodes(int nd_pts[], bool shared[], int fc, std::vector<MeshNode>& nd_ar, std::vector<MeshFace>& fc_ar);

	void print_info(std::vector<MeshNode>& nd_ar);
};

#endif