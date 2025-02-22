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
	SpatialGrid node_grid;
	std::vector<MeshElement> elements;
	int el_ct;
	int el_cap;
	SpatialGrid element_grid;
	std::vector<MeshFace> faces;
	int fc_ct;
	int fc_cap;
	SpatialGrid face_grid;
	
	std::vector<int> grid_out1;
	std::vector<int> grid_out2;
	int g_olen;

	int num_bound_nds;
	double avg_proj;
	double max_proj;
	double max_edge_len;
	double glob_proj_wt;
	int max_num_els;

	int new_el;
	std::vector<MeshFace> new_el_fcs;
	int new_nd;

public:
	Mesher();

	void read_input(std::string file_name);

	void init_boundary_normals();

	void prep();

	bool check_new_el(MeshElement& new_el, std::vector<MeshFace>& new_faces);

	bool add_face_if_absent(int new_el);

	bool adopt_connected_nd(int fc_i, double tgt_pt[], double srch_rad);

	bool adopt_any_nd(int fc_i, double tgt_pt[], double srch_rad);

	bool create_new_nd(int fc_i, double tgt_pt[]);

	bool generate_mesh();

	void distribute_nodes();

	void write_output(std::string file_name);

	void print_current_mesh();
};

#endif
