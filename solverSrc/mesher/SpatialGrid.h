#ifndef spatialgrid
#define spatialgrid
#include <list>
#include <vector>
#include "MeshNode.h"
#include "MeshElement.h"
#include "MeshFace.h"

//class mesh_ent {
//private:
//	MeshNode* nd;
//	MeshElement* el;
//	MeshFace* fc;
//	mesh_ent* next;
//
//public:
//	mesh_ent();
//
//	void set_pt(MeshNode* new_nd);
//
//	void set_pt(MeshElement* new_el);
//
//	void set_pt(MeshFace* new_fc);
//
//	void set_next(mesh_ent* new_next);
//
//	MeshNode* get_pt(MeshNode* pt);
//
//	MeshElement* get_pt(MeshElement* pt);
//
//	MeshFace* get_pt(MeshFace* pt);
//
//	mesh_ent* get_next();
//
//};
//
//class entity_list {
//private:
//	mesh_ent* first;
//	mesh_ent* last;
//	int length;
//
//public:
//	entity_list();
//
//	void add_entry(MeshNode* new_nd);
//
//	void add_entry(MeshElement* new_el);
//
//	void add_entry(MeshFace* new_fc);
//
//	int copy_to_array(mesh_ent* out_lst[], int max_len);
//
//	~entity_list();
//};

class IntList {
public:
	std::list<int> i_lst;

	IntList();

	int copy_to_vector(std::vector<int>& in_vec, int st_i, int max_len);
};

class SpatialGrid {
private:
	double x_min;
	double x_sp;
	int x_bins;
	double y_min;
	double y_sp;
	int y_bins;
	double z_min;
	double z_sp;
	int z_bins;
	std::vector<IntList> list_ar;

public:
	SpatialGrid();

	void initialize(double x_range[], double x_spacing, double y_range[], double y_spacing, double z_range[], double z_spacing);

	void add_ent(int label, double crd[]);

	int get_in_xyzrange(std::vector<int>& out_lst, int max_len, double x_range[], double y_range[], double z_range[]);

	int get_in_radius(std::vector<int>& out_list, int max_len, double pt[], double rad);
};

#endif
