#include "Mesher.h"
#include "constants.h"
#include "utilities.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

Mesher::Mesher() {
	nd_ct = 0;
	fc_ct = 0;
	el_ct = 0;
	grid_out1.clear();
	grid_out2.clear();
	glob_proj_wt = 0.75;
	max_num_els = max_int;
	new_el = max_int;
	new_el_fcs = vector<MeshFace>(4);
	new_nd = max_int;
	return;
}

void Mesher::read_input(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4];
	int hd_ld_space[4];
	string data[3];
	int data_len;
	int new_nd;
	double nd_crd[3];
	int new_fc;
	int fc_nd[3];

	nd_ct = 0;
	fc_ct = 0;

	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "nodes" && data_len == 3) {
				nd_ct++;
			}
			else if (headings[0] == "faces" && data_len == 3) {
				fc_ct++;
			}
			else if (headings[0] == "globProjWt" && data_len == 1) {
				glob_proj_wt = stod(data[0]);
			}
			else if (headings[0] == "maxNumEls" && data_len == 1) {
				max_num_els = stoi(data[0]);
			}
		}
		in_file.close();
	}
	else {
		string er_st = "Error: Could not open input file " + file_name + " for unstructured 3D mesh generation.";
		throw invalid_argument(er_st);
	}

	int rad = 1;
	while (12*rad*rad < nd_ct) {
		rad++;
	}

	nd_cap = 15 * rad * rad * rad;
	el_cap = 6 * nd_cap;
	if (max_num_els == max_int) {
		max_num_els = el_cap;
	}
	else if (el_cap > max_num_els) {
		el_cap = max_num_els;
	}
	fc_cap = 2 * el_cap;

	nodes = vector<MeshNode>(nd_cap);
	elements = vector<MeshElement>(el_cap);
	faces = vector<MeshFace>(fc_cap);

	nd_ct = 0;
	fc_ct = 0;

	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "nodes" && data_len == 3) {
				nodes[nd_ct].coord[0] = stod(data[0]);
				nodes[nd_ct].coord[1] = stod(data[1]);
				nodes[nd_ct].coord[2] = stod(data[2]);
				nd_ct++;
			}
			else if (headings[0] == "faces" && data_len == 3) {
				faces[fc_ct].nodes[0] = stoi(data[0]);
				faces[fc_ct].nodes[1] = stoi(data[1]);
				faces[fc_ct].nodes[2] = stoi(data[2]);
				fc_ct++;
			}
		}
		in_file.close();
	}
	
	return;
}

void Mesher::init_boundary_normals() {
	int i1;
	int lst_len;
	double cent[3];
	double* norm_dir;
	int srch_dir;
	double proj;
	double x_range[2];
	double y_range[2];
	double z_range[2];
	double vec[3];
	double int_out[3];
	double dp;
	bool intersects;
	int int_ct;
	int fci;
	int fci2;
	for (fci = 0; fci < fc_ct; fci++) {
		MeshFace& this_fc = faces[fci];
		this_fc.get_centroid(cent,nodes);
		norm_dir = &this_fc.norm_dir[0];
		srch_dir = 0;
		if (abs(norm_dir[1]) > abs(norm_dir[0])) {
			srch_dir = 1;
		}
		if (abs(norm_dir[2]) > abs(norm_dir[1])) {
			srch_dir = 2;
		}
		if (srch_dir == 0) {
			vec[0] = 1.0;
			vec[1] = 7.52893558402e-4;
			vec[2] = 5.39890329009e-4;
			x_range[0] = 1;
			x_range[1] = 0;
			y_range[0] = cent[1] - avg_proj;
			y_range[1] = cent[1] + avg_proj;
			z_range[0] = cent[2] - avg_proj;
			z_range[1] = cent[2] + avg_proj;
		}
		else if (srch_dir == 1) {
			vec[0] = 7.52893558402e-4;
			vec[1] = 1.0;
			vec[2] = 5.39890329009e-4;
			x_range[0] = cent[0] - avg_proj;
			x_range[1] = cent[0] + avg_proj;
			y_range[0] = 1;
			y_range[1] = 0;
			z_range[0] = cent[2] - avg_proj;
			z_range[1] = cent[2] + avg_proj;
		}
		else {
			vec[0] = 7.52893558402e-4;
			vec[1] = 5.39890329009e-4;
			vec[2] = 1.0;
			x_range[0] = cent[0] - avg_proj;
			x_range[1] = cent[0] + avg_proj;
			y_range[0] = cent[1] - avg_proj;
			y_range[1] = cent[1] + avg_proj;
			z_range[0] = 1;
			z_range[1] = 0;
		}
		lst_len = face_grid.get_in_xyzrange(grid_out1, g_olen, x_range, y_range, z_range);
		int_ct = 1;
		for (i1 = 0; i1 < lst_len; i1++) {
			fci2 = grid_out1[i1];
			MeshFace& this_fc2 = faces[fci2];
			if (fci != fci2) {
				intersects = this_fc2.get_intersection(int_out, cent, vec, nodes);
				if (intersects && int_out[0] > 0.0) {
					int_ct *= -1;
				}
			}
		}
		// if in_ct > 0, number of intersections is even, vec points out of the surface
		dp = vec[0] * norm_dir[0] + vec[1] * norm_dir[1] + vec[2] * norm_dir[2];
		if ((int_ct > 0 && dp > 0.0) || (int_ct < 0 && dp < 0.0)) {
			norm_dir[0] *= -1.0;
			norm_dir[1] *= -1.0;
			norm_dir[2] *= -1.0;
		}
		//
		if (cent[0] > 9.99) {
			dp = dp;
		}
		//
	}
	return;
}

void Mesher::prep() {
	int i1;
	int i2;
	num_bound_nds = nd_ct;

	//allocate the grid output list
	int num_faces = fc_ct;
	if (num_faces > num_bound_nds) {
		g_olen = 2 * num_faces;
	}
	else {
		g_olen = 2 * num_bound_nds;
	}
	grid_out1 = vector<int>(g_olen);
	grid_out2 = vector<int>(g_olen);

	//initialize grids
	double x_range[2] = { 1.0e+100,-1.0e+100 };
	double y_range[2] = { 1.0e+100,-1.0e+100 };
	double z_range[2] = { 1.0e+100,-1.0e+100 };
	double* crd;
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& this_nd = nodes[i1];
		crd = &this_nd.coord[0];
		if (crd[0] < x_range[0]) {
			x_range[0] = crd[0];
		}
		if (crd[0] > x_range[1]) {
			x_range[1] = crd[0];
		}
		if (crd[1] < y_range[0]) {
			y_range[0] = crd[1];
		}
		if (crd[1] > y_range[1]) {
			y_range[1] = crd[1];
		}
		if (crd[2] < z_range[0]) {
			z_range[0] = crd[2];
		}
		if (crd[2] > z_range[1]) {
			z_range[1] = crd[2];
		}
	}

	double spacing = 0.0;
	avg_proj = 0.0;
	max_proj = 0.0;
	max_edge_len = 0.0;
	for (i1 = 0; i1 < fc_ct; i1++) {
		MeshFace& this_fc = faces[i1];
		spacing = this_fc.proj_dist;
		avg_proj += spacing;
		if (spacing > max_proj) {
			max_proj = spacing;
		}
		spacing = this_fc.get_longest_edge_len(nodes);
		if (spacing > max_edge_len) {
			max_edge_len = spacing;
		}
	}
	avg_proj /= num_faces;
	spacing = avg_proj;
	node_grid.initialize(x_range, spacing, y_range, spacing, z_range, spacing);
	element_grid.initialize(x_range, spacing, y_range, spacing, z_range, spacing);
	face_grid.initialize(x_range, spacing, y_range, spacing, z_range, spacing);

	//add boundary nodes and faces to their grids
	double cent[3];
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& this_nd = nodes[i1];
		node_grid.add_ent(i1, this_nd.coord);
	}

	for (i1 = 0; i1 < fc_ct; i1++) {
		MeshFace& this_fc = faces[i1];
		this_fc.get_centroid(cent,nodes);
		face_grid.add_ent(i1, cent);
	}

	// check boundary completeness

	int lst_len;
	int lst_nd[3];
	bool shared[3];
	int num_shared;
	int neighb_ct;
	for (i2 = 0; i2 < fc_ct; i2++) {
		MeshFace& this_fc = faces[i2];
		this_fc.get_centroid(cent,nodes);
		lst_len = face_grid.get_in_radius(grid_out1, g_olen, cent, 1.01*max_edge_len);
		neighb_ct = 0;
		for (i1 = 0; i1 < lst_len; i1++) {
			MeshFace& lst_fc = faces[grid_out1[i1]];
			//fc_nds = lst_fc->get_nd_pt();
			//cout << fc_nds[0]->get_label() << ", " << fc_nds[1]->get_label() << ", " << fc_nds[2]->get_label() << endl;
			num_shared = lst_fc.get_shared_nodes(lst_nd, shared, i2, nodes, faces);
			if (num_shared == 2) {
				neighb_ct++;
			}
		}
		if (neighb_ct != 3) {
			//cout << "neighbCt " << neighb_ct;
			//fc_nds = this_fc->get_nd_pt();
			//cout << fc_nds[0]->get_label() << ", " << fc_nds[1]->get_label() << ", " << fc_nds[2]->get_label() << endl;
			string er_st = "Error: Invalid boundary surface input to 3D unstructured mesh generator.\n";
			er_st = er_st + "Boundary must be a completely closed surface of triangular faces.";
			throw invalid_argument(er_st);
		}
	}

	// check face normal directions

	init_boundary_normals();

	return;
}

bool Mesher::check_new_el(MeshElement& new_el, std::vector<MeshFace>& new_faces) {
	int i1;
	int i2;
	int i3;
	int i4;
	int n1;
	int n2;
	int lst_len;
	double cent[3];
	int el_edges[12] = { 0,1,0,2,0,3,1,2,1,3,2,3 };
	int face_edges[6] = { 0,1,0,2,1,2 };
	double* pt;
	double* pt2;
	double vec[3];
	double out_p[3];
	bool intersects;
	int* el_nds = &new_el.nodes[0];
	int* fc_nds = nullptr;
	int* fc_els = nullptr;
	int fc_nd_out[3];
	bool shared[3];

	new_el.get_centroid(cent,nodes);
	lst_len = face_grid.get_in_radius(grid_out2, g_olen, cent, 1.01 * max_edge_len);
	for (i1 = 0; i1 < lst_len; i1++) {
		MeshFace& this_fc = faces[grid_out2[i1]];

		for (i2 = 0; i2 < 4; i2++) {
			// any face of new element already exists and is closed
			i3 = this_fc.get_shared_nodes(fc_nd_out, shared, i2, nodes, new_faces);
			fc_els = &this_fc.elements[0];
			if (i3 == 3 && fc_els[1]) {
				return false;
			}
			// any edge of new element intersects existing edge
			intersects = this_fc.edges_intersect(i2, 1.0e-6 * avg_proj, nodes, new_faces);
			if (intersects) {
				return false;
			}
		}

		// any edge of new element intersects existing face
		i3 = 0;
		for (i2 = 0; i2 < 6; i2++) {
			n1 = el_edges[i3];
			n2 = el_edges[i3 + 1];
			pt = &nodes[el_nds[n1]].coord[0];
			pt2 = &nodes[el_nds[n2]].coord[0];
			vec[0] = pt2[0] - pt[0];
			vec[1] = pt2[1] - pt[1];
			vec[2] = pt2[2] - pt[2];
			intersects = this_fc.get_intersection(out_p, pt, vec, nodes);
			if (intersects && out_p[0] < 0.9999999999 && out_p[0] > 0.0000000001) {
				return false;
			}
			i3 += 2;
		}

		// any edge of existing face intersects face of new element

		fc_nds = &this_fc.nodes[0];
		i3 = 0;
		for (i2 = 0; i2 < 3; i2++) {
			n1 = face_edges[i3];
			n2 = face_edges[i3 + 1];
			pt = &nodes[fc_nds[n1]].coord[0];
			pt2 = &nodes[fc_nds[n2]].coord[0];
			vec[0] = pt2[0] - pt[0];
			vec[1] = pt2[1] - pt[1];
			vec[2] = pt2[2] - pt[2];
			for (i4 = 0; i4 < 4; i4++) {
				intersects = new_faces[i4].get_intersection(out_p, pt, vec, nodes);
				if (intersects && out_p[0] < 0.9999999999 && out_p[0] > 0.0000000001) {
					return false;
				}
			}
			i3 += 2;
		}
	}

	// any node of new element in an existing element

	lst_len = element_grid.get_in_radius(grid_out2, g_olen, cent, 1.01 * max_edge_len);
	for (i1 = 0; i1 < lst_len; i1++) {
		MeshElement& this_el = elements[grid_out2[i1]];
		for (i2 = 0; i2 < 4; i2++) {
			pt = &nodes[el_nds[i2]].coord[0];
			intersects = this_el.point_in(pt,nodes);
			if (intersects) {
				return false;
			}
		}
	}

	// any existing node in the new element

	lst_len = node_grid.get_in_radius(grid_out2, g_olen, cent, 1.01 * max_edge_len);
	for (i1 = 0; i1 < lst_len; i1++) {
		MeshNode& this_nd = nodes[grid_out2[i1]];
		pt = &this_nd.coord[0];
		intersects = new_el.point_in(pt, nodes);
		if (intersects) {
			return false;
		}
	}

	return true;
}

bool Mesher::add_face_if_absent(int new_el) {
	int i1;
	double cent[3];
	MeshFace& new_face = faces[fc_ct];
	new_face.get_centroid(cent, nodes);
	MeshFace* this_fc = nullptr;
	int this_nds[3];
	int* fc_els;
	int new_fc_els[2];
	bool shared[3];
	int num_shared;

	int lst_len = face_grid.get_in_radius(grid_out2, g_olen, cent, 0.1*avg_proj);
	for (i1 = 0; i1 < lst_len; i1++) {
		MeshFace& this_fc = faces[grid_out2[i1]];
		num_shared = this_fc.get_shared_nodes(this_nds, shared, fc_ct, nodes, faces);
		if (num_shared == 3) {
			this_fc.elements[1] = new_el;
			return false;
		}
	}

	face_grid.add_ent(fc_ct, cent);
	fc_ct++;

	return true;
}

bool Mesher::adopt_connected_nd(int fc_i, double tgt_pt[], double srch_rad) {
	int i1;
	int i2;
	MeshFace& this_fc = faces[fc_i];
	int face_nds[3];
	bool shared[3];
	int num_shared;
	int un_shared;
	double* crd;
	double d_vec[3];
	double dist;
	double dp;
	double cent[3];
	int* new_el_nds = &elements[el_ct].nodes[0];
	bool el_check;
	bool fc_added;
	MeshFace* new_fc = nullptr;

	new_el_fcs[0].copy_data(this_fc);
	
	double* norm_dir = &this_fc.norm_dir[0];
	int* this_fc_nds = &this_fc.nodes[0];
	int* this_fc_els = &this_fc.elements[0];
	this_fc.get_centroid(cent, nodes);
	int lst_len = face_grid.get_in_radius(grid_out1, g_olen, tgt_pt, 1.01 * max_edge_len);
	for (i1 = 0; i1 < lst_len; i1++) {
		MeshFace& list_fc = faces[grid_out1[i1]];
		num_shared = list_fc.get_shared_nodes(face_nds, shared, fc_i, nodes, faces);
		if (num_shared == 2) {
			un_shared = max_int;
			for (i2 = 0; i2 < 3; i2++) {
				if (!shared[i2]) {
					un_shared = face_nds[i2];
				}
			}
			crd = &nodes[un_shared].coord[0];
			d_vec[0] = crd[0] - tgt_pt[0];
			d_vec[1] = crd[1] - tgt_pt[1];
			d_vec[2] = crd[2] - tgt_pt[2];
			dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
			if (dist < srch_rad) {
				d_vec[0] = crd[0] - cent[0];
				d_vec[1] = crd[1] - cent[1];
				d_vec[2] = crd[2] - cent[2];
				dp = d_vec[0] * norm_dir[0] + d_vec[1] * norm_dir[1] + d_vec[2] * norm_dir[2];
				if (dp > 1.0e-3*avg_proj) {
					new_el_nds[0] = this_fc_nds[0];
					new_el_nds[1] = this_fc_nds[1];
					new_el_nds[2] = this_fc_nds[2];
					new_el_nds[3] = un_shared;

					new_fc = &new_el_fcs[1];
					new_fc->nodes[0] = this_fc_nds[0];
					new_fc->nodes[1] = this_fc_nds[1];
					new_fc->nodes[2] = un_shared;
					new_fc->elements[0] = el_ct;
					new_fc->elements[1] = max_int;
					new_fc->init_norm_dir(nodes);

					new_fc = &new_el_fcs[2];
					new_fc->nodes[0] = this_fc_nds[1];
					new_fc->nodes[1] = this_fc_nds[2];
					new_fc->nodes[2] = un_shared;
					new_fc->elements[0] = el_ct;
					new_fc->elements[1] = max_int;
					new_fc->init_norm_dir(nodes);

					new_fc = &new_el_fcs[3];
					new_fc->nodes[0] = this_fc_nds[2];
					new_fc->nodes[1] = this_fc_nds[0];
					new_fc->nodes[2] = un_shared;
					new_fc->elements[0] = el_ct;
					new_fc->elements[1] = max_int;
					new_fc->init_norm_dir(nodes);

					el_check = check_new_el(elements[el_ct], new_el_fcs);
					if (el_check) {
						elements[el_ct].get_centroid(cent, nodes);
						element_grid.add_ent(el_ct, cent);

						this_fc.elements[1] = el_ct;

						new_el_fcs[1].norm_dir_from_el_cent(cent, nodes);
						faces[fc_ct].copy_data(new_el_fcs[1]);
						fc_added = add_face_if_absent(el_ct);

						new_el_fcs[2].norm_dir_from_el_cent(cent, nodes);
						faces[fc_ct].copy_data(new_el_fcs[2]);
						fc_added = add_face_if_absent(el_ct);

						new_el_fcs[3].norm_dir_from_el_cent(cent, nodes);
						faces[fc_ct].copy_data(new_el_fcs[3]);
						fc_added = add_face_if_absent(el_ct);

						el_ct++;
						return true;
					}
				}
			}
		}
	}

	return false;
}

bool Mesher::adopt_any_nd(int fc_i, double tgt_pt[], double srch_rad) {
	int i1;
	int i2;
	MeshFace& this_fc = faces[fc_i];
	MeshFace* new_fc = nullptr;
	MeshNode* list_nd = nullptr;
	int ndi;
	MeshNode* face_nds[3];
	double* crd;
	double d_vec[3];
	double dist;
	double dp;
	double cent[3];
	int* new_el_nds = &elements[el_ct].nodes[0];
	bool el_check;
	bool fc_added;

	new_el_fcs[0] = this_fc;

	double* norm_dir = &this_fc.norm_dir[0];
	int* this_fc_nds = &this_fc.nodes[0];
	int* this_fc_els = &this_fc.elements[0];
	this_fc.get_centroid(cent, nodes);
	int lst_len = node_grid.get_in_radius(grid_out1,g_olen,tgt_pt,srch_rad);
	for (i1 = 0; i1 < lst_len; i1++) {
		//list_nd = grid_out1[i1]->get_pt(list_nd);
		ndi = grid_out1[i1];
		MeshNode& list_nd = nodes[ndi];
		if (ndi != this_fc_nds[0] && ndi != this_fc_nds[1] && ndi != this_fc_nds[2]) {
			crd = &list_nd.coord[0];
			d_vec[0] = crd[0] - tgt_pt[0];
			d_vec[1] = crd[1] - tgt_pt[1];
			d_vec[2] = crd[2] - tgt_pt[2];
			dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
			if (dist < srch_rad) {
				d_vec[0] = crd[0] - cent[0];
				d_vec[1] = crd[1] - cent[1];
				d_vec[2] = crd[2] - cent[2];
				dp = d_vec[0] * norm_dir[0] + d_vec[1] * norm_dir[1] + d_vec[2] * norm_dir[2];
				if (dp > 1.0e-3*avg_proj) {
					new_el_nds[0] = this_fc_nds[0];
					new_el_nds[1] = this_fc_nds[1];
					new_el_nds[2] = this_fc_nds[2];
					new_el_nds[3] = ndi;

					new_fc = &new_el_fcs[1];
					new_fc->nodes[0] = this_fc_nds[0];
					new_fc->nodes[1] = this_fc_nds[1];
					new_fc->nodes[2] = ndi;
					new_fc->elements[0] = el_ct;
					new_fc->elements[1] = max_int;
					new_fc->init_norm_dir(nodes);

					new_fc = &new_el_fcs[2];
					new_fc->nodes[0] = this_fc_nds[1];
					new_fc->nodes[1] = this_fc_nds[2];
					new_fc->nodes[2] = ndi;
					new_fc->elements[0] = el_ct;
					new_fc->elements[1] = max_int;
					new_fc->init_norm_dir(nodes);

					new_fc = &new_el_fcs[3];
					new_fc->nodes[0] = this_fc_nds[2];
					new_fc->nodes[1] = this_fc_nds[0];
					new_fc->nodes[2] = ndi;
					new_fc->elements[0] = el_ct;
					new_fc->elements[1] = max_int;
					new_fc->init_norm_dir(nodes);

					el_check = check_new_el(elements[el_ct], new_el_fcs);
					if (el_check) {
						//elements.add_ent(new_el);
						elements[el_ct].get_centroid(cent, nodes);
						element_grid.add_ent(el_ct, cent);

						this_fc.elements[1] = el_ct;

						new_el_fcs[1].norm_dir_from_el_cent(cent, nodes);
						faces[fc_ct].copy_data(new_el_fcs[1]);
						fc_added = add_face_if_absent(el_ct);

						new_el_fcs[2].norm_dir_from_el_cent(cent, nodes);
						faces[fc_ct].copy_data(new_el_fcs[2]);
						fc_added = add_face_if_absent(el_ct);

						new_el_fcs[3].norm_dir_from_el_cent(cent, nodes);
						faces[fc_ct].copy_data(new_el_fcs[3]);
						fc_added = add_face_if_absent(el_ct);

						el_ct++;

						return true;
					}
				}
			}
		}
	}

	return false;
}

bool Mesher::create_new_nd(int fc_i, double tgt_pt[]) {
	nodes[nd_ct].coord[0] = tgt_pt[0];
	nodes[nd_ct].coord[1] = tgt_pt[1];
	nodes[nd_ct].coord[2] = tgt_pt[2];
	MeshFace& this_fc = faces[fc_i];
	bool el_check;
	bool fc_added;
	double cent[3];

	int* this_fc_nds = &this_fc.nodes[0];
	int* this_fc_els = &this_fc.elements[0];
	int* new_el_nds = &elements[el_ct].nodes[0];
	new_el_nds[0] = this_fc_nds[0];
	new_el_nds[1] = this_fc_nds[1];
	new_el_nds[2] = this_fc_nds[2];
	new_el_nds[3] = nd_ct;
	elements[el_ct].get_centroid(cent, nodes);

	new_el_fcs[0] = this_fc;

	MeshFace* new_fc = &new_el_fcs[1];
	new_fc->nodes[0] = this_fc_nds[0];
	new_fc->nodes[1] = this_fc_nds[1];
	new_fc->nodes[2] = nd_ct;
	new_fc->elements[0] = el_ct;
	new_fc->elements[1] = max_int;
	new_fc->init_norm_dir(nodes);
	new_fc->norm_dir_from_el_cent(cent, nodes);

	new_fc = &new_el_fcs[2];
	new_fc->nodes[0] = this_fc_nds[1];
	new_fc->nodes[1] = this_fc_nds[2];
	new_fc->nodes[2] = nd_ct;
	new_fc->elements[0] = el_ct;
	new_fc->elements[1] = max_int;
	new_fc->init_norm_dir(nodes);
	new_fc->norm_dir_from_el_cent(cent, nodes);

	new_fc = &new_el_fcs[3];
	new_fc->nodes[0] = this_fc_nds[2];
	new_fc->nodes[1] = this_fc_nds[0];
	new_fc->nodes[2] = nd_ct;
	new_fc->elements[0] = el_ct;
	new_fc->elements[1] = max_int;
	new_fc->init_norm_dir(nodes);
	new_fc->norm_dir_from_el_cent(cent, nodes);

	el_check = check_new_el(elements[el_ct], new_el_fcs);
	if (el_check) {
		node_grid.add_ent(nd_ct, tgt_pt);

		element_grid.add_ent(el_ct, cent);

		this_fc.elements[1] = el_ct;

		faces[fc_ct].copy_data(new_el_fcs[1]);
		fc_added = add_face_if_absent(el_ct);

		faces[fc_ct].copy_data(new_el_fcs[2]);
		fc_added = add_face_if_absent(el_ct);

		faces[fc_ct].copy_data(new_el_fcs[3]);
		fc_added = add_face_if_absent(el_ct);

		nd_ct++;
		el_ct++;

		return true;
	}

	return false;
}

bool Mesher::generate_mesh() {
	int fc_i;
	int* fc_els;
	double fc_cent[3];
	double fc_proj;
	double* fc_norm;
	double proj;
	double tgt_pt[3];
	double srch_rad;
	bool edge_closed;
	
	bool el_added = true;
	while (el_added) {
		el_added = false;
		fc_i = 0;
		while (fc_i < fc_ct && el_ct < max_num_els) {
			MeshFace& this_fc = faces[fc_i];
			fc_els = &this_fc.elements[0];
			if (fc_els[1] == max_int) {
				this_fc.get_centroid(fc_cent, nodes);
				fc_proj = this_fc.proj_dist;
				proj = glob_proj_wt * avg_proj + (1.0 - glob_proj_wt) * fc_proj;
				fc_norm = &this_fc.norm_dir[0];
				tgt_pt[0] = fc_cent[0] + 0.5 * proj * fc_norm[0];
				tgt_pt[1] = fc_cent[1] + 0.5 * proj * fc_norm[1];
				tgt_pt[2] = fc_cent[2] + 0.5 * proj * fc_norm[2];
				srch_rad = 0.75 * proj;
				edge_closed = adopt_connected_nd(fc_i, tgt_pt, srch_rad);
				if (!edge_closed) {
					edge_closed = adopt_any_nd(fc_i, tgt_pt, srch_rad);
				}
				if (!edge_closed) {
					tgt_pt[0] = fc_cent[0] + proj * fc_norm[0];
					tgt_pt[1] = fc_cent[1] + proj * fc_norm[1];
					tgt_pt[2] = fc_cent[2] + proj * fc_norm[2];
					edge_closed = create_new_nd(fc_i, tgt_pt);
				}
				if (!edge_closed) {
					tgt_pt[0] = fc_cent[0] + 0.5 * proj * fc_norm[0];
					tgt_pt[1] = fc_cent[1] + 0.5 * proj * fc_norm[1];
					tgt_pt[2] = fc_cent[2] + 0.5 * proj * fc_norm[2];
					edge_closed = create_new_nd(fc_i, tgt_pt);
				}
				if (edge_closed) {
					el_added = true;
				}
			}
			fc_i++;
		}
	}
	
	if (el_ct >= max_num_els) {
		return false;
	}

	return true;
}

void Mesher::distribute_nodes() {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int i9;
	int dim = nd_ct * 3;
	int* el_nds;
	double e_vol;
	double avg_wt;
	double* crd;
	double res;
	double alpha;
	double beta;
	double r_next;
	double dp;

	vector<double> el_mat(144);
	vector<double> x_vec(dim);
	vector<double> h_vec(dim);
	vector<double> z_vec(dim);
	vector<double> g_vec(dim);
	vector<double> w_vec(dim);
	vector<double> d_mat(dim);
	vector<double> p_mat(dim);
	vector<double> p_inv(dim);
	
	i1 = el_ct;
	vector<double> el_wt(i1);

	// make static element matrix
	for (i1 = 0; i1 < 144; i1++) {
		el_mat[i1] = 0.0;
	}

	i3 = 0;
	for (i1 = 0; i1 < 12; i1++) {
		if (i3 < 9) {
			i4 = 0;
		}
		else {
			i4 = i3 - 9;
		}
		i5 = i1 * 12;
		i6 = (i1 + 1) * 12;
		while (i4 <= (i3 + 9)) {
			if (i4 >= i5 && i4 < i6) {
				el_mat[i4] = -1.0;
			}
			i4 += 3;
		}
		el_mat[i3] = 3.0;
		i3 += 13;
	}

	// determine element weights
	i1 = 0;
	avg_wt = 0.0;
	for (i2 = 0; i2 < el_ct; i2++) {
		MeshElement& this_el = elements[i2];
		e_vol = this_el.get_volume(nodes);
		el_wt[i1] = e_vol;
		avg_wt += e_vol;
		i1++;
	}
	avg_wt /= i1;
	for (i2 = 0; i2 < i1; i2++) {
		el_wt[i2] /= avg_wt;
	}

	// initialize matrices
	for (i1 = 0; i1 < dim; i1++) {
		d_mat[i1] = 0.0;
		p_mat[i1] = 30.0;
		g_vec[i1] = 0.0;
	}

	//this_nd = nodes.get_first();
	i1 = 0;
	for (i1 = 0; i1 < num_bound_nds; i1++) {
		MeshNode& this_nd = nodes[i1];
		crd = &this_nd.coord[0];
		for (i3 = 0; i3 < 3; i3++) {
			i4 = i1 * 3 + i3;
			d_mat[i4] = 100000.0;
			p_mat[i4] += 100000.0;
			g_vec[i4] = -100000.0 * crd[i3];
		}
	}

	for (i1 = 0; i1 < dim; i1++) {
		p_inv[i1] = 1.0 / p_mat[i1];
	}

	// intialize vectors

	res = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		x_vec[i1] = 0.0;
		w_vec[i1] = p_inv[i1] * g_vec[i1];
		h_vec[i1] = -w_vec[i1];
		res += w_vec[i1] * g_vec[i1];
	}

	i1 = 0;
	while (i1 < dim && res > 1.0e-12) {
		for (i2 = 0; i2 < dim; i2++) {
			z_vec[i2] = 0.0;
		}
		i2 = 0;
		for (i2 = 0; i2 < el_ct; i2++) {
			MeshElement& this_el = elements[i2];
			el_nds = &this_el.nodes[0];
			i7 = 0; //index in el_mat
			for (i3 = 0; i3 < 4; i3++) {
				for (i4 = 0; i4 < 3; i4++) {
					i8 = 3*el_nds[i3] + i4; // global row
					for (i5 = 0; i5 < 4; i5++) {
						for (i6 = 0; i6 < 3; i6++) {
							i9 = 3 * el_nds[i5] + i6; //global col
							z_vec[i8] += (el_wt[i2] * el_mat[i7] * h_vec[i9]);
							i7++;
						}
					}
				}
			}
		}
		for (i2 = 0; i2 < dim; i2++) {
			z_vec[i2] += (d_mat[i2] * h_vec[i2]);
		}
		dp = 0.0;
		for (i2 = 0; i2 < dim; i2++) {
			dp += (z_vec[i2] * h_vec[i2]);
		}
		alpha = res / dp;
		r_next = 0.0;
		for (i2 = 0; i2 < dim; i2++) {
			x_vec[i2] += (alpha * h_vec[i2]);
			g_vec[i2] += (alpha * z_vec[i2]);
			w_vec[i2] = p_inv[i2] * g_vec[i2];
			r_next += (g_vec[i2] * w_vec[i2]);
		}
		beta = r_next / res;
		for (i2 = 0; i2 < dim; i2++) {
			h_vec[i2] = -w_vec[i2] + beta * h_vec[i2];
		}
		res = r_next;
		i1++;
	}

	i1 = 0;
	for (i1 = num_bound_nds; i1 < nd_ct; i1++) {
		MeshNode& this_nd = nodes[i1];
		crd = &this_nd.coord[0];
		for (i3 = 0; i3 < 3; i3++) {
			i4 = 3 * i1 + i3;
			crd[i3] = x_vec[i4];
		}
	}

	return;
}

void Mesher::write_output(string file_name) {
	int i1;
	ofstream out_file;
	int* el_nds;
	double* crd;

	out_file.open(file_name);
	
	out_file << "nodes:\n";
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& this_nd = nodes[i1];
		crd = &this_nd.coord[0];
		out_file << "  - [" << crd[0] << ", " << crd[1] << ", " << crd[2] << "]\n";
	}

	out_file << "elements:\n";
	for (i1 = 0; i1 < el_ct; i1++) {
		MeshElement& this_el = elements[i1];
		el_nds = &this_el.nodes[0];
		out_file << "  - [" << el_nds[0];
		for (i1 = 1; i1 < 4; i1++) {
			out_file << ", " << el_nds[i1];
		}
		out_file << "]\n";
	}

	out_file.close();
	return;
}

void Mesher::print_current_mesh() {
	int i1;
	double* crd;
	int* el_nds;
	int* fc_els;
	
	cout << "Current Mesh:" << endl;
	cout << "Nodes:" << endl;
	
	//this_nd = nodes.get_first();
	for (i1 = 0; i1 < nd_ct; i1++) {
		MeshNode& this_nd = nodes[i1];
		crd = &this_nd.coord[0];
		cout << crd[0] << ", " << crd[1] << ", " << crd[2] << endl;
	}
	cout << "Elements:" << endl;
	//this_el = elements.get_first();
	for (i1 = 0; i1 < el_ct; i1++) {
		MeshElement& this_el = elements[i1];
		el_nds = &this_el.nodes[0];
		cout << el_nds[0] << ", ";
		cout << el_nds[1] << ", ";
		cout << el_nds[2] << ", ";
		cout << el_nds[3] << endl;
	}

	cout << "Faces:" << endl;
	//this_fc = faces.get_first();
	for (i1 = 0; i1 < fc_ct; i1++) {
		MeshFace& this_fc = faces[i1];
		cout << "  nodes: ";
		el_nds = &this_fc.nodes[0];
		cout << el_nds[0] << ", ";
		cout << el_nds[1] << ", ";
		cout << el_nds[2] << ", ";
		fc_els = &this_fc.elements[0];
		cout << "element pt: ";
		cout << fc_els[0] << ", " << fc_els[1] << endl;
	}

	return;
}