#include <vector>
#include "FaceClass.h"
#include "DiffDoubClass.h"
#include "NodeClass.h"
#include "DesignVariableClass.h"
#include "matrixFunctions.h"

using namespace std;

Face::Face() {
	num_nds = 0;
	on_surf = true;
	return;
}

void Face::set_node(int place, int loc_nd, int glob_nd) {
	loc_nodes[place] = loc_nd;
	glob_nodes[place] = glob_nd;
	return;
}

void Face::sorted_nodes(int srt_nds[]) {
	int i1;
	int i2;
	int i3;
	int i4;
	int swap;
	for (i1 = 0; i1 < num_nds; i1++) {
		srt_nds[i1] = glob_nodes[i1];
	}
	i3 = num_nds - 1;
	for (i1 = 0; i1 < i3; i1++) {
		for (i2 = 0; i2 < i3; i2++) {
			i4 = i2 + 1;
			if(srt_nds[i4] < srt_nds[i2]) {
				swap = srt_nds[i2];
				srt_nds[i2] = srt_nds[i4];
				srt_nds[i4] = swap;
			}
		}
	}
	return;
}

int Face::get_low_nd() {
	int i1;
	int low_nd = glob_nodes[0];
	for (i1 = 1; i1 < num_nds; i1++) {
		if(glob_nodes[i1] < low_nd) {
			low_nd = glob_nodes[i1];
		}
	}
	return low_nd;
}

//dup1

void Face::get_area_normal_dfd0(DiffDoub0& area, DiffDoub0 norm[], vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub0 v1[3];
	DiffDoub0 v2[3];
	DiffDoub0 tmp_v[3];
	DiffDoub0 tmp;

	if (num_nds == 4) {
		nd_ar[glob_nodes[2]].get_crd_dfd0(v1, dv_ar);
		nd_ar[glob_nodes[0]].get_crd_dfd0(tmp_v, dv_ar);
		v1[0].sub(tmp_v[0]);
		v1[1].sub(tmp_v[1]);
		v1[2].sub(tmp_v[2]);
		nd_ar[glob_nodes[3]].get_crd_dfd0(v2, dv_ar);
		nd_ar[glob_nodes[1]].get_crd_dfd0(tmp_v, dv_ar);
		v2[0].sub(tmp_v[0]);
		v2[1].sub(tmp_v[1]);
		v2[2].sub(tmp_v[2]);
	}
	else {
		nd_ar[glob_nodes[1]].get_crd_dfd0(v1, dv_ar);
		nd_ar[glob_nodes[0]].get_crd_dfd0(tmp_v, dv_ar);
		v1[0].sub(tmp_v[0]);
		v1[1].sub(tmp_v[1]);
		v1[2].sub(tmp_v[2]);
		nd_ar[glob_nodes[2]].get_crd_dfd0(v2, dv_ar);
		nd_ar[glob_nodes[0]].get_crd_dfd0(tmp_v, dv_ar);
		v2[0].sub(tmp_v[0]);
		v2[1].sub(tmp_v[1]);
		v2[2].sub(tmp_v[2]);
	}

	cross_prod_dfd0(norm, v1, v2);
	area.set_val_dfd0(norm[0]);
	area.sqr();
	tmp.set_val_dfd0(norm[1]);
	tmp.sqr();
	area.add(tmp);
	tmp.set_val_dfd0(norm[2]);
	tmp.sqr();
	area.add(tmp);
	area.sqt();

	tmp.set_val(1.0);
	tmp.dvd(area);
	norm[0].mult(tmp);
	norm[1].mult(tmp);
	norm[2].mult(tmp);

	tmp.set_val(0.5);
	area.mult(tmp);

	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

void Face::get_area_normal_dfd1(DiffDoub1& area, DiffDoub1 norm[], vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub1 v1[3];
	DiffDoub1 v2[3];
	DiffDoub1 tmp_v[3];
	DiffDoub1 tmp;

	if (num_nds == 4) {
		nd_ar[glob_nodes[2]].get_crd_dfd1(v1, dv_ar);
		nd_ar[glob_nodes[0]].get_crd_dfd1(tmp_v, dv_ar);
		v1[0].sub(tmp_v[0]);
		v1[1].sub(tmp_v[1]);
		v1[2].sub(tmp_v[2]);
		nd_ar[glob_nodes[3]].get_crd_dfd1(v2, dv_ar);
		nd_ar[glob_nodes[1]].get_crd_dfd1(tmp_v, dv_ar);
		v2[0].sub(tmp_v[0]);
		v2[1].sub(tmp_v[1]);
		v2[2].sub(tmp_v[2]);
	}
	else {
		nd_ar[glob_nodes[1]].get_crd_dfd1(v1, dv_ar);
		nd_ar[glob_nodes[0]].get_crd_dfd1(tmp_v, dv_ar);
		v1[0].sub(tmp_v[0]);
		v1[1].sub(tmp_v[1]);
		v1[2].sub(tmp_v[2]);
		nd_ar[glob_nodes[2]].get_crd_dfd1(v2, dv_ar);
		nd_ar[glob_nodes[0]].get_crd_dfd1(tmp_v, dv_ar);
		v2[0].sub(tmp_v[0]);
		v2[1].sub(tmp_v[1]);
		v2[2].sub(tmp_v[2]);
	}

	cross_prod_dfd1(norm, v1, v2);
	area.set_val_dfd1(norm[0]);
	area.sqr();
	tmp.set_val_dfd1(norm[1]);
	tmp.sqr();
	area.add(tmp);
	tmp.set_val_dfd1(norm[2]);
	tmp.sqr();
	area.add(tmp);
	area.sqt();

	tmp.set_val(1.0);
	tmp.dvd(area);
	norm[0].mult(tmp);
	norm[1].mult(tmp);
	norm[2].mult(tmp);

	tmp.set_val(0.5);
	area.mult(tmp);

	return;
}

//end dup
 
//end skip 
 
 
 
//face_pt_list

FacePtList::FacePtList() {
	return;
}

void FacePtList::add_face(int new_i) {
	fc_list.push_back(new_i);
	return;
}

bool FacePtList::add_if_absent(int new_i, vector<Face>& glob_faces) {
	int i1;
	int new_num_nds;
	int new_srtd[8];
	int this_num_nds;
	int this_srtd[8];
	bool all_match;

	Face& new_fc = glob_faces[new_i];
	new_num_nds = new_fc.num_nds;
	new_fc.sorted_nodes(new_srtd);
	for (auto& fi : fc_list) {
		Face& this_fc = glob_faces[fi];
		this_num_nds = this_fc.num_nds;
		if (this_num_nds == new_num_nds) {
			this_fc.sorted_nodes(this_srtd);
			all_match = true;
			for (i1 = 0; i1 < this_num_nds; i1++) {
				if (this_srtd[i1] != new_srtd[i1]) {
					all_match = false;
				}
			}
			if (all_match) {
				new_fc.on_surf = false;
				this_fc.on_surf = false;
				return false;
			}
		}
	}

	fc_list.push_back(new_i);
	return true;
}