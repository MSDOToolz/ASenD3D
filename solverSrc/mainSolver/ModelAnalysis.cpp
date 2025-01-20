#include <cmath>
#include <iostream>
#include <iomanip>
#include "ModelClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "LowerTriMatClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"
#include "FaceClass.h"
#include "DiffDoubClass.h"
#include "LoadClass.h"
#include "matrixFunctions.h"

using namespace std;

const int max_int = 2000000000;

void Model::reorder_nodes(int block_dim) {
	int i1;
	int i2;
	int i3;
	int nd1;
	int nd2;
	int min_ct;
	int min_nd;
	double min_dist;
	int since_restart;
	double dist;
	int el_num_nds;
	int el_dof_per_nd;
	int num_nodes = nodes.size();
	
	vector<Set> nodal_conn(num_nodes);
	vector<int> node_inserted(num_nodes);
	
	// build nodal connectivity
	for (auto& el : elements) {
		el_num_nds = el.num_nds;
		el_dof_per_nd = el.dof_per_nd;
		for (i1 = 0; i1 < el_num_nds; i1++) {
			nd1 = el.nodes[i1];
			if (el_dof_per_nd > 3) {
				nodes[nd1].num_dof = el_dof_per_nd;
			}
			for (i2 = i1+1; i2 < el_num_nds; i2++) {
				nd2 = el.nodes[i2];
				nodal_conn[nd1].add_if_absent(nd2);
				nodal_conn[nd2].add_if_absent(nd1);
			}
		}
	}

	for (auto& this_const : elastic_const.const_vec) {
		for (auto& term1 : this_const.terms) {
			list<int>& t1_labs = node_sets[term1.ns_ptr].labels;
			i1 = t1_labs.size();
			for(auto& term2 : this_const.terms) {
				list<int>& t2_labs = node_sets[term2.ns_ptr].labels;
				i2 = t2_labs.size();
				if (i1 > 1 && i2 > 1) {
					auto iter1 = t1_labs.begin();
					auto end1 = t1_labs.end();
					auto iter2 = t2_labs.begin();
					auto end2 = t2_labs.end();
					while (iter1 != end1 && iter2 != end2) {
						nodal_conn[*iter1].add_if_absent(*iter2);
						nodal_conn[*iter2].add_if_absent(*iter1);
						++iter1;
						++iter2;
					}
				}
				else {
					for (auto& nd1 : t1_labs) {
						for (auto& nd2 : t2_labs) {
							nodal_conn[nd1].add_if_absent(nd2);
							nodal_conn[nd2].add_if_absent(nd1);
						}
					}
				}
			}
		}
	}

	for (auto& this_const : thermal_const.const_vec) {
		for (auto& term1 : this_const.terms) {
			list<int>& t1_labs = node_sets[term1.ns_ptr].labels;
			for (auto& nd1 : t1_labs) {
				for (auto& term2 : this_const.terms) {
					list<int>& t2_labs = node_sets[term2.ns_ptr].labels;
					for (auto& nd2 : t2_labs) {
						nodal_conn[nd1].add_if_absent(nd2);
						nodal_conn[nd2].add_if_absent(nd1);
					}
				}
			}
		}
	}
	
	// find Node with least connectivity
	min_ct = num_nodes;
	min_nd = 0;
	for (i1 = 0; i1 < num_nodes; i1++) {
		i2 = nodal_conn[i1].labels.size();
		if(i2 < min_ct) {
			min_nd = i1;
			min_ct = i2;
		}
		node_inserted[i1] = 0;
	}
	
	// put nodes into an integer list in level order.
	list<int> ordered_nds;
	
	ordered_nds.push_back(min_nd);
	auto this_nd = ordered_nds.begin();
	node_inserted[min_nd] = 1;
	since_restart = 0;
	
	while(ordered_nds.size() < num_nodes) {
		nd1 = *this_nd;
		auto neighbor_nd = nodal_conn[nd1].labels.begin();
		while(neighbor_nd != nodal_conn[nd1].labels.end() && since_restart < block_dim) {
			nd2 = *neighbor_nd;
			if(node_inserted[nd2] == 0) {
				ordered_nds.push_back(nd2);
				node_inserted[nd2] = 1;
				since_restart++;
			}
			++neighbor_nd;
		}
		if(since_restart >= block_dim || this_nd == ordered_nds.end()) {
			this_nd = ordered_nds.end();
			--this_nd;
			nd1 = *this_nd;
			min_dist = 1.0e+100;
			for (i1 = 0; i1 < num_nodes; i1++) {
				if(node_inserted[i1] == 0) {
					dist = get_dist(nodes[nd1].coord, nodes[i1].coord);
					if(dist < min_dist) {
						min_dist = dist;
						min_nd = i1;
					}
				}
			}
			ordered_nds.push_back(min_nd);
			node_inserted[min_nd] = 1;
			since_restart++;
			if(since_restart >= block_dim) {
				since_restart = 0;
			}
		}
		++this_nd;
	}
	
	// update the global degree of freedom indexes for the nodes
	
    i2 = 0; // index in elastic matrix
	i3 = 0; // sorted rank for solid nodes
	for (auto& ndi : ordered_nds) {
		Node& this_nd = nodes[ndi];
		if (!this_nd.fluid) {
			this_nd.sorted_rank = i3;
			i3++;
			this_nd.dof_index[0] = i2;
			i2++;
			this_nd.dof_index[1] = i2;
			i2++;
			this_nd.dof_index[2] = i2;
			i2++;
			if (this_nd.num_dof == 6) {
				this_nd.dof_index[3] = i2;
				i2++;
				this_nd.dof_index[4] = i2;
				i2++;
				this_nd.dof_index[5] = i2;
				i2++;
			}
		}
	}
	el_mat_dim = i2;
	elastic_mat.set_dim(el_mat_dim);
	non_frc_el_mat.set_dim(el_mat_dim);
	therm_mat.set_dim(nodes.size());

	JobCommand& scmd = job[solve_cmd];
	if (scmd.solver_method == "iterative" && scmd.max_it == 0) {
		scmd.max_it = el_mat_dim;
	}

	for (auto& this_el : elements) {
		i3 = this_el.num_int_dof;
		if (i3 > 0) {
			this_el.int_dof_index = i2;
			i2 += i3;
		}
	}
	tot_glob_dof = i2;
	
	temp_v1 = vector<double>(el_mat_dim);
	temp_v2 = vector<double>(el_mat_dim);
	temp_v3 = vector<double>(el_mat_dim);
	temp_v4 = vector<double>(el_mat_dim);
	temp_v5 = vector<double>(el_mat_dim);
	temp_v6 = vector<double>(el_mat_dim);
	temp_d1 = vector<DiffDoub0>(el_mat_dim);

	i3 = nodes.size();
	d_ld_u = vector<double>(tot_glob_dof);
	d_ld_v = vector<double>(el_mat_dim);
	d_ld_a = vector<double>(el_mat_dim);
	d_ld_t = vector<double>(i3);
	d_ld_tdot = vector<double>(i3);
	u_adj = vector<double>(el_mat_dim);
	v_adj = vector<double>(el_mat_dim);
	a_adj = vector<double>(el_mat_dim);
	t_adj = vector<double>(i3);
	tdot_adj = vector<double>(i3);

	d_rud_d = vector<DiffDoub1>(el_mat_dim);
	d_rtd_d = vector<DiffDoub1>(i3);

	el_in_d = vector<int>(elements.size());

	i3 = design_vars.size();
	if (i3 > 0) {
		d_ld_d = vector<double>(i3);
	}

	return;
}

void Model::build_constraint_mats() {
	for (auto& this_const : elastic_const.const_vec) {
		this_const.build_mat(nodes,node_sets);
	}
	for (auto& this_const : thermal_const.const_vec) {
		this_const.build_mat(nodes,node_sets);
	}
	return;
}

void Model::update_reference() {
	int i1;
	int i2;
	int i3;
	// Set the Section pointers for all elements and Material pointers for all sections
	string el_set;
	Set *es_ptr;
	string mat_name;
	i2 = 0;
	for (auto& this_sec : sections) {
		el_set = this_sec.el_set_name;
		i1 = es_map.at(el_set);
		for (auto& eli : element_sets[i1].labels) {
			elements[eli].sect_ptr = i2;
		}
		i3 = 0;
		for (auto& this_mat : materials) {
			mat_name = this_mat.name;
			if(mat_name == this_sec.mat_name) {
				this_sec.mat_ptr = i3;
			}
			for (auto& this_lay : this_sec.layers) {
				if(mat_name == this_lay.mat_name) {
					this_lay.mat_ptr = i3;
				}
			}
			i3++;
		}
		i2++;
	}
	
	// Set Node & Element Set pointers in loads, constraints, design variables and objectives
	string nd_set;
	for (auto& this_load : elastic_loads) {
		nd_set = this_load.node_set;
		if (key_in_map(ns_map,nd_set)) {
			i1 = ns_map.at(nd_set);
			this_load.nd_set_ptr = i1;
		}
		else {
			el_set = this_load.element_set;
			i1 = es_map.at(el_set);
			this_load.el_set_ptr = i1;
		}
	}
	for (auto& this_load : thermal_loads) {
		nd_set = this_load.node_set;
		if (key_in_map(ns_map,nd_set)) {
			i1 = ns_map.at(nd_set);
			this_load.nd_set_ptr = i1;
		}
		else {
			el_set = this_load.element_set;
			i1 = es_map.at(el_set);
			this_load.el_set_ptr = i1;
		}
	}

	for (auto& this_const : elastic_const.const_vec) {
		for (auto& this_cterm : this_const.terms) {
			nd_set = this_cterm.node_set;
			i1 = ns_map.at(nd_set);
			this_cterm.ns_ptr = i1;
		}
	}
	for (auto& this_const : thermal_const.const_vec) {
		for (auto& this_cterm : this_const.terms) {
			nd_set = this_cterm.node_set;
			i1 = ns_map.at(nd_set);
			this_cterm.ns_ptr = i1;
		}
	}

	for (auto& this_dv : design_vars) {
		nd_set = this_dv.nd_set_name;
		if (key_in_map(ns_map,nd_set)) {
			i1 = ns_map.at(nd_set);
			this_dv.nd_set_ptr = i1;
		}
		else {
			el_set = this_dv.el_set_name;
			i1 = es_map.at(el_set);
			this_dv.el_set_ptr = i1;
		}
	}

	for (auto& this_term : obj.terms) {
		el_set = this_term.el_set_name;
		if (key_in_map(es_map,el_set)) {
			i1 = es_map.at(el_set);
			this_term.el_set_ptr = i1;
		}
		else {
			nd_set = this_term.nd_set_name;
			i1 = ns_map.at(nd_set);
			this_term.nd_set_ptr = i1;
		}
	}
	
	// build dv reference list for nodes and elements
	int coef_len;
	double const_coef;
	int dvi = 0;
	for (auto& this_dv : design_vars) {
		coef_len = this_dv.coefs.size();
		if(this_dv.el_set_ptr < max_int) {
			if(coef_len < 2) {
				if(coef_len == 0) {
					const_coef = 1.0;
				} else {
					const_coef = *this_dv.coefs.begin();
				}
				for (auto& eli : element_sets[this_dv.el_set_ptr].labels) {
					elements[eli].add_design_variable(dvi,const_coef);
				}
			}
			else {
				list<int>& set_labs = element_sets[this_dv.el_set_ptr].labels;
				auto set_iter = set_labs.begin();
				auto set_end = set_labs.end();
				auto coef_iter = this_dv.coefs.begin();
				auto coef_end = this_dv.coefs.end();
				while (set_iter != set_end && coef_iter != coef_end) {
					elements[*set_iter].add_design_variable(dvi, *coef_iter);
					++set_iter;
					++coef_iter;
				}
		    }
		}
		
		if (this_dv.nd_set_ptr < max_int) {
			if (coef_len < 2) {
				if (coef_len == 0) {
					const_coef = 1.0;
				}
				else {
					const_coef = *this_dv.coefs.begin();
				}
				for (auto& ndi : node_sets[this_dv.nd_set_ptr].labels) {
					nodes[ndi].add_design_variable(dvi, const_coef);
				}
			}
			else {
				list<int>& set_labs = node_sets[this_dv.nd_set_ptr].labels;
				auto set_iter = set_labs.begin();
				auto set_end = set_labs.end();
				auto coef_iter = this_dv.coefs.begin();
				auto coef_end = this_dv.coefs.end();
				while (set_iter != set_end && coef_iter != coef_end) {
					nodes[*set_iter].add_design_variable(dvi, *coef_iter);
					++set_iter;
					++coef_iter;
				}
			}
		}
		dvi++;
	}
	
	// build comprehensive Element list for each design variable
	
	int el_label;
	int el_num_nds;
	for (auto& this_el : elements) {
		el_label = this_el.label;
		for (auto& dv : this_el.design_vars) {
			design_vars[dv.int_dat].add_comp_el(el_label);
			this_el.add_comp_dvar(dv.int_dat);
		}
		el_num_nds = this_el.num_nds;
		for (i1 = 0; i1 < el_num_nds; i1++) {
			Node& this_nd = nodes[this_el.nodes[i1]];
			for (auto& dv : this_nd.d_var_lst) {
				DesignVariable& this_dv = design_vars[dv.int_dat];
				if (this_dv.category == "nodeCoord") {
					this_dv.add_comp_el(el_label);
					this_el.add_comp_dvar(dv.int_dat);
				}
			}
		}
	}
	
	return;
}

void Model::find_surface_faces() {
	int i1;
	int low_nd;
	bool added;

	int fc_ct = 0;
	for (auto& el : elements) {
		fc_ct += el.num_faces;
	}

	faces = vector<Face>(fc_ct);

	i1 = 0;
	for (auto& el : elements) {
		el.initialize_faces(faces,i1);
	}
	
	int num_nodes = nodes.size();
	vector<FacePtList> f_larray(num_nodes);
	for (auto& this_el : elements) {
		if (this_el.dof_per_nd == 3) {
			for (auto& this_fc : this_el.faces) {
				low_nd = faces[this_fc].get_low_nd();
				added = f_larray[low_nd].add_if_absent(this_fc,faces);
			}
		}
	}
	
	return;
}

void Model::prep_matrix_factorizations() {
	double zero_ar[9] = { 0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };

	JobCommand& scmd = job[solve_cmd];
	if (scmd.thermal) {
		for (auto& this_nd : nodes) {
			this_nd.initialize_temp();
			if (scmd.dynamic) {
				this_nd.update_tdot(scmd.newmark_gamma, scmd.time_step);
			}
		}
		if (!therm_lt.is_allocated()) {
			build_thermal_soln_load(temp_v1, true);
			therm_lt.allocate_from_sparse_mat(therm_mat, thermal_const, scmd.solver_block_dim);
		}
		if (!therm_scaled) {
			scale_thermal_const();
		}
	}

	if (scmd.elastic) {
		for (auto& this_nd : nodes) {
			this_nd.initialize_disp();
			if (scmd.dynamic) {
				this_nd.update_vel_acc(scmd.newmark_beta, scmd.newmark_gamma, scmd.time_step);
			}
		}
		for (auto& this_el : elements) {
			this_el.set_int_disp(zero_ar);
			this_el.set_int_prev_disp(zero_ar);
		}
		if (!elastic_lt.is_allocated()) {
			build_elastic_soln_load(temp_v1, true, true);
			scale_elastic_const();
			elastic_lt.allocate_from_sparse_mat(elastic_mat, elastic_const, 6 * scmd.solver_block_dim);
		}
	}

	return;
}

void Model::analysis_prep() {
	int i1;
	int i2;
	int num_nds;
	int block_dim;
	double time;
	double t_inc;

	if (solve_cmd < max_int) {
		JobCommand& scmd = job[solve_cmd];
		if (scmd.solver_method == "direct") {
			scmd.solver_block_dim = 2000000000;
		}
		else {
			i1 = 6;
			num_nds = nodes.size();
			while ((i1 * i1) < num_nds) {
				i1 += 6;
			}
			if (scmd.solver_block_dim == 2000000000) {
				scmd.solver_block_dim = i1;
			}
		}
		block_dim = scmd.solver_block_dim;

		if (scmd.static_load_time.size() == 0) {
			scmd.static_load_time.push_back(0.0);
		}
	}
	else {
		block_dim = 2000000000;
	}

	i1 = 0;
	for (auto& sec : sections) {
		i2 = sec.layers.size();
		if (i2 > i1) {
			i1 = i2;
		}
	}
	if (i1 > 0) {
		d0_pre.allocate_layers(i1);
		d1_pre.allocate_layers(i1);
	}

	update_reference();
	reorder_nodes(block_dim);
	build_constraint_mats();
	prep_matrix_factorizations();
	find_surface_faces();

	an_prep_run = true;
	return;
}

void Model::build_elastic_app_load(vector<double>& app_ld, double time) {
	// construct loads from input file
	int i1;
	int num_dof;
	int dof_ind;
	string ld_type;
	DiffDoub0 nd_dvld[6];

	for (i1 = 0; i1 < el_mat_dim; i1++) {
		temp_d1[i1].set_val(0.0);
	}

	//loads from Model input file
	for (auto& this_load : elastic_loads) {
		ld_type = this_load.type;
		if(time >= this_load.active_time[0] && time <= this_load.active_time[1]) {
			if(ld_type == "nodalForce") {
				for (auto& ndi : node_sets[this_load.nd_set_ptr].labels) {
					Node& this_nd = nodes[ndi];
					num_dof = this_nd.num_dof;
					for (i1 = 0; i1 < num_dof; i1++) {
						dof_ind = this_nd.dof_index[i1];
						app_ld[dof_ind] += this_load.load[i1];
					}
				}
			}
			else {
				for (auto& eli : element_sets[this_load.el_set_ptr].labels) {
					Element& this_el = elements[eli];
					this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
					JobCommand& scmd = job[solve_cmd];
					this_el.get_app_load(temp_d1, this_load, scmd.nonlinear_geom, d0_pre, sections, faces, nodes, design_vars);
				}
			}
		}
	}

	for (i1 = 0; i1 < el_mat_dim; i1++) {
		app_ld[i1] += temp_d1[i1].val;
	}
	
	// design variable dependent loads.
	for (auto& this_nd : nodes) {
		this_nd.get_elastic_dvload(nd_dvld,design_vars);
		num_dof = this_nd.num_dof;
		for (i1 = 0; i1 < num_dof; i1++) {
			dof_ind = this_nd.dof_index[i1];
			app_ld[dof_ind] += nd_dvld[i1].val;
		}
	}
	
	return;
}

void Model::build_thermal_app_load(vector<double>& app_ld, double time) {
	// construct loads from input file
	int i1;
	int tot_nodes;
	int num_dof;
	int dof_ind;
	Node* this_nd;
	Element* this_el;
	string ld_type;
	double act_time[2];
	double nd_load[6];
	DiffDoub0 nd_dvld;
	Set* this_set;

	tot_nodes = nodes.size();
	for (i1 = 0; i1 < tot_nodes; i1++) {
		temp_d1[i1].set_val(0.0);
	}

	//loads from Model input file
	for (auto& this_load : thermal_loads) {
		ld_type = this_load.type;
		if (time >= this_load.active_time[0] && time <= this_load.active_time[1]) {
			if (ld_type == "nodalHeatGen") {
				for (auto& ndi : node_sets[this_load.nd_set_ptr].labels) {
					Node& this_nd = nodes[ndi];
					dof_ind = this_nd.sorted_rank;
					app_ld[dof_ind] += nd_load[0];
				}
			}
			else {
				for (auto& eli : element_sets[this_load.el_set_ptr].labels) {
					Element& this_el = elements[eli];
					this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
					this_el.get_app_therm_load(temp_d1, this_load, d0_pre, sections, faces, nodes, design_vars);
				}
			}
		}
	}

	for (i1 = 0; i1 < tot_nodes; i1++) {
		app_ld[i1] += temp_d1[i1].val;
	}

	// design variable dependent loads.
	for (auto& this_nd : nodes) {
		this_nd.get_thermal_dvload(nd_dvld, design_vars);
		dof_ind = this_nd.sorted_rank;
		app_ld[dof_ind] += nd_dvld.val;
	}

	return;
}

void Model::build_elastic_soln_load(vector<double>& soln_ld, bool build_mat, bool full_ref) {
	int i1;
	
	for (i1 = 0; i1 < el_mat_dim; i1++) {
		temp_d1[i1].set_val(0.0);
	}
	if(build_mat) {
		elastic_mat.zero_all();
	}
	if (full_ref) {
		non_frc_el_mat.zero_all();
	}
	
	JobCommand& scmd = job[solve_cmd];
	for (auto& this_el : elements) {
		this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
		if (this_el.type == 21) {
			this_el.get_ru(temp_d1, elastic_mat, build_mat, scmd, d0_pre, nodes, design_vars);
		}
		else {
			if (full_ref) {
				this_el.get_ru(temp_d1, non_frc_el_mat, true, scmd, d0_pre, nodes, design_vars);
			}
			else {
				this_el.get_ru(temp_d1, non_frc_el_mat, false, scmd, d0_pre, nodes, design_vars);
			}
		}
	}
	
	if (build_mat) {
		elastic_mat.add_matrix(non_frc_el_mat);
	}
	
	for (i1 = 0; i1 < el_mat_dim; i1++) {
		soln_ld[i1]-= temp_d1[i1].val;
	}
	
	return;
}

void Model::build_thermal_soln_load(vector<double>& soln_ld, bool build_mat) {
	int i1;
	int num_nodes;

	num_nodes = nodes.size();
	for (i1 = 0; i1 < num_nodes; i1++) {
		temp_d1[i1].set_val(0.0);
	}

	if (build_mat) {
		therm_mat.zero_all();
	}

	JobCommand& scmd = job[solve_cmd];
	for (auto& this_el : elements) {
		this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
		this_el.get_rt(temp_d1, therm_mat, build_mat, scmd, d0_pre, nodes);
	}

	for (i1 = 0; i1 < num_nodes; i1++) {
		soln_ld[i1] -= temp_d1[i1].val;
	}

	return;
}

void Model::scale_elastic_const() {
	double scale_fact = 100000.0*elastic_mat.get_max_abs_val();
	elastic_const.set_scale_fact(scale_fact);
	elastic_scaled = true;
	return;
}

void Model::scale_thermal_const() {
	double scale_fact = 100000.0 * therm_mat.get_max_abs_val();
	thermal_const.set_scale_fact(scale_fact);
	therm_scaled = true;
	return;
}

void Model::build_elastic_const_load(vector<double>& const_ld) {
	int i1;
	int nd_dof;
	int glob_ind;
	for (auto& this_nd : nodes) {
		nd_dof = this_nd.num_dof;
		for (i1 = 0; i1 < nd_dof; i1++) {
			glob_ind = this_nd.dof_index[i1];
			temp_v2[glob_ind] = this_nd.displacement[i1];
		}
	}
	elastic_const.get_total_load(const_ld, temp_v2, temp_v3, el_mat_dim);
	return;
}

void Model::build_thermal_const_load(vector<double>& const_ld) {
	int i1;
	int num_nodes;
	double nd_temp;
	int glob_ind;
	for (auto& this_nd : nodes) {
		nd_temp = this_nd.temperature;
		glob_ind = this_nd.sorted_rank;
		temp_v2[glob_ind] = nd_temp;
	}
	num_nodes = nodes.size();
	thermal_const.get_total_load(const_ld, temp_v2, temp_v3, num_nodes);
	return;
}

void Model::solve_step(double time, double app_ld_fact, bool full_ref) {
	int i1;
	int i2;
	int num_nodes;
	int ideb;
	int max_nlit;
	double d_unorm;
	double d_utol;
	double absd_u;
	int ndof;
	int dof_ind;
	double nd_del_disp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	Node *this_nd;
	Element *this_el;

	JobCommand& cmd = job[solve_cmd];
	if(cmd.thermal) {
		num_nodes = nodes.size();
		for (i1 = 0; i1 < num_nodes; i1++) {
			temp_v1[i1] = 0.0;
			temp_v2[11] = 0.0;
		}
		build_thermal_app_load(temp_v1, time);
		build_thermal_soln_load(temp_v1, false);
		if (!therm_scaled) {
			scale_thermal_const();
		}
		build_thermal_const_load(temp_v1);
		if (cmd.solver_method == "direct") {
			therm_lt.ldl_solve(temp_v2, temp_v1);
		}
		else {
			conj_grad_sparse(temp_v2, therm_mat, thermal_const, therm_lt, temp_v1, cmd.conv_tol, cmd.max_it);
			//g_mres_sparse(temp_v2, *therm_mat, *thermal_const, *therm_lt, temp_v1, cmd->conv_tol, cmd->max_it, cmd->solver_block_dim);
		}
		for (auto& this_nd : nodes) {
			dof_ind = this_nd.sorted_rank;
			this_nd.temperature += temp_v2[dof_ind];
			if (cmd.dynamic) {
				this_nd.update_tdot(cmd.newmark_gamma, cmd.time_step);
			}
		}
	}
	
	if(cmd.elastic) {
		if(cmd.nonlinear_geom) {
			max_nlit = 50;
		} else {
			max_nlit = 1;
		}
		d_utol = 1.0e-12;
		d_unorm = 1.0;
		i2 = 0;
		while(i2 < max_nlit && d_unorm > d_utol) {
			for (i1 = 0; i1 < el_mat_dim; i1++) {
				temp_v1[i1] = 0.0;
				temp_v2[i1] = 0.0;
			}
			build_elastic_app_load(temp_v1,time);
			for (i1 = 0; i1 < el_mat_dim; i1++) {
				temp_v1[i1] *= app_ld_fact;
			}
			build_elastic_soln_load(temp_v1,cmd.nonlinear_geom,full_ref);
			if (!elastic_scaled) {
				scale_elastic_const();
			}
			build_elastic_const_load(temp_v1);
			for (auto& this_el : elements) {
				if(this_el.num_int_dof > 0) {
				    this_el.update_external(temp_v1,1,nodes,d0_pre.scr_mat1,d0_pre.scr_mat2);
				}
			}
			if (cmd.nonlinear_geom) {
				elastic_lt.populate_from_sparse_mat(elastic_mat, elastic_const);
				elastic_lt.ldl_factor();
			}
			if(cmd.solver_method == "direct") {
				elastic_lt.ldl_solve(temp_v2,temp_v1);
			}
			else {
				conj_grad_sparse(temp_v2, elastic_mat, elastic_const, elastic_lt, temp_v1, cmd.conv_tol, cmd.max_it);
				//g_mres_sparse(temp_v2, *elastic_mat, *elastic_const, *elastic_lt, temp_v1, cmd->conv_tol, cmd->max_it, 6*cmd->solver_block_dim);
			}
			for (auto& this_nd : nodes) {
				ndof = this_nd.num_dof;
				for (i1 = 0; i1 < ndof; i1++) {
					dof_ind = this_nd.dof_index[i1];
					nd_del_disp[i1] = temp_v2[dof_ind];
				}
				this_nd.add_to_displacement(nd_del_disp);
				if(cmd.dynamic) {
					this_nd.update_vel_acc(cmd.newmark_beta,cmd.newmark_gamma,cmd.time_step);
				}
			}
			for (auto& this_el : elements) {
				if(this_el.num_int_dof > 0) {
				    this_el.update_internal(temp_v2,1,nodes,d0_pre.scr_mat1,d0_pre.scr_mat2);
				}
			}
			d_unorm = 0.0;
			for (i1 = 0; i1 < el_mat_dim; i1++) {
				absd_u = abs(temp_v2[i1]);
				if(absd_u > d_unorm) {
					d_unorm = absd_u;
				}
			}
			if (cmd.nonlinear_geom) {
				cout << "Nonlinear iteration: " << i2 << ", max solution step: " << d_unorm << endl;
			}
			i2++;
		}
	}

	
	return;
}

void Model::solve() {
	JobCommand& cmd = job[solve_cmd];
	int i1;
	int i2;
	int i3;
	int since_ref;
	double app_ld_fact;
	double time;
	int ld_steps = cmd.load_ramp_steps;
	double c1 = 1.0/(ld_steps*ld_steps);
	double zero_ar[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	if(!an_prep_run) {
		analysis_prep();
	}

	if(cmd.thermal) {
		therm_lt.populate_from_sparse_mat(therm_mat, thermal_const);
		therm_lt.ldl_factor();
	}
	
	if(cmd.elastic && !cmd.nonlinear_geom) {
		elastic_lt.populate_from_sparse_mat(elastic_mat, elastic_const);
		cout << "factoring stiffness matrix" << endl;
		elastic_lt.ldl_factor();
		cout << "finished factoring" << endl;
	}
	
	if(cmd.dynamic) {
		if (cmd.save_soln_hist) {
			write_time_step_soln(0);
		}
		time = 0.0;
		i1 = 1;
		since_ref = 0;
		while (time < cmd.sim_period) {
			since_ref++;
			if (since_ref == cmd.full_reform) {
				solve_step(time, 1.0, true);
				since_ref = 0;
			}
			else {
				solve_step(time, 1.0, false);
			}
			for (auto& this_nd : nodes) {
				if (cmd.thermal) {
					this_nd.advance_temp();
					this_nd.update_tdot(cmd.newmark_gamma,cmd.time_step);
				}
				if (cmd.elastic) {
					this_nd.advance_disp();
					this_nd.update_vel_acc(cmd.newmark_beta,cmd.newmark_gamma,cmd.time_step);
				}
			}
			if (cmd.elastic) {
				for (auto& this_el : elements) {
					this_el.advance_int_disp();
				}
			}
			if (cmd.save_soln_hist) {
				write_time_step_soln(i1);
				time_steps_saved = i1;
				i1++;
			}
			time += cmd.time_step;
		}
	} else {
		i3 = 0;
		for (auto& this_ld : cmd.static_load_time) {
			if (cmd.nonlinear_geom) {
				for (i1 = 0; i1 < ld_steps; i1++) {
					i2 = ld_steps - i1 - 1;
					app_ld_fact = 1.0 - c1 * i2 * i2;
					solve_step(this_ld, app_ld_fact, true);
				}
			}
			else {
				app_ld_fact = 1.0;
				solve_step(this_ld, app_ld_fact, true);
			}
			if (cmd.save_soln_hist) {
				for (auto& this_nd : nodes) {
					if (cmd.thermal) {
						this_nd.advance_temp();
					}
					if (cmd.elastic) {
						this_nd.advance_disp();
					}
				}
				if (cmd.elastic) {
					for (auto& this_el : elements) {
						this_el.advance_int_disp();
					}
				}
				write_time_step_soln(i3);
				time_steps_saved = i3 + 1;
			}
			i3++;
		}
	}
	
	return;
}

void Model::zero_solution(list<string>& fields) {
	double zero_ar[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	int n_int_dof;

	bool disp_inc = false;
	for (auto& this_field : fields) {
		if (this_field == "displacement") {
			disp_inc = true;
		}
	}

	for (auto& this_nd : nodes) {
		for (auto& this_field : fields) {
			if (this_field == "temperature") {
				this_nd.temperature = 0.0;
			}
			else if (this_field == "tdot") {
				this_nd.temp_change_rate = 0.0;
			}
			else if (this_field == "displacement") {
				this_nd.set_displacement(zero_ar);
			}
			else if (this_field == "velocity") {
				this_nd.set_velocity(zero_ar);
			}
			else if (this_field == "acceleration") {
				this_nd.set_acceleration(zero_ar);
			}
		}
	}

	if (disp_inc) {
		for (auto& this_el : elements) {
			n_int_dof = this_el.num_int_dof;
			if (n_int_dof > 0) {
				this_el.set_int_disp(zero_ar);
			}
		}
	}

	return;
}

void Model::eigen_solve() {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int n_dof;
	int glob_ind;
	int step_inc;
	double nd_disp[6];
	vector<DiffDoub0> rvec(33);
	vector<double> rv2(33);
	vector<double> d_rd_a(1089);
	double zeros[9] = {0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};
	double shft;
	double v_kv;
	double dp;
	double m_avg;
	double m_min;

	JobCommand& cmd = job[modal_cmd];
	JobCommand& scmd = job[solve_cmd];

	if (!an_prep_run) {
		analysis_prep();
	}

	if (eig_vecs.size() == 0) {
		i1 = cmd.num_modes;
		eig_vals = vector<double>(i1);
		load_fact = vector<double>(i1);
		i1 *= el_mat_dim;
		eig_vecs = vector<double>(i1);
		diag_mass = vector<double>(el_mat_dim);
	}

	if (cmd.type == "buckling") {
		for (i1 = 0; i1 < el_mat_dim; i1++) {
			diag_mass[i1] = 1.0;
		}
	}
	else {
		for (i1 = 0; i1 < el_mat_dim; i1++) {
			diag_mass[i1] = 0.0;
		}
		for (auto& this_el : elements) {
			this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
			i2 = this_el.num_nds * this_el.dof_per_nd;
			for (i1 = 0; i1 < i2; i1++) {
				d0_pre.glob_acc[i1].set_val(1.0);
			}
			this_el.get_rum(rvec, d_rd_a, false, true, scmd.nonlinear_geom, d0_pre, nodes, design_vars);
			for (i1 = 0; i1 < i2; i1++) {
				rv2[i1] = rvec[i1].val;
			}
			this_el.add_to_glob_vec(rv2, diag_mass, false, false, nodes);
		}
	}
	
	m_avg = 0.0;
	for (i1 = 0; i1 < el_mat_dim; i1++) {
		m_avg += abs(diag_mass[i1]);
	}
	m_avg /= el_mat_dim;
	m_min = 1.0e-6 * m_avg;
	for (i1 = 0; i1 < el_mat_dim; i1++) {
		if (diag_mass[i1] < m_min) {
			cout << "Warning: small or negative value in diagonal mass matrix: " << diag_mass[i1] << endl;
			cout << "Adjusting to minimum value" << m_min << endl;
			diag_mass[i1] = m_min;
		}
	}

	cmd.nonlinear_geom = scmd.nonlinear_geom;
	cmd.dynamic = scmd.dynamic;
	scmd.nonlinear_geom = true;
	scmd.dynamic = false;
	cout << "building stiffness matrix" << endl;
	build_elastic_soln_load(temp_v1, true, true);
	cout << "finished building matrix" << endl;
	if (cmd.type == "buckling" || cmd.type == "frequency") {
		for (i1 = 0; i1 < el_mat_dim; i1++) {
			shft = -cmd.tgt_eval * diag_mass[i1];
			elastic_mat.add_entry(i1, i1, shft);
		}
		if (cmd.solver_method == "direct") {
			elastic_lt.populate_from_sparse_mat(elastic_mat, elastic_const);
			cout << "factoring stiffness matrix" << endl;
			elastic_lt.ldl_factor();
			cout << "finished factoring.  beginning eigensolve" << endl;
			eigen_sparse_direct(eig_vals, eig_vecs, cmd.num_modes, elastic_lt, diag_mass, el_mat_dim);
			cout << "finished eigensolve" << endl;
		}
		else {

		}
		for (i1 = 0; i1 < cmd.num_modes; i1++) {
			eig_vals[i1] += cmd.tgt_eval;
		}
	}
	else {
		//step_inc = time_steps_saved / cmd.num_modes;
		//if (step_inc == 0) {
		//	string er_str = "Error: not enough solution time steps have been saved to find the requested number of active eigenmodes.\n";
		//	er_str += "Be sure to Set the saveSolnHist option to true in the dynamic solve command.\n";
		//	throw invalid_argument(er_str);
		//}
		//for (i1 = 0; i1 < cmd.num_modes; i1++) {
		//	i2 = (i1 + 1) * step_inc;
		//	i3 = i1 * el_mat_dim;
		//	read_time_step_soln(i2);
		//	this_nd = nodes.get_first();
		//	while (this_nd) {
		//		this_nd->get_prev_disp(nd_disp);
		//		n_dof = this_nd->get_num_dof();
		//		for (i4 = 0; i4 < n_dof; i4++) {
		//			glob_ind = this_nd->get_dof_index(i4);
		//			eig_vecs[i3 + glob_ind] = nd_disp[i4];
		//		}
		//		this_nd = this_nd->get_next();
		//	}
		//}
		////get_nearest_evec_rq(*elastic_mat, *elastic_const, diag_mass, eig_vecs, eig_vals, cmd->num_modes, cmd->max_it);
		//get_nearest_evec_subspace(elastic_mat, elastic_const, diag_mass, eig_vecs, eig_vals, cmd->num_modes);
	}

	if (cmd.type == "buckling") {
		write_time_step_soln(-1);
		for (auto& this_nd : nodes) {
			this_nd.advance_disp();
			this_nd.set_displacement(zeros);
		}
		for (auto& this_el : elements) {
			this_el.advance_int_disp();
			this_el.set_int_disp(zeros);
		}
		build_elastic_soln_load(temp_v1, true, true);
		i2 = 0;
		for (i1 = 0; i1 < cmd.num_modes; i1++) {
			for (i3 = 0; i3 < el_mat_dim; i3++) {
				temp_v1[i3] = 0.0;
			}
			sub_vec(temp_v2, eig_vecs, i2, i2 + el_mat_dim);
			elastic_mat.vector_multiply(temp_v1, temp_v2, false);
			v_kv = 0.0;
			for (i3 = 0; i3 < el_mat_dim; i3++) {
				v_kv += eig_vecs[i2 + i3] * temp_v1[i3];
			}
			if ((v_kv - eig_vals[i1]) != 0.0) {
				load_fact[i1] = v_kv / (v_kv - eig_vals[i1]);
			}
			else {
				load_fact[i1] = 1.0e+100;
			}
			i2 += el_mat_dim;
		}
		for (auto& this_nd : nodes) {
			this_nd.backstep_disp();
		}
		for (auto& this_el : elements) {
			this_el.backstep_int_disp();
		}
		read_time_step_soln(-1);
	}
	else {
		for (i1 = 0; i1 < cmd.num_modes; i1++) {
			if (eig_vals[i1] >= 0.0) {
				load_fact[i1] = 0.159154943091895335 * sqrt(eig_vals[i1]);
			}
			else {
				load_fact[i1] = -1.0;
			}
		}
	}

	scmd.nonlinear_geom = cmd.nonlinear_geom;
	scmd.dynamic = cmd.dynamic;

	return;
}


void Model::set_soln_to_mode(string field, int mode, double max_val) {
	int i1;
	int i2;
	int mode_st;
	int n_dof;
	int glob_ind;
	string disp_fields = "displacement velocity acceleration";
	string er_str;
	double nd_dat[6];
	double scale_fact;
	double ab_val;
	
	mode_st = mode * el_mat_dim;
	i2 = mode_st;
	scale_fact = 0.0;
	for (i1 = 0; i1 < el_mat_dim; i1++) {
		ab_val = abs(eig_vecs[i2]);
		if (ab_val > scale_fact) {
			scale_fact = ab_val;
		}
		i2++;
	}
	scale_fact = max_val / scale_fact;

	i1 = disp_fields.find(field);
	if (i1 > -1) {
	    for (auto& this_nd : nodes) {
			n_dof = this_nd.num_dof;
			for (i2 = 0; i2 < n_dof; i2++) {
				glob_ind = this_nd.dof_index[i2];
				nd_dat[i2] = scale_fact * eig_vecs[mode_st + glob_ind];
			}
			if (field == "displacement") {
				this_nd.set_displacement(nd_dat);
			}
			else if (field == "velocity") {
				this_nd.set_velocity(nd_dat);
			}
			else {
				this_nd.set_acceleration(nd_dat);
			}
		}
	}
	else {
		er_str = "Warning: " + field + " is not a valid field to Set to a mode. Aborting setSolnToMode()";
		cout << er_str << endl;
	}

	return;
}

void Model::augmentd_ld_u() {
	int i1;
	int i2;
	int i3;
	int num_nds;
	int nd_dof;
	JobCommand& scmd = job[solve_cmd];
	double dt = scmd.time_step;
	double bet = scmd.newmark_beta;
	double gam = scmd.newmark_gamma;
	double d_tdotd_tp = -1.0 / (gam * dt);
	double d_tdotd_tdotp = -(1.0 - gam) / gam;
	double d_ad_up = 1.0 / (dt * dt * (bet - gam));
	double d_ad_vp = d_ad_up * dt;
	double d_ad_ap = d_ad_up * dt * dt * (0.5 + bet - gam);
	double d_vd_up = dt * gam * d_ad_up;
	double d_vd_vp = 1.0 + dt * gam * d_ad_vp;
	double d_vd_ap = dt * (1.0 - gam) + dt * gam * d_ad_ap;

	double c1 = d_ad_up;
	double c2 = d_ad_ap;
	double c3 = dt * (1.0 - gam);
	double c4 = 1.0 / (dt*(bet - gam));

	vector<DiffDoub0> rvec(30);
	vector<double> mmat(900);
	vector<double> dmat(900);
	vector<double> el_adj(30);
	vector<double> eld_ld_t(10);
	vector<double> eld_ld_tdot(10);
	vector<double> eld_ld_u(30);
	vector<double> eld_ld_v(30);
	vector<double> eld_ld_a(30);

	Element* this_el;

	if (scmd.thermal) {
		for (auto& this_el : elements) {
			this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
			this_el.get_rtm(rvec, mmat, true, true, d0_pre);
			this_el.get_el_vec(el_adj, t_adj, true, false, nodes);
			num_nds = this_el.num_nds;
			i3 = 0;
			for (i1 = 0; i1 < num_nds; i1++) {
				eld_ld_t[i1] = 0.0;
				eld_ld_tdot[i1] = 0.0;
				for (i2 = 0; i2 < num_nds; i2++) {
					eld_ld_t[i1] -= d_tdotd_tp * mmat[i3] * el_adj[i2];
					eld_ld_tdot[i1] -= d_tdotd_tdotp * mmat[i3] * el_adj[i2];
					i3++;
				}
			}
			this_el.add_to_glob_vec(eld_ld_t, d_ld_t, true, false, nodes);
			this_el.add_to_glob_vec(eld_ld_tdot, d_ld_tdot, true, false, nodes);
		}
		i2 = nodes.size();
		for (i1 = 0; i1 < i2; i1++) {
			d_ld_t[i1] += tdot_adj[i1];
			d_ld_tdot[i1] += c3 * tdot_adj[i1];
		}
	}

	if (scmd.elastic) {
		for (auto& this_el : elements) {
			this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
			this_el.get_rum(rvec, mmat, true, true, scmd.nonlinear_geom, d0_pre, nodes, design_vars);
			this_el.get_rud(rvec, dmat, true, scmd, d0_pre, nodes, design_vars);
			this_el.get_el_vec(el_adj, u_adj, false, false, nodes);
			nd_dof = this_el.num_nds * this_el.dof_per_nd;
			i3 = 0;
			for (i1 = 0; i1 < nd_dof; i1++) {
				eld_ld_u[i1] = 0.0;
				eld_ld_v[i1] = 0.0;
				eld_ld_a[i1] = 0.0;
				for (i2 = 0; i2 < nd_dof; i2++) {
					eld_ld_u[i1] -= (d_ad_up * mmat[i3] + d_vd_up * dmat[i3]) * el_adj[i2];
					eld_ld_a[i1] -= (d_ad_ap * mmat[i3] + d_vd_ap * dmat[i3]) * el_adj[i2];
					eld_ld_v[i1] -= (d_ad_vp * mmat[i3] + d_vd_vp * dmat[i3]) * el_adj[i2];
					i3++;
				}
			}
			this_el.add_to_glob_vec(eld_ld_u, d_ld_u, false, false, nodes);
			this_el.add_to_glob_vec(eld_ld_a, d_ld_a, false, false, nodes);
			this_el.add_to_glob_vec(eld_ld_v, d_ld_v, false, false, nodes);
		}
		for (i1 = 0; i1 < el_mat_dim; i1++) {
			d_ld_u[i1] += c1 * a_adj[i1];
			d_ld_a[i1] += (c2 * a_adj[i1] + c3 * v_adj[i1]);
			d_ld_v[i1] += c4 * a_adj[i1] + v_adj[i1];
		}
	}

	return;
}

void Model::solve_for_adjoint(double time, bool full_ref) {
	int i1;
	int i2;
	int i3;
	int tot_nodes;
	int el_num_nds;
	int el_dof_per_nd;
	int el_nd_dof;
	int el_int_dof;
	double c1;
	double c2;
	double c3;
	Element* this_el;
	vector<DiffDoub0> rvec(33);
	vector<double> d_rd_u(1089);
	vector<double> d_rd_t(330);
	vector<double> eld_ld_t(10);
	vector<double> el_adj(33);
	double* int_adj;
	vector<double> scr1(33);
	vector<double> scr2(33);
	JobCommand& scmd = job[solve_cmd];

	obj.calculated_ld_u(d_ld_u, d_ld_v, d_ld_a, d_ld_t, d_ld_tdot, time, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d0_pre);
    
	if (scmd.elastic) {
		if (scmd.dynamic) {
			c1 = scmd.time_step * scmd.newmark_gamma;
			c2 = scmd.time_step;
			c2 = 1.0 / (c2 * c2 * (scmd.newmark_beta - scmd.newmark_gamma));
			for (i1 = 0; i1 < el_mat_dim; i1++) {
				v_adj[i1] = d_ld_v[i1];
				a_adj[i1] = (d_ld_a[i1] + c1 * v_adj[i1]);
				d_ld_u[i1] -= c2 * a_adj[i1];
			}
		}
		if (scmd.nonlinear_geom) {
			build_elastic_soln_load(temp_v1, true, full_ref);
			elastic_lt.populate_from_sparse_mat(elastic_mat, elastic_const);
			elastic_lt.ldl_factor();
		}
		for (auto& this_el : elements) {
			this_el.set_intd_ld_u(d_ld_u);
			this_el.update_external(d_ld_u, 0, nodes, scr1, scr2);
		}
		if (scmd.solver_method == "direct") {
			elastic_lt.ldl_solve(u_adj, d_ld_u);
		}
		else {
			conj_grad_sparse(u_adj, elastic_mat, elastic_const, elastic_lt, d_ld_u, scmd.conv_tol, scmd.max_it);
			//g_mres_sparse(u_adj, *elastic_mat, *elastic_const, *elastic_lt, d_ld_u, solve_cmd->conv_tol, solve_cmd->max_it, 6*solve_cmd->solver_block_dim);
		}
		for (auto& this_el : elements) {
			this_el.update_internal(u_adj, 0, nodes, scr1, scr2);
		}
	}

	if (scmd.thermal) {
		if (scmd.elastic) {
			for (auto& this_el : elements) {
				this_el.get_stress_prereq(d0_pre, sections, materials, nodes, design_vars);
				this_el.get_ruk(rvec, d_rd_u, d_rd_t, true, scmd.nonlinear_geom, d0_pre, nodes, design_vars);
				this_el.get_el_vec(el_adj, u_adj, false, false, nodes);
				el_num_nds = this_el.num_nds;
				el_dof_per_nd = this_el.dof_per_nd;
				el_int_dof = this_el.num_int_dof;
				el_nd_dof = el_num_nds * el_dof_per_nd;
				for (i1 = 0; i1 < el_num_nds; i1++) {
					eld_ld_t[i1] = 0.0;
					i3 = i1;
					for (i2 = 0; i2 < el_nd_dof; i2++) {
						//i3 = i2 * el_num_nds + i1;
						eld_ld_t[i1] -= d_rd_t[i3] * el_adj[i2];
						i3 += el_num_nds;
					}
					for (i2 = 0; i2 < el_int_dof; i2++) {
						eld_ld_t[i1] -= d_rd_t[i3] * el_adj[i2];
						i3 += el_num_nds;
					}
				}
				this_el.add_to_glob_vec(eld_ld_t, d_ld_t, true, false, nodes);
			}
		}
		tot_nodes = nodes.size();
		if (scmd.dynamic) {
			c3 = -1.0 / (scmd.time_step * scmd.newmark_gamma);
			for (i1 = 0; i1 < tot_nodes; i1++) {
				tdot_adj[i1] = c3 * d_ld_tdot[i1];
				d_ld_t[i1] -= tdot_adj[i1];
			}
		}
		if (scmd.solver_method == "direct") {
			therm_lt.ldl_solve(t_adj, d_ld_t);
		}
		else {
			conj_grad_sparse(t_adj, therm_mat, thermal_const, therm_lt, d_ld_t, scmd.conv_tol, scmd.max_it);
			//g_mres_sparse(t_adj, *therm_mat, *thermal_const, *therm_lt, d_ld_t, solve_cmd->conv_tol, solve_cmd->max_it, solve_cmd->solver_block_dim);
		}
	}
	
	return;
}

void Model::d_rthermald_d(int d_var_num) {
	int i1;
	int tot_nodes;
	int num_dof;
	int glob_ind;
	DiffDoub0 dv_val;

	JobCommand& scmd = job[solve_cmd];
	DesignVariable& this_dv = design_vars[d_var_num];
	this_dv.get_value(dv_val);
	this_dv.diff_val.set_val(dv_val.val,1.0);

	tot_nodes = nodes.size();
	for (i1 = 0; i1 < tot_nodes; i1++) {
		d_rtd_d[i1].set_val(0.0);
	}

	for (i1 = 0; i1 < elements.size(); i1++) {
		el_in_d[i1] = 0;
	}

	for (auto& eli : this_dv.comp_el_list) {
		el_in_d[eli] = 1;
	}
	
	// applied contribution
	for (auto& this_ld : thermal_loads) {
		if (this_ld.type != "nodalHeatGen") {
			for (auto& eli : element_sets[this_ld.el_set_ptr].labels) {
				if (el_in_d[eli]) {
					Element& this_el = elements[eli];
					this_el.get_stress_prereq(d1_pre, sections, materials, nodes, design_vars);
					this_el.get_app_therm_load(d_rtd_d, this_ld, d1_pre, sections, faces, nodes, design_vars);
				}
			}
		}
	}

	for (i1 = 0; i1 < tot_nodes; i1++) {
		d_rtd_d[i1].neg();
	}

	// solution-dependent contribution of load
	for (auto& eli : this_dv.comp_el_list) {
		Element& this_el = elements[eli];
		this_el.get_stress_prereq(d1_pre, sections, materials, nodes, design_vars);
		this_el.get_rt(d_rtd_d,therm_mat,false,scmd,d1_pre,nodes);
	}


	// design variable dependent contribution
	if (this_dv.nd_set_ptr < max_int) {
		DiffDoub1 nd_ld;
		for (auto& ndi : node_sets[this_dv.nd_set_ptr].labels) {
			Node& this_nd = nodes[ndi];
			this_nd.get_thermal_dvload(nd_ld, design_vars);
			glob_ind = this_nd.sorted_rank;
			d_rtd_d[glob_ind].sub(nd_ld);
		}
	}

	this_dv.diff_val.set_val(dv_val.val, 0.0);
	
	return;
}

void Model::d_relasticd_d(int d_var_num) {
	int i1;
	int num_dof;
	int glob_ind;
	DiffDoub0 dv_val;

	JobCommand& scmd = job[solve_cmd];
	DesignVariable& this_dv = design_vars[d_var_num];
	this_dv.get_value(dv_val);
	this_dv.diff_val.set_val(dv_val.val, 1.0);

	for (i1 = 0; i1 < el_mat_dim; i1++) {
		d_rud_d[i1].set_val(0.0);
	}

	for (i1 = 0; i1 < elements.size(); i1++) {
		el_in_d[i1] = 0;
	}

	for (auto& eli : this_dv.comp_el_list) {
		el_in_d[eli] = 1;
	}
	
	// applied contribution

    for (auto& this_ld : elastic_loads) {
		if (this_ld.type != "nodalForce") {
			for (auto eli : element_sets[this_ld.el_set_ptr].labels) {
				if (el_in_d[eli] > 0) {
					Element& this_el = elements[eli];
					this_el.get_stress_prereq(d1_pre, sections, materials, nodes, design_vars);
					this_el.get_app_load(d_rud_d, this_ld, scmd.nonlinear_geom, d1_pre, sections, faces, nodes, design_vars);
				}
			}
		}
	}

	for (i1 = 0; i1 < el_mat_dim; i1++) {
		d_rud_d[i1].neg();
	}

	// solution-dependent contribution of load
	for (auto& eli : this_dv.comp_el_list) {
		Element& this_el = elements[eli];
		this_el.get_stress_prereq(d1_pre, sections, materials, nodes, design_vars);
		this_el.get_ru(d_rud_d, elastic_mat, false, scmd, d1_pre, nodes, design_vars);
	}


	// design variable dependent contribution

	if (this_dv.nd_set_ptr < max_int) {
		DiffDoub1 nd_ld[6];
		for (auto& ndi : node_sets[this_dv.nd_set_ptr].labels) {
			Node& this_nd = nodes[ndi];
			this_nd.get_elastic_dvload(nd_ld, design_vars);
			num_dof = this_nd.num_dof;
			for (i1 = 0; i1 < num_dof; i1++) {
				glob_ind = this_nd.dof_index[i1];
				nd_ld[i1].neg();
				d_rud_d[glob_ind].add(nd_ld[i1]);
			}
		}
	}

	this_dv.diff_val.set_val(dv_val.val, 0.0);

	return;
}

void Model::get_objective() {
	int i1;
	double time;
	Node* this_nd;
	Element* this_el;
	string er_str;

	JobCommand& scmd = job[solve_cmd];
	if (!scmd.save_soln_hist) {
		er_str = "Error: The solution history must be saved to calculate the Objective gradient.\n";
		er_str += "Save history by setting saveSolnHist=True in the solve command.\n";
		throw invalid_argument(er_str);
	}

	if (scmd.dynamic) {
		obj.clear_values();
		read_time_step_soln(time_steps_saved);
		i1 = time_steps_saved - 1;
		time = scmd.time_step * time_steps_saved;
		while (i1 >= 0) {
			for (auto& this_nd : nodes) {
				if (scmd.elastic) {
					this_nd.backstep_disp();
				}
				if (scmd.thermal) {
					this_nd.backstep_temp();
				}
			}
			if (scmd.elastic) {
				for (auto& this_el : elements) {
					this_el.backstep_int_disp();
				}
			}
			read_time_step_soln(i1);
			obj.calculate_terms(time, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d0_pre);
			i1--;
			time -= scmd.time_step;
		}
	}
	else {
		obj.clear_values();
		i1 = 0;
		for (auto& this_ld_tm : scmd.static_load_time) {
			read_time_step_soln(i1);
			for (auto& this_nd : nodes) {
				if (scmd.elastic) {
					this_nd.backstep_disp();
				}
				if (scmd.thermal) {
					this_nd.backstep_temp();
				}
			}
			if (scmd.elastic) {
				for (auto& this_el : elements) {
					this_el.backstep_int_disp();
				}
			}
			obj.calculate_terms(this_ld_tm, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d0_pre);
			i1++;
		}
	}
	return;
}

void Model::get_obj_gradient() {
	int i1;
	int i2;
	int i3;
	int i4;
	double time;
	int since_ref;
	int num_dv = design_vars.size();
	Element* this_el;
	Node* this_nd;
	string er_str;

	JobCommand& scmd = job[solve_cmd];
	if (!scmd.save_soln_hist) {
		er_str = "Error: The solution history must be saved to calculate the Objective gradient.\n";
		er_str += "Save history by setting saveSolnHist=True in the solve command.\n";
		throw invalid_argument(er_str);
	}

	obj.clear_values();
	for (i1 = 0; i1 < num_dv; i1++) {
		d_ld_d[i1] = 0.0;
	}

	if (scmd.dynamic) {
		time = scmd.time_step * time_steps_saved;
		i1 = time_steps_saved;
		since_ref = 0;
		read_time_step_soln(i1);
		while (i1 > 0) {
			for (i2 = 0; i2 < nodes.size(); i2++) {
				d_ld_t[i2] = 0.0;
				d_ld_tdot[i2] = 0.0;
			}
			for (i2 = 0; i2 < el_mat_dim; i2++) {
				d_ld_a[i2] = 0.0;
				d_ld_v[i2] = 0.0;
			}
			for (i2 = 0; i2 < tot_glob_dof; i2++) {
				d_ld_u[i2] = 0.0;
			}
			if (i1 < time_steps_saved) {
				augmentd_ld_u();
			}
			for (auto& this_nd : nodes) {
				if (scmd.elastic) {
					this_nd.backstep_disp();
				}
				if (scmd.thermal) {
					this_nd.backstep_temp();
				}
			}
			for (auto& this_el : elements) {
				this_el.backstep_int_disp();
			}
			read_time_step_soln(i1 - 1);
			obj.calculate_terms(time, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d0_pre);
			since_ref++;
			if (since_ref == scmd.full_reform) {
				solve_for_adjoint(time, true);
				since_ref = 0;
			}
			else {
				solve_for_adjoint(time, false);
			}
			obj.calculated_ld_d(d_ld_d, time, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d1_pre);
			for (i2 = 0; i2 < num_dv; i2++) {
				if (scmd.thermal) {
					d_rthermald_d(i2);
					for (i3 = 0; i3 < nodes.size(); i3++) {
						d_ld_d[i2] -= t_adj[i3] * d_rtd_d[i3].dval;
					}
				}
				if (scmd.elastic) {
					d_relasticd_d(i2);
					for (i3 = 0; i3 < el_mat_dim; i3++) {
						d_ld_d[i2] -= u_adj[i3] * d_rud_d[i3].dval;
					}
					for (auto& this_el : elements) {
						d_ld_d[i2] -= this_el.get_int_adjd_rd_d();
					}
				}
			}
			i1--;
		}
	}
	else {
		i4 = 0;
		for (auto& this_ld : scmd.static_load_time) {
			read_time_step_soln(i4);
			for (auto& this_nd : nodes) {
				if (scmd.elastic) {
					this_nd.backstep_disp();
				}
				if (scmd.thermal) {
					this_nd.backstep_temp();
				}
			}
			for (auto& this_el : elements) {
				this_el.backstep_int_disp();
			}
			obj.calculate_terms(this_ld, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d0_pre);
			for (i2 = 0; i2 < nodes.size(); i2++) {
				d_ld_t[i2] = 0.0;
				d_ld_tdot[i2] = 0.0;
			}
			for (i2 = 0; i2 < el_mat_dim; i2++) {
				d_ld_a[i2] = 0.0;
				d_ld_v[i2] = 0.0;
			}
			for (i2 = 0; i2 < tot_glob_dof; i2++) {
				d_ld_u[i2] = 0.0;
			}
			solve_for_adjoint(this_ld,true);
			obj.calculated_ld_d(d_ld_d, this_ld, scmd.nonlinear_geom, nodes, elements, node_sets, element_sets, sections, materials, design_vars, d1_pre);
			for (i1 = 0; i1 < num_dv; i1++) {
				if (scmd.thermal) {
					d_rthermald_d(i1);
					for (i3 = 0; i3 < nodes.size(); i3++) {
						d_ld_d[i1] -= t_adj[i3] * d_rtd_d[i3].dval;
					}
				}
				if (scmd.elastic) {
					d_relasticd_d(i1);
					for (i2 = 0; i2 < el_mat_dim; i2++) {
						d_ld_d[i1] -= u_adj[i2] * d_rud_d[i2].dval;
					}
					for (auto& this_el : elements) {
						d_ld_d[i1] -= this_el.get_int_adjd_rd_d();
					}
				}
			}
			i4++;
		}
	}
	return;
}