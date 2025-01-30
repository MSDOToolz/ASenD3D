#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "ModelClass.h"
#include "constants.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"
#include "JobClass.h"

using namespace std;


void Model::write_time_step_soln(int t_step) {
	int i1;
	int dof_per_nd;
	int num_int_dof;
	double nd_dat[6];
	double el_dat[9];
	char* buf;
	JobCommand& scmd = job[solve_cmd];
	string full_file = scmd.file_name + "/solnTStep" + to_string(t_step) + ".out";
	ofstream out_file;
	out_file.open(full_file, std::ofstream::binary);

	for (auto& this_nd : nodes) {
		if (scmd.thermal) {
			buf = reinterpret_cast<char*>(&this_nd.prev_temp);
			out_file.write(buf, 8);
			buf = reinterpret_cast<char*>(&this_nd.prev_tdot);
			out_file.write(buf, 8);
		}
		if (scmd.elastic) {
			dof_per_nd = this_nd.num_dof;
			i1 = 8 * dof_per_nd;
			buf = reinterpret_cast<char*>(&this_nd.prev_disp[0]);
			out_file.write(buf, i1);
			buf = reinterpret_cast<char*>(&this_nd.prev_vel[0]);
			out_file.write(buf, i1);
			buf = reinterpret_cast<char*>(&this_nd.prev_acc[0]);
			out_file.write(buf, i1);
		}
	}

	if (scmd.elastic) {
		buf = reinterpret_cast<char*>(&el_dat[0]);
		for (auto& this_el : elements) {
			num_int_dof = this_el.num_int_dof;
			if (num_int_dof > 0) {
				i1 = 8 * num_int_dof;
				vec_to_ar(el_dat, this_el.int_prev_disp, 0, num_int_dof);
				out_file.write(buf, i1);
			}
		}
	}

	out_file.close();
	return;
}

void Model::write_node_results(string file_name, string node_set, list<string>& fields, int time_step) {
	int i1;
	int glob_ind;
	int set_pt;
	string err_str;
	string this_field;
	double nd_dat[6];
	int ndof;
	double time;

	JobCommand& scmd = job[solve_cmd];
	
	if(time_step < max_int) {
		// read the results from the time step file and store them in nodes
		read_time_step_soln(time_step);
		for (auto& nd_pt : nodes) {
			if (scmd.thermal) {
				nd_pt.backstep_temp();
			}
			if (scmd.elastic) {
				nd_pt.backstep_disp();
			}
		}
		if (scmd.dynamic) {
			time = scmd.time_step * time_step;
		}
		else {
			auto this_ld_tm = scmd.static_load_time.begin();
			i1 = 0;
			while (this_ld_tm != scmd.static_load_time.end() && i1 < time_step) {
				++this_ld_tm;
				i1++;
			}
			time = *this_ld_tm;
		}
	}
	else {
		time = -1.0;
	}
	
	if (key_in_map(ns_map,node_set)) {
		set_pt = ns_map.at(node_set);
	}
	else {
		err_str = "Warning: there is no Node Set named " + node_set + ". Defaulting to all nodes in writeNodeResults";
		cout << err_str << endl;
		set_pt = ns_map.at("all");
	}
	Set& this_set = node_sets[set_pt];
	
	ofstream out_file;
	out_file.open(file_name);
	out_file << setprecision(12);
	
	out_file << "nodeResults:\n";
	out_file << "    time: " << time << "\n";
	out_file << "    nodeSet: " << node_set << "\n";
	
	for (auto& this_field : fields) {
		out_file << "    " << this_field << ":\n";
		// -----------------
		// calculate reaction force if necessary
		if (this_field == "reactionForce") {
			for (i1 = 0; i1 < el_mat_dim; i1++) {
				temp_v1[i1] = 0.0;
			}
			build_elastic_soln_load(temp_v1, false, false);
			for (auto& nd_label : this_set.labels) {
				Node& nd_pt = nodes[nd_label];
				out_file << "        - [" << nd_label;
				ndof = nd_pt.num_dof;
				for (i1 = 0; i1 < ndof; i1++) {
					glob_ind = nd_pt.dof_index[i1];
					out_file << ", " << -temp_v1[glob_ind];
				}
				out_file << "]\n";
			}
		}
		else if (this_field == "reactionHeatGen") {
			for (i1 = 0; i1 < el_mat_dim; i1++) {
				temp_v1[i1] = 0.0;
			}
			build_thermal_soln_load(temp_v1, false);
			for (auto& nd_label : this_set.labels) {
				Node& nd_pt = nodes[nd_label];
				out_file << "        - [" << nd_label;
				glob_ind = nd_pt.sorted_rank;
				out_file << ", " << -temp_v1[glob_ind] << "\n";
			}
		}
		else {
			for (auto& nd_label : this_set.labels) {
				Node& nd_pt = nodes[nd_label];
				out_file << "        - [" << nd_label;
				if (this_field == "displacement") {
					ndof = nd_pt.num_dof;
					for (i1 = 0; i1 < ndof; i1++) {
						out_file << ", " << nd_pt.displacement[i1];
					}
					out_file << "]\n";
				}
				else if (this_field == "velocity") {
					ndof = nd_pt.num_dof;
					for (i1 = 0; i1 < ndof; i1++) {
						out_file << ", " << nd_pt.velocity[i1];
					}
					out_file << "]\n";
				}
				else if (this_field == "acceleration") {
					ndof = nd_pt.num_dof;
					for (i1 = 0; i1 < ndof; i1++) {
						out_file << ", " << nd_pt.acceleration[i1];
					}
					out_file << "]\n";
				}
				else if (this_field == "temperature") {
					nd_dat[0] = nd_pt.temperature;
					out_file << ", " << nd_dat[0] << "]\n";
				}
				else if (this_field == "tdot") {
					nd_dat[0] = nd_pt.temp_change_rate;
					out_file << ", " << nd_dat[0] << "]\n";
				}
			}
		}
	}
	
	out_file.close();
	
	return;
}

void Model::write_element_results(string file_name, string el_set, list<string>& fields, string position, int time_step) {
	int i1;
	int i2;
	int i3;
	string err_str;
	string this_field;
	int set_pt;
	int type;
	int num_ip;
	int num_lay;
	double int_pts[24];
	double cent[3] = { 0.0,0.0,0.0 };
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	double seden;
	DiffDoub0 def[9];
	DiffDoub0 frc_mom[9];
	DiffDoub0 t_grad[3];
	DiffDoub0 flux[3];
	string field_list;
	double time;

	JobCommand& scmd = job[solve_cmd];

	if (time_step < max_int) {
		// read the results from the time step file and store them in nodes
		read_time_step_soln(time_step);
		for (auto& nd_pt : nodes) {
			if (scmd.thermal) {
				nd_pt.backstep_temp();
			}
			if (scmd.elastic) {
				nd_pt.backstep_disp();
			}
		}
		if (scmd.elastic) {
			for (auto& el_pt : elements) {
				el_pt.backstep_int_disp();
			}
		}
		if (scmd.dynamic) {
			time = scmd.time_step * time_step;
		}
		else {
			auto this_ld_tm = scmd.static_load_time.begin();
			i1 = 0;
			while (this_ld_tm != scmd.static_load_time.end() && i1 < time_step) {
				++this_ld_tm;
				i1++;
			}
			time = *this_ld_tm;
		}
	}
	else {
		time = -1.0;
	}

	if (key_in_map(es_map,el_set)) {
		set_pt = es_map.at(el_set);
	}
	else {
		err_str = "Warning: there is no Element Set named " + el_set + ". Defaulting to all elements in writeNodeResults";
		cout << err_str << endl;
		set_pt = es_map.at("all");
	}
	Set& this_set = element_sets[set_pt];

	ofstream out_file;
	out_file.open(file_name);
	out_file << setprecision(12);

	out_file << "elementResults:\n";
	out_file << "    time: " << time << "\n";
	out_file << "    elSet: " << el_set << "\n";

	for (auto& this_field : fields) {
		out_file << "    " << this_field << ":\n";
		field_list = "stress strain";
		i2 = field_list.find(this_field);
		if (i2 > -1) {
			out_file << "    ##  - [element label, integration pt, layer, S11, S22, S33, S12, S13, S23]\n";
		}
		field_list = "strainEnergyDen";
		i2 = field_list.find(this_field);
		if (i2 > -1) {
			out_file << "    ##  - [element label, integration pt, layer, strain energy]\n";
		}
		field_list = "sectionDef sectionFrcMom";
		i2 = field_list.find(this_field);
		if (i2 > -1) {
			out_file << "    ## for shells:  - [element label, integration pt, S11, S22, S12, K11, K22, K12]\n";
			out_file << "    ## for beams:  - [element label, integration pt, S11, S12, S13, K11, K12, K13]\n";
		}
		if (this_field == "heatFlux") {
			out_file << "    ##  - [element label, integration pt, layer, f1, f2, f3]\n";
		}
		if (this_field == "tempGradient") {
			out_file << "    ##  - [element label, integration pt, layer, dT/dx, dT/dy, dT/dz]\n";
		}

		for (auto& el_label : this_set.labels) {
			Element& el_pt = elements[el_label];
			type = el_pt.type;
			
			field_list = "stress strain strainEnergyDen";
			i2 = field_list.find(this_field);
			if (i2 > -1) {
				el_pt.get_stress_prereq_dfd0(d0_pre, sections, materials, nodes, design_vars);
				if (position == "intPts") {
					num_ip = el_pt.num_ip;
					vec_to_ar(int_pts, el_pt.int_pts, 0, 3 * num_ip);
				}
				else {
					num_ip = 1;
					int_pts[0] = el_pt.s_cent[0];
					int_pts[1] = el_pt.s_cent[1];
					int_pts[2] = el_pt.s_cent[2];
				}
				for (i1 = 0; i1 < num_ip; i1++) {
					if (type == 3 || type == 41) {
						num_lay = sections[el_pt.sect_ptr].layers.size();
					} else {
						num_lay = 1;
					}
					for (i2 = 0; i2 < num_lay; i2++) {
						el_pt.get_stress_strain_dfd0(stress, strain, &int_pts[3 * i1], i2, scmd.nonlinear_geom, d0_pre);
						out_file << "        - [" << el_label << ", ";
						if (this_field == "strain") {
							out_file << i1 << ", " << i2 << ", " << strain[0].val;
							for (i3 = 1; i3 < 6; i3++) {
								out_file << ", " << strain[i3].val;
							}
							out_file << "]\n";
						} else if(this_field == "stress") {
							out_file << i1 << ", " << i2 << ", " << stress[0].val;
							for (i3 = 1; i3 < 6; i3++) {
								out_file << ", " << stress[i3].val;
							}
							out_file << "]\n";
						} else {
							seden = 0.0;
							for (i3 = 0; i3 < 6; i3++) {
								seden += stress[i3].val * strain[i3].val;
							}
							seden *= 0.5;
							out_file << i1 << ", " << i2 << ", " << seden << "]\n";
						}
					}
				}
			}

			field_list = "sectionDef sectionFrcMom";
			i2 = field_list.find(this_field);
			if (i2 > -1 && el_pt.dof_per_nd == 6) {
				el_pt.get_stress_prereq_dfd0(d0_pre, sections, materials, nodes, design_vars);
				if (position == "intPts") {
					num_ip = el_pt.num_ip;
					vec_to_ar(int_pts, el_pt.int_pts, 0, 3 * num_ip);
				}
				else {
					num_ip = 1;
					int_pts[0] = el_pt.s_cent[0];
					int_pts[1] = el_pt.s_cent[1];
					int_pts[2] = el_pt.s_cent[2];
				}
				for (i1 = 0; i1 < num_ip; i1++) {
					type = el_pt.type;
					//el_pt->get_stress_strain_dfd0(stress, strain, &int_pts[3 * i1], i2, solve_cmd->nonlinear_geom, st_pre);
					el_pt.get_def_frc_mom_dfd0(def, frc_mom, &int_pts[3 * i1], scmd.nonlinear_geom, d0_pre);
					out_file << "        - [" << el_label << ", ";
					if (this_field == "sectionDef") {
						out_file << i1 << ", " << def[0].val;
						for (i3 = 1; i3 < 6; i3++) {
							out_file << ", " << def[i3].val;
						}
						out_file << "]\n";
					}
					else if (this_field == "sectionFrcMom") {
						out_file << i1 << ", " << frc_mom[0].val;
						for (i3 = 1; i3 < 6; i3++) {
							out_file << ", " << frc_mom[i3].val;
						}
						out_file << "]\n";
					}
				}
			}

			field_list = "tempGradient heatFlux";
			i2 = field_list.find(this_field);
			if (i2 > -1) {
				el_pt.get_stress_prereq_dfd0(d0_pre, sections, materials, nodes, design_vars);
				if (position == "intPts") {
					num_ip = el_pt.num_ip;
					vec_to_ar(int_pts, el_pt.int_pts, 0, 3 * num_ip);
				}
				else {
					num_ip = 1;
					int_pts[0] = el_pt.s_cent[0];
					int_pts[1] = el_pt.s_cent[1];
					int_pts[2] = el_pt.s_cent[2];
				}
				for (i1 = 0; i1 < num_ip; i1++) {
					type = el_pt.type;
					if (type == 3 || type == 41) {
						num_lay = sections[el_pt.sect_ptr].layers.size();
					}
					else {
						num_lay = 1;
					}
					for (i2 = 0; i2 < num_lay; i2++) {
						//el_pt->get_stress_strain_dfd0(stress, strain, &int_pts[3 * i1], i2, solve_cmd->nonlinear_geom, st_pre);
						el_pt.get_flux_tgrad_dfd0(flux, t_grad, &int_pts[3 * i1], i2, d0_pre);
						out_file << "        - [" << el_label << ", " << i1 << ", " << i2 << ", ";
						if (this_field == "tempGradient") {
							out_file << t_grad[0].val << ", " << t_grad[1].val << ", " << t_grad[2].val << "]\n";
						}
						else if (this_field == "heatFlux") {
							out_file << flux[0].val << ", " << flux[1].val << ", " << flux[2].val << "]\n";
						}
					}
				}
			}

		}
	}

	out_file.close();

	return;
}

void Model::write_modal_results(string file_name, bool write_modes) {
	int i1;
	int i2;
	int i3;
	int nd;
	int dof_per_nd;
	int glob_ind;
	JobCommand& mcmd = job[modal_cmd];
	int n_mds = mcmd.num_modes;
	Node* this_nd;

	ofstream out_file;
	out_file.open(file_name);
	out_file << setprecision(12);

	out_file << "modalResults:\n";
	out_file << "    eigenValues:\n";
	for (i1 = 0; i1 < n_mds; i1++) {
		out_file << "      - " << eig_vals[i1] << "\n";
	}
	if (mcmd.type == "buckling") {
		out_file << "    loadFactors:\n";
	}
	else {
		out_file << "    frequencies:\n";
	}
	for (i1 = 0; i1 < n_mds; i1++) {
		out_file << "      - " << load_fact[i1] << "\n";
	}

	if (mcmd.write_modes) {
		out_file << "    modes:\n";
		for (i1 = 0; i1 < n_mds; i1++) {
			out_file << "      - mode: " << i1 << "\n";
			out_file << "        displacement:\n";
			for (auto& this_nd : nodes) {
				nd = this_nd.label;
				out_file << "          - [" << nd;
				dof_per_nd = this_nd.num_dof;
				for (i2 = 0; i2 < dof_per_nd; i2++) {
					glob_ind = this_nd.dof_index[i2];
					i3 = i1 * el_mat_dim + glob_ind;
					out_file << ", " << eig_vecs[i3];
				}
				out_file << "]\n";
			}
		}
	}

	out_file.close();

	return;
}

void Model::write_objective(string file_name, list<string>& include_fields, bool write_grad) {
	int i1;
	int num_dv;
	double tot_obj;
	double this_val;
	string fld_str;

	ofstream out_file;
	out_file.open(file_name);
	out_file << setprecision(12);

	out_file << "objective:\n";
	out_file << "    terms:\n";
	tot_obj = 0.0;
	for (auto& this_term : obj.terms) {
		this_val = this_term.value;
		out_file << "        - value: " << this_val << "\n";
		for (auto& fld_str : include_fields) {
			if (fld_str == "category") {
				out_file << "          category: " << this_term.category << "\n";
			}
			else if (fld_str == "operator") {
				out_file << "          operator: " << this_term.optr << "\n";
			}
			else if (fld_str == "component") {
				out_file << "          component: " << this_term.component << "\n";
			}
			else if (fld_str == "layer") {
				out_file << "          Layer: " << this_term.layer << "\n";
			}
			else if (fld_str == "coefficient") {
				out_file << "          coefficient: " << this_term.coef << "\n";
			}
			else if (fld_str == "exponent") {
				out_file << "          exponent: " << this_term.expnt << "\n";
			}
			else if (fld_str == "elementSet") {
				out_file << "          elementSet: " << this_term.el_set_name << "\n";
			}
			else if (fld_str == "nodeSet") {
				out_file << "          nodeSet: " << this_term.nd_set_name << "\n";
			}
			else if (fld_str == "activeTime") {
				out_file << "          activeTime: [" << this_term.active_time[0] << ", " << this_term.active_time[1] << "]\n";
			}
		}
		tot_obj += this_val;
	}
	out_file << "    totalValue: " << tot_obj << "\n";
	if (write_grad) {
		out_file << "objectiveGradient:\n";
		num_dv = design_vars.size();
		for (i1 = 0; i1 < num_dv; i1++) {
			out_file << "    - [" << i1 << ", " << d_ld_d[i1] << "]\n";
		}
	}

	out_file.close();

	return;
}