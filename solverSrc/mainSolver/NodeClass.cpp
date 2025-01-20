#include "NodeClass.h"
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;

Node::Node() {
	label = 0;
	fluid = false;
	num_dof = 3;
	sorted_rank = 0;
	coord[0] = 0.0;
	coord[1] = 0.0;
	coord[2] = 0.0;
	int i1;
	for (i1=0; i1 < 6; i1++) {
		displacement[i1] = 0.0;
		velocity[i1] = 0.0;
		acceleration[i1] = 0.0;
		prev_disp[i1] = 0.0;
		prev_vel[i1] = 0.0;
		prev_acc[i1] = 0.0;
		initial_disp[i1] = 0.0;
		initial_vel[i1] = 0.0;
		initial_acc[i1] = 0.0;
		dof_index[i1] = 0;
	}
	for (i1 = 0; i1 < 3; i1++) {
		fl_vel[i1] = 0.0;
		fl_vel_dot[i1] = 0.0;
		prev_fl_vel[i1] = 0.0;
		p_fl_vel_lf[i1] = 0.0;
		prev_fl_vel_dot[i1] = 0.0;
		initial_fl_vel[i1] = 0.0;
		initial_fl_vel_dot[i1] = 0.0;
	}
	temperature = 0.0;
	temp_change_rate = 0.0;
	fl_den = 0.0;
	fl_den_dot = 0.0;
	turb_e = 0.0;
	turb_edot = 0.0;
	prev_temp = 0.0;
	p_temp_lf = 0.0;
	prev_tdot = 0.0;
	prev_fl_den = 0.0;
	p_fl_den_lf = 0.0;
	prev_fl_den_dot = 0.0;
	prev_turb_e = 0.0;
	prev_turb_edot = 0.0;
	p_turb_elf = 0.0;
	initial_temp = 0.0;
	initial_tdot = 0.0;
	initial_fl_den = 0.0;
	initial_fl_den_dot = 0.0;
	initial_turb_e = 0.0;
	initial_turb_edot = 0.0;
	d_var_lst.clear();
}

void Node::set_crd(double new_crd[]) {
	coord[0] = new_crd[0];
	coord[1] = new_crd[1];
	coord[2] = new_crd[2];
	return;
}

void Node::set_displacement(double new_disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] = new_disp[i1];
	}
	return;
}

void Node::set_velocity(double new_vel[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		velocity[i1] = new_vel[i1];
	}
	return;
}

void Node::set_acceleration(double new_acc[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		acceleration[i1] = new_acc[i1];
	}
	return;
}

void Node::set_fl_vel(double new_vel[]) {
	fl_vel[0] = new_vel[0];
	fl_vel[1] = new_vel[1];
	fl_vel[2] = new_vel[2];
	return;
}

void Node::set_fl_vel_dot(double new_fl_vdot[]) {
	fl_vel_dot[0] = new_fl_vdot[0];
	fl_vel_dot[1] = new_fl_vdot[1];
	fl_vel_dot[2] = new_fl_vdot[2];
	return;
}

void Node::add_to_displacement(double del_disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] += del_disp[i1];
	}	
	return;
}

void Node::add_to_fl_vel(double del_fl_vel[]) {
	fl_vel[0] += del_fl_vel[0];
	fl_vel[1] += del_fl_vel[1];
	fl_vel[2] += del_fl_vel[2];
	return;
}

void Node::set_initial_disp(double new_disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		initial_disp[i1] = new_disp[i1];
	}
	return;
}

void Node::set_initial_vel(double new_vel[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		initial_vel[i1] = new_vel[i1];
	}
	return;
}

void Node::set_initial_acc(double new_acc[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		initial_acc[i1] = new_acc[i1];
	}
	return;
}

void Node::set_initial_fl_vel(double new_fl_vel[]) {
	initial_fl_vel[0] = new_fl_vel[0];
	initial_fl_vel[1] = new_fl_vel[1];
	initial_fl_vel[2] = new_fl_vel[2];
	return;
}

void Node::set_initial_fl_vdot(double new_fl_vdot[]) {
	initial_fl_vel_dot[0] = new_fl_vdot[0];
	initial_fl_vel_dot[1] = new_fl_vdot[1];
	initial_fl_vel_dot[2] = new_fl_vdot[2];
	return;
}

void Node::set_prev_disp(double new_disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prev_disp[i1] = new_disp[i1];
	}
	return;
}

void Node::set_prev_vel(double new_vel[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prev_vel[i1] = new_vel[i1];
	}
	return;
}

void Node::set_prev_acc(double new_acc[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prev_acc[i1] = new_acc[i1];
	}
	return;
}

void Node::set_prev_fl_vel(double new_fl_vel[]) {
	prev_fl_vel[0] = new_fl_vel[0];
	prev_fl_vel[1] = new_fl_vel[1];
	prev_fl_vel[2] = new_fl_vel[2];
	return;
}

void Node::set_pfl_vel_lf(double new_vel[]) {
	p_fl_vel_lf[0] = new_vel[0];
	p_fl_vel_lf[1] = new_vel[1];
	p_fl_vel_lf[2] = new_vel[2];
	return;
}

void Node::set_prev_fl_vdot(double new_fl_vdot[]) {
	prev_fl_vel_dot[0] = new_fl_vdot[0];
	prev_fl_vel_dot[1] = new_fl_vdot[1];
	prev_fl_vel_dot[2] = new_fl_vdot[2];
	return;
}

void Node::initialize_disp() {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prev_disp[i1] = initial_disp[i1];
		prev_vel[i1] = initial_vel[i1];
		prev_acc[i1] = initial_acc[i1];
		displacement[i1] = initial_disp[i1];
	}
	return;
}

void Node::initialize_fl_vel() {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		prev_fl_vel[i1] = initial_fl_vel[i1];
		prev_fl_vel_dot[i1] = initial_fl_vel_dot[i1];
		fl_vel[i1] = initial_fl_vel[i1];
	}
	return;
}

void Node::initialize_temp() {
	prev_temp = initial_temp;
	prev_tdot = initial_tdot;
	temperature = initial_temp;
	return;
}

void Node::initialize_fl_den() {
	prev_fl_den = initial_fl_den;
	prev_fl_den_dot = initial_fl_den_dot;
	return;
}

void Node::initialize_turb_e() {
	prev_turb_e = initial_turb_e;
	prev_turb_edot = initial_turb_edot;
	return;
}

void Node::update_vel_acc(double nm_beta, double nm_gamma, double del_t) {
	double c1;
	double c2;
	int i1;
	c1 = 1.0 / (del_t * del_t * (nm_beta - nm_gamma));
	c2 = del_t * del_t * (0.5 + nm_beta - nm_gamma);
	for (i1 = 0; i1 < 6; i1++) {
		acceleration[i1] = c1 * (prev_disp[i1] - displacement[i1] + del_t * prev_vel[i1] + c2 * prev_acc[i1]);
		velocity[i1] = prev_vel[i1] + del_t * ((1.0 - nm_gamma) * prev_acc[i1] + nm_gamma * acceleration[i1]);
	}
	
	return;
}

void Node::update_fl_vel_dot(double nm_gamma, double del_t) {
	int i1;
	double c1;
	double c2;
	c1 = 1.0 / nm_gamma;
	c2 = 1.0 / del_t;
	for (i1 = 0; i1 < 3; i1++) {
		fl_vel_dot[i1] = c1 * (c2 * (fl_vel[i1] - p_fl_vel_lf[i1]) - (1.0 - nm_gamma) * prev_fl_vel_dot[i1]);
	}
	return;
}

void Node::update_tdot(double nm_gamma, double del_t) {
	double c1;
	double c2;
	c1 = 1.0 / nm_gamma;
	c2 = 1.0 / del_t;
	if (fluid) {
		temp_change_rate = c1 * (c2 * (temperature - p_temp_lf) - (1.0 - nm_gamma) * prev_tdot);
	}
	else {
		temp_change_rate = c1 * (c2 * (temperature - prev_temp) - (1.0 - nm_gamma) * prev_tdot);
	}
	return;
}

void Node::update_fl_den_dot(double nm_gamma, double del_t) {
	double c1;
	double c2;
	c1 = 1.0 / nm_gamma;
	c2 = 1.0 / del_t;
	if (fluid) {
		fl_den_dot = c1 * (c2 * (fl_den - p_fl_den_lf) - (1.0 - nm_gamma) * prev_fl_den_dot);
	}
	else {
	    fl_den_dot = c1 * (c2 * (fl_den - prev_fl_den) - (1.0 - nm_gamma) * prev_fl_den_dot);
	}
	return;
}

void Node::update_turb_edot(double nm_gamma, double del_t) {
	double c1;
	double c2;
	c1 = 1.0 / nm_gamma;
	c2 = 1.0 / del_t;
	if (fluid) {
		turb_edot = c1 * (c2 * (turb_e - p_turb_elf) - (1.0 - nm_gamma) * prev_turb_edot);
	}
	return;
}

void Node::advance_disp() {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		prev_disp[i1] = displacement[i1];
		prev_vel[i1] = velocity[i1];
		prev_acc[i1] = acceleration[i1];
	}
	return;	
}

void Node::advance_fl_vel() {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		prev_fl_vel[i1] = fl_vel[i1];
		prev_fl_vel_dot[i1] = fl_vel_dot[i1];
	}
	return;
}

void Node::advance_temp() {
	prev_temp = temperature;
	prev_tdot = temp_change_rate;
	return;
}

void Node::advance_fl_den() {
	prev_fl_den = fl_den;
	prev_fl_den_dot = fl_den_dot;
	return;
}

void Node::advance_turb_e() {
	prev_turb_e = turb_e;
	prev_turb_edot = turb_edot;
	return;
}

void Node::backstep_disp() {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
		displacement[i1] = prev_disp[i1];
		velocity[i1] = prev_vel[i1];
		acceleration[i1] = prev_acc[i1];
	}
	return;
}

void Node::backstep_fl_vel() {
	int i1;
	for (i1 = 0; i1 < 3; i1++) {
		fl_vel[i1] = prev_fl_vel[i1];
		fl_vel_dot[i1] = prev_fl_vel_dot[i1];
	}
	return;
}

void Node::backstep_temp() {
	temperature = prev_temp;
	temp_change_rate = prev_tdot;
	return;
}

void Node::backstep_fl_den() {
	fl_den = prev_fl_den;
	fl_den_dot = prev_fl_den_dot;
	return;
}

void Node::backstep_turb_e() {
	turb_e = prev_turb_e;
	turb_edot = prev_turb_edot;
	return;
}

void Node::add_design_variable(int d_index, double coef) {
	IDCapsule dv;
	dv.int_dat = d_index;
	dv.doub_dat = coef;
	d_var_lst.push_back(dv);
	return;
}

//dup1
void Node::get_crd(DiffDoub0 crd_out[], vector<DesignVariable>& dv_ar) {
	crd_out[0].set_val(coord[0]);
	crd_out[1].set_val(coord[1]);
	crd_out[2].set_val(coord[2]);
	int d_index;
	DiffDoub0 d_val;
	int comp;
	string cat;
	DiffDoub0 coef;
	for (auto& dv : d_var_lst) {
		d_index = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_index];
		this_dv.get_value(d_val);
		cat = this_dv.category;
		comp = this_dv.component - 1;
		if(cat == "nodeCoord") {
			coef.set_val(dv.doub_dat);
			coef.mult(d_val);
			crd_out[comp].add(coef);
		}
	}
	return;
}

void Node::get_disp(DiffDoub0 disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].set_val(displacement[i1]);
	}
	return;
}

void Node::get_elastic_dvload(DiffDoub0 ld[], vector<DesignVariable>& dv_ar) {
	int i1;
	int d_index;
	DiffDoub0 d_val;
	DiffDoub0 coef;
	string cat;
	int comp;
	
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].set_val(0.0);
	}
	for (auto& dv : d_var_lst) {
		d_index = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_index];
		this_dv.get_value(d_val);
		cat = this_dv.category;
		comp = this_dv.component - 1;
		if(cat == "elasticLoad") {
			coef.set_val(dv.doub_dat);
			coef.mult(d_val);
			ld[comp].add(coef);
		}
	}
	return;
}

void Node::get_thermal_dvload(DiffDoub0& ld, vector<DesignVariable>& dv_ar) {
	int d_index;
	DiffDoub0 d_val;
	DiffDoub0 coef;
	string cat;
	
	ld.set_val(0.0);
	for (auto& dv : d_var_lst) {
		d_index = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_index];
		this_dv.get_value(d_val);
		cat = this_dv.category;
		if (cat == "thermalLoad") {
			coef.set_val(dv.doub_dat);
			coef.mult(d_val);
			ld.add(coef);
		}
	}
	return;
}

//end dup
 
//skip 
 
//diff_doub1 versions: 
//dup1
void Node::get_crd(DiffDoub1 crd_out[], vector<DesignVariable>& dv_ar) {
	crd_out[0].set_val(coord[0]);
	crd_out[1].set_val(coord[1]);
	crd_out[2].set_val(coord[2]);
	int d_index;
	DiffDoub1 d_val;
	int comp;
	string cat;
	DiffDoub1 coef;
	for (auto& dv : d_var_lst) {
		d_index = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_index];
		this_dv.get_value(d_val);
		cat = this_dv.category;
		comp = this_dv.component - 1;
		if(cat == "nodeCoord") {
			coef.set_val(dv.doub_dat);
			coef.mult(d_val);
			crd_out[comp].add(coef);
		}
	}
	return;
}

void Node::get_disp(DiffDoub1 disp[]) {
	int i1;
	for (i1 = 0; i1 < 6; i1++) {
	    disp[i1].set_val(displacement[i1]);
	}
	return;
}

void Node::get_elastic_dvload(DiffDoub1 ld[], vector<DesignVariable>& dv_ar) {
	int i1;
	int d_index;
	DiffDoub1 d_val;
	DiffDoub1 coef;
	string cat;
	int comp;
	
	for (i1 = 0; i1 < 6; i1++) {
		ld[i1].set_val(0.0);
	}
	for (auto& dv : d_var_lst) {
		d_index = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_index];
		this_dv.get_value(d_val);
		cat = this_dv.category;
		comp = this_dv.component - 1;
		if(cat == "elasticLoad") {
			coef.set_val(dv.doub_dat);
			coef.mult(d_val);
			ld[comp].add(coef);
		}
	}
	return;
}

void Node::get_thermal_dvload(DiffDoub1& ld, vector<DesignVariable>& dv_ar) {
	int d_index;
	DiffDoub1 d_val;
	DiffDoub1 coef;
	string cat;
	
	ld.set_val(0.0);
	for (auto& dv : d_var_lst) {
		d_index = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_index];
		this_dv.get_value(d_val);
		cat = this_dv.category;
		if (cat == "thermalLoad") {
			coef.set_val(dv.doub_dat);
			coef.mult(d_val);
			ld.add(coef);
		}
	}
	return;
}

//end dup
 
//end skip 
 
 
