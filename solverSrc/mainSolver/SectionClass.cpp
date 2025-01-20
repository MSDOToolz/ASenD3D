#include "SectionClass.h"
#include "matrixFunctions.h"
#include <string>
#include <cmath>
#include <vector>

using namespace std;

const int max_int = 2000000000;

Material::Material() {
	name = "";
	int i1;
	for(i1 = 0; i1 < 36; i1++) {
		stiffness[i1] = 0.0;
		damping[i1] = 0.0;
	}
	return;
}

void Material::set_stiffness(int row, int col, double val) {
	stiffness[(6 * row + col)] = val;
	stiffness[(6 * col + row)] = val;
	return;
}

void Material::set_damping(int row, int col, double val) {
	damping[(6 * row + col)] = val;
	damping[(6 * col + row)] = val;
	return;
}

Fluid::Fluid() {
	name = "";
	return;
}

Layer::Layer() {
	mat_name = "";
	mat_ptr = max_int;
	return;
}

Section::Section() {
	type = "";
	int i1;
	for (i1 = 0; i1 < 36; i1++) {
		stiffness[i1] = 0.0;
		mass[i1] = 0.0;
		damping[i1] = 0.0;
	}
	for (i1 = 1; i1 < 8; i1++) {
		orientation[i1] = 0.0;
	}
	orientation[0] = 1.0;
	orientation[4] = 1.0;
	orientation[8] = 1.0;
	z_offset = 0.0;
	area = 0.0;
	for (i1 = 0; i1 < 5; i1++) {
		area_moment[i1] = 0.0;
	}
	for (i1 = 0; i1 < 6; i1++) {
		exp_load_coef[i1] = 0.0;
	}
	conductivity = 0.0;
	spec_heat = 0.0;
	pot_coef = 0.0;
	pot_exp = 1.0;
	damp_coef = 0.0;
	damp_exp = 1.0;
	cond_coef = 0.0;
	rad_coef = 0.0;
	den_vis_coef = 0.0;
	temp_vis_coef = 0.0;
	turb_vis_coef = 0.0;
	grad_vturb_coef = 0.0;
	diss_turb_coef = 0.0;
	enth_coef = 0.0;
	enth_exp = 1.0;
	pres_coef = 0.0;
	pres_exp = 1.0;
	ref_temp = 0.0;
	ref_den = 0.0;
	ref_turb_e = 0.0;
	ref_enth = 0.0;
	mat_ptr = max_int;
	fl_ptr = max_int;
	return;
}

void Section::set_orientation(double new_ori[]) {
	orientation[0] = new_ori[0];
	orientation[1] = new_ori[1];
	orientation[2] = new_ori[2];
	orientation[3] = new_ori[3];
	orientation[4] = new_ori[4];
	orientation[5] = new_ori[5];
	
	double mag = orientation[0]*orientation[0] + orientation[1]*orientation[1] + orientation[2]*orientation[2];
	mag = 1.0/sqrt(mag);
	orientation[0] = mag*orientation[0];
	orientation[1] = mag*orientation[1];
	orientation[2] = mag*orientation[2];
	
	cross_prod(&orientation[6],&orientation[0],&orientation[3]);
	
	mag = orientation[6]*orientation[6] + orientation[7]*orientation[7] + orientation[8]*orientation[8];
	mag = 1.0/sqrt(mag);
	orientation[6] = mag*orientation[6];
	orientation[7] = mag*orientation[7];
	orientation[8] = mag*orientation[8];
	
	cross_prod(&orientation[3],&orientation[6],&orientation[0]);
	
	return;
}

void Section::set_area_moment(double new_i[]) {
    area_moment[0] = new_i[0];
	area_moment[1] = new_i[1];
	area_moment[2] = new_i[2];
	area_moment[3] = new_i[3];
	area_moment[4] = new_i[4];
	return;
}

void Section::set_stiffness(int row, int col, double val) {
	stiffness[(6*row+col)] = val;
	stiffness[(6*col+row)] = val;
	return;
}

void Section::set_mass(int row, int col, double val) {
	mass[(6*row+col)] = val;
	mass[(6*col+row)] = val;
	return;
}

void Section::set_damping(int row, int col, double val) {
	damping[(6 * row + col)] = val;
	damping[(6 * col + row)] = val;
	return;
}

void Section::set_exp_ld(double new_exp_ld[]) {
	exp_load_coef[0] = new_exp_ld[0];
	exp_load_coef[1] = new_exp_ld[1];
	exp_load_coef[2] = new_exp_ld[2];
	exp_load_coef[3] = new_exp_ld[3];
	exp_load_coef[4] = new_exp_ld[4];
	exp_load_coef[5] = new_exp_ld[5];
	return;
}

int Section::get_layer_mat_ptr(int layi) {
	int i1 = 0;
	for (auto& lay : layers) {
		if (i1 == layi) {
			return lay.mat_ptr;
		}
		i1++;
	}
}