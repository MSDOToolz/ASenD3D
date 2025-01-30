#include <string>
#include <cmath>
#include "LoadClass.h"
#include "constants.h"
#include "SetClass.h"

using namespace std;

Load::Load() {
	type = "";
	active_time[0] = 0.0;
	active_time[1] = 1.0e+100;
	nd_set_ptr = max_int;
	el_set_ptr = max_int;
	return;
}

void Load::set_act_time(double new_at[]) {
	active_time[0] = new_at[0];
	active_time[1] = new_at[1];
	return;
}

void Load::set_load(double new_ld[]) {
	load[0] = new_ld[0];
	load[1] = new_ld[1];
	load[2] = new_ld[2];
	load[3] = new_ld[3];
	load[4] = new_ld[4];
	load[5] = new_ld[5];
	return;
}

void Load::set_norm_dir(double new_ndir[]) {
	normal_dir[0] = new_ndir[0];
	normal_dir[1] = new_ndir[1];
	normal_dir[2] = new_ndir[2];
	double mag = normal_dir[0]*normal_dir[0] + normal_dir[1]*normal_dir[1] + normal_dir[2]*normal_dir[2];
	mag = 1.0/sqrt(mag);
	normal_dir[0] = mag*normal_dir[0];
	normal_dir[1] = mag*normal_dir[1];
	normal_dir[2] = mag*normal_dir[2];
	return;
}

void Load::set_center(double new_cent[]) {
	center[0] = new_cent[0];
	center[1] = new_cent[1];
	center[2] = new_cent[2];
	return;
}

void Load::set_axis(double new_axis[]) {
	axis[0] = new_axis[0];
	axis[1] = new_axis[1];
	axis[2] = new_axis[2];
	double mag = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
	mag = 1.0/sqrt(mag);
	axis[0] = mag*axis[0];
	axis[1] = mag*axis[1];
	axis[2] = mag*axis[2];
	return;
}