#include "DesignVariableClass.h"
#include <string>
#include <vector>
#include "constants.h"
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "SetClass.h"

using namespace std;

DesignVariable::DesignVariable() {
	category = "";
	component = 1;
	layer = 0;
	el_set_name = "";
	el_set_ptr = max_int;
	nd_set_name = "";
	nd_set_ptr = max_int;
	coefs.clear();
	active_time[0] = 0.0;
	active_time[1] = 1.0e+100;
	value.set_val(0.0);
	diff_val.set_val(0.0);
	comp_el_list.clear();
}

void DesignVariable::set_active_time(double new_at[]) {
	active_time[0] = new_at[0];
	active_time[1] = new_at[1];
}

void DesignVariable::get_value_dfd0(DiffDoub0& inp) {
	inp.set_val_dfd0(value);
	return;
}

void DesignVariable::get_value_dfd1(DiffDoub1& inp) {
	inp.set_val_dfd1(diff_val);
	return;
}

void DesignVariable::add_comp_el(int eli) {
	for (auto& el : comp_el_list) {
		if (el == eli) {
			return;
		}
	}
	comp_el_list.push_back(eli);
	return;
}
