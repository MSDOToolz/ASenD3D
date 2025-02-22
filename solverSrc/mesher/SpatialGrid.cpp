#include "SpatialGrid.h"
#include <list>
#include <vector>

using namespace std;
IntList::IntList() {
	i_lst.clear();
	return;
}

int IntList::copy_to_vector(vector<int>& in_vec, int st_i, int max_len) {
	int i1 = st_i;
	for (auto& i2 : i_lst) {
		if (i1 >= max_len) {
			return i1;
		}
		in_vec[i1] = i2;
		i1++;
	}
	return i1;
}

SpatialGrid::SpatialGrid() {
	list_ar.clear();
	return;
}

void SpatialGrid::initialize(double x_range[], double x_spacing, double y_range[], double y_spacing, double z_range[], double z_spacing) {
	x_min = x_range[0];
	x_sp = x_spacing;
	y_min = y_range[0];
	y_sp = y_spacing;
	z_min = z_range[0];
	z_sp = z_spacing;
	x_bins = (x_range[1] - x_min) / x_sp + 1;
	y_bins = (y_range[1] - y_min) / y_sp + 1;
	z_bins = (z_range[1] - z_min) / z_sp + 1;
	int tot_bins = x_bins * y_bins * z_bins;
	list_ar = vector<IntList>(tot_bins);
	return;
}

void SpatialGrid::add_ent(int label, double crd[]) {
	int x_b = (crd[0] - x_min) / x_sp;
	if (x_b < 0) {
		x_b = 0;
	}
	if (x_b >= x_bins) {
		x_b = x_bins - 1;
	}
	int y_b = (crd[1] - y_min) / y_sp;
	if (y_b < 0) {
		y_b = 0;
	}
	if (y_b >= y_bins) {
		y_b = y_bins - 1;
	}
	int z_b = (crd[2] - z_min) / z_sp;
	if (z_b < 0) {
		z_b = 0;
	}
	if (z_b >= z_bins) {
		z_b = z_bins - 1;
	}
	int ind = (z_b*y_bins + y_b)*x_bins + x_b;
	list_ar[ind].i_lst.push_back(label);
}

int SpatialGrid::get_in_xyzrange(vector<int>& out_lst, int max_len, double x_range[], double y_range[], double z_range[]) {
	int i_min;
	int i_max;
	int j_min;
	int j_max;
	int k_min;
	int k_max;

	if (x_range[0] < x_range[1]) {
		i_min = (x_range[0] - x_min) / x_sp;
		if (i_min < 0) {
			i_min = 0;
		}
		i_max = (x_range[1] - x_min) / x_sp;
		if (i_max >= x_bins) {
			i_max = x_bins - 1;
		}
	}
	else {
		i_min = 0;
		i_max = x_bins - 1;
	}

	if (y_range[0] < y_range[1]) {
		j_min = (y_range[0] - y_min) / y_sp;
		if (j_min < 0) {
			j_min = 0;
		}
		j_max = (y_range[1] - y_min) / y_sp;
		if (j_max >= y_bins) {
			j_max = y_bins - 1;
		}
	}
	else {
		j_min = 0;
		j_max = y_bins - 1;
	}

	if (z_range[0] < z_range[1]) {
		k_min = (z_range[0] - z_min) / z_sp;
		if (k_min < 0) {
			k_min = 0;
		}
		k_max = (z_range[1] - z_min) / z_sp;
		if (k_max >= z_bins) {
			k_max = z_bins - 1;
		}
	}
	else {
		k_min = 0;
		k_max = z_bins - 1;
	}

	int i;
	int j;
	int k;
	int ind;
	int lst_len = 0;
	int num_added;
	for (k = k_min; k <= k_max; k++) {
		for (j = j_min; j <= j_max; j++) {
			for (i = i_min; i <= i_max; i++) {
				ind = (k * y_bins + j)*x_bins + i;
				lst_len = list_ar[ind].copy_to_vector(out_lst, lst_len, max_len);
			}
		}
	}

	return lst_len;
}

int SpatialGrid::get_in_radius(vector<int>& out_list, int max_len, double pt[], double rad) {
	double range[6];
	range[0] = pt[0] - rad;
	range[1] = pt[0] + rad;
	range[2] = pt[1] - rad;
	range[3] = pt[1] + rad;
	range[4] = pt[2] - rad;
	range[5] = pt[2] + rad;

	return get_in_xyzrange(out_list, max_len, &range[0], &range[2], &range[4]);
}