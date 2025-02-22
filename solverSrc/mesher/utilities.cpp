#include "utilities.h"
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

const double tol = 1.0e-12;

void cross_prod(double prod[], double v1[], double v2[]) {
	prod[0] = v1[1] * v2[2] - v1[2] * v2[1];
	prod[1] = v1[2] * v2[0] - v1[0] * v2[2];
	prod[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return;
}

void q_rfactor(double mat[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i2_min;
	int i2_max;
	int i3_min;
	int i3_max;
	int k11;
	int k12;
	int k22;
	int k23;
	double theta;
	double sth;
	double cth;
	double p1;
	double p2;

	for (i1 = st_col; i1 <= end_col; i1++) {
		i2_min = st_row + (i1 - st_col) + 1;
		if (tri_diag == 0) {
			i2_max = end_row;
		}
		else {
			i2_max = st_row + (i1 - st_col) + 1;
			if (i2_max > end_row) {
				i2_max = end_row;
			}
		}
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			k11 = (i2_min - 1) * col_dim + i1;
			k12 = i2 * col_dim + i1;
			if (abs(mat[k11]) < tol) {
				mat[k11] = tol;
			}
			theta = atan(mat[k12] / mat[k11]);
			sth = sin(theta);
			cth = cos(theta);
			i3_min = i1;
			if (tri_diag == 2) {
				i3_max = i1 + 2;
				if (i3_max > end_col) {
					i3_max = end_col;
				}
			}
			else {
				i3_max = end_col;
			}
			for (i3 = i3_min; i3 <= i3_max; i3++) {
				k22 = (i2_min - 1) * col_dim + i3;
				k23 = i2 * col_dim + i3;
				p1 = cth * mat[k22] + sth * mat[k23];
				p2 = -sth * mat[k22] + cth * mat[k23];
				mat[k22] = p1;
				mat[k23] = p2;
			}
			mat[k12] = theta;
		}
	}
	return;
}

void solveq_rx_eqb(double x_vec[], double mat[], double b_vec[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i2_min;
	int i2_max;
	int k11;
	int k12;
	double theta;
	double sth;
	double cth;
	double p1;
	double p2;
	double row_sum;

	for (i1 = st_col; i1 <= end_col; i1++) {
		i2_min = st_row + (i1 - st_col) + 1;
		if (tri_diag == 0) {
			i2_max = end_row;
		}
		else {
			i2_max = st_row + (i1 - st_col) + 1;
			if (i2_max > end_row) {
				i2_max = end_row;
			}
		}
		i3 = i2_min - 1;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			k12 = i2 * col_dim + i1;
			theta = mat[k12];
			sth = sin(theta);
			cth = cos(theta);
			p1 = cth * b_vec[i3] + sth * b_vec[i2];
			p2 = -sth * b_vec[i3] + cth * b_vec[i2];
			b_vec[i3] = p1;
			b_vec[i2] = p2;
		}
		x_vec[i1] = 0.0;
	}

	for (i1 = end_col; i1 >= st_col; i1--) {
		i3 = st_row + (i1 - st_col);
		i2_min = i1 + 1;
		if (tri_diag == 2) {
			i2_max = i1 + 2;
			if (i2_max > end_col) {
				i2_max = end_col;
			}
		}
		else {
			i2_max = end_col;
		}
		row_sum = 0.0;
		k11 = i3 * col_dim + i2_min;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			row_sum += mat[k11] * x_vec[i2];
			k11++;
		}
		k11 = i3 * col_dim + i1;
		x_vec[i1] = (b_vec[i3] - row_sum) / mat[k11];
	}

	return;
}

void read_input_line(string& file_line, string headings[], int hd_ld_space[], string data[], int& data_len) {
	int i1;
	int i2;
	int ln_len;
	int wrd_len;
	//getline(in_file, file_line);
	i1 = file_line.find("#");
	if (i1 > -1) {
		file_line = file_line.substr(0, i1);
	}
	file_line = file_line + " ";
	ln_len = file_line.length();
	i1 = file_line.find(":");
	data_len = 0;
	if (i1 > -1) {
		i2 = file_line.find_first_not_of(" -\n\t");
		wrd_len = i1 - i2;
		if (headings[0] == "" || hd_ld_space[0] == i2) {
			headings[0] = file_line.substr(i2, wrd_len);
			hd_ld_space[0] = i2;
			headings[1] = "";
			hd_ld_space[1] = 0;
			headings[2] = "";
			hd_ld_space[2] = 0;
			headings[3] = "";
			hd_ld_space[3] = 0;
		}
		else if (headings[1] == "" || hd_ld_space[1] == i2) {
			headings[1] = file_line.substr(i2, wrd_len);
			hd_ld_space[1] = i2;
			headings[2] = "";
			hd_ld_space[2] = 0;
			headings[3] = "";
			hd_ld_space[3] = 0;
		}
		else if (headings[2] == "" || hd_ld_space[2] == i2) {
			headings[2] = file_line.substr(i2, wrd_len);
			hd_ld_space[2] = i2;
			headings[3] = "";
			hd_ld_space[3] = 0;
		}
		else {
			headings[3] = file_line.substr(i2, wrd_len);
			hd_ld_space[3] = i2;
		}
		i1++;
		while (i1 < ln_len) {
			file_line = file_line.substr(i1);
			i1 = file_line.find_first_not_of(" ,[]\t\n");
			if (i1 > -1) {
				file_line = file_line.substr(i1);
				ln_len = file_line.length();
				i1 = file_line.find_first_of(" ,[]\t\n");
				if (i1 > -1) {
					data[data_len] = file_line.substr(0, i1);
					data_len++;
				}
				else {
					i1 = ln_len;
				}
			}
			else {
				i1 = ln_len;
			}
		}
	}
	else {
		i1 = file_line.find("- ");
		if (i1 > -1) {
			i1++;
			while (i1 < ln_len) {
				file_line = file_line.substr(i1);
				i1 = file_line.find_first_not_of(" ,[]\t\n");
				if (i1 > -1) {
					file_line = file_line.substr(i1);
					ln_len = file_line.length();
					i1 = file_line.find_first_of(" ,[]\t\n");
					if (i1 > -1) {
						data[data_len] = file_line.substr(0, i1);
						data_len++;
					}
					else {
						i1 = ln_len;
					}
				}
				else {
					i1 = ln_len;
				}
			}
		}
	}

	return;
}