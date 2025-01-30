#include "matrixFunctions.h"
#include "constants.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "LUMatClass.h"
#include "LowerTriMatClass.h"
#include "ConstraintClass.h"

using namespace std;

void sub_vec(vector<double>& sub_v, vector<double>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		sub_v[i2] = v_in[i1];
		i2++;
	}
}

void return_sv(vector<double>& sub_v, vector<double>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		v_in[i1] = sub_v[i2];
		i2++;
	}
}

void vec_to_ar(double ar[], std::vector<double>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		ar[i2] = vc[i1];
		i2++;
	}
}

double get_dist(double p1[], double p2[]) {
	double dist;
	double vec[3];
	vec[0] = p1[0] - p2[0];
	vec[1] = p1[1] - p2[1];
	vec[2] = p1[2] - p2[2];
	dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	return dist;
}

void cross_prod(double prod[], double v1[], double v2[]) {
	prod[0] = v1[1]*v2[2] - v1[2]*v2[1];
	prod[1] = v1[2]*v2[0] - v1[0]*v2[2];
	prod[2] = v1[0]*v2[1] - v1[1]*v2[0];
	return;
}

void q_rfactor(vector<double>& mat, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
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
		if(tri_diag == 0) {
			i2_max = end_row;
		} else {
			i2_max = st_row + (i1 - st_col) + 1;
			if(i2_max > end_row) {
				i2_max = end_row;
			}
		}
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			k11 = (i2_min-1)*col_dim + i1;
			k12 = i2*col_dim + i1;
			if(abs(mat[k11]) < tol) {
				mat[k11] = tol;
			}
			theta = atan(mat[k12]/mat[k11]);
			sth = sin(theta);
			cth = cos(theta);
			i3_min = i1;
			if(tri_diag == 2) {
				i3_max = i1 + 2;
				if(i3_max > end_col) {
					i3_max = end_col;
				}
			} else {
				i3_max = end_col;
			}
			for (i3 = i3_min; i3 <= i3_max; i3++) {
				k22 = (i2_min-1)*col_dim + i3;
				k23 = i2*col_dim + i3;
				p1 = cth*mat[k22] + sth*mat[k23];
				p2 = -sth*mat[k22] + cth*mat[k23];
				mat[k22] = p1;
				mat[k23] = p2;
			}
			mat[k12] = theta;
		}
	}
	return;
}

void solveq_rx_eqb(vector<double>& x_vec, vector<double>& mat, vector<double>& b_vec, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
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
		if(tri_diag == 0) {
			i2_max = end_row;
		} else {
			i2_max = st_row + (i1 - st_col) + 1;
			if(i2_max > end_row) {
				i2_max = end_row;
			}
		}
		i3 = i2_min - 1;
		for (i2 = i2_min; i2 <=i2_max; i2++) {
			k12 = i2*col_dim + i1;
			theta = mat[k12];
			sth = sin(theta);
			cth = cos(theta);
			p1 = cth*b_vec[i3] + sth*b_vec[i2];
			p2 = -sth*b_vec[i3] + cth*b_vec[i2];
		    b_vec[i3] = p1;
			b_vec[i2] = p2;
		}
		x_vec[i1] = 0.0;
	}
	
	for (i1 = end_col; i1 >= st_col; i1--) {
		i3 = st_row + (i1 - st_col);
		i2_min = i1 + 1;
		if(tri_diag == 2) {
			i2_max = i1 + 2;
			if(i2_max > end_col) {
				i2_max = end_col;
			}
		} else {
			i2_max = end_col;
		}
		row_sum = 0.0;
		k11 = i3*col_dim + i2_min;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			row_sum+= mat[k11]*x_vec[i2];
			k11++;
		}
		k11 = i3*col_dim + i1;
		x_vec[i1] = (b_vec[i3] - row_sum)/mat[k11];
	}
	
	return;
}

void conj_grad_sparse(vector<double>& soln, SparseMat& mat, ConstraintList& cnst, LowerTriMat& pc_mat, vector<double>& rhs, double conv_tol, int max_it) {
	int i1;
	int i2;
	int i3;
	int dim = mat.dim;
	int status_freq = dim / 24;
	double res;
	double r_next;
	double alpha;
	double dp;
	double beta;
	vector<double> h_vec(dim);
	vector<double> g_vec(dim);
	vector<double> z_vec(dim);
	vector<double> w_vec(dim);
	vector<double> t_vec(dim);

	if (!pc_mat.pos_def()) {
		cout << "Warning: preconditioning matrix for conjugate gradient solver is not positive definite." << endl;
	}

	for (i1 = 0; i1 < dim; i1++) {
		soln[i1] = 0.0;
		g_vec[i1] = -rhs[i1];
	}

	pc_mat.ldl_solve(w_vec, g_vec);

	res = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		h_vec[i1] = -w_vec[i1];
		res += g_vec[i1] * w_vec[i1];
	}

	i1 = 0;
	i3 = 0;
	while (i1 < max_it && res > conv_tol) {
		for (i2 = 0; i2 < dim; i2++) {
			z_vec[i2] = 0.0;
		}
		mat.vector_multiply(z_vec, h_vec, false);
		cnst.get_total_vec_mult(z_vec, h_vec, t_vec);
		dp = 0.0;
		for (i2 = 0; i2 < dim; i2++) {
			dp += h_vec[i2] * z_vec[i2];
		}
		alpha = res / dp;
		for (i2 = 0; i2 < dim; i2++) {
			soln[i2] += alpha * h_vec[i2];
			g_vec[i2] += alpha * z_vec[i2];
		}
		pc_mat.ldl_solve(w_vec, g_vec);
		r_next = 0.0;
		for (i2 = 0; i2 < dim; i2++) {
			r_next += g_vec[i2] * w_vec[i2];
		}
		beta = r_next / res;
		for (i2 = 0; i2 < dim; i2++) {
			h_vec[i2] = -w_vec[i2] + beta * h_vec[i2];
		}
		res = r_next;

		if (i3 == status_freq) {
			cout << "Conjugate gradient iteration: " << i1 << ", residual: " << res << endl;
			i3 = 0;
		}
		i1++;
		i3++;
	}

	if (i1 == max_it) {
		cout << "Warning: Conjugate gradient solver did not converge to specified tolerance, " << conv_tol << " within tne maximum number of iterations, " << max_it << endl;
	}

	cout << "Total CG iterations: " << i1 << endl;
	cout << "Final residual norm: " << res << endl; cout << "Final residual norm: " << res << endl;

	return;
}

void g_mres_sparse(vector<double>& soln, SparseMat& mat, ConstraintList& cnst, LUMat& pc_mat, vector<double>& rhs, double conv_tol, int max_it, int restart) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int it_ct;
	double res_nrm;
	double tmp;
	int dim = mat.dim;

	vector<double> res_vec(dim);
	vector<double> tmp_v(dim);
	vector<double> tmp_v2(dim);
	vector<double> tmp_v3(dim);
	i1 = dim * (restart + 1);
	vector<double> h_mat(i1);
	i1 = restart * (restart + 1);
	vector<double> phi_mat(i1);

	for (i1 = 0; i1 < dim; i1++) {
		soln[i1] = 0.0;
		tmp_v[i1] = -rhs[i1];
	}

	pc_mat.lu_solve(res_vec, tmp_v, false);
	//pc_mat.ldl_solve(res_vec, tmp_v);
	res_nrm = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		res_nrm += res_vec[i1] * res_vec[i1];
	}
	res_nrm = sqrt(res_nrm);
	cout << "Initial residual norm: " << res_nrm << endl;

	it_ct = 0;
	while (res_nrm > conv_tol && it_ct < max_it) {
		tmp = 1.0 / res_nrm;
		for (i1 = 0; i1 < dim; i1++) {
			h_mat[i1] = tmp * res_vec[i1];
		}
		i2 = restart * (restart + 1);
		for (i1 = 0; i1 < i2; i1++) {
			phi_mat[i1] = 0.0;
		}
		//generate basis vectors
		for (i1 = 1; i1 < (restart + 1); i1++) {
			//multiply by previous vector
			for (i2 = 0; i2 < dim; i2++) {
				tmp_v[i2] = 0.0;
			}
			i2 = dim * (i1 - 1);
			sub_vec(tmp_v3, h_mat, i2, i2 + dim);
			mat.vector_multiply(tmp_v, tmp_v3, false);
			cnst.get_total_vec_mult(tmp_v, tmp_v3, tmp_v2);
			//pc_mat.ldl_solve(tmp_v2, tmp_v);
			pc_mat.lu_solve(tmp_v2, tmp_v, false);
			//orthogonalize with all previous vectors
			for (i2 = 0; i2 < i1; i2++) {
				i4 = dim * i2;
				tmp = 0.0;
				for (i3 = 0; i3 < dim; i3++) {
					tmp += tmp_v2[i3] * h_mat[i4];
					i4++;
				}
				i5 = i2 * restart + (i1 - 1);
				phi_mat[i5] = tmp;
				i4 = dim * i2;
				for (i3 = 0; i3 < dim; i3++) {
					tmp_v2[i3] -= tmp * h_mat[i4];
					i4++;
				}
			}
			//calculate magnitude and Set new unit basis std::vector
			tmp = 0.0;
			for (i3 = 0; i3 < dim; i3++) {
				tmp += tmp_v2[i3] * tmp_v2[i3];
			}
			tmp = sqrt(tmp);
			i5 = i1 * restart + (i1 - 1);
			phi_mat[i5] = tmp;
			tmp = 1.0 / tmp;
			i4 = dim * i1;
			for (i2 = 0; i2 < dim; i2++) {
				h_mat[i4] = tmp * tmp_v2[i2];
				i4++;
			}
		}
		//find the least squares solution
		i3 = 0;
		for (i1 = 0; i1 <= restart; i1++) {
			tmp_v[i1] = 0.0;
			for (i2 = 0; i2 < dim; i2++) {
				tmp_v[i1] -= h_mat[i3] * res_vec[i2];
				i3++;
			}
		}
		q_rfactor(phi_mat, restart, 0, restart, 0, (restart-1), 1);
		solveq_rx_eqb(tmp_v2, phi_mat, tmp_v, restart, 0, restart, 0, (restart-1), 1);
		//update the solution std::vector
		i3 = 0;
		for (i1 = 0; i1 < restart; i1++) {
			for (i2 = 0; i2 < dim; i2++) {
				soln[i2] += h_mat[i3] * tmp_v2[i1];
				i3++;
			}
		}
		//update residual std::vector
		for (i1 = 0; i1 < dim; i1++) {
			tmp_v[i1] = -rhs[i1];
		}
		mat.vector_multiply(tmp_v, soln, false);
		cnst.get_total_vec_mult(tmp_v, soln, tmp_v2);
		//pc_mat.ldl_solve(res_vec, tmp_v);
		pc_mat.lu_solve(res_vec, tmp_v, false);
		res_nrm = 0.0;
		for (i1 = 0; i1 < dim; i1++) {
			res_nrm += res_vec[i1] * res_vec[i1];
		}
		res_nrm = sqrt(res_nrm);
		it_ct += restart;
		cout << "Iteration: " << it_ct << ",  Residual Norm: " << res_nrm << endl;
	}

	if (res_nrm > conv_tol) {
		cout << "Warning: GMRES solver did not converge to the requested tolerance of " << conv_tol << endl;
		cout << "Residual norm after " << it_ct << " iterations: " << res_nrm << endl;
	}

	return;
}

void sym_factor(vector<double>& mat, vector<double>& q_mat, int mat_dim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double theta;
	double sth;
	double cth;
	double a1;
	double a2;
	
	i3 = 0;
	for (i1 = 0; i1 < mat_dim; i1++) {
		for (i2 = 0; i2 < mat_dim; i2++) {
			if(i2 == i1) {
				q_mat[i3] = 1.0;
			} else {
				q_mat[i3] = 0.0;
			}
			i3++;
		}
	}
	
	for (i1 = 0; i1 < mat_dim; i1++) {
		for (i2 = i1+2; i2 < mat_dim; i2++) {
			i3 = (i1+1)*mat_dim + i1;
			if(abs(mat[i3]) < tol) {
				theta = pi_2;
			} else {
				i4 = i2*mat_dim + i1;
				theta = atan(mat[i4]/mat[i3]);
			}
			sth = sin(theta);
			cth = cos(theta);
			i4 = (i1+1)*mat_dim + i1;
			i5 = i2*mat_dim + i1;
			for (i3 = i1; i3 < mat_dim; i3++) {
				a1 = cth*mat[i4] + sth*mat[i5];
				a2 = -sth*mat[i4] + cth*mat[i5];
				mat[i4] = a1;
				mat[i5] = a2;
				i4++;
				i5++;
 			}
			i4 = i1 + 1;
			i5 = i2;
			for (i3 = 0; i3 < mat_dim; i3++) {
				a1 = cth*mat[i4] + sth*mat[i5];
				a2 = -sth*mat[i4] + cth*mat[i5];
				mat[i4] = a1;
				mat[i5] = a2;
				a1 = cth*q_mat[i4] + sth*q_mat[i5];
				a2 = -sth*q_mat[i4] + cth*q_mat[i5];
				q_mat[i4] = a1;
				q_mat[i5] = a2;
				i4+= mat_dim;
				i5+= mat_dim;
			}
		}
	}		
	
	return;
}

void get_char_fun(DiffDoub1& c_fun, vector<DiffDoub1>& mat, int mat_dim, vector<double>& e_vals, double lam, int tri_diag) {
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
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;

	for (i1 = 0; i1 <= (mat_dim - 1); i1++) {
		i2_min = i1 + 1;
		if(tri_diag == 0) {
			i2_max = mat_dim - 1;
		} else {
			i2_max = i1 + 1;
			if(i2_max > (mat_dim - 1)) {
				i2_max = mat_dim - 1;
			}
		}
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			k11 = (i2_min-1)*mat_dim + i1;
			k12 = i2*mat_dim + i1;
			if(abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.set_val_dfd1(mat[k12]);
			tmp.set_val_dfd1(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.set_val_dfd1(theta);
			sth.sn();
			cth.set_val_dfd1(theta);
			cth.cs();
			i3_min = i1;
			if(tri_diag == 2) {
				i3_max = i1 + 2;
				if(i3_max > (mat_dim - 1)) {
					i3_max = mat_dim - 1;
				}
			} else {
				i3_max = mat_dim - 1;
			}
			for (i3 = i3_min; i3 <= i3_max; i3++) {
				k22 = (i2_min-1)*mat_dim + i3;
				k23 = i2*mat_dim + i3;
				p1.set_val_dfd1(cth);
				p1.mult(mat[k22]);
				tmp.set_val_dfd1(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.set_val_dfd1(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.set_val_dfd1(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].set_val_dfd1(p1);
				mat[k23].set_val_dfd1(p2);
			}
			mat[k12].set_val_dfd1(theta);
		}
	}
	
	c_fun.set_val_dfd1(mat[0]);
	i2 = mat_dim + 1;
	for (i1 = 1; i1 < mat_dim; i1++) {
		c_fun.mult(mat[i2]);
		while(abs(c_fun.val) > max_mag) {
			c_fun.val*= min_mag;
			c_fun.dval*= min_mag;
		}
		while(abs(c_fun.val) < min_mag) {
			c_fun.val*= max_mag;
			c_fun.dval*= max_mag;
		}
		i2+= (mat_dim + 1);
	}
	
	double term;
	for (auto& ev : e_vals) {
		term = ev - lam;
		tmp.set_val_2(term,-1.0);
		c_fun.dvd(tmp);
		while(abs(c_fun.val) > max_mag) {
			c_fun.val*= min_mag;
			c_fun.dval*= min_mag;
		}
		while(abs(c_fun.val) < min_mag) {
			c_fun.val*= max_mag;
			c_fun.dval*= max_mag;
		}
	}
	
	return;
}

void get_evals(vector<double>& e_vals, vector<double>& mat, int mat_dim, double lam_init, double conv_tol, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int mat_size = mat_dim*mat_dim;
	vector<double> e_val_list;
	vector<DiffDoub1> mat_copy(mat_size);
	DiffDoub1 c_fun;
	double lam = lam_init;
	double d_lam;
	double d_trm;
	double back_step = 1.0e+6*conv_tol;
	int num_found = 0;
	int max_it = 100*mat_dim;
	int it = 0;
	while(num_found < mat_dim && it < max_it) {
		i3 = 0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			for (i2 = 0; i2 < mat_dim; i2++) {
				if(i1 == i2) {
					d_trm = mat[i3] - lam;
					mat_copy[i3].set_val_2(d_trm,-1.0);
				} else {
					mat_copy[i3].set_val(mat[i3]);
				}
				i3++;
			}
		}
		get_char_fun(c_fun,mat_copy,mat_dim,e_val_list,lam,tri_diag);
		d_lam = -c_fun.val/c_fun.dval;
		lam+= d_lam;
		if(abs(d_lam) < conv_tol) {
			//e_val_list.add_entry(lam);
			e_val_list.push_back(lam);
			num_found++;
			lam-= back_step;
		}
		it++;
	}
	//doub_list_ent *this_ent = e_val_list.get_first();
	i1 = 0;
	for (auto& ev : e_val_list) {
		e_vals[i1] = ev;
		i1++;
	}
	return;
}

double get_low_eval(vector<double>& mat, int mat_dim) {
	vector<double> vec(mat_dim);
	vector<double> prev_vec(mat_dim);
	double high;
	double tmp;
	double mag;
	double dp;
	int i1;
	int i2;
	int i3;
	int it;
	int max_it;
	
	mag = 0.0;
	for (i1 = 0; i1 < mat_dim; i1++) {
		tmp = sin(1.0*i1);
		prev_vec[i1] = tmp;
		mag+= tmp*tmp;
	}
	mag = 1.0/sqrt(mag);
	for (i1 = 0; i1 < mat_dim; i1++) {
		prev_vec[i1]*= mag;
	}
	
	bool conv = false;
	max_it = 100*mat_dim;
	it = 0;
	while(!conv && it < max_it) {
		i3 = 0;
		mag = 0.0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			vec[i1] = 0.0;
			for (i2 = 0; i2 < mat_dim; i2++) {
				vec[i1]+= mat[i3]*prev_vec[i2];
				i3++;
			}
			mag+= vec[i1]*vec[i1];
		}
		mag = sqrt(mag);
		dp = 0.0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			dp+= vec[i1]*prev_vec[i1];
		}
		if(abs(dp) > 0.99999999*mag) {
			conv = true;
		} else {
			mag = 1.0/mag;
			for (i1 = 0; i1 < mat_dim; i1++) {
				prev_vec[i1] = mag*vec[i1];
			}
		}
		it++;
	}
	
	if(dp < 0.0) {
		return dp;
	}
	
	high = dp;
    
    i3 = 0;	
	for (i1 = 0; i1 < mat_dim; i1++) {
		mat[i3]-= high;
		i3+= (mat_dim + 1);
	}
	
	mag = 0.0;
	for (i1 = 0; i1 < mat_dim; i1++) {
		tmp = sin(1.0*i1);
		prev_vec[i1] = tmp;
		mag+= tmp*tmp;
	}
	mag = 1.0/sqrt(mag);
	for (i1 = 0; i1 < mat_dim; i1++) {
		prev_vec[i1]*= mag;
	}
	
	conv = false;
	it = 0;
	while(!conv && it < max_it) {
		i3 = 0;
		mag = 0.0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			vec[i1] = 0.0;
			for (i2 = 0; i2 < mat_dim; i2++) {
				vec[i1]+= mat[i3]*prev_vec[i2];
				i3++;
			}
			mag+= vec[i1]*vec[i1];
		}
		mag = sqrt(mag);
		dp = 0.0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			dp+= vec[i1]*prev_vec[i1];
		}
		if(abs(dp) > 0.99999999*mag) {
			conv = true;
		} else {
			mag = 1.0/mag;
			for (i1 = 0; i1 < mat_dim; i1++) {
				prev_vec[i1] = mag*vec[i1];
			}
		}
		it++;
	}
	
	i3 = 0;	
	for (i1 = 0; i1 < mat_dim; i1++) {
		mat[i3]+= high;
		i3+= (mat_dim + 1);
	}
	
	dp+= high;
	return dp;
}

void eigen_solve(vector<double>& e_vals, vector<double>& e_vecs, vector<double>& mat, int mat_dim, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double tmp;
	int mat_size = mat_dim*mat_dim;
	int new_td = tri_diag;
	vector<double> q_mat;
	if(tri_diag == 0) {
		q_mat = vector<double>(mat_size);
		sym_factor(mat,q_mat,mat_dim);
		new_td = 1;
	}
	
    double low = get_low_eval(mat,mat_dim);
	
	double mat_mag = 0.0;
	for (i1 = 0; i1 < mat_size; i1++) {
		tmp = abs(mat[i1]);
		if(tmp > mat_mag) {
			mat_mag = tmp;
		}
	}
	
	low-= 0.001*mat_mag;
	double conv_tol = 1.0e-12*mat_mag;
	get_evals(e_vals,mat,mat_dim,low,conv_tol,new_td);
	
	//doub_list_ent *this_val = e_vals.get_first();
	vector<double> mat_copy(mat_size);
	vector<double> vec(mat_dim);
	vector<double> prev_vec(mat_dim);
	vector<double> b_vec(mat_dim);
	double shift;
	bool conv;
	double mag;
	double dp;
	int it;
	int max_it = 100*mat_dim;
	for (i5 = 0; i5 < mat_dim; i5++) {
		shift = e_vals[i5] + 1.0e-8*mat_mag;
		i3 = 0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			for (i2 = 0; i2 < mat_dim; i2++) {
				if (i1 == i2) {
					mat_copy[i3] = mat[i3] - shift;
				}
				else {
					mat_copy[i3] = mat[i3];
				}
				i3++;
			}
		}
		q_rfactor(mat_copy, mat_dim, 0, (mat_dim - 1), 0, (mat_dim - 1), new_td);
		mag = 0.0;
		for (i1 = 0; i1 < mat_dim; i1++) {
			tmp = sin(1.0 * i1 * (i5 + 1));
			prev_vec[i1] = tmp;
			mag += tmp * tmp;
		}
		mag = 1.0 / sqrt(mag);
		for (i1 = 0; i1 < mat_dim; i1++) {
			prev_vec[i1] = mag * prev_vec[i1];
		}
		conv = false;
		it = 0;
		while (!conv && it < max_it) {
			for (i1 = 0; i1 < mat_dim; i1++) {
				b_vec[i1] = prev_vec[i1];
			}
			solveq_rx_eqb(vec, mat_copy, b_vec, mat_dim, 0, (mat_dim - 1), 0, (mat_dim - 1), new_td);
			mag = 0.0;
			dp = 0.0;
			for (i1 = 0; i1 < mat_dim; i1++) {
				mag += vec[i1] * vec[i1];
				dp += vec[i1] * prev_vec[i1];
			}
			mag = sqrt(mag);
			if (abs(dp) > 0.999999 * mag) {
				conv = true;
				if (tri_diag == 0) {
					i3 = 0;
					//i4 = i5;
					i4 = i5 * mat_dim;
					for (i1 = 0; i1 < mat_dim; i1++) {
						//i4 = i5*mat_dim + i1
						e_vecs[i4] = 0.0;
						for (i2 = 0; i2 < mat_dim; i2++) {
							e_vecs[i4] += q_mat[i3] * prev_vec[i2];
							i3++;
						}
						//i4 += mat_dim;
						i4++;
					}
				}
				else {
					//i4 = i5;
					i4 = i5 * mat_dim;
					for (i1 = 0; i1 < mat_dim; i1++) {
						//i4 = i5*mat_dim + i1;
						e_vecs[i4] = prev_vec[i1];
						//i4 += mat_dim;
						i4++;
					}
				}
			}
			else {
				mag = 1.0 / mag;
				for (i1 = 0; i1 < mat_dim; i1++) {
					prev_vec[i1] = mag * vec[i1];
				}
			}
			it++;
		}
	}
	return;
}

void sym_eigen_solve(vector<double>& e_vals, vector<double>& e_vecs, vector<double>& mat, int mat_dim, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i3_min;
	int i3_max;
	int i4;
	int i5;
	int loop_ct;
	double tmp;
	double mag;
	double dp;

	vector<double> temp_v1(mat_dim);
	vector<double> temp_v2(mat_dim);
	vector<double> e_vtmp(mat_dim*mat_dim);
	vector<double> q_mat(mat_dim * mat_dim);
	vector<int> srt_order(mat_dim);

	if (tri_diag == 0) {
		sym_factor(mat, q_mat, mat_dim);
	}

	for (i1 = 0; i1 < mat_dim; i1++) {
		for (i2 = 0; i2 < mat_dim; i2++) {
			tmp = sin(1.0 * i2 * (i1+1));
			temp_v1[i2] = tmp;
		}
		for (i2 = 0; i2 < i1; i2++) {
			dp = 0.0;
			i4 = i2 * mat_dim;
			for (i3 = 0; i3 < mat_dim; i3++) {
				dp += temp_v1[i3] * e_vtmp[i4];
				i4++;
			}
			i4 = i2 * mat_dim;
			for (i3 = 0; i3 < mat_dim; i3++) {
				temp_v1[i3] -= dp*e_vtmp[i4];
				i4++;
			}
		}
		mag = 0.0;
		for (i2 = 0; i2 < mat_dim; i2++) {
			mag += temp_v1[i2] * temp_v1[i2];
		}
		mag = 1.0 / sqrt(mag);
		for (i2 = 0; i2 < mat_dim; i2++) {
			temp_v1[i2] *= mag;
		}
		dp = 0.0;
		loop_ct = 0;
		while (abs(dp) < 0.99999999*mag && loop_ct < 10000) {
			// multiply by mat
			for (i2 = 0; i2 < mat_dim; i2++) {
				temp_v2[i2] = 0.0;
				if (tri_diag == 1) {
					i3_min = i2 - 1;
					if (i3_min < 0) {
						i3_min = 0;
					}
					i3_max = i2 + 2;
					if (i3_max > mat_dim) {
						i3_max = mat_dim;
					}
				}
				else {
					i3_min = 0;
					i3_max = mat_dim;
				}
				i4 = i2 * mat_dim + i3_min;
				for (i3 = i3_min; i3 < i3_max; i3++) {
					//i4 = i2 * mat_dim + i3;
					temp_v2[i2] += mat[i4] * temp_v1[i3];
					i4++;
				}
			}
			// orthogonalize with previous std::vectors
			for (i2 = 0; i2 < i1; i2++) {
				dp = 0.0;
				i4 = i2 * mat_dim;
				for (i3 = 0; i3 < mat_dim; i3++) {
					dp += temp_v2[i3] * e_vtmp[i4];
					i4++;
				}
				i4 = i2 * mat_dim;
				for (i3 = 0; i3 < mat_dim; i3++) {
					temp_v2[i3] -= dp * e_vtmp[i4];
					i4++;
				}
			}
			// get mag, dp
			mag = 0.0;
			dp = 0.0;
			for (i2 = 0; i2 < mat_dim; i2++) {
				mag += temp_v2[i2] * temp_v2[i2];
				dp += temp_v1[i2] * temp_v2[i2];
			}
			mag = sqrt(mag);
			// update temp_v1
			tmp = 1.0 / mag;
			for (i2 = 0; i2 < mat_dim; i2++) {
				temp_v1[i2] = tmp * temp_v2[i2];
			}
			loop_ct++;
		}
		e_vals[i1] = dp;
		if (tri_diag == 0) {
			i4 = 0;
			for (i2 = 0; i2 < mat_dim; i2++) {
				i5 = i1 * mat_dim + i2;
				e_vtmp[i5] = 0.0;
				for (i3 = 0; i3 < mat_dim; i3++) {
					e_vtmp[i5] += q_mat[i4] * temp_v1[i3];
					i4++;
				}
			}
		}
		else {
			i3 = i1 * mat_dim;
			for (i2 = 0; i2 < mat_dim; i2++) {
				e_vtmp[i3] = temp_v1[i2];
				i3++;
			}
		}
	}

	// sort pairs

	for (i1 = 0; i1 < mat_dim; i1++) {
		srt_order[i1] = i1;
	}

	for (i1 = 0; i1 < mat_dim; i1++) {
		for (i2 = 0; i2 < (mat_dim - 1); i2++) {
			if (e_vals[i2 + 1] < e_vals[i2]) {
				tmp = e_vals[i2];
				e_vals[i2] = e_vals[i2 + 1];
				e_vals[i2 + 1] = tmp;
				i3 = srt_order[i2];
				srt_order[i2] = srt_order[i2 + 1];
				srt_order[i2 + 1] = i3;
			}
		}
	}

	for (i1 = 0; i1 < mat_dim; i1++) {
		i2 = mat_dim * i1;
		i3 = mat_dim * srt_order[i1];
		for (i4 = 0; i4 < mat_dim; i4++) {
			e_vecs[i2] = e_vtmp[i3];
			i2++;
			i3++;
		}
	}

	return;
}

void eigen_sparse_direct(vector<double>& e_vals, vector<double>& e_vecs, int num_pairs, LowerTriMat& mat, vector<double>& mass_mat, int mat_dim) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	
	int h_cols = 5*num_pairs;
	i1 = h_cols*mat_dim;
	vector<double> h_mat(i1);
	vector<double> t_vec1(mat_dim);
	vector<double> t_vec2(mat_dim);
	i1 = h_cols*(h_cols - 1);
	vector<double> coef_mat(i1);
	
	for (i2 = 0; i2 < i1; i2++) {
		coef_mat[i2] = 0.0;
	}
	
	double mag = 0.0;
	double tmp;
	for (i2 = 0; i2 < mat_dim; i2++) {
		mass_mat[i2] = sqrt(mass_mat[i2]);
	}

	for (i2 = mat_dim; i2 < (2 * mat_dim); i2++) {
		h_mat[i2] = sin(1.0 * i2);
	}

	sub_vec(t_vec1, h_mat, 0, mat_dim);
	sub_vec(t_vec2, h_mat, mat_dim, mat_dim + mat_dim);
	mat.ldl_solve(t_vec1, t_vec2);
	return_sv(t_vec1, h_mat, 0, mat_dim);

	mag = 0.0;
	for (i2 = 0; i2 < mat_dim; i2++) {
		mag += h_mat[i2] * h_mat[i2];
	}
	mag = 1.0 / sqrt(mag);

	for (i2 = 0; i2 < mat_dim; i2++) {
		h_mat[i2] *= mag;
	}
	
	double dp;
	double st_vec;
	for (i1 = 1; i1 < h_cols; i1++) {
		i3 = mat_dim*(i1-1);
		for (i2 = 0; i2 < mat_dim; i2++) {
			t_vec1[i2] = mass_mat[i2]*h_mat[i3];
			i3++;
		}
		mat.ldl_solve(t_vec2,t_vec1);
		for (i2 = 0; i2 < mat_dim; i2++) {
			t_vec1[i2] = t_vec2[i2]*mass_mat[i2];
		}
		if(i1 == 1) {
			st_vec = 0;
		} else {
			st_vec = i1 - 2;
		}
		for (i2 = st_vec; i2 < i1; i2++) {
			dp = 0.0;
			i3 = mat_dim*i2;
			for (i4 = 0; i4 < mat_dim; i4++) {
				dp+= t_vec1[i4]*h_mat[i3];
				i3++;
			}
			i3 = mat_dim*i2;
			for (i4 = 0; i4 < mat_dim; i4++) {
				t_vec1[i4]-= dp*h_mat[i3];
				i3++;
			}
			i3 = (h_cols-1)*i2 + (i1-1);
			coef_mat[i3] = dp;
		}
		mag = 0.0;
		for (i2 = 0; i2 < mat_dim; i2++) {
			mag+= t_vec1[i2]*t_vec1[i2];
		}
		mag = sqrt(mag);
		i3 = h_cols*i1 - 1; // (h_cols-1)*i1 + (i1 - 1)
		coef_mat[i3] = mag;
		mag = 1.0/mag;
		i3 = mat_dim*i1;
		for (i2 = 0; i2 < mat_dim; i2++) {
			h_mat[i3] = mag*t_vec1[i2];
			i3++;
		}
	}
	
	i1 = h_cols - 1;
	vector<double> coef_vals(i1);
	vector<double> coef_vecs(i1*i1);
	//eigen_solve(coef_vals, coef_vecs, coef_mat, i1, 2);
	sym_eigen_solve(coef_vals, coef_vecs, coef_mat, i1, 1);

	for (i3 = 0; i3 < mat_dim; i3++) {
		mass_mat[i3] = 1.0/mass_mat[i3];
	}
	
	i1 = 0;
	i2 = h_cols - 2;
	while (i1 < num_pairs && i2 >= 0) {
		for (i3 = 0; i3 < mat_dim; i3++) {
			t_vec1[i3] = 0.0;
		}
		i5 = 0;
		for (i3 = 0; i3 < (h_cols-1); i3++) {
			//i6 = (h_cols-1)*i3 + i2;
			i6 = i2 * (h_cols - 1) + i3;
			for (i4 = 0; i4 < mat_dim; i4++) {
				t_vec1[i4] += h_mat[i5] * coef_vecs[i6];
				i5++;
			}
		}
		mag = 0.0;
		for (i3 = 0; i3 < mat_dim; i3++) {
			tmp = t_vec1[i3]*mass_mat[i3];
			t_vec2[i3] = tmp;
			mag+= tmp*tmp;
		}
		mag = 1.0/sqrt(mag);
		for (i3 = 0; i3 < mat_dim; i3++) {
			t_vec2[i3] *= mag;
		}
		if (i1 > 0) {
			tmp = 0.0;
			i4 = mat_dim * (i1 - 1);
			for (i3 = 0; i3 < mat_dim; i3++) {
				tmp += e_vecs[i4] * t_vec2[i3];
				i4++;
			}
			if (abs(tmp) < 0.9) {
				for (i3 = 0; i3 < mat_dim; i3++) {
					e_vecs[i4] = t_vec2[i3];
					i4++;
				}
				e_vals[i1] = 1.0 / coef_vals[i2];
				i1++;
				i2--;
			}
			else {
				i2--;
			}
		}
		else {
			i4 = 0;
			for (i3 = 0; i3 < mat_dim; i3++) {
				//i4 = i1 * mat_dim + i3;
				e_vecs[i4] = t_vec2[i3];
				i4++;
			}
			e_vals[i1] = 1.0 / coef_vals[i2];
			i1++;
			i2--;
		}
	}
	
	for (i1 = 0; i1 < mat_dim; i1++) {
		tmp = 1.0/mass_mat[i1];
		mass_mat[i1] = tmp*tmp;
	}
	
	return;
}

double ray_quot(vector<double>& grad, vector<double>& kv, vector<double>& mv, SparseMat& mat, ConstraintList& cnst, vector<double>& mass_mat, vector<double>& in_vec) {
	int i1;
	double r_c;
	int dim = mat.dim;
	double v_kv;
	double v_mv = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		kv[i1] = 0.0;
		mv[i1] = mass_mat[i1] * in_vec[i1];
		v_mv += in_vec[i1] * mv[i1];
	}
	mat.vector_multiply(kv, in_vec, false);
	cnst.get_total_vec_mult(kv, in_vec, grad);
	v_kv = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		v_kv += in_vec[i1] * kv[i1];
	}
	r_c = v_kv / v_mv;
	for (i1 = 0; i1 < dim; i1++) {
		grad[i1] = 2.0 * ((kv[i1] / v_mv) - r_c * (mv[i1] / v_mv));
	}
	return r_c;
}

double unitize_vec(vector<double>& vec, int dim) {
	int i1;
	double mag;
	double minv;
	mag = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		mag += vec[i1] * vec[i1];
	}
	mag = sqrt(mag);
	minv = 1.0 / mag;
	for (i1 = 0; i1 < dim; i1++) {
		vec[i1] *= minv;
	}
	return mag;
}

void get_nearest_evec_rq(SparseMat& mat, ConstraintList& cnst, vector<double>& mass_mat, vector<double>& in_vecs, vector<double>& e_vals, int num_vecs, int max_it) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dim;
	double dp;
	double mag;
	double minv;
	double res;
	double r_q0;
	double r_q1;
	double step_len;
	bool reduced;
	
	dim = mat.dim;

	vector<double> grad0(dim);
	vector<double> grad1(dim);
	vector<double> d2_rq(dim);
	vector<double> v_step(dim);
	vector<double> kv(dim);
	vector<double> mv(dim);
	vector<double> t_vec1(dim);

	for (i1 = 0; i1 < num_vecs; i1++) {
		i2 = i1 * dim;
		sub_vec(t_vec1, in_vecs, i2, i2 + dim);
		unitize_vec(t_vec1, dim);
		r_q0 = ray_quot(grad0, kv, mv, mat, cnst, mass_mat, t_vec1);
		return_sv(t_vec1, in_vecs, i2, i2 + dim);
		res = 1.0;
		i3 = 0;
		while (i3 < max_it && res > 1.0e-6) {
			for (i4 = 0; i4 < dim; i4++) {
				v_step[i4] = -grad0[i4];
				if (v_step[i4] = 0.0) {
					v_step[i4] = 1.0e-12;
				}
			}
			unitize_vec(v_step, dim);
			for (i4 = 0; i4 < dim; i4++) {
				in_vecs[i2 + i4] += 0.01*v_step[i4];
			}
			sub_vec(t_vec1, in_vecs, i2, i2 + dim);
			r_q1 = ray_quot(grad1, kv, mv, mat, cnst, mass_mat, t_vec1);
			for (i4 = 0; i4 < dim; i4++) {
				in_vecs[i2 + i4] -= 0.01 * v_step[i4];
			}
			for (i4 = 0; i4 < dim; i4++) {
				d2_rq[i4] = (grad1[i4] - grad0[i4]) / (0.01 * v_step[i4]);
				if (d2_rq[i4] != 0.0) {
					v_step[i4] = -grad0[i4] / d2_rq[i4];
				}
				else {
					v_step[i4] = 0.0;
				}
				in_vecs[i2 + i4] += v_step[i4];
			}
			sub_vec(t_vec1, in_vecs, i2, i2 + dim);
			unitize_vec(t_vec1, dim);
			r_q0 = ray_quot(grad0, kv, mv, mat, cnst, mass_mat, t_vec1);
			return_sv(t_vec1, in_vecs, i2, i2 + dim);
			unitize_vec(kv, dim);
			unitize_vec(mv, dim);
			dp = 0.0;
			for (i4 = 0; i4 < dim; i4++) {
				dp += kv[i4] * mv[i4];
			}
			res = 1.0 - abs(dp);
			i3++;
		}
		cout << "Warning: eigenstd::vector " << i1 << " did not converge within the max iterations." << endl;
		e_vals[i1] = r_q0;
	}

}

//dup1
void q_rfactor_dfd0(vector<DiffDoub0>& mat, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
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
	DiffDoub0 theta;
	DiffDoub0 sth;
	DiffDoub0 cth;
	DiffDoub0 p1;
	DiffDoub0 p2;
	DiffDoub0 tmp;

	for (i1 = st_col; i1 <= end_col; i1++) {
		i2_min = st_row + (i1 - st_col) + 1;
		if(tri_diag == 0) {
			i2_max = end_row;
		} else {
			i2_max = st_row + (i1 - st_col) + 1;
			if(i2_max > end_row) {
				i2_max = end_row;
			}
		}
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			k11 = (i2_min-1)*col_dim + i1;
			k12 = i2*col_dim + i1;
			if(abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.set_val_dfd0(mat[k12]);
			tmp.set_val_dfd0(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.set_val_dfd0(theta);
			sth.sn();
			cth.set_val_dfd0(theta);
			cth.cs();
			i3_min = i1;
			if(tri_diag == 2) {
				i3_max = i1 + 2;
				if(i3_max > end_col) {
					i3_max = end_col;
				}
			} else {
				i3_max = end_col;
			}
			for (i3 = i3_min; i3 <= i3_max; i3++) {
				k22 = (i2_min-1)*col_dim + i3;
				k23 = i2*col_dim + i3;
				p1.set_val_dfd0(cth);
				p1.mult(mat[k22]);
				tmp.set_val_dfd0(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.set_val_dfd0(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.set_val_dfd0(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].set_val_dfd0(p1);
				mat[k23].set_val_dfd0(p2);
			}
			mat[k12].set_val_dfd0(theta);
		}
	}
	return;
}

void q_rfactor_ar_dfd0(DiffDoub0 mat[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
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
	DiffDoub0 theta;
	DiffDoub0 sth;
	DiffDoub0 cth;
	DiffDoub0 p1;
	DiffDoub0 p2;
	DiffDoub0 tmp;

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
			if (abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.set_val_dfd0(mat[k12]);
			tmp.set_val_dfd0(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.set_val_dfd0(theta);
			sth.sn();
			cth.set_val_dfd0(theta);
			cth.cs();
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
				p1.set_val_dfd0(cth);
				p1.mult(mat[k22]);
				tmp.set_val_dfd0(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.set_val_dfd0(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.set_val_dfd0(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].set_val_dfd0(p1);
				mat[k23].set_val_dfd0(p2);
			}
			mat[k12].set_val_dfd0(theta);
		}
	}
	return;
}

void solveq_rx_eqb_dfd0(vector<DiffDoub0>& x_vec, vector<DiffDoub0>& mat, vector<DiffDoub0>& b_vec, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i2_min;
	int i2_max;
	int k11;
	int k12;
	DiffDoub0 theta;
	DiffDoub0 sth;
	DiffDoub0 cth;
	DiffDoub0 p1;
	DiffDoub0 p2;
	DiffDoub0 tmp;
	DiffDoub0 row_sum;

    for (i1 = st_col; i1 <= end_col; i1++) {
		i2_min = st_row + (i1 - st_col) + 1;
		if(tri_diag == 0) {
			i2_max = end_row;
		} else {
			i2_max = st_row + (i1 - st_col) + 1;
			if(i2_max > end_row) {
				i2_max = end_row;
			}
		}
		i3 = i2_min - 1;
		for (i2 = i2_min; i2 <=i2_max; i2++) {
			k12 = i2*col_dim + i1;
			theta.set_val_dfd0(mat[k12]);
			sth.set_val_dfd0(theta);
			sth.sn();
			cth.set_val_dfd0(theta);
			cth.cs();
			p1.set_val_dfd0(cth);
			p1.mult(b_vec[i3]);
			tmp.set_val_dfd0(sth);
			tmp.mult(b_vec[i2]);
			p1.add(tmp);
			p2.set_val_dfd0(sth);
			p2.neg();
			p2.mult(b_vec[i3]);
			tmp.set_val_dfd0(cth);
			tmp.mult(b_vec[i2]);
			p2.add(tmp);
		    b_vec[i3].set_val_dfd0(p1);
			b_vec[i2].set_val_dfd0(p2);
		}
		x_vec[i1].set_val(0.0);
	}
	
	
	for (i1 = end_col; i1 >= st_col; i1--) {
		i3 = st_row + (i1 - st_col);
		i2_min = i1 + 1;
		if(tri_diag == 2) {
			i2_max = i1 + 2;
			if(i2_max > end_col) {
				i2_max = end_col;
			}
		} else {
			i2_max = end_col;
		}
		row_sum.set_val(0.0);
		k11 = i3*col_dim + i2_min;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			tmp.set_val_dfd0(mat[k11]);
			tmp.mult(x_vec[i2]);
			row_sum.add(tmp);
			k11++;
		}
		tmp.set_val_dfd0(b_vec[i3]);
		tmp.sub(row_sum);
		k11 = i3*col_dim + i1;
		tmp.dvd(mat[k11]);
		x_vec[i1].set_val_dfd0(tmp);
	}
	
	return;
}

void solveq_rx_eqb_ar_dfd0(DiffDoub0 x_vec[], DiffDoub0 mat[], DiffDoub0 b_vec[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i2_min;
	int i2_max;
	int k11;
	int k12;
	DiffDoub0 theta;
	DiffDoub0 sth;
	DiffDoub0 cth;
	DiffDoub0 p1;
	DiffDoub0 p2;
	DiffDoub0 tmp;
	DiffDoub0 row_sum;

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
			theta.set_val_dfd0(mat[k12]);
			sth.set_val_dfd0(theta);
			sth.sn();
			cth.set_val_dfd0(theta);
			cth.cs();
			p1.set_val_dfd0(cth);
			p1.mult(b_vec[i3]);
			tmp.set_val_dfd0(sth);
			tmp.mult(b_vec[i2]);
			p1.add(tmp);
			p2.set_val_dfd0(sth);
			p2.neg();
			p2.mult(b_vec[i3]);
			tmp.set_val_dfd0(cth);
			tmp.mult(b_vec[i2]);
			p2.add(tmp);
			b_vec[i3].set_val_dfd0(p1);
			b_vec[i2].set_val_dfd0(p2);
		}
		x_vec[i1].set_val(0.0);
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
		row_sum.set_val(0.0);
		k11 = i3 * col_dim + i2_min;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			tmp.set_val_dfd0(mat[k11]);
			tmp.mult(x_vec[i2]);
			row_sum.add(tmp);
			k11++;
		}
		tmp.set_val_dfd0(b_vec[i3]);
		tmp.sub(row_sum);
		k11 = i3 * col_dim + i1;
		tmp.dvd(mat[k11]);
		x_vec[i1].set_val_dfd0(tmp);
	}

	return;
}

void get_det_inv_dfd0(DiffDoub0& det, vector<DiffDoub0>& inv, vector<DiffDoub0>& mat, int col_dim, int tri_diag, vector<DiffDoub0>& x_vec, vector<DiffDoub0>& b_vec) {
	q_rfactor_dfd0(mat, col_dim, 0, (col_dim-1), 0, (col_dim-1), tri_diag);
	int i1;
    int i2;
	int i3;
	det.set_val(1.0);
	for (i1 = 0; i1 < col_dim; i1++) {
		for (i2 = 0; i2 < col_dim; i2++) {
			if(i1 == i2) {
				b_vec[i2].set_val(1.0);
			} else {
			    b_vec[i2].set_val(0.0);
			}
		}
		solveq_rx_eqb_dfd0(x_vec,mat,b_vec,col_dim,0,(col_dim-1),0,(col_dim-1),tri_diag);
		for (i2 = 0; i2 < col_dim; i2++) {
			i3 = i2*col_dim + i1;
			inv[i3].set_val_dfd0(x_vec[i2]);
		}
		i3 = i1*col_dim + i1;
		det.mult(mat[i3]);
	}
	return;
}

void get_det_inv_ar_dfd0(DiffDoub0& det, DiffDoub0 inv[], DiffDoub0 mat[], int col_dim, int tri_diag, DiffDoub0 x_vec[], DiffDoub0 b_vec[]) {
	q_rfactor_ar_dfd0(mat, col_dim, 0, (col_dim - 1), 0, (col_dim - 1), tri_diag);
	int i1;
	int i2;
	int i3;
	det.set_val(1.0);
	for (i1 = 0; i1 < col_dim; i1++) {
		for (i2 = 0; i2 < col_dim; i2++) {
			if (i1 == i2) {
				b_vec[i2].set_val(1.0);
			}
			else {
				b_vec[i2].set_val(0.0);
			}
		}
		solveq_rx_eqb_ar_dfd0(x_vec, mat, b_vec, col_dim, 0, (col_dim - 1), 0, (col_dim - 1), tri_diag);
		for (i2 = 0; i2 < col_dim; i2++) {
			i3 = i2 * col_dim + i1;
			inv[i3].set_val_dfd0(x_vec[i2]);
		}
		i3 = i1 * col_dim + i1;
		det.mult(mat[i3]);
	}
	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void q_rfactor_dfd1(vector<DiffDoub1>& mat, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
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
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;

	for (i1 = st_col; i1 <= end_col; i1++) {
		i2_min = st_row + (i1 - st_col) + 1;
		if(tri_diag == 0) {
			i2_max = end_row;
		} else {
			i2_max = st_row + (i1 - st_col) + 1;
			if(i2_max > end_row) {
				i2_max = end_row;
			}
		}
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			k11 = (i2_min-1)*col_dim + i1;
			k12 = i2*col_dim + i1;
			if(abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.set_val_dfd1(mat[k12]);
			tmp.set_val_dfd1(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.set_val_dfd1(theta);
			sth.sn();
			cth.set_val_dfd1(theta);
			cth.cs();
			i3_min = i1;
			if(tri_diag == 2) {
				i3_max = i1 + 2;
				if(i3_max > end_col) {
					i3_max = end_col;
				}
			} else {
				i3_max = end_col;
			}
			for (i3 = i3_min; i3 <= i3_max; i3++) {
				k22 = (i2_min-1)*col_dim + i3;
				k23 = i2*col_dim + i3;
				p1.set_val_dfd1(cth);
				p1.mult(mat[k22]);
				tmp.set_val_dfd1(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.set_val_dfd1(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.set_val_dfd1(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].set_val_dfd1(p1);
				mat[k23].set_val_dfd1(p2);
			}
			mat[k12].set_val_dfd1(theta);
		}
	}
	return;
}

void q_rfactor_ar_dfd1(DiffDoub1 mat[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
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
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;

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
			if (abs(mat[k11].val) < tol) {
				mat[k11].val = tol;
			}
			theta.set_val_dfd1(mat[k12]);
			tmp.set_val_dfd1(mat[k11]);
			theta.dvd(tmp);
			theta.atn();
			sth.set_val_dfd1(theta);
			sth.sn();
			cth.set_val_dfd1(theta);
			cth.cs();
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
				p1.set_val_dfd1(cth);
				p1.mult(mat[k22]);
				tmp.set_val_dfd1(sth);
				tmp.mult(mat[k23]);
				p1.add(tmp);
				p2.set_val_dfd1(sth);
				p2.neg();
				p2.mult(mat[k22]);
				tmp.set_val_dfd1(cth);
				tmp.mult(mat[k23]);
				p2.add(tmp);
				mat[k22].set_val_dfd1(p1);
				mat[k23].set_val_dfd1(p2);
			}
			mat[k12].set_val_dfd1(theta);
		}
	}
	return;
}

void solveq_rx_eqb_dfd1(vector<DiffDoub1>& x_vec, vector<DiffDoub1>& mat, vector<DiffDoub1>& b_vec, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i2_min;
	int i2_max;
	int k11;
	int k12;
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;
	DiffDoub1 row_sum;

    for (i1 = st_col; i1 <= end_col; i1++) {
		i2_min = st_row + (i1 - st_col) + 1;
		if(tri_diag == 0) {
			i2_max = end_row;
		} else {
			i2_max = st_row + (i1 - st_col) + 1;
			if(i2_max > end_row) {
				i2_max = end_row;
			}
		}
		i3 = i2_min - 1;
		for (i2 = i2_min; i2 <=i2_max; i2++) {
			k12 = i2*col_dim + i1;
			theta.set_val_dfd1(mat[k12]);
			sth.set_val_dfd1(theta);
			sth.sn();
			cth.set_val_dfd1(theta);
			cth.cs();
			p1.set_val_dfd1(cth);
			p1.mult(b_vec[i3]);
			tmp.set_val_dfd1(sth);
			tmp.mult(b_vec[i2]);
			p1.add(tmp);
			p2.set_val_dfd1(sth);
			p2.neg();
			p2.mult(b_vec[i3]);
			tmp.set_val_dfd1(cth);
			tmp.mult(b_vec[i2]);
			p2.add(tmp);
		    b_vec[i3].set_val_dfd1(p1);
			b_vec[i2].set_val_dfd1(p2);
		}
		x_vec[i1].set_val(0.0);
	}
	
	
	for (i1 = end_col; i1 >= st_col; i1--) {
		i3 = st_row + (i1 - st_col);
		i2_min = i1 + 1;
		if(tri_diag == 2) {
			i2_max = i1 + 2;
			if(i2_max > end_col) {
				i2_max = end_col;
			}
		} else {
			i2_max = end_col;
		}
		row_sum.set_val(0.0);
		k11 = i3*col_dim + i2_min;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			tmp.set_val_dfd1(mat[k11]);
			tmp.mult(x_vec[i2]);
			row_sum.add(tmp);
			k11++;
		}
		tmp.set_val_dfd1(b_vec[i3]);
		tmp.sub(row_sum);
		k11 = i3*col_dim + i1;
		tmp.dvd(mat[k11]);
		x_vec[i1].set_val_dfd1(tmp);
	}
	
	return;
}

void solveq_rx_eqb_ar_dfd1(DiffDoub1 x_vec[], DiffDoub1 mat[], DiffDoub1 b_vec[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag) {
	int i1;
	int i2;
	int i3;
	int i2_min;
	int i2_max;
	int k11;
	int k12;
	DiffDoub1 theta;
	DiffDoub1 sth;
	DiffDoub1 cth;
	DiffDoub1 p1;
	DiffDoub1 p2;
	DiffDoub1 tmp;
	DiffDoub1 row_sum;

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
			theta.set_val_dfd1(mat[k12]);
			sth.set_val_dfd1(theta);
			sth.sn();
			cth.set_val_dfd1(theta);
			cth.cs();
			p1.set_val_dfd1(cth);
			p1.mult(b_vec[i3]);
			tmp.set_val_dfd1(sth);
			tmp.mult(b_vec[i2]);
			p1.add(tmp);
			p2.set_val_dfd1(sth);
			p2.neg();
			p2.mult(b_vec[i3]);
			tmp.set_val_dfd1(cth);
			tmp.mult(b_vec[i2]);
			p2.add(tmp);
			b_vec[i3].set_val_dfd1(p1);
			b_vec[i2].set_val_dfd1(p2);
		}
		x_vec[i1].set_val(0.0);
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
		row_sum.set_val(0.0);
		k11 = i3 * col_dim + i2_min;
		for (i2 = i2_min; i2 <= i2_max; i2++) {
			tmp.set_val_dfd1(mat[k11]);
			tmp.mult(x_vec[i2]);
			row_sum.add(tmp);
			k11++;
		}
		tmp.set_val_dfd1(b_vec[i3]);
		tmp.sub(row_sum);
		k11 = i3 * col_dim + i1;
		tmp.dvd(mat[k11]);
		x_vec[i1].set_val_dfd1(tmp);
	}

	return;
}

void get_det_inv_dfd1(DiffDoub1& det, vector<DiffDoub1>& inv, vector<DiffDoub1>& mat, int col_dim, int tri_diag, vector<DiffDoub1>& x_vec, vector<DiffDoub1>& b_vec) {
	q_rfactor_dfd1(mat, col_dim, 0, (col_dim-1), 0, (col_dim-1), tri_diag);
	int i1;
    int i2;
	int i3;
	det.set_val(1.0);
	for (i1 = 0; i1 < col_dim; i1++) {
		for (i2 = 0; i2 < col_dim; i2++) {
			if(i1 == i2) {
				b_vec[i2].set_val(1.0);
			} else {
			    b_vec[i2].set_val(0.0);
			}
		}
		solveq_rx_eqb_dfd1(x_vec,mat,b_vec,col_dim,0,(col_dim-1),0,(col_dim-1),tri_diag);
		for (i2 = 0; i2 < col_dim; i2++) {
			i3 = i2*col_dim + i1;
			inv[i3].set_val_dfd1(x_vec[i2]);
		}
		i3 = i1*col_dim + i1;
		det.mult(mat[i3]);
	}
	return;
}

void get_det_inv_ar_dfd1(DiffDoub1& det, DiffDoub1 inv[], DiffDoub1 mat[], int col_dim, int tri_diag, DiffDoub1 x_vec[], DiffDoub1 b_vec[]) {
	q_rfactor_ar_dfd1(mat, col_dim, 0, (col_dim - 1), 0, (col_dim - 1), tri_diag);
	int i1;
	int i2;
	int i3;
	det.set_val(1.0);
	for (i1 = 0; i1 < col_dim; i1++) {
		for (i2 = 0; i2 < col_dim; i2++) {
			if (i1 == i2) {
				b_vec[i2].set_val(1.0);
			}
			else {
				b_vec[i2].set_val(0.0);
			}
		}
		solveq_rx_eqb_ar_dfd1(x_vec, mat, b_vec, col_dim, 0, (col_dim - 1), 0, (col_dim - 1), tri_diag);
		for (i2 = 0; i2 < col_dim; i2++) {
			i3 = i2 * col_dim + i1;
			inv[i3].set_val_dfd1(x_vec[i2]);
		}
		i3 = i1 * col_dim + i1;
		det.mult(mat[i3]);
	}
	return;
}

//end dup
 
//end skip 
 
 
 
//dup2
void sub_vec_dfd0(vector<DiffDoub0>& sub_v, vector<DiffDoub0>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		sub_v[i2].set_val_dfd0(v_in[i1]);
		i2++;
	}
}

void return_sv_dfd0(vector<DiffDoub0>& sub_v, vector<DiffDoub0>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		v_in[i1].set_val_dfd0(sub_v[i2]);
		i2++;
	}
}

void vec_to_ar_dfd0(DiffDoub0 ar[], std::vector<DiffDoub0>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		ar[i2].set_val_dfd0(vc[i1]);
		i2++;
	}
}

void ar_to_vec_dfd0(DiffDoub0 ar[], std::vector<DiffDoub0>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		vc[i1].set_val_dfd0(ar[i2]);
		i2++;
	}
}

void mat_mul_dfd0(vector<DiffDoub0>& prod, vector<DiffDoub0>& mat1, vector<DiffDoub0>& mat2, int m1_rows, int m1_cols, int m2_cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub0 tmp;
	
	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1_rows; i1++) {
		for (i2 = 0; i2 < m2_cols; i2++) {
			i6 = i2;
			prod[i4].set_val(0.0);
			for (i3 = 0; i3 < m1_cols; i3++) {
				tmp.set_val_dfd0(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2_cols;
			}
			i5 -= m1_cols;
			i4++;
		}
		i5 += m1_cols;
	}
	return;
}

void mat_mul_ar_dfd0(DiffDoub0 prod[], DiffDoub0 mat1[], DiffDoub0 mat2[], int m1_rows, int m1_cols, int m2_cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub0 tmp;

	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1_rows; i1++) {
		for (i2 = 0; i2 < m2_cols; i2++) {
			i6 = i2;
			prod[i4].set_val(0.0);
			for (i3 = 0; i3 < m1_cols; i3++) {
				tmp.set_val_dfd0(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6 += m2_cols;
			}
			i5 -= m1_cols;
			i4++;
		}
		i5 += m1_cols;
	}
	return;
}

void transpose_dfd0(vector<DiffDoub0>& mat_t, vector<DiffDoub0>& mat, int row_dim, int col_dim) {
	DiffDoub0 tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < row_dim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < col_dim; i2++) {
			mat_t[i4].set_val_dfd0(mat[i3]);
			i3++;
			i4+= row_dim;
		}
	}
	return;
}

void transpose_ar_dfd0(DiffDoub0 mat_t[], DiffDoub0 mat[], int row_dim, int col_dim) {
	DiffDoub0 tmp;
	int i1;
	int i2;
	int i3;
	int i4;

	i3 = 0;
	for (i1 = 0; i1 < row_dim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < col_dim; i2++) {
			mat_t[i4].set_val_dfd0(mat[i3]);
			i3++;
			i4 += row_dim;
		}
	}
	return;
}

void cross_prod_dfd0(DiffDoub0 prod[], DiffDoub0 v1[], DiffDoub0 v2[]) {
	DiffDoub0 tmp;
	
	prod[0].set_val_dfd0(v1[1]);
	prod[0].mult(v2[2]);
	tmp.set_val_dfd0(v1[2]);
	tmp.mult(v2[1]);
	prod[0].sub(tmp);
	
	prod[1].set_val_dfd0(v1[2]);
	prod[1].mult(v2[0]);
	tmp.set_val_dfd0(v1[0]);
	tmp.mult(v2[2]);
	prod[1].sub(tmp);
	
	prod[2].set_val_dfd0(v1[0]);
	prod[2].mult(v2[1]);
	tmp.set_val_dfd0(v1[1]);
	tmp.mult(v2[0]);
	prod[2].sub(tmp);
	return;
}

void rotate_orient_dfd0(DiffDoub0 inst_ori[], DiffDoub0 loc_ori[], DiffDoub0 rot[]) {
	DiffDoub0 mag;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 one_half;
	DiffDoub0 a1[9];
	DiffDoub0 a2[9];
	DiffDoub0 a3[9];
	int i1;
	int i2;
	int i3;
	
	one_half.set_val(0.5);
	
	mag.set_val_dfd0(rot[0]);
	mag.sqr();
	tmp.set_val_dfd0(rot[1]);
	tmp.sqr();
	tmp2.set_val_dfd0(rot[2]);
	tmp2.sqr();
	mag.add(tmp);
	mag.add(tmp2);
	mag.sqt();
	if(mag.val < magtol) {
		DiffDoub0 loc_rot[3];
		i3 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			loc_rot[i1].set_val(0.0);
			for (i2 = 0; i2 < 3; i2++) {
				tmp.set_val_dfd0(loc_ori[i3]);
				tmp.mult(rot[i2]);
				loc_rot[i1].add(tmp);
				i3++;
			}
		}
		a1[0].set_val(1.0);
		tmp.set_val_dfd0(loc_rot[1]);
		tmp.sqr();
		tmp2.set_val_dfd0(loc_rot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[0].sub(tmp);
		
		a1[1].set_val(0.5);
		a1[1].mult(loc_rot[0]);
		a1[1].mult(loc_rot[1]);
		a1[1].add(loc_rot[2]);
		
		a1[2].set_val(0.5);
		a1[2].mult(loc_rot[0]);
		a1[2].mult(loc_rot[2]);
		a1[2].sub(loc_rot[1]);
		
		a1[3].set_val(0.5);
		a1[3].mult(loc_rot[0]);
		a1[3].mult(loc_rot[1]);
		a1[3].sub(loc_rot[2]);
		
		a1[4].set_val(1.0);
		tmp.set_val_dfd0(loc_rot[0]);
		tmp.sqr();
		tmp2.set_val_dfd0(loc_rot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[4].sub(tmp);
		
		a1[5].set_val(0.5);
		a1[5].mult(loc_rot[1]);
		a1[5].mult(loc_rot[2]);
		a1[5].add(loc_rot[0]);
		
		a1[6].set_val(0.5);
		a1[6].mult(loc_rot[0]);
		a1[6].mult(loc_rot[2]);
		a1[6].add(loc_rot[1]);
		
		a1[7].set_val(0.5);
		a1[7].mult(loc_rot[1]);
		a1[7].mult(loc_rot[2]);
		a1[7].sub(loc_rot[0]);
		
		a1[8].set_val(1.0);
		tmp.set_val_dfd0(loc_rot[0]);
		tmp.sqr();
		tmp2.set_val_dfd0(loc_rot[1]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[8].sub(tmp);
		
		mat_mul_ar_dfd0(inst_ori, a1, loc_ori, 3, 3, 3);
	} else {
		DiffDoub0 sth;
		DiffDoub0 cth;
		
		sth.set_val_dfd0(mag);
		sth.sn();
		cth.set_val_dfd0(mag);
		cth.cs();
		
		DiffDoub0 unit_rot[3];
		tmp.set_val(1.0);
		tmp.dvd(mag);
		unit_rot[0].set_val_dfd0(rot[0]);
		unit_rot[0].mult(tmp);
		unit_rot[1].set_val_dfd0(rot[1]);
		unit_rot[1].mult(tmp);
		unit_rot[2].set_val_dfd0(rot[2]);
		unit_rot[2].mult(tmp);
		
		a1[0].set_val_dfd0(unit_rot[0]);
		a1[1].set_val_dfd0(unit_rot[1]);
		a1[2].set_val_dfd0(unit_rot[2]);
		
		i1 = 0;
		if(abs(unit_rot[1].val) < abs(unit_rot[0].val)) {
			i1 = 1;
		}
		if(abs(unit_rot[2].val) < abs(unit_rot[i1].val)) {
			i1 = 2;
		}
		tmp.set_val(1.0);
		tmp2.set_val_dfd0(unit_rot[i1]);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp.sqt(); // tmp = sqrt(1 - unit_rot[i1]^2)
		a1[3+i1].set_val_dfd0(tmp);
		for (i2 = 0; i2 < 3; i2++) {
			if(i2 != i1) {
				i3 = 3 + i2;
				a1[i3].set_val_dfd0(a1[i1]);
				a1[i3].neg();
				a1[i3].mult(a1[i2]);
				a1[i3].dvd(tmp);
			}
		}
		
		cross_prod_dfd0(&a1[6], &a1[0], &a1[3]);
		
		a2[0].set_val_dfd0(a1[0]);
		a2[1].set_val_dfd0(a1[1]);
		a2[2].set_val_dfd0(a1[2]);
		
		a2[3].set_val_dfd0(cth);
		a2[3].mult(a1[3]);
		tmp.set_val_dfd0(sth);
		tmp.mult(a1[6]);
		a2[3].sub(tmp);
		
		a2[4].set_val_dfd0(cth);
		a2[4].mult(a1[4]);
		tmp.set_val_dfd0(sth);
		tmp.mult(a1[7]);
		a2[4].sub(tmp);
		
		a2[5].set_val_dfd0(cth);
		a2[5].mult(a1[5]);
		tmp.set_val_dfd0(sth);
		tmp.mult(a1[8]);
		a2[5].sub(tmp);
		
		a2[6].set_val_dfd0(sth);
		a2[6].mult(a1[3]);
		tmp.set_val_dfd0(cth);
		tmp.mult(a1[6]);
		a2[6].add(tmp);
		
		a2[7].set_val_dfd0(sth);
		a2[7].mult(a1[4]);
		tmp.set_val_dfd0(cth);
		tmp.mult(a1[7]);
		a2[7].add(tmp);
		
		a2[8].set_val_dfd0(sth);
		a2[8].mult(a1[5]);
		tmp.set_val_dfd0(cth);
		tmp.mult(a1[8]);
		a2[8].add(tmp);
		
		transpose_ar_dfd0(a3,a2,3,3); // a3 = a2^t
		mat_mul_ar_dfd0(a2,a3,a1,3,3,3); //a2 = a2^t*a1
		mat_mul_ar_dfd0(inst_ori,loc_ori,a2,3,3,3);
	}
	return;
}
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup2
void sub_vec_dfd1(vector<DiffDoub1>& sub_v, vector<DiffDoub1>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		sub_v[i2].set_val_dfd1(v_in[i1]);
		i2++;
	}
}

void return_sv_dfd1(vector<DiffDoub1>& sub_v, vector<DiffDoub1>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		v_in[i1].set_val_dfd1(sub_v[i2]);
		i2++;
	}
}

void vec_to_ar_dfd1(DiffDoub1 ar[], std::vector<DiffDoub1>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		ar[i2].set_val_dfd1(vc[i1]);
		i2++;
	}
}

void ar_to_vec_dfd1(DiffDoub1 ar[], std::vector<DiffDoub1>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		vc[i1].set_val_dfd1(ar[i2]);
		i2++;
	}
}

void mat_mul_dfd1(vector<DiffDoub1>& prod, vector<DiffDoub1>& mat1, vector<DiffDoub1>& mat2, int m1_rows, int m1_cols, int m2_cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub1 tmp;
	
	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1_rows; i1++) {
		for (i2 = 0; i2 < m2_cols; i2++) {
			i6 = i2;
			prod[i4].set_val(0.0);
			for (i3 = 0; i3 < m1_cols; i3++) {
				tmp.set_val_dfd1(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2_cols;
			}
			i5 -= m1_cols;
			i4++;
		}
		i5 += m1_cols;
	}
	return;
}

void mat_mul_ar_dfd1(DiffDoub1 prod[], DiffDoub1 mat1[], DiffDoub1 mat2[], int m1_rows, int m1_cols, int m2_cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub1 tmp;

	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1_rows; i1++) {
		for (i2 = 0; i2 < m2_cols; i2++) {
			i6 = i2;
			prod[i4].set_val(0.0);
			for (i3 = 0; i3 < m1_cols; i3++) {
				tmp.set_val_dfd1(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6 += m2_cols;
			}
			i5 -= m1_cols;
			i4++;
		}
		i5 += m1_cols;
	}
	return;
}

void transpose_dfd1(vector<DiffDoub1>& mat_t, vector<DiffDoub1>& mat, int row_dim, int col_dim) {
	DiffDoub1 tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < row_dim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < col_dim; i2++) {
			mat_t[i4].set_val_dfd1(mat[i3]);
			i3++;
			i4+= row_dim;
		}
	}
	return;
}

void transpose_ar_dfd1(DiffDoub1 mat_t[], DiffDoub1 mat[], int row_dim, int col_dim) {
	DiffDoub1 tmp;
	int i1;
	int i2;
	int i3;
	int i4;

	i3 = 0;
	for (i1 = 0; i1 < row_dim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < col_dim; i2++) {
			mat_t[i4].set_val_dfd1(mat[i3]);
			i3++;
			i4 += row_dim;
		}
	}
	return;
}

void cross_prod_dfd1(DiffDoub1 prod[], DiffDoub1 v1[], DiffDoub1 v2[]) {
	DiffDoub1 tmp;
	
	prod[0].set_val_dfd1(v1[1]);
	prod[0].mult(v2[2]);
	tmp.set_val_dfd1(v1[2]);
	tmp.mult(v2[1]);
	prod[0].sub(tmp);
	
	prod[1].set_val_dfd1(v1[2]);
	prod[1].mult(v2[0]);
	tmp.set_val_dfd1(v1[0]);
	tmp.mult(v2[2]);
	prod[1].sub(tmp);
	
	prod[2].set_val_dfd1(v1[0]);
	prod[2].mult(v2[1]);
	tmp.set_val_dfd1(v1[1]);
	tmp.mult(v2[0]);
	prod[2].sub(tmp);
	return;
}

void rotate_orient_dfd1(DiffDoub1 inst_ori[], DiffDoub1 loc_ori[], DiffDoub1 rot[]) {
	DiffDoub1 mag;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 one_half;
	DiffDoub1 a1[9];
	DiffDoub1 a2[9];
	DiffDoub1 a3[9];
	int i1;
	int i2;
	int i3;
	
	one_half.set_val(0.5);
	
	mag.set_val_dfd1(rot[0]);
	mag.sqr();
	tmp.set_val_dfd1(rot[1]);
	tmp.sqr();
	tmp2.set_val_dfd1(rot[2]);
	tmp2.sqr();
	mag.add(tmp);
	mag.add(tmp2);
	mag.sqt();
	if(mag.val < magtol) {
		DiffDoub1 loc_rot[3];
		i3 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			loc_rot[i1].set_val(0.0);
			for (i2 = 0; i2 < 3; i2++) {
				tmp.set_val_dfd1(loc_ori[i3]);
				tmp.mult(rot[i2]);
				loc_rot[i1].add(tmp);
				i3++;
			}
		}
		a1[0].set_val(1.0);
		tmp.set_val_dfd1(loc_rot[1]);
		tmp.sqr();
		tmp2.set_val_dfd1(loc_rot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[0].sub(tmp);
		
		a1[1].set_val(0.5);
		a1[1].mult(loc_rot[0]);
		a1[1].mult(loc_rot[1]);
		a1[1].add(loc_rot[2]);
		
		a1[2].set_val(0.5);
		a1[2].mult(loc_rot[0]);
		a1[2].mult(loc_rot[2]);
		a1[2].sub(loc_rot[1]);
		
		a1[3].set_val(0.5);
		a1[3].mult(loc_rot[0]);
		a1[3].mult(loc_rot[1]);
		a1[3].sub(loc_rot[2]);
		
		a1[4].set_val(1.0);
		tmp.set_val_dfd1(loc_rot[0]);
		tmp.sqr();
		tmp2.set_val_dfd1(loc_rot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[4].sub(tmp);
		
		a1[5].set_val(0.5);
		a1[5].mult(loc_rot[1]);
		a1[5].mult(loc_rot[2]);
		a1[5].add(loc_rot[0]);
		
		a1[6].set_val(0.5);
		a1[6].mult(loc_rot[0]);
		a1[6].mult(loc_rot[2]);
		a1[6].add(loc_rot[1]);
		
		a1[7].set_val(0.5);
		a1[7].mult(loc_rot[1]);
		a1[7].mult(loc_rot[2]);
		a1[7].sub(loc_rot[0]);
		
		a1[8].set_val(1.0);
		tmp.set_val_dfd1(loc_rot[0]);
		tmp.sqr();
		tmp2.set_val_dfd1(loc_rot[1]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[8].sub(tmp);
		
		mat_mul_ar_dfd1(inst_ori, a1, loc_ori, 3, 3, 3);
	} else {
		DiffDoub1 sth;
		DiffDoub1 cth;
		
		sth.set_val_dfd1(mag);
		sth.sn();
		cth.set_val_dfd1(mag);
		cth.cs();
		
		DiffDoub1 unit_rot[3];
		tmp.set_val(1.0);
		tmp.dvd(mag);
		unit_rot[0].set_val_dfd1(rot[0]);
		unit_rot[0].mult(tmp);
		unit_rot[1].set_val_dfd1(rot[1]);
		unit_rot[1].mult(tmp);
		unit_rot[2].set_val_dfd1(rot[2]);
		unit_rot[2].mult(tmp);
		
		a1[0].set_val_dfd1(unit_rot[0]);
		a1[1].set_val_dfd1(unit_rot[1]);
		a1[2].set_val_dfd1(unit_rot[2]);
		
		i1 = 0;
		if(abs(unit_rot[1].val) < abs(unit_rot[0].val)) {
			i1 = 1;
		}
		if(abs(unit_rot[2].val) < abs(unit_rot[i1].val)) {
			i1 = 2;
		}
		tmp.set_val(1.0);
		tmp2.set_val_dfd1(unit_rot[i1]);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp.sqt(); // tmp = sqrt(1 - unit_rot[i1]^2)
		a1[3+i1].set_val_dfd1(tmp);
		for (i2 = 0; i2 < 3; i2++) {
			if(i2 != i1) {
				i3 = 3 + i2;
				a1[i3].set_val_dfd1(a1[i1]);
				a1[i3].neg();
				a1[i3].mult(a1[i2]);
				a1[i3].dvd(tmp);
			}
		}
		
		cross_prod_dfd1(&a1[6], &a1[0], &a1[3]);
		
		a2[0].set_val_dfd1(a1[0]);
		a2[1].set_val_dfd1(a1[1]);
		a2[2].set_val_dfd1(a1[2]);
		
		a2[3].set_val_dfd1(cth);
		a2[3].mult(a1[3]);
		tmp.set_val_dfd1(sth);
		tmp.mult(a1[6]);
		a2[3].sub(tmp);
		
		a2[4].set_val_dfd1(cth);
		a2[4].mult(a1[4]);
		tmp.set_val_dfd1(sth);
		tmp.mult(a1[7]);
		a2[4].sub(tmp);
		
		a2[5].set_val_dfd1(cth);
		a2[5].mult(a1[5]);
		tmp.set_val_dfd1(sth);
		tmp.mult(a1[8]);
		a2[5].sub(tmp);
		
		a2[6].set_val_dfd1(sth);
		a2[6].mult(a1[3]);
		tmp.set_val_dfd1(cth);
		tmp.mult(a1[6]);
		a2[6].add(tmp);
		
		a2[7].set_val_dfd1(sth);
		a2[7].mult(a1[4]);
		tmp.set_val_dfd1(cth);
		tmp.mult(a1[7]);
		a2[7].add(tmp);
		
		a2[8].set_val_dfd1(sth);
		a2[8].mult(a1[5]);
		tmp.set_val_dfd1(cth);
		tmp.mult(a1[8]);
		a2[8].add(tmp);
		
		transpose_ar_dfd1(a3,a2,3,3); // a3 = a2^t
		mat_mul_ar_dfd1(a2,a3,a1,3,3,3); //a2 = a2^t*a1
		mat_mul_ar_dfd1(inst_ori,loc_ori,a2,3,3,3);
	}
	return;
}
//end dup
 
//DiffDoub2 versions: 
//dup2
void sub_vec_dfd2(vector<DiffDoub2>& sub_v, vector<DiffDoub2>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		sub_v[i2].set_val_dfd2(v_in[i1]);
		i2++;
	}
}

void return_sv_dfd2(vector<DiffDoub2>& sub_v, vector<DiffDoub2>& v_in, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		v_in[i1].set_val_dfd2(sub_v[i2]);
		i2++;
	}
}

void vec_to_ar_dfd2(DiffDoub2 ar[], std::vector<DiffDoub2>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		ar[i2].set_val_dfd2(vc[i1]);
		i2++;
	}
}

void ar_to_vec_dfd2(DiffDoub2 ar[], std::vector<DiffDoub2>& vc, int st, int end) {
	int i1;
	int i2 = 0;
	for (i1 = st; i1 < end; i1++) {
		vc[i1].set_val_dfd2(ar[i2]);
		i2++;
	}
}

void mat_mul_dfd2(vector<DiffDoub2>& prod, vector<DiffDoub2>& mat1, vector<DiffDoub2>& mat2, int m1_rows, int m1_cols, int m2_cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub2 tmp;
	
	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1_rows; i1++) {
		for (i2 = 0; i2 < m2_cols; i2++) {
			i6 = i2;
			prod[i4].set_val(0.0);
			for (i3 = 0; i3 < m1_cols; i3++) {
				tmp.set_val_dfd2(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6+= m2_cols;
			}
			i5 -= m1_cols;
			i4++;
		}
		i5 += m1_cols;
	}
	return;
}

void mat_mul_ar_dfd2(DiffDoub2 prod[], DiffDoub2 mat1[], DiffDoub2 mat2[], int m1_rows, int m1_cols, int m2_cols) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	DiffDoub2 tmp;

	i4 = 0;
	i5 = 0;
	for (i1 = 0; i1 < m1_rows; i1++) {
		for (i2 = 0; i2 < m2_cols; i2++) {
			i6 = i2;
			prod[i4].set_val(0.0);
			for (i3 = 0; i3 < m1_cols; i3++) {
				tmp.set_val_dfd2(mat1[i5]);
				tmp.mult(mat2[i6]);
				prod[i4].add(tmp);
				i5++;
				i6 += m2_cols;
			}
			i5 -= m1_cols;
			i4++;
		}
		i5 += m1_cols;
	}
	return;
}

void transpose_dfd2(vector<DiffDoub2>& mat_t, vector<DiffDoub2>& mat, int row_dim, int col_dim) {
	DiffDoub2 tmp;
	int i1;
	int i2;
	int i3;
	int i4;
	
	i3 = 0;
	for (i1 = 0; i1 < row_dim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < col_dim; i2++) {
			mat_t[i4].set_val_dfd2(mat[i3]);
			i3++;
			i4+= row_dim;
		}
	}
	return;
}

void transpose_ar_dfd2(DiffDoub2 mat_t[], DiffDoub2 mat[], int row_dim, int col_dim) {
	DiffDoub2 tmp;
	int i1;
	int i2;
	int i3;
	int i4;

	i3 = 0;
	for (i1 = 0; i1 < row_dim; i1++) {
		i4 = i1;
		for (i2 = 0; i2 < col_dim; i2++) {
			mat_t[i4].set_val_dfd2(mat[i3]);
			i3++;
			i4 += row_dim;
		}
	}
	return;
}

void cross_prod_dfd2(DiffDoub2 prod[], DiffDoub2 v1[], DiffDoub2 v2[]) {
	DiffDoub2 tmp;
	
	prod[0].set_val_dfd2(v1[1]);
	prod[0].mult(v2[2]);
	tmp.set_val_dfd2(v1[2]);
	tmp.mult(v2[1]);
	prod[0].sub(tmp);
	
	prod[1].set_val_dfd2(v1[2]);
	prod[1].mult(v2[0]);
	tmp.set_val_dfd2(v1[0]);
	tmp.mult(v2[2]);
	prod[1].sub(tmp);
	
	prod[2].set_val_dfd2(v1[0]);
	prod[2].mult(v2[1]);
	tmp.set_val_dfd2(v1[1]);
	tmp.mult(v2[0]);
	prod[2].sub(tmp);
	return;
}

void rotate_orient_dfd2(DiffDoub2 inst_ori[], DiffDoub2 loc_ori[], DiffDoub2 rot[]) {
	DiffDoub2 mag;
	DiffDoub2 tmp;
	DiffDoub2 tmp2;
	DiffDoub2 one_half;
	DiffDoub2 a1[9];
	DiffDoub2 a2[9];
	DiffDoub2 a3[9];
	int i1;
	int i2;
	int i3;
	
	one_half.set_val(0.5);
	
	mag.set_val_dfd2(rot[0]);
	mag.sqr();
	tmp.set_val_dfd2(rot[1]);
	tmp.sqr();
	tmp2.set_val_dfd2(rot[2]);
	tmp2.sqr();
	mag.add(tmp);
	mag.add(tmp2);
	mag.sqt();
	if(mag.val < magtol) {
		DiffDoub2 loc_rot[3];
		i3 = 0;
		for (i1 = 0; i1 < 3; i1++) {
			loc_rot[i1].set_val(0.0);
			for (i2 = 0; i2 < 3; i2++) {
				tmp.set_val_dfd2(loc_ori[i3]);
				tmp.mult(rot[i2]);
				loc_rot[i1].add(tmp);
				i3++;
			}
		}
		a1[0].set_val(1.0);
		tmp.set_val_dfd2(loc_rot[1]);
		tmp.sqr();
		tmp2.set_val_dfd2(loc_rot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[0].sub(tmp);
		
		a1[1].set_val(0.5);
		a1[1].mult(loc_rot[0]);
		a1[1].mult(loc_rot[1]);
		a1[1].add(loc_rot[2]);
		
		a1[2].set_val(0.5);
		a1[2].mult(loc_rot[0]);
		a1[2].mult(loc_rot[2]);
		a1[2].sub(loc_rot[1]);
		
		a1[3].set_val(0.5);
		a1[3].mult(loc_rot[0]);
		a1[3].mult(loc_rot[1]);
		a1[3].sub(loc_rot[2]);
		
		a1[4].set_val(1.0);
		tmp.set_val_dfd2(loc_rot[0]);
		tmp.sqr();
		tmp2.set_val_dfd2(loc_rot[2]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[4].sub(tmp);
		
		a1[5].set_val(0.5);
		a1[5].mult(loc_rot[1]);
		a1[5].mult(loc_rot[2]);
		a1[5].add(loc_rot[0]);
		
		a1[6].set_val(0.5);
		a1[6].mult(loc_rot[0]);
		a1[6].mult(loc_rot[2]);
		a1[6].add(loc_rot[1]);
		
		a1[7].set_val(0.5);
		a1[7].mult(loc_rot[1]);
		a1[7].mult(loc_rot[2]);
		a1[7].sub(loc_rot[0]);
		
		a1[8].set_val(1.0);
		tmp.set_val_dfd2(loc_rot[0]);
		tmp.sqr();
		tmp2.set_val_dfd2(loc_rot[1]);
		tmp2.sqr();
		tmp.add(tmp2);
		tmp.mult(one_half);
		a1[8].sub(tmp);
		
		mat_mul_ar_dfd2(inst_ori, a1, loc_ori, 3, 3, 3);
	} else {
		DiffDoub2 sth;
		DiffDoub2 cth;
		
		sth.set_val_dfd2(mag);
		sth.sn();
		cth.set_val_dfd2(mag);
		cth.cs();
		
		DiffDoub2 unit_rot[3];
		tmp.set_val(1.0);
		tmp.dvd(mag);
		unit_rot[0].set_val_dfd2(rot[0]);
		unit_rot[0].mult(tmp);
		unit_rot[1].set_val_dfd2(rot[1]);
		unit_rot[1].mult(tmp);
		unit_rot[2].set_val_dfd2(rot[2]);
		unit_rot[2].mult(tmp);
		
		a1[0].set_val_dfd2(unit_rot[0]);
		a1[1].set_val_dfd2(unit_rot[1]);
		a1[2].set_val_dfd2(unit_rot[2]);
		
		i1 = 0;
		if(abs(unit_rot[1].val) < abs(unit_rot[0].val)) {
			i1 = 1;
		}
		if(abs(unit_rot[2].val) < abs(unit_rot[i1].val)) {
			i1 = 2;
		}
		tmp.set_val(1.0);
		tmp2.set_val_dfd2(unit_rot[i1]);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp.sqt(); // tmp = sqrt(1 - unit_rot[i1]^2)
		a1[3+i1].set_val_dfd2(tmp);
		for (i2 = 0; i2 < 3; i2++) {
			if(i2 != i1) {
				i3 = 3 + i2;
				a1[i3].set_val_dfd2(a1[i1]);
				a1[i3].neg();
				a1[i3].mult(a1[i2]);
				a1[i3].dvd(tmp);
			}
		}
		
		cross_prod_dfd2(&a1[6], &a1[0], &a1[3]);
		
		a2[0].set_val_dfd2(a1[0]);
		a2[1].set_val_dfd2(a1[1]);
		a2[2].set_val_dfd2(a1[2]);
		
		a2[3].set_val_dfd2(cth);
		a2[3].mult(a1[3]);
		tmp.set_val_dfd2(sth);
		tmp.mult(a1[6]);
		a2[3].sub(tmp);
		
		a2[4].set_val_dfd2(cth);
		a2[4].mult(a1[4]);
		tmp.set_val_dfd2(sth);
		tmp.mult(a1[7]);
		a2[4].sub(tmp);
		
		a2[5].set_val_dfd2(cth);
		a2[5].mult(a1[5]);
		tmp.set_val_dfd2(sth);
		tmp.mult(a1[8]);
		a2[5].sub(tmp);
		
		a2[6].set_val_dfd2(sth);
		a2[6].mult(a1[3]);
		tmp.set_val_dfd2(cth);
		tmp.mult(a1[6]);
		a2[6].add(tmp);
		
		a2[7].set_val_dfd2(sth);
		a2[7].mult(a1[4]);
		tmp.set_val_dfd2(cth);
		tmp.mult(a1[7]);
		a2[7].add(tmp);
		
		a2[8].set_val_dfd2(sth);
		a2[8].mult(a1[5]);
		tmp.set_val_dfd2(cth);
		tmp.mult(a1[8]);
		a2[8].add(tmp);
		
		transpose_ar_dfd2(a3,a2,3,3); // a3 = a2^t
		mat_mul_ar_dfd2(a2,a3,a1,3,3,3); //a2 = a2^t*a1
		mat_mul_ar_dfd2(inst_ori,loc_ori,a2,3,3,3);
	}
	return;
}
//end dup
 
//end skip 
 
 
 
//dup1
void d_orid_thet_dfd0(DiffDoub0 inst_ori[], DiffDoub0 loc_ori[], DiffDoub0 rot[], int v1, int v2) {
	if(v1 + v2 == 0) {
		rotate_orient_dfd0(inst_ori, loc_ori, rot);
		return;
	}
//preserve
    DiffDoub1 d_rot[3];
    DiffDoub2 d2_rot[3];
	DiffDoub1 d_loc_ori[9];
	DiffDoub2 d2_loc_ori[9];
	DiffDoub1 d_inst_ori[9];
	DiffDoub2 d2_inst_ori[9];
//end preserve
    int i1;
	bool is_diff = loc_ori[0].diff_type();
	double param;
	double param2;
	if(!is_diff) {
		if(v1*v2 == 0) {
			for (i1 = 0; i1 < 9; i1++) {
				param = loc_ori[i1].val;
			    d_loc_ori[i1].set_val(param);
		    }
			v1 = v1 + v2;
			for (i1 = 1; i1 < 4; i1++) {
				param = rot[i1-1].val;
				if(i1 == v1) {
					d_rot[i1-1].set_val_2(param,1.0);
				} else {
					d_rot[i1-1].set_val_2(param,0.0);
				}
			}
			rotate_orient_dfd1(d_inst_ori, d_loc_ori, d_rot);
			for (i1 = 0; i1 < 9; i1++) {
				param = d_inst_ori[i1].dval;
				inst_ori[i1].set_val(param);
			}
		} else {
			for (i1 = 0; i1 < 9; i1++) {
				param = loc_ori[i1].val;
				param2 = loc_ori[i1].dval;
			    d2_loc_ori[i1].set_val_2(param,param2);
		    }
			if(v1 == v2) {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2_rot[i1-1].set_val_3(param,1.0,1.0);
					} else {
						d2_rot[i1-1].set_val(param);
					}
				}	
			} else {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2_rot[i1-1].set_val_2(param,1.0);
					} else if(i1 == v2) {
						d2_rot[i1-1].set_val_3(param,0.0,1.0);
					} else {
						d2_rot[i1-1].set_val_2(param,0.0);
					}
				}
			}
			rotate_orient_dfd2(d2_inst_ori, d2_loc_ori, d2_rot);
			for (i1 = 0; i1 < 9; i1++) {
				param = d2_inst_ori[i1].dv12;
				inst_ori[i1].set_val(param);
			}
		}
	} else {
		for (i1 = 0; i1 < 9; i1++) {
			param = loc_ori[i1].val;
			param2 = loc_ori[i1].dval;
			d2_loc_ori[i1].set_val_2(param,param2);
		}
		v1 = v1 + v2;
		for (i1 = 1; i1 < 4; i1++) {
			param = rot[i1-1].val;
			if(i1 == v1) {
				d2_rot[i1-1].set_val_3(param,0.0,1.0);
			} else {
				d2_rot[i1-1].set_val_2(param,0.0);
			}
		}
		rotate_orient_dfd2(d2_inst_ori, d2_loc_ori, d2_rot);
		for (i1 = 0; i1 < 9; i1++) {
			param = d2_inst_ori[i1].dv2;
			param2 = d2_inst_ori[i1].dv12;
			inst_ori[i1].set_val_2(param,param2);
		}
	}
	return;
}

//end dup  
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void d_orid_thet_dfd1(DiffDoub1 inst_ori[], DiffDoub1 loc_ori[], DiffDoub1 rot[], int v1, int v2) {
	if(v1 + v2 == 0) {
		rotate_orient_dfd1(inst_ori, loc_ori, rot);
		return;
	}
    DiffDoub1 d_rot[3];
    DiffDoub2 d2_rot[3];
	DiffDoub1 d_loc_ori[9];
	DiffDoub2 d2_loc_ori[9];
	DiffDoub1 d_inst_ori[9];
	DiffDoub2 d2_inst_ori[9];
    int i1;
	bool is_diff = loc_ori[0].diff_type();
	double param;
	double param2;
	if(!is_diff) {
		if(v1*v2 == 0) {
			for (i1 = 0; i1 < 9; i1++) {
				param = loc_ori[i1].val;
			    d_loc_ori[i1].set_val(param);
		    }
			v1 = v1 + v2;
			for (i1 = 1; i1 < 4; i1++) {
				param = rot[i1-1].val;
				if(i1 == v1) {
					d_rot[i1-1].set_val_2(param,1.0);
				} else {
					d_rot[i1-1].set_val_2(param,0.0);
				}
			}
			rotate_orient_dfd1(d_inst_ori, d_loc_ori, d_rot);
			for (i1 = 0; i1 < 9; i1++) {
				param = d_inst_ori[i1].dval;
				inst_ori[i1].set_val(param);
			}
		} else {
			for (i1 = 0; i1 < 9; i1++) {
				param = loc_ori[i1].val;
				param2 = loc_ori[i1].dval;
			    d2_loc_ori[i1].set_val_2(param,param2);
		    }
			if(v1 == v2) {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2_rot[i1-1].set_val_3(param,1.0,1.0);
					} else {
						d2_rot[i1-1].set_val(param);
					}
				}	
			} else {
				for (i1 = 1; i1 < 4; i1++) {
					param = rot[i1-1].val;
					if(i1 == v1) {
						d2_rot[i1-1].set_val_2(param,1.0);
					} else if(i1 == v2) {
						d2_rot[i1-1].set_val_3(param,0.0,1.0);
					} else {
						d2_rot[i1-1].set_val_2(param,0.0);
					}
				}
			}
			rotate_orient_dfd2(d2_inst_ori, d2_loc_ori, d2_rot);
			for (i1 = 0; i1 < 9; i1++) {
				param = d2_inst_ori[i1].dv12;
				inst_ori[i1].set_val(param);
			}
		}
	} else {
		for (i1 = 0; i1 < 9; i1++) {
			param = loc_ori[i1].val;
			param2 = loc_ori[i1].dval;
			d2_loc_ori[i1].set_val_2(param,param2);
		}
		v1 = v1 + v2;
		for (i1 = 1; i1 < 4; i1++) {
			param = rot[i1-1].val;
			if(i1 == v1) {
				d2_rot[i1-1].set_val_3(param,0.0,1.0);
			} else {
				d2_rot[i1-1].set_val_2(param,0.0);
			}
		}
		rotate_orient_dfd2(d2_inst_ori, d2_loc_ori, d2_rot);
		for (i1 = 0; i1 < 9; i1++) {
			param = d2_inst_ori[i1].dv2;
			param2 = d2_inst_ori[i1].dv12;
			inst_ori[i1].set_val_2(param,param2);
		}
	}
	return;
}

//end dup  
 
//end skip 
 
 
 
