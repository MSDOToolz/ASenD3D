#ifndef matrixfun
#define matrixfun
#include <vector>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "LUMatClass.h"
#include "LowerTriMatClass.h"
#include "ConstraintClass.h"

void sub_vec(std::vector<double>& sub_v, std::vector<double>& v_in, int st, int end);

void return_sv(std::vector<double>& sub_v, std::vector<double>& v_in, int st, int end);

void vec_to_ar(double ar[], std::vector<double>& vc, int st, int end);

double get_dist(double p1[], double p2[]);

void cross_prod(double prod[], double v1[], double v2[]);

void q_rfactor(std::vector<double>& mat, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void solveq_rx_eqb(std::vector<double>& x_vec, std::vector<double>& mat, std::vector<double>& b_vec, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void conj_grad_sparse(std::vector<double>& soln, SparseMat& mat, ConstraintList& cnst, LowerTriMat& pc_mat, std::vector<double>& rhs, double conv_tol, int max_it);

void g_mres_sparse(std::vector<double>& soln, SparseMat& mat, ConstraintList& cnst, LUMat& pc_mat, std::vector<double>& rhs, double conv_tol, int max_it, int restart);

void sym_factor(std::vector<double>& mat, std::vector<double>& q_mat, int mat_dim);

void get_char_fun(DiffDoub1& c_fun, std::vector<DiffDoub1>& mat, int mat_dim, std::vector<double>& e_vals, double lam, int tri_diag);

void get_evals(std::vector<double>& e_vals, std::vector<double>& mat, int mat_dim, double lam_init, double conv_tol, int tri_diag);

double get_low_eval(std::vector<double>& mat, int mat_dim);

void eigen_solve(std::vector<double>& e_vals, std::vector<double>& e_vecs, std::vector<double>& mat, int mat_dim, int tri_diag);

void sym_eigen_solve(std::vector<double>& e_vals, std::vector<double>& e_vecs, std::vector<double>& mat, int mat_dim, int tri_diag);

void eigen_sparse_direct(std::vector<double>& e_vals, std::vector<double>& e_vecs, int num_pairs, LowerTriMat& mat, std::vector<double>& mass_mat, int mat_dim);

double ray_quot(std::vector<double>& grad, std::vector<double>& kv, std::vector<double>& mv, SparseMat& mat, ConstraintList& cnst, std::vector<double>& mass_mat, std::vector<double>& in_vec);

double unitize_vec(std::vector<double>& vec, int dim);

void get_nearest_evec_rq(SparseMat& mat, ConstraintList& cnst, std::vector<double>& mass_mat, std::vector<double>& in_vecs, std::vector<double>& e_vals, int num_vecs, int max_it);


//dup1
void q_rfactor_dfd0(std::vector<DiffDoub0>& mat, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void q_rfactor_ar_dfd0(DiffDoub0 mat[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void solveq_rx_eqb_dfd0(std::vector<DiffDoub0>& x_vec, std::vector<DiffDoub0>& mat, std::vector<DiffDoub0>& b_vec, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void solveq_rx_eqb_ar_dfd0(DiffDoub0 x_vec[], DiffDoub0 mat[], DiffDoub0 b_vec[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void get_det_inv_dfd0(DiffDoub0& det, std::vector<DiffDoub0>& inv, std::vector<DiffDoub0>& mat, int col_dim, int tri_diag, std::vector<DiffDoub0>& x_vec, std::vector<DiffDoub0>& b_vec);

void get_det_inv_ar_dfd0(DiffDoub0& det, DiffDoub0 inv[], DiffDoub0 mat[], int col_dim, int tri_diag, DiffDoub0 x_vec[], DiffDoub0 b_vec[]);
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void q_rfactor_dfd1(std::vector<DiffDoub1>& mat, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void q_rfactor_ar_dfd1(DiffDoub1 mat[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void solveq_rx_eqb_dfd1(std::vector<DiffDoub1>& x_vec, std::vector<DiffDoub1>& mat, std::vector<DiffDoub1>& b_vec, int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void solveq_rx_eqb_ar_dfd1(DiffDoub1 x_vec[], DiffDoub1 mat[], DiffDoub1 b_vec[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void get_det_inv_dfd1(DiffDoub1& det, std::vector<DiffDoub1>& inv, std::vector<DiffDoub1>& mat, int col_dim, int tri_diag, std::vector<DiffDoub1>& x_vec, std::vector<DiffDoub1>& b_vec);

void get_det_inv_ar_dfd1(DiffDoub1& det, DiffDoub1 inv[], DiffDoub1 mat[], int col_dim, int tri_diag, DiffDoub1 x_vec[], DiffDoub1 b_vec[]);
//end dup
 
//end skip 
 
 
 
//dup2
void sub_vec_dfd0(std::vector<DiffDoub0>& sub_v, std::vector<DiffDoub0>& in_v, int st, int end);

void return_sv_dfd0(std::vector<DiffDoub0>& sub_v, std::vector<DiffDoub0>& in_v, int st, int end);

void vec_to_ar_dfd0(DiffDoub0 ar[], std::vector<DiffDoub0>& vc, int st, int end);

void ar_to_vec_dfd0(DiffDoub0 ar[], std::vector<DiffDoub0>& vc, int st, int end);

void mat_mul_dfd0(std::vector<DiffDoub0>& prod, std::vector<DiffDoub0>& mat1, std::vector<DiffDoub0>& mat2, int m1_rows, int m1_cols, int m2_cols);

void mat_mul_ar_dfd0(DiffDoub0 prod[], DiffDoub0 mat1[], DiffDoub0 mat2[], int m1_rows, int m1_cols, int m2_cols);

void transpose_dfd0(std::vector<DiffDoub0>& mat_t, std::vector<DiffDoub0>& mat, int row_dim, int col_dim);

void transpose_ar_dfd0(DiffDoub0 mat_t[], DiffDoub0 mat[], int row_dim, int col_dim);

void cross_prod_dfd0(DiffDoub0 prod[], DiffDoub0 v1[], DiffDoub0 v2[]);

void rotate_orient_dfd0(DiffDoub0 inst_ori[], DiffDoub0 loc_ori[], DiffDoub0 rot[]);
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup2
void sub_vec_dfd1(std::vector<DiffDoub1>& sub_v, std::vector<DiffDoub1>& in_v, int st, int end);

void return_sv_dfd1(std::vector<DiffDoub1>& sub_v, std::vector<DiffDoub1>& in_v, int st, int end);

void vec_to_ar_dfd1(DiffDoub1 ar[], std::vector<DiffDoub1>& vc, int st, int end);

void ar_to_vec_dfd1(DiffDoub1 ar[], std::vector<DiffDoub1>& vc, int st, int end);

void mat_mul_dfd1(std::vector<DiffDoub1>& prod, std::vector<DiffDoub1>& mat1, std::vector<DiffDoub1>& mat2, int m1_rows, int m1_cols, int m2_cols);

void mat_mul_ar_dfd1(DiffDoub1 prod[], DiffDoub1 mat1[], DiffDoub1 mat2[], int m1_rows, int m1_cols, int m2_cols);

void transpose_dfd1(std::vector<DiffDoub1>& mat_t, std::vector<DiffDoub1>& mat, int row_dim, int col_dim);

void transpose_ar_dfd1(DiffDoub1 mat_t[], DiffDoub1 mat[], int row_dim, int col_dim);

void cross_prod_dfd1(DiffDoub1 prod[], DiffDoub1 v1[], DiffDoub1 v2[]);

void rotate_orient_dfd1(DiffDoub1 inst_ori[], DiffDoub1 loc_ori[], DiffDoub1 rot[]);
//end dup
 
//DiffDoub2 versions: 
//dup2
void sub_vec_dfd2(std::vector<DiffDoub2>& sub_v, std::vector<DiffDoub2>& in_v, int st, int end);

void return_sv_dfd2(std::vector<DiffDoub2>& sub_v, std::vector<DiffDoub2>& in_v, int st, int end);

void vec_to_ar_dfd2(DiffDoub2 ar[], std::vector<DiffDoub2>& vc, int st, int end);

void ar_to_vec_dfd2(DiffDoub2 ar[], std::vector<DiffDoub2>& vc, int st, int end);

void mat_mul_dfd2(std::vector<DiffDoub2>& prod, std::vector<DiffDoub2>& mat1, std::vector<DiffDoub2>& mat2, int m1_rows, int m1_cols, int m2_cols);

void mat_mul_ar_dfd2(DiffDoub2 prod[], DiffDoub2 mat1[], DiffDoub2 mat2[], int m1_rows, int m1_cols, int m2_cols);

void transpose_dfd2(std::vector<DiffDoub2>& mat_t, std::vector<DiffDoub2>& mat, int row_dim, int col_dim);

void transpose_ar_dfd2(DiffDoub2 mat_t[], DiffDoub2 mat[], int row_dim, int col_dim);

void cross_prod_dfd2(DiffDoub2 prod[], DiffDoub2 v1[], DiffDoub2 v2[]);

void rotate_orient_dfd2(DiffDoub2 inst_ori[], DiffDoub2 loc_ori[], DiffDoub2 rot[]);
//end dup
 
//end skip 
 
 
 
//dup1

void d_orid_thet_dfd0(DiffDoub0 inst_ori[], DiffDoub0 loc_ori[], DiffDoub0 rot[], int v1, int v2);

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

void d_orid_thet_dfd1(DiffDoub1 inst_ori[], DiffDoub1 loc_ori[], DiffDoub1 rot[], int v1, int v2);

//end dup
 
//end skip 
 
 
 
#endif