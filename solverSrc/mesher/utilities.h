#ifndef utilities
#define utilities
#include <string>
#include <fstream>

void cross_prod(double prod[], double v1[], double v2[]);

void q_rfactor(double mat[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void solveq_rx_eqb(double x_vec[], double mat[], double b_vec[], int col_dim, int st_row, int end_row, int st_col, int end_col, int tri_diag);

void read_input_line(std::string& file_line, std::string headings[], int hd_ld_space[], std::string data[], int& data_len);

#endif
