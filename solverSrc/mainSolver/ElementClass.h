#ifndef element
#define element
#include <vector>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"
#include "FaceClass.h"
#include "LoadClass.h"
#include "JobClass.h"

//dup1
class DiffDoub0StressPrereq {
public:
	std::vector<DiffDoub0> glob_nds;  //glob_nds[30];
	std::vector<DiffDoub0> loc_nds;// [30] ;
	std::vector<DiffDoub0> loc_ori;// [9] ;
	std::vector<DiffDoub0> inst_ori;// [720] ;
	std::vector<DiffDoub0> glob_disp;// [60] ;
	std::vector<DiffDoub0> glob_vel;
	std::vector<DiffDoub0> glob_acc;// [30];
	std::vector<DiffDoub0> glob_temp;
	std::vector<DiffDoub0> glob_tdot;
	std::vector<DiffDoub0> cmat;// [81] ;
	std::vector<DiffDoub0> mmat;// [36];
	std::vector<DiffDoub0> dmat;
	std::vector<DiffDoub0> therm_exp;
	std::vector<DiffDoub0> einit;
	std::vector<DiffDoub0> tcmat;
	DiffDoub0 spec_heat;
	std::vector<DiffDoub0> bmat;
	std::vector<DiffDoub0> cbmat;
	std::vector<DiffDoub0> layer_z;
	std::vector<DiffDoub0> layer_thk;
	std::vector<DiffDoub0> layer_ang;
	std::vector<DiffDoub0> layer_q;
	std::vector<DiffDoub0> layer_d;
	std::vector<DiffDoub0> layer_te;
	std::vector<DiffDoub0> layer_e0;
	std::vector<DiffDoub0> layer_den;
	std::vector<DiffDoub0> layer_tc;
	std::vector<DiffDoub0> layer_sh;
	std::vector<DiffDoub0> frc_fld_coef;
	std::vector<DiffDoub0> frc_fld_exp;
	std::vector<DiffDoub0> thrm_fld_coef;
	DiffDoub0 ref_temp;
	DiffDoub0 mass_per_el;
	std::vector<double> scr_mat1;
	std::vector<double> scr_mat2;
	std::vector<double> scr_mat3;
	std::vector<double> scr_mat4;
	std::vector<double> scr_mat5;
	std::vector<DiffDoub0> scr_vec1;
	std::vector<DiffDoub0> scr_vec2;
	std::vector<DiffDoub0> scr_vec3;
	std::vector<DiffDoub0> scr_vec4;
	std::vector<DiffDoub0> scr_vec5;

	int current_lay_len;

	DiffDoub0StressPrereq();

	void allocate_layers_dfd0(int num_layers);

};

class DiffDoub0FlPrereq {
public:
	std::vector<DiffDoub0> glob_nds;
	std::vector<DiffDoub0> glob_disp;
	std::vector<DiffDoub0> glob_vel;
	std::vector<DiffDoub0> fl_den;
	std::vector<DiffDoub0> fl_vel;
	std::vector<DiffDoub0> fl_temp;
	std::vector<DiffDoub0> fl_turb_e;
	std::vector<DiffDoub0> fl_den_dot;
	std::vector<DiffDoub0> fl_vel_dot;
	std::vector<DiffDoub0> fl_tdot;
	std::vector<DiffDoub0> fl_turb_edot;
	DiffDoub0 ref_turb_e;
	DiffDoub0 grad_vturb_coef;
	DiffDoub0 diss_turb_coef;
	DiffDoub0 ref_visc;
	DiffDoub0 den_vis_coef;
	DiffDoub0 temp_vis_coef;
	DiffDoub0 turb_vis_coef;
	DiffDoub0 ref_enth;
	DiffDoub0 den_enth_coef;
	DiffDoub0 den_enth_exp;
	DiffDoub0 den_pres_coef;
	DiffDoub0 den_pres_exp;
	DiffDoub0 ref_den;
	DiffDoub0 ref_temp;
	DiffDoub0 therm_cond;
	DiffDoub0 spec_heat;
	DiffDoub0 i_gconst;
	std::vector<double> scratch;

	DiffDoub0FlPrereq();

};

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
class DiffDoub1StressPrereq {
public:
	std::vector<DiffDoub1> glob_nds;  //glob_nds[30];
	std::vector<DiffDoub1> loc_nds;// [30] ;
	std::vector<DiffDoub1> loc_ori;// [9] ;
	std::vector<DiffDoub1> inst_ori;// [720] ;
	std::vector<DiffDoub1> glob_disp;// [60] ;
	std::vector<DiffDoub1> glob_vel;
	std::vector<DiffDoub1> glob_acc;// [30];
	std::vector<DiffDoub1> glob_temp;
	std::vector<DiffDoub1> glob_tdot;
	std::vector<DiffDoub1> cmat;// [81] ;
	std::vector<DiffDoub1> mmat;// [36];
	std::vector<DiffDoub1> dmat;
	std::vector<DiffDoub1> therm_exp;
	std::vector<DiffDoub1> einit;
	std::vector<DiffDoub1> tcmat;
	DiffDoub1 spec_heat;
	std::vector<DiffDoub1> bmat;
	std::vector<DiffDoub1> cbmat;
	std::vector<DiffDoub1> layer_z;
	std::vector<DiffDoub1> layer_thk;
	std::vector<DiffDoub1> layer_ang;
	std::vector<DiffDoub1> layer_q;
	std::vector<DiffDoub1> layer_d;
	std::vector<DiffDoub1> layer_te;
	std::vector<DiffDoub1> layer_e0;
	std::vector<DiffDoub1> layer_den;
	std::vector<DiffDoub1> layer_tc;
	std::vector<DiffDoub1> layer_sh;
	std::vector<DiffDoub1> frc_fld_coef;
	std::vector<DiffDoub1> frc_fld_exp;
	std::vector<DiffDoub1> thrm_fld_coef;
	DiffDoub1 ref_temp;
	DiffDoub1 mass_per_el;
	std::vector<double> scr_mat1;
	std::vector<double> scr_mat2;
	std::vector<double> scr_mat3;
	std::vector<double> scr_mat4;
	std::vector<double> scr_mat5;
	std::vector<DiffDoub1> scr_vec1;
	std::vector<DiffDoub1> scr_vec2;
	std::vector<DiffDoub1> scr_vec3;
	std::vector<DiffDoub1> scr_vec4;
	std::vector<DiffDoub1> scr_vec5;

	int current_lay_len;

	DiffDoub1StressPrereq();

	void allocate_layers_dfd1(int num_layers);

};

class DiffDoub1FlPrereq {
public:
	std::vector<DiffDoub1> glob_nds;
	std::vector<DiffDoub1> glob_disp;
	std::vector<DiffDoub1> glob_vel;
	std::vector<DiffDoub1> fl_den;
	std::vector<DiffDoub1> fl_vel;
	std::vector<DiffDoub1> fl_temp;
	std::vector<DiffDoub1> fl_turb_e;
	std::vector<DiffDoub1> fl_den_dot;
	std::vector<DiffDoub1> fl_vel_dot;
	std::vector<DiffDoub1> fl_tdot;
	std::vector<DiffDoub1> fl_turb_edot;
	DiffDoub1 ref_turb_e;
	DiffDoub1 grad_vturb_coef;
	DiffDoub1 diss_turb_coef;
	DiffDoub1 ref_visc;
	DiffDoub1 den_vis_coef;
	DiffDoub1 temp_vis_coef;
	DiffDoub1 turb_vis_coef;
	DiffDoub1 ref_enth;
	DiffDoub1 den_enth_coef;
	DiffDoub1 den_enth_exp;
	DiffDoub1 den_pres_coef;
	DiffDoub1 den_pres_exp;
	DiffDoub1 ref_den;
	DiffDoub1 ref_temp;
	DiffDoub1 therm_cond;
	DiffDoub1 spec_heat;
	DiffDoub1 i_gconst;
	std::vector<double> scratch;

	DiffDoub1FlPrereq();

};

//end dup
 
//end skip 
 
 
 
class Element {
	public:
	    int type;
		int label;
	    std::vector<int> nodes;
		int num_nds;
		int dof_per_nd;
		int num_int_dof;
		int n_dim;
		int def_dim;
		std::vector<int> dof_table;  //dof_table[2*i] = nd/basis function of dof i, dof_table[2*i+1] = component of dof i
		int int_dof_index;
		std::vector<double> int_pts;
		std::vector<double> ip_wt;
		std::vector<double> nd_spts;
		double s_cent[3];
		int num_ip;
		int num_faces;
		std::list<int> faces;
		std::vector<double> internal_disp;
		std::vector<double> int_prev_disp;
		std::vector<double> internald_ldu;
		std::vector<double> internal_adj;
		std::vector<DiffDoub1> internal_ru;
		std::vector<double> internal_mat;
		std::list<IDCapsule> design_vars;
		std::list<int> comp_dvars;
		int sect_ptr;
		
		Element();

	    void initialize_type(int new_type);

        void set_nodes(int new_nds[]);

        void initialize_faces(std::vector<Face>& glob_fc_lst, int& fi);	

		void set_int_disp(double new_disp[]);

		void set_int_prev_disp(double new_disp[]);

		void advance_int_disp();

		void backstep_int_disp();

		void set_intd_ld_u(std::vector<double>& globd_ld_u);

		int get_num_layers(std::vector<Section>& sec_lst);
		
		void add_design_variable(int d_index, double coef);

		void add_comp_dvar(int d_index);

//dup1

// properties
		void get_gen_prop_dfd0(DiffDoub0& prop, std::string prop_key, std::vector<DesignVariable>& dv_ar);

		void get_layer_thk_z_dfd0(std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_z, DiffDoub0& z_offset, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_q_dfd0(std::vector<DiffDoub0>& lay_q, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_d_dfd0(std::vector<DiffDoub0>& lay_d, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_angle_dfd0(std::vector<DiffDoub0>& lay_ang, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_th_exp_dfd0(std::vector<DiffDoub0>& lay_th_exp, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_einit_dfd0(std::vector<DiffDoub0>& lay_einit, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_den_dfd0(std::vector<DiffDoub0>& layer_den, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_cond_dfd0(std::vector<DiffDoub0>& lay_cond, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_spec_heat_dfd0(std::vector<DiffDoub0>& lay_sh, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void transform_strain_dfd0(DiffDoub0 stn_new[], DiffDoub0 stn_orig[], DiffDoub0& angle);

		void transform_q_dfd0(DiffDoub0 q_new[], DiffDoub0 q_orig[], DiffDoub0& angle);

		void get_solid_stiff_dfd0(std::vector<DiffDoub0>& cmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_abd_dfd0(std::vector<DiffDoub0>& cmat, std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_z, std::vector<DiffDoub0>& lay_q, std::vector<DiffDoub0>& lay_ang, std::vector<Section>& sec_ar);

		void get_beam_stiff_dfd0(std::vector<DiffDoub0>& cmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_thermal_exp_dfd0(std::vector<DiffDoub0>& th_exp, std::vector<DiffDoub0>& einit, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_exp_load_dfd0(std::vector<DiffDoub0>& exp_ld, std::vector<DiffDoub0>& e0_ld, std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_z, std::vector<DiffDoub0>& lay_q, std::vector<DiffDoub0>& lay_th_exp, std::vector<DiffDoub0>& lay_einit, std::vector<DiffDoub0>& lay_ang, std::vector<Section>& sec_ar);

		void get_beam_exp_load_dfd0(std::vector<DiffDoub0>& exp_ld, std::vector<DiffDoub0>& e0_ld, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_density_dfd0(DiffDoub0& den, int layer, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_mass_dfd0(std::vector<DiffDoub0>& mmat, std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_z, std::vector<DiffDoub0>& lay_den, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_beam_mass_dfd0(std::vector<DiffDoub0>& mmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_solid_damp_dfd0(std::vector<DiffDoub0>& dmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_damp_dfd0(std::vector<DiffDoub0>& dmat, std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_z, std::vector<DiffDoub0>& lay_d, std::vector<DiffDoub0>& lay_ang, std::vector<Section>& sec_ar);

		void get_beam_damp_dfd0(std::vector<DiffDoub0>& dmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_conductivity_dfd0(std::vector<DiffDoub0>& t_cond, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_cond_dfd0(std::vector<DiffDoub0>& t_cond, std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_ang, std::vector<DiffDoub0>& lay_cond, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_beam_cond_dfd0(std::vector<DiffDoub0>& t_cond, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_specific_heat_dfd0(DiffDoub0& spec_heat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_spec_heat_dfd0(DiffDoub0& spec_heat, std::vector<DiffDoub0>& lay_thk, std::vector<DiffDoub0>& lay_sh, std::vector<DiffDoub0>& lay_den, std::vector<Section>& sec_ar);

		void get_beam_spec_heat_dfd0(DiffDoub0& spec_heat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

        void get_nd_crds_dfd0(std::vector<DiffDoub0>& x_glob, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_loc_ori_dfd0(std::vector<DiffDoub0>& loc_ori, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);
		
		void correct_orient_dfd0(std::vector<DiffDoub0>& loc_ori, std::vector<DiffDoub0>& x_glob);

		void get_frc_fld_const_dfd0(std::vector<DiffDoub0>& coef, std::vector<DiffDoub0>& exp, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_thrm_fld_const_dfd0(std::vector<DiffDoub0>& coef, DiffDoub0& ref_t, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_mass_per_el_dfd0(DiffDoub0& mass_per_el, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

// solution fields
		void get_nd_disp_dfd0(std::vector<DiffDoub0>& glob_disp, std::vector<Node>& nd_ar);

		void get_nd_vel_dfd0(std::vector<DiffDoub0>& glob_vel, std::vector<Node>& nd_ar);

		void get_nd_acc_dfd0(std::vector<DiffDoub0>& glob_acc, std::vector<Node>& nd_ar);

		void get_nd_fl_vel_dfd0(std::vector<DiffDoub0>& fl_vel, std::vector<Node>& nd_ar);

		void get_nd_fl_vdot_dfd0(std::vector<DiffDoub0>& fl_vdot, std::vector<Node>& nd_ar);

		void get_nd_temp_dfd0(std::vector<DiffDoub0>& glob_temp, std::vector<Node>& nd_ar);

		void get_nd_tdot_dfd0(std::vector<DiffDoub0>& glob_tdot, std::vector<Node>& nd_ar);

		void get_nd_fl_den_dfd0(std::vector<DiffDoub0>& fl_den, std::vector<Node>& nd_ar);

		void get_nd_fl_den_dot_dfd0(std::vector<DiffDoub0>& fl_den_dot, std::vector<Node>& nd_ar);

		void get_nd_turb_e_dfd0(std::vector<DiffDoub0>& turb_e, std::vector<Node>& nd_ar);

		void get_nd_turb_edot_dfd0(std::vector<DiffDoub0>& turb_edot, std::vector<Node>& nd_ar);

		void eval_n_dfd0(DiffDoub0 n_vec[], DiffDoub0 d_nds[], double spt[]);
		
		void get_ip_data_dfd0(DiffDoub0 n_vec[], DiffDoub0 d_ndx[], DiffDoub0& det_j, std::vector<DiffDoub0>& loc_nds, double spt[]);
		
		void get_inst_ori_dfd0(std::vector<DiffDoub0>& inst_ori_mat, std::vector<DiffDoub0>& loc_ori, std::vector<DiffDoub0>& glob_disp, int stat);
		
		void get_inst_disp_dfd0(DiffDoub0 inst_disp[], std::vector<DiffDoub0>& glob_disp, std::vector<DiffDoub0>& inst_ori_mat, std::vector<DiffDoub0>& loc_ori, std::vector<DiffDoub0>& x_glob, bool n_lgeom, int dv1, int dv2);

		void get_stress_prereq_dfd0(DiffDoub0StressPrereq& pre, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_fluid_prereq_dfd0(DiffDoub0FlPrereq& pre, std::vector<Section>& sec_ar, std::vector<Fluid>& fl_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_volume_dfd0(DiffDoub0& vol, DiffDoub0StressPrereq& pre, int layer, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_section_def_dfd0(DiffDoub0 sec_def[], std::vector<DiffDoub0>& glob_disp, std::vector<DiffDoub0>& inst_ori_mat, std::vector<DiffDoub0>& loc_ori, std::vector<DiffDoub0>& x_glob, DiffDoub0 d_ndx[], DiffDoub0 n_vec[], bool n_lgeom, int dv1, int dv2);
		
		void get_solid_strain_dfd0(DiffDoub0 strain[], DiffDoub0 ux[], DiffDoub0 d_ndx[], std::vector<DiffDoub0>& loc_ori, int dv1, int dv2, bool n_lgeom);

		void get_stress_strain_dfd0(DiffDoub0 stress[], DiffDoub0 strain[], double spt[], int layer, bool n_lgeom, DiffDoub0StressPrereq& pre);

		void d_stress_straind_u_dfd0(std::vector<DiffDoub0>& dsd_u, std::vector<DiffDoub0>& ded_u, std::vector<DiffDoub0>& dsd_t, double spt[], int layer, bool n_lgeom, DiffDoub0StressPrereq& pre);

		void get_def_frc_mom_dfd0(DiffDoub0 def[], DiffDoub0 frc_mom[], double spt[], bool n_lgeom, DiffDoub0StressPrereq& pre);

		void d_def_frc_momd_u_dfd0(std::vector<DiffDoub0>& d_defd_u, std::vector<DiffDoub0>& d_frc_momd_u, std::vector<DiffDoub0>& d_frc_momd_t, double spt[], bool n_lgeom, DiffDoub0StressPrereq& pre);

		void get_flux_tgrad_dfd0(DiffDoub0 flux[], DiffDoub0 t_grad[], double spt[], int layer, DiffDoub0StressPrereq& pre);

		void d_flux_tgradd_t_dfd0(std::vector<DiffDoub0>& d_fd_t, std::vector<DiffDoub0>& d_tg, double spt[], int layer, DiffDoub0StressPrereq& pre);

		void put_vec_to_glob_mat_dfd0(SparseMat& q_mat, std::vector<DiffDoub0>& el_qvec, bool for_therm, int mat_row, std::vector<Node>& nd_ar);
		
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

// properties
		void get_gen_prop_dfd1(DiffDoub1& prop, std::string prop_key, std::vector<DesignVariable>& dv_ar);

		void get_layer_thk_z_dfd1(std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_z, DiffDoub1& z_offset, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_q_dfd1(std::vector<DiffDoub1>& lay_q, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_d_dfd1(std::vector<DiffDoub1>& lay_d, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_angle_dfd1(std::vector<DiffDoub1>& lay_ang, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_th_exp_dfd1(std::vector<DiffDoub1>& lay_th_exp, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_einit_dfd1(std::vector<DiffDoub1>& lay_einit, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_den_dfd1(std::vector<DiffDoub1>& layer_den, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_cond_dfd1(std::vector<DiffDoub1>& lay_cond, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_layer_spec_heat_dfd1(std::vector<DiffDoub1>& lay_sh, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void transform_strain_dfd1(DiffDoub1 stn_new[], DiffDoub1 stn_orig[], DiffDoub1& angle);

		void transform_q_dfd1(DiffDoub1 q_new[], DiffDoub1 q_orig[], DiffDoub1& angle);

		void get_solid_stiff_dfd1(std::vector<DiffDoub1>& cmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_abd_dfd1(std::vector<DiffDoub1>& cmat, std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_z, std::vector<DiffDoub1>& lay_q, std::vector<DiffDoub1>& lay_ang, std::vector<Section>& sec_ar);

		void get_beam_stiff_dfd1(std::vector<DiffDoub1>& cmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_thermal_exp_dfd1(std::vector<DiffDoub1>& th_exp, std::vector<DiffDoub1>& einit, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_exp_load_dfd1(std::vector<DiffDoub1>& exp_ld, std::vector<DiffDoub1>& e0_ld, std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_z, std::vector<DiffDoub1>& lay_q, std::vector<DiffDoub1>& lay_th_exp, std::vector<DiffDoub1>& lay_einit, std::vector<DiffDoub1>& lay_ang, std::vector<Section>& sec_ar);

		void get_beam_exp_load_dfd1(std::vector<DiffDoub1>& exp_ld, std::vector<DiffDoub1>& e0_ld, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_density_dfd1(DiffDoub1& den, int layer, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_mass_dfd1(std::vector<DiffDoub1>& mmat, std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_z, std::vector<DiffDoub1>& lay_den, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_beam_mass_dfd1(std::vector<DiffDoub1>& mmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_solid_damp_dfd1(std::vector<DiffDoub1>& dmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_damp_dfd1(std::vector<DiffDoub1>& dmat, std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_z, std::vector<DiffDoub1>& lay_d, std::vector<DiffDoub1>& lay_ang, std::vector<Section>& sec_ar);

		void get_beam_damp_dfd1(std::vector<DiffDoub1>& dmat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_conductivity_dfd1(std::vector<DiffDoub1>& t_cond, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_cond_dfd1(std::vector<DiffDoub1>& t_cond, std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_ang, std::vector<DiffDoub1>& lay_cond, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_beam_cond_dfd1(std::vector<DiffDoub1>& t_cond, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_specific_heat_dfd1(DiffDoub1& spec_heat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

		void get_shell_spec_heat_dfd1(DiffDoub1& spec_heat, std::vector<DiffDoub1>& lay_thk, std::vector<DiffDoub1>& lay_sh, std::vector<DiffDoub1>& lay_den, std::vector<Section>& sec_ar);

		void get_beam_spec_heat_dfd1(DiffDoub1& spec_heat, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar);

        void get_nd_crds_dfd1(std::vector<DiffDoub1>& x_glob, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_loc_ori_dfd1(std::vector<DiffDoub1>& loc_ori, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);
		
		void correct_orient_dfd1(std::vector<DiffDoub1>& loc_ori, std::vector<DiffDoub1>& x_glob);

		void get_frc_fld_const_dfd1(std::vector<DiffDoub1>& coef, std::vector<DiffDoub1>& exp, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_thrm_fld_const_dfd1(std::vector<DiffDoub1>& coef, DiffDoub1& ref_t, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

		void get_mass_per_el_dfd1(DiffDoub1& mass_per_el, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);

// solution fields
		void get_nd_disp_dfd1(std::vector<DiffDoub1>& glob_disp, std::vector<Node>& nd_ar);

		void get_nd_vel_dfd1(std::vector<DiffDoub1>& glob_vel, std::vector<Node>& nd_ar);

		void get_nd_acc_dfd1(std::vector<DiffDoub1>& glob_acc, std::vector<Node>& nd_ar);

		void get_nd_fl_vel_dfd1(std::vector<DiffDoub1>& fl_vel, std::vector<Node>& nd_ar);

		void get_nd_fl_vdot_dfd1(std::vector<DiffDoub1>& fl_vdot, std::vector<Node>& nd_ar);

		void get_nd_temp_dfd1(std::vector<DiffDoub1>& glob_temp, std::vector<Node>& nd_ar);

		void get_nd_tdot_dfd1(std::vector<DiffDoub1>& glob_tdot, std::vector<Node>& nd_ar);

		void get_nd_fl_den_dfd1(std::vector<DiffDoub1>& fl_den, std::vector<Node>& nd_ar);

		void get_nd_fl_den_dot_dfd1(std::vector<DiffDoub1>& fl_den_dot, std::vector<Node>& nd_ar);

		void get_nd_turb_e_dfd1(std::vector<DiffDoub1>& turb_e, std::vector<Node>& nd_ar);

		void get_nd_turb_edot_dfd1(std::vector<DiffDoub1>& turb_edot, std::vector<Node>& nd_ar);

		void eval_n_dfd1(DiffDoub1 n_vec[], DiffDoub1 d_nds[], double spt[]);
		
		void get_ip_data_dfd1(DiffDoub1 n_vec[], DiffDoub1 d_ndx[], DiffDoub1& det_j, std::vector<DiffDoub1>& loc_nds, double spt[]);
		
		void get_inst_ori_dfd1(std::vector<DiffDoub1>& inst_ori_mat, std::vector<DiffDoub1>& loc_ori, std::vector<DiffDoub1>& glob_disp, int stat);
		
		void get_inst_disp_dfd1(DiffDoub1 inst_disp[], std::vector<DiffDoub1>& glob_disp, std::vector<DiffDoub1>& inst_ori_mat, std::vector<DiffDoub1>& loc_ori, std::vector<DiffDoub1>& x_glob, bool n_lgeom, int dv1, int dv2);

		void get_stress_prereq_dfd1(DiffDoub1StressPrereq& pre, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_fluid_prereq_dfd1(DiffDoub1FlPrereq& pre, std::vector<Section>& sec_ar, std::vector<Fluid>& fl_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_volume_dfd1(DiffDoub1& vol, DiffDoub1StressPrereq& pre, int layer, std::vector<Section>& sec_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_section_def_dfd1(DiffDoub1 sec_def[], std::vector<DiffDoub1>& glob_disp, std::vector<DiffDoub1>& inst_ori_mat, std::vector<DiffDoub1>& loc_ori, std::vector<DiffDoub1>& x_glob, DiffDoub1 d_ndx[], DiffDoub1 n_vec[], bool n_lgeom, int dv1, int dv2);
		
		void get_solid_strain_dfd1(DiffDoub1 strain[], DiffDoub1 ux[], DiffDoub1 d_ndx[], std::vector<DiffDoub1>& loc_ori, int dv1, int dv2, bool n_lgeom);

		void get_stress_strain_dfd1(DiffDoub1 stress[], DiffDoub1 strain[], double spt[], int layer, bool n_lgeom, DiffDoub1StressPrereq& pre);

		void d_stress_straind_u_dfd1(std::vector<DiffDoub1>& dsd_u, std::vector<DiffDoub1>& ded_u, std::vector<DiffDoub1>& dsd_t, double spt[], int layer, bool n_lgeom, DiffDoub1StressPrereq& pre);

		void get_def_frc_mom_dfd1(DiffDoub1 def[], DiffDoub1 frc_mom[], double spt[], bool n_lgeom, DiffDoub1StressPrereq& pre);

		void d_def_frc_momd_u_dfd1(std::vector<DiffDoub1>& d_defd_u, std::vector<DiffDoub1>& d_frc_momd_u, std::vector<DiffDoub1>& d_frc_momd_t, double spt[], bool n_lgeom, DiffDoub1StressPrereq& pre);

		void get_flux_tgrad_dfd1(DiffDoub1 flux[], DiffDoub1 t_grad[], double spt[], int layer, DiffDoub1StressPrereq& pre);

		void d_flux_tgradd_t_dfd1(std::vector<DiffDoub1>& d_fd_t, std::vector<DiffDoub1>& d_tg, double spt[], int layer, DiffDoub1StressPrereq& pre);

		void put_vec_to_glob_mat_dfd1(SparseMat& q_mat, std::vector<DiffDoub1>& el_qvec, bool for_therm, int mat_row, std::vector<Node>& nd_ar);
		
//end dup
 
//end skip 
 
 
 
		void get_el_vec(std::vector<double>& el_vec, std::vector<double>& glob_vec, bool for_therm, bool intnl, std::vector<Node>& nd_ar);

		void add_to_glob_vec(std::vector<double>& el_vec, std::vector<double>& glob_vec, bool for_therm, bool intnl, std::vector<Node>& nd_ar);
 
 
// equations
		void condense_mat(std::vector<double>& mat, std::vector<double>& scr1, std::vector<double>& scr2);
		
		void update_external(std::vector<double>& ext_vec, int for_soln, std::vector<Node>& nd_ar, std::vector<double>& scr1, std::vector<double>& scr2);
		
		void update_internal(std::vector<double>& ext_vec, int for_soln, std::vector<Node>& nd_ar, std::vector<double>& scr1, std::vector<double>& scr2);

		double get_int_adjd_rd_d();

//dup1
        void get_ruk_dfd0(std::vector<DiffDoub0>& rvec, std::vector<double>& d_rdu, std::vector<double>& d_rd_t, bool get_matrix, bool n_lgeom, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_rum_dfd0(std::vector<DiffDoub0>& rvec, std::vector<double>& d_rd_a, bool get_matrix, bool actual_props, bool n_lgeom, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_rud_dfd0(std::vector<DiffDoub0>& rvec, std::vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_ru_dfd0(std::vector<DiffDoub0>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_rtk_dfd0(std::vector<DiffDoub0>& rvec, std::vector<double>& d_rd_t, bool get_matrix, DiffDoub0StressPrereq& pre);

		void get_rtm_dfd0(std::vector<DiffDoub0>& rvec, std::vector<double>& d_rd_tdot, bool get_matrix, bool actual_props, DiffDoub0StressPrereq& pre);

		void get_rt_dfd0(std::vector<DiffDoub0>& glob_r, SparseMat& globd_rd_t, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar);

		void get_ru_frc_fld_dfd0(std::vector<DiffDoub0>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar);

		void get_rt_frc_fld_dfd0(std::vector<DiffDoub0>& glob_r, SparseMat& globd_rd_t, std::vector<double>& d_rd_u, std::vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& nd_ar);

		void get_app_load_dfd0(std::vector<DiffDoub0>& app_ld, Load& ld_pt, bool n_lgeom, DiffDoub0StressPrereq& pre, std::vector<Section>& sec_ar, std::vector<Face>& fc_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_app_therm_load_dfd0(std::vector<DiffDoub0>& app_ld, Load& ld_pt, DiffDoub0StressPrereq& pre, std::vector<Section>& sec_ar, std::vector<Face>& fc_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_rf_unst_dfd0(std::vector<DiffDoub0>& rvec, DiffDoub0FlPrereq& pre);

		void get_rf_conv_dfd0(std::vector<DiffDoub0>& rvec, DiffDoub0FlPrereq& pre);

		void get_rf_visc_dfd0(std::vector<DiffDoub0>& rvec, DiffDoub0FlPrereq& pre);

		void get_rf_pres_dfd0(std::vector<DiffDoub0>& rvec, DiffDoub0FlPrereq& pre);

		void get_rf_load_dfd0(std::vector<DiffDoub0>& rvec, DiffDoub0FlPrereq& pre, std::vector<DiffDoub0>& body_f, DiffDoub0& body_q);

		void get_rf_dfd0(std::vector<DiffDoub0>& glob_r, SparseMat& globd_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub0FlPrereq& pre);

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
        void get_ruk_dfd1(std::vector<DiffDoub1>& rvec, std::vector<double>& d_rdu, std::vector<double>& d_rd_t, bool get_matrix, bool n_lgeom, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_rum_dfd1(std::vector<DiffDoub1>& rvec, std::vector<double>& d_rd_a, bool get_matrix, bool actual_props, bool n_lgeom, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_rud_dfd1(std::vector<DiffDoub1>& rvec, std::vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_ru_dfd1(std::vector<DiffDoub1>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);
		
		void get_rtk_dfd1(std::vector<DiffDoub1>& rvec, std::vector<double>& d_rd_t, bool get_matrix, DiffDoub1StressPrereq& pre);

		void get_rtm_dfd1(std::vector<DiffDoub1>& rvec, std::vector<double>& d_rd_tdot, bool get_matrix, bool actual_props, DiffDoub1StressPrereq& pre);

		void get_rt_dfd1(std::vector<DiffDoub1>& glob_r, SparseMat& globd_rd_t, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar);

		void get_ru_frc_fld_dfd1(std::vector<DiffDoub1>& glob_r, SparseMat& globd_rdu, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar);

		void get_rt_frc_fld_dfd1(std::vector<DiffDoub1>& glob_r, SparseMat& globd_rd_t, std::vector<double>& d_rd_u, std::vector<double>& d_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& nd_ar);

		void get_app_load_dfd1(std::vector<DiffDoub1>& app_ld, Load& ld_pt, bool n_lgeom, DiffDoub1StressPrereq& pre, std::vector<Section>& sec_ar, std::vector<Face>& fc_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_app_therm_load_dfd1(std::vector<DiffDoub1>& app_ld, Load& ld_pt, DiffDoub1StressPrereq& pre, std::vector<Section>& sec_ar, std::vector<Face>& fc_ar, std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		void get_rf_unst_dfd1(std::vector<DiffDoub1>& rvec, DiffDoub1FlPrereq& pre);

		void get_rf_conv_dfd1(std::vector<DiffDoub1>& rvec, DiffDoub1FlPrereq& pre);

		void get_rf_visc_dfd1(std::vector<DiffDoub1>& rvec, DiffDoub1FlPrereq& pre);

		void get_rf_pres_dfd1(std::vector<DiffDoub1>& rvec, DiffDoub1FlPrereq& pre);

		void get_rf_load_dfd1(std::vector<DiffDoub1>& rvec, DiffDoub1FlPrereq& pre, std::vector<DiffDoub1>& body_f, DiffDoub1& body_q);

		void get_rf_dfd1(std::vector<DiffDoub1>& glob_r, SparseMat& globd_rd_v, bool get_matrix, JobCommand& cmd, DiffDoub1FlPrereq& pre);

//end dup
 
//end skip 
 
 
 
};


#endif