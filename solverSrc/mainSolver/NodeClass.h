#ifndef nodeclass
#define nodeclass
#include <string>
#include <iostream>
#include <vector>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "DesignVariableClass.h"

class Node {
	public:
	    int label;
		bool fluid;
		int num_dof;
		int dof_index[6];
		int sorted_rank;
	    double coord[3];
		double displacement[6];
		double velocity[6];
		double fl_vel[3];
		double fl_vel_dot[3];
		double acceleration[6];
		double temperature;
		double temp_change_rate;
		double fl_den;
		double fl_den_dot;
		double turb_e;
		double turb_edot;
		double prev_disp[6];
		double prev_vel[6];
		double prev_fl_vel[3];
		double p_fl_vel_lf[3];
		double prev_fl_vel_dot[3];
		double prev_acc[6];
		double prev_temp;
		double p_temp_lf;
		double prev_tdot;
		double prev_fl_den;
		double p_fl_den_lf;
		double prev_fl_den_dot;
		double prev_turb_e;
		double prev_turb_edot;
		double p_turb_elf;
		double initial_disp[6];
		double initial_vel[6];
		double initial_fl_vel[3];
		double initial_fl_vel_dot[3];
		double initial_acc[6];
		double initial_temp;
		double initial_tdot;
		double initial_fl_den;
		double initial_fl_den_dot;
		double initial_turb_e;
		double initial_turb_edot;
		std::list<IDCapsule> d_var_lst;
		
	    Node();
		
		void set_crd(double new_crd[]);
		
		void set_displacement(double new_disp[]);

		void set_velocity(double new_vel[]);

		void set_acceleration(double new_acc[]);

		void set_fl_vel(double new_vel[]);

		void set_fl_vel_dot(double new_fl_vdot[]);
		
		void add_to_displacement(double del_disp[]);

		void add_to_fl_vel(double del_fl_vel[]);
		
		void set_initial_disp(double new_disp[]);
		
		void set_initial_vel(double new_vel[]);
		
		void set_initial_acc(double new_acc[]);

		void set_initial_fl_vel(double new_fl_vel[]);

		void set_initial_fl_vdot(double new_fl_vdot[]);

		void set_prev_disp(double new_disp[]);

		void set_prev_vel(double new_vel[]);

		void set_prev_acc(double new_acc[]);

		void set_prev_fl_vel(double new_fl_vel[]);

		void set_pfl_vel_lf(double new_vel[]);

		void set_prev_fl_vdot(double new_fl_vdot[]);
		
		void initialize_disp();

		void initialize_fl_vel();

		void initialize_temp();

		void initialize_fl_den();

		void initialize_turb_e();
		
		void update_vel_acc(double nm_beta, double nm_gamma, double del_t);

		void update_fl_vel_dot(double nm_gamma, double del_t);
		
		void update_tdot(double nm_gamma, double del_t);

		void update_fl_den_dot(double nm_gamma, double del_t);

		void update_turb_edot(double nm_gamma, double del_t);
		
		void advance_disp();

		void advance_fl_vel();

		void advance_temp();

		void advance_fl_den();

		void advance_turb_e();

		void backstep_disp();

		void backstep_fl_vel();

		void backstep_temp();

		void backstep_fl_den();

		void backstep_turb_e();
		
		void add_design_variable(int d_index, double coef);
		
//dup1
		void get_crd_dfd0(DiffDoub0 crd_out[], std::vector<DesignVariable>& dv_ar);
		
		void get_disp_dfd0(DiffDoub0 disp[]);
		
		void get_elastic_dvload_dfd0(DiffDoub0 ld[], std::vector<DesignVariable>& dv_ar);

		void get_thermal_dvload_dfd0(DiffDoub0& ld, std::vector<DesignVariable>& dv_ar);
//end dup	
 
//skip 
 
//DiffDoub1 versions: 
//dup1
		void get_crd_dfd1(DiffDoub1 crd_out[], std::vector<DesignVariable>& dv_ar);
		
		void get_disp_dfd1(DiffDoub1 disp[]);
		
		void get_elastic_dvload_dfd1(DiffDoub1 ld[], std::vector<DesignVariable>& dv_ar);

		void get_thermal_dvload_dfd1(DiffDoub1& ld, std::vector<DesignVariable>& dv_ar);
//end dup	
 
//end skip 
 
 
 
};

#endif