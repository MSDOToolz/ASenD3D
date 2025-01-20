#ifndef model
#define model
#include <fstream>
#include <iostream>
#include <map>
#include "JobClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "FaceClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "LoadClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"
#include "DiffDoubClass.h"

class Model {
	public:
	    std::vector<Node> nodes;
		std::vector<Element> elements;
		std::vector<Face> faces;
		std::vector<Set> node_sets;
		std::map<std::string, int> ns_map;
		std::vector<Set> element_sets;
		std::map<std::string, int> es_map;
		std::vector<Section> sections;
		std::vector<Material> materials;
		std::vector<Fluid> fluids;
		ConstraintList elastic_const;
		ConstraintList thermal_const;
		std::vector<Load> elastic_loads;
		std::vector<Load> thermal_loads;
		std::vector<DesignVariable> design_vars;
		Objective obj;
		std::vector<JobCommand> job;
		
		int el_mat_dim;
		int tot_glob_dof;
		bool an_prep_run;
		int time_steps_saved;
		int solve_cmd;
		int modal_cmd;
		DiffDoub0StressPrereq d0_pre;
		DiffDoub1StressPrereq d1_pre;
		
		SparseMat elastic_mat;
		LowerTriMat elastic_lt;
		SparseMat non_frc_el_mat;
		bool elastic_scaled;
		SparseMat therm_mat;
		LowerTriMat therm_lt;
		bool therm_scaled;

		std::vector<double> eig_vecs;
		std::vector<double> eig_vals;
		std::vector<double> diag_mass;
		std::vector<double> load_fact;
		
		std::vector<double> temp_v1;
		std::vector<double> temp_v2;
		std::vector<double> temp_v3;
		std::vector<double> temp_v4;
		std::vector<double> temp_v5;
		std::vector<double> temp_v6;

		std::vector<DiffDoub0> temp_d1;

		std::vector<double> d_ld_u;
		std::vector<double> d_ld_v;
		std::vector<double> d_ld_a;
		std::vector<double> d_ld_t;
		std::vector<double> d_ld_tdot;
		std::vector<double> u_adj;
		std::vector<double> v_adj;
		std::vector<double> a_adj;
		std::vector<double> t_adj;
		std::vector<double> tdot_adj;

		std::vector<DiffDoub1> d_rud_d;
		std::vector<DiffDoub1> d_rtd_d;

		std::vector<int> el_in_d;

		std::vector<double> d_ld_d;
	
	    Model();

		bool key_in_map(std::map<std::string, int>& in_map, std::string& key);

		bool is_int(std::string& in_str);

		bool is_doub(std::string& in_str);
		
		void execute_job();
		
		// input
		
		void read_input_line(std::string& file_line,std::string headings[],int hd_ld_space[], std::string data[], int& data_len);
		
		void read_job(std::string file_name);
		
		void read_model_input(std::string file_name);
		
		void read_constraint_input(std::string file_name);
		
		void read_load_input(std::string file_name);
		
		void read_initial_state(std::string file_name);
		
		void read_des_var_input(std::string file_name);
		
		void read_objective_input(std::string file_name);
		
		void read_des_var_values(std::string file_name);

		void read_node_results(std::string file_name);

		void read_time_step_soln(int t_step);
		
		// analysis
		
		void reorder_nodes(int block_dim);

		void build_constraint_mats();
		
		void update_reference();
		
		void find_surface_faces();

		void prep_matrix_factorizations();
		
		void analysis_prep();
		
		void build_elastic_app_load(std::vector<double>& app_ld, double time);

		void build_thermal_app_load(std::vector<double>& app_ld, double time);
		
		void build_elastic_soln_load(std::vector<double>& soln_ld, bool build_mat, bool full_ref);

		void build_thermal_soln_load(std::vector<double>& soln_ld, bool build_mat);

		void scale_elastic_const();

		void scale_thermal_const();

		void build_elastic_const_load(std::vector<double>& const_ld);

		void build_thermal_const_load(std::vector<double>& const_ld);
		
		void solve_step(double time, double app_ld_fact, bool full_ref);
		
		void solve();

		void zero_solution(std::list<std::string>& fields);
		
		void eigen_solve();

	/*	void backup_elastic();

		void restore_elastic();*/

		void set_soln_to_mode(std::string field, int mode, double max_val);

		void augmentd_ld_u();

		void solve_for_adjoint(double time, bool full_ref);

		void d_rthermald_d(int d_var_num);

		void d_relasticd_d(int d_var_num);

		void get_objective();
		
		void get_obj_gradient();
		
		// output

		void write_time_step_soln(int t_step);
		
		void write_node_results(std::string file_name, std::string node_set, std::list<std::string>& fields, int time_step);

		void write_element_results(std::string file_name, std::string el_set, std::list<std::string>& fields, std::string position, int time_step);

		void write_modal_results(std::string file_name, bool write_modes);

		void write_objective(std::string file_name, std::list<std::string>& include_fields, bool write_grad);
		
		//
};

#endif