#ifndef jobclass
#define jobclass
#include "ListEntClass.h"
#include <vector>
#include <string>

class JobCommand {
	public:
	    std::string cmd_string;
		std::string file_name;
		
		// for solve
		bool dynamic;
		bool elastic;
		int load_ramp_steps;
		double newmark_beta;
		double newmark_gamma;
		double ray_damp_ck;
		double ray_damp_cm;
		bool nonlinear_geom;
		bool save_soln_hist;
		double sim_period;
		bool lump_mass;
		int full_reform;
		int solver_bandwidth;
		int solver_block_dim;
		std::string solver_method;
		int max_it;
		double conv_tol;
		std::list<double> static_load_time;
		bool thermal;
		double time_step;
		
		// modal
		std::string type;
		int num_modes;
		double tgt_eval;
		
		// set_soln_to_mode
		std::string soln_field;
		int mode;
		double max_amplitude;
		
		// write_node_results
		std::string node_set;
		std::list<std::string> fields;
		std::list<int> time_steps;
		std::string time_step_tag;
		
		// write_element_results
		std::string element_set;
		std::string position;
		
		// write_modal_results
		bool write_modes;
		
		// write_element_properties
		std::list<std::string> properties;
		
		// write_objective
		std::list<std::string> obj_include;
		bool write_gradient;
		
		JobCommand();
};

#endif