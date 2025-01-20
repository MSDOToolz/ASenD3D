#include "JobClass.h"
#include <string>

using namespace std;

JobCommand::JobCommand() {
	cmd_string = "";
	file_name = "";
	
	// for solve
	dynamic = false;
	elastic = true;
	load_ramp_steps = 1;
	newmark_beta = 0.25;
	newmark_gamma = 0.5;
	ray_damp_ck = 0.0;
	ray_damp_cm = 0.0;
	nonlinear_geom = false;
	save_soln_hist = false;
	lump_mass = false;
	full_reform = 1;
	sim_period = 1.0;
	solver_bandwidth = 2000000000;
	solver_block_dim = 2000000000;
	solver_method = "direct";
	max_it = 0;
	conv_tol = 1.0e-12;
	thermal = false;
	time_step = 1.0;
	
	// modal
	type = "buckling";
	num_modes = 10;
	tgt_eval = 0.0;
	
	// set_soln_to_mode
	soln_field = "displacement";
	mode = 1;
	max_amplitude = 1.0;
	
	// write_node_results
	node_set = "all";
	time_step_tag = "";
	
	// write_element_results
	element_set = "all";
	position = "centroid";
	
	// write_modal_results
	write_modes = true;
	
	// write_objective
	write_gradient = true;	
	return;
}

