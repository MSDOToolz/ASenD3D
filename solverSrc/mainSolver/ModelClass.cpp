#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include "ModelClass.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"

using namespace std;

const int max_int = 2000000000;

Model::Model() {
	an_prep_run = false;
	time_steps_saved = 0;
	el_mat_dim = 0;
	tot_glob_dof = 0;
	elastic_scaled = false;
	therm_scaled = false;
	solve_cmd = max_int;
	modal_cmd = max_int;

	return;
}

bool Model::key_in_map(map<string, int>& in_map, string& key) {
	int i1;
	try {
		i1 = in_map.at(key);
		return true;
	}
	catch (...) {
		return false;
	}
}

bool Model::is_int(string& in_str) {
	int i1;
	try {
		i1 = stoi(in_str);
		return true;
	}
	catch (...) {
		return false;
	}
}

bool Model::is_doub(string& in_str) {
	double d1;
	try {
		d1 = stod(in_str);
		return true;
	}
	catch (...) {
		return false;
	}
}

void Model::execute_job() {
	int i1;
	int i2;
	int num_tsteps;
	double time;
	string cmd_str;
	string file_name;
	string exten;
	string full_fname;

	i1 = 0;
	for (auto& this_cmd : job) {
		cmd_str = this_cmd.cmd_string;
		file_name = this_cmd.file_name;
		if(cmd_str == "readModelInput") {
			cout << "reading Model input: " << file_name << endl;
			read_model_input(file_name);
			read_constraint_input(file_name);
			read_load_input(file_name);
			read_initial_state(file_name);
		} else if(cmd_str == "readConstraints") {
			cout << "reading constraints: " << file_name << endl;
			read_constraint_input(file_name);
		} else if(cmd_str == "readLoads") {
			cout << "reading loads: " << file_name << endl;
			read_load_input(file_name);
		} else if(cmd_str == "readInitialState") {
			cout << "reading initial state: " << file_name << endl;
			read_initial_state(file_name);
		} else if(cmd_str == "readDesignVarInput") {
			cout << "reading design variable input: " << file_name << endl;
			read_des_var_input(file_name);
		} else if(cmd_str == "readObjectiveInput") {
			cout << "reading Objective input: " << file_name << endl;
			read_objective_input(file_name);
		}
		else if (cmd_str == "readDesignVarValues") {
			read_des_var_values(file_name);
		}
		else if (cmd_str == "readNodeResults") {
			read_node_results(this_cmd.file_name);
		}
		else if (cmd_str == "solvePrep") {
			solve_cmd = i1;
			analysis_prep();
		}
		else if (cmd_str == "solve") {
			cout << "running main analysis " << endl;
			solve_cmd = i1;
			solve();
		}
		else if (cmd_str == "zeroSolution") {
			zero_solution(this_cmd.fields);
		}
		else if (cmd_str == "modalAnalysis") {
			cout << "running modal analysis " << file_name << endl;
			modal_cmd = i1;
			eigen_solve();
		}
		else if (cmd_str == "setSolnToMode") {
			set_soln_to_mode(this_cmd.soln_field, this_cmd.mode, this_cmd.max_amplitude);
		}
		else if (cmd_str == "calcObjective") {
			cout << "calculating Objective function " << file_name << endl;
			get_objective();
		}
		else if (cmd_str == "calcObjGradient") {
			cout << "calculating Objective gradient " << file_name << endl;
			get_obj_gradient();
		} else if (cmd_str == "writeNodeResults" || cmd_str == "writeElementResults") {
			if (cmd_str == "writeNodeResults") {
				cout << "writing Node results" << endl;
			}
			else {
				cout << "writing Element results" << endl;
			}
			num_tsteps = this_cmd.time_steps.size();
			if (num_tsteps > 0) {
				i2 = this_cmd.file_name.find(".");
				if (i2 > -1) {
					exten = this_cmd.file_name.substr(i2);
					file_name = this_cmd.file_name.substr(0, i2);
				}
				else {
					exten = "";
					file_name = this_cmd.file_name;
				}
				for (auto& ts : this_cmd.time_steps) {
					full_fname = file_name + "_timestep" + to_string(ts) + exten;
					if (cmd_str == "writeNodeResults") {
						write_node_results(full_fname, this_cmd.node_set, this_cmd.fields, ts);
					}
					else {
						write_element_results(full_fname, this_cmd.element_set, this_cmd.fields, this_cmd.position, ts);
					}
				}
			}
			else {
				if (cmd_str == "writeNodeResults") {
					write_node_results(this_cmd.file_name, this_cmd.node_set, this_cmd.fields, max_int);
				}
				else {
					write_element_results(this_cmd.file_name, this_cmd.element_set, this_cmd.fields, this_cmd.position, max_int);
				}
			}
		}
		else if (cmd_str == "writeModalResults") {
			cout << "writing modal results" << endl;
			write_modal_results(this_cmd.file_name, this_cmd.write_modes);
		}
		else if (cmd_str == "writeObjective") {
			cout << "writing Objective results" << endl;
			write_objective(this_cmd.file_name, this_cmd.obj_include, this_cmd.write_gradient);
		}
		
		i1++;
	}
	
	return;
}