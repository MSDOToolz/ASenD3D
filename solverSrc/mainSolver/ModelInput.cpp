#include <fstream>
#include <iostream>
#include <string>
#include "ModelClass.h"
#include "constants.h"
#include "JobClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"

using namespace std;

// input

void Model::read_input_line(string& file_line, string headings[], int hd_ld_space[], string data[], int& data_len) {
	int i1;
	int i2;
	int ln_len;
	int wrd_len;
	//getline(in_file,file_line);
	i1 = file_line.find("#");
	if(i1 > -1) {
		file_line = file_line.substr(0,i1);
	}
	file_line = file_line + " ";
	ln_len = file_line.length();
	i1 = file_line.find(":");
	data_len = 0;
	if(i1 > -1) {
		i2 = file_line.find_first_not_of(" -\n\t");
		wrd_len = i1 - i2;
		if(headings[0] == "" || hd_ld_space[0] == i2) {
			headings[0] = file_line.substr(i2,wrd_len);
			hd_ld_space[0] = i2;
			headings[1] = "";
			hd_ld_space[1] = 0;
			headings[2] = "";
			hd_ld_space[2] = 0;
			headings[3] = "";
			hd_ld_space[3] = 0;
		} else if(headings[1] == "" || hd_ld_space[1] == i2) {
			headings[1] = file_line.substr(i2,wrd_len);
			hd_ld_space[1] = i2;
			headings[2] = "";
			hd_ld_space[2] = 0;
			headings[3] = "";
			hd_ld_space[3] = 0;
		} else if(headings[2] == "" || hd_ld_space[2] == i2) {
			headings[2] = file_line.substr(i2,wrd_len);
			hd_ld_space[2] = i2;
			headings[3] = "";
			hd_ld_space[3] = 0;
		} else {
			headings[3] = file_line.substr(i2,wrd_len);
			hd_ld_space[3] = i2;
		}
		i1++;
		while(i1 < ln_len) {
			file_line = file_line.substr(i1);
			i1 = file_line.find_first_not_of(" ,[]\t\n");
			if(i1 > -1) {
				file_line = file_line.substr(i1);
				ln_len = file_line.length();
				i1 = file_line.find_first_of(" ,[]\t\n");
				if(i1 > -1) {
				    data[data_len] = file_line.substr(0,i1);
				    data_len++;
				} else {
					i1 = ln_len;
				}
			} else {
				i1 = ln_len;
			}
		}
	} else {
		i1 = file_line.find("- ");
		if(i1 > -1) {
			i1++;
			while(i1 < ln_len) {
				file_line = file_line.substr(i1);
				i1 = file_line.find_first_not_of(" ,[]\t\n");
				if(i1 > -1) {
					file_line = file_line.substr(i1);
					ln_len = file_line.length();
					i1 = file_line.find_first_of(" ,[]\t\n");
					if(i1 > -1) {
						data[data_len] = file_line.substr(0,i1);
						data_len++;
					} else {
						i1 = ln_len;
					}
				} else {
					i1 = ln_len;
				}
			}
		}
	}
	
	return;
}

void Model::read_job(string file_name) {
	int i1;
	int i2;
	int i3;
	int i4;
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	string prev_ld_hd = "";
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	string err_st;
	int cmd_ct;
	
	cmd_ct = 0;
	in_file.open(file_name);
	while (!in_file.eof()) {
		getline(in_file, file_line);
		read_input_line(file_line, headings, hd_ld_space, data, data_len);
		if (headings[1] == "command" && data_len == 1) {
			cmd_ct++;
		}
	}

	job = vector<JobCommand>(cmd_ct);
	//in_file.seekg(0, std::ios::beg);
	in_file.close();
	in_file.open(file_name);
	
	if(in_file) {
		cmd_ct = max_int;
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[1] == "command" && data_len == 1) {
				if (cmd_ct == max_int) {
					cmd_ct = 0;
				}
				else {
					cmd_ct++;
				}
				job[cmd_ct].cmd_string = data[0];
			} else if(headings[1] == "fileName" && data_len == 1) {
				job[cmd_ct].file_name = data[0];
			} else if(headings[1] == "dynamic" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].dynamic = true;
				} else {
					job[cmd_ct].dynamic = false;
				}
			}
			else if (headings[1] == "elastic" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].elastic = true;
				} else {
					job[cmd_ct].elastic = false;
				}
			} else if(headings[1] == "loadRampSteps" && data_len == 1) {
				job[cmd_ct].load_ramp_steps = stoi(data[0]);
			} else if(headings[1] == "newmarkBeta" && data_len == 1) {
				job[cmd_ct].newmark_beta = stod(data[0]);
			} else if(headings[1] == "newmarkGamma" && data_len == 1) {
				job[cmd_ct].newmark_gamma = stod(data[0]);
			} else if(headings[1] == "nonlinearGeom" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].nonlinear_geom = true;
				} else {
					job[cmd_ct].nonlinear_geom = false;
				}
			} else if(headings[1] == "saveSolnHist" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].save_soln_hist = true;
				} else {
					job[cmd_ct].save_soln_hist = false;
				}
			}
			else if (headings[1] == "solnHistDir" && data_len == 1) {
				job[cmd_ct].file_name = data[0];
			}
			else if (headings[1] == "lumpMass" && data_len == 1) {
				if (data[0] == "yes") {
					job[cmd_ct].lump_mass = true;
				}
				else {
					job[cmd_ct].lump_mass = false;
				}
			}
			else if (headings[1] == "fullReformFreq") {
				job[cmd_ct].full_reform = stoi(data[0]);
			}
			else if (headings[1] == "simPeriod" && data_len == 1) {
				job[cmd_ct].sim_period = stod(data[0]);
			} else if(headings[1] == "solverBandwidth" && data_len == 1) {
				job[cmd_ct].solver_bandwidth = stoi(data[0]);
			} else if(headings[1] == "solverBlockDim" && data_len == 1) {
				job[cmd_ct].solver_block_dim = stoi(data[0]);
			} else if(headings[1] == "solverMethod" && data_len == 1) {
				job[cmd_ct].solver_method = data[0];
			}
			else if (headings[1] == "maxIterations" && data_len == 1) {
				job[cmd_ct].max_it = stoi(data[0]);
			}
			else if (headings[1] == "convergenceTol" && data_len == 1) {
				job[cmd_ct].conv_tol = stod(data[0]);
			}
			else if (headings[1] == "staticLoadTime" && data_len == 1) {
				job[cmd_ct].static_load_time.push_back(stod(data[0]));
				//new_cmd->static_load_time = stod(data[0]);
			} else if(headings[1] == "thermal" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].thermal = true;
				} else {
					job[cmd_ct].thermal = false;
				}
			} else if(headings[1] == "timeStep" && data_len == 1) {
				job[cmd_ct].time_step = stod(data[0]);
			} else if(headings[1] == "type" && data_len == 1) {
				job[cmd_ct].type = data[0];
			} else if(headings[1] == "numModes" && data_len == 1) {
				job[cmd_ct].num_modes = stoi(data[0]);
			}
			else if (headings[1] == "targetEigenvalue" && data_len == 1) {
				job[cmd_ct].tgt_eval = stod(data[0]);
			}
			else if (headings[1] == "solnField" && data_len == 1) {
				job[cmd_ct].soln_field = data[0];
			} else if(headings[1] == "mode" && data_len == 1) {
				job[cmd_ct].mode = stoi(data[0]);
			} else if(headings[1] == "maxAmplitude" && data_len == 1) {
				job[cmd_ct].max_amplitude = stod(data[0]);
			} else if(headings[1] == "nodeSet" && data_len == 1) {
				job[cmd_ct].node_set = data[0];
			} else if(headings[1] == "fields" && data_len == 1) {
				job[cmd_ct].fields.push_back(data[0]);
			} else if(headings[1] == "timeSteps") {
				if(data_len == 1) {
					if(data[0] == "all") {
						job[cmd_ct].time_step_tag = "all";
					} else if(data[0] == "last") {
						job[cmd_ct].time_step_tag = "last";
					} else {
						if (is_int(data[0])) {
							i1 = stoi(data[0]);
							job[cmd_ct].time_steps.push_back(i1);
						} 
						else {
							err_st = "Warning: possible invalid entry for " + headings[0] + " timeSteps in job file " + file_name;
							cout << err_st << endl;
						}
					}
				} else if(data_len > 1) {
					i2 = stoi(data[0]);
					i3 = stoi(data[1]);
					if(data_len == 3) {
						i4 = stoi(data[2]);
					} else {
						i4 = 1;
					}
					for (i1 = i2; i1 < i3; i1+= i4) {
						job[cmd_ct].time_steps.push_back(i1);
					}
				}
			} else if(headings[1] == "elementSet" && data_len == 1) {
				job[cmd_ct].element_set = data[0];
			}
			else if (headings[1] == "position" && data_len == 1) {
				job[cmd_ct].position = data[0];
			}
			else if (headings[1] == "writeModes" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].write_modes = true;
				} else {
					job[cmd_ct].write_modes = false;
				}
			} else if(headings[1] == "properties" && data_len == 1) {
				job[cmd_ct].properties.push_back(data[0]);
			} else if(headings[1] == "include" && data_len == 1) {
				job[cmd_ct].obj_include.push_back(data[0]);
			} else if(headings[1] == "writeGradient" && data_len == 1) {
				if(data[0] == "yes") {
					job[cmd_ct].write_gradient = true;
				} else {
					job[cmd_ct].write_gradient = false;
				}
			}
			
			prev_ld_hd = headings[0];
		}
	} else {
		err_st = "Error: could not open job file: " + file_name;
		throw invalid_argument(err_st);
	}

	in_file.close();
	
	return;
}
		
void Model::read_model_input(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	
	int i1;
	int i2;
	int i3;
	int i4;
	double doub_inp[10] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	int int_inp[10] = {0,0,0,0,0,0,0,0,0,0};
	int el_type;
	
	int nd_ct = 0;
	int el_ct = 0;
	int ns_ct = 0;
	int es_ct = 0;
	int sec_ct = 0;
	int mat_ct = 0;
	int fl_ct = 0;
	int lay_ct = 0;
	int set_ind = 0;

	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "nodes" && data_len == 4) {
				nd_ct++;
			}
			else if (headings[0] == "elements") {
				if (headings[1] == "connectivity" && data_len > 1) {
					el_ct++;
				}
			}
			else if (headings[0] == "sets") {
				if (headings[1] == "node") {
					if (headings[2] == "name" && data_len == 1) {
						ns_ct++;
					}
				}
				else if (headings[1] == "element") {
					if (headings[2] == "name" && data_len == 1) {
						es_ct++;
					}
				}
			}
			else if (headings[0] == "sections") {
				if (headings[1] == "type" && data_len == 1) {
					sec_ct++;
				}
			}
			else if (headings[0] == "materials") {
				if (headings[1] == "name" && data_len == 1) {
					mat_ct++;
				}
			}
			else if (headings[0] == "fluids") {
				if (headings[1] == "name" && data_len == 1) {
					fl_ct++;
				}
			}
		}
	}

	nodes = vector<Node>(nd_ct);
	elements = vector<Element>(el_ct);
	node_sets = vector<Set>(ns_ct + nd_ct + 1);
	element_sets = vector<Set>(es_ct + el_ct + 1);
	sections = vector<Section>(sec_ct);
	materials = vector<Material>(mat_ct);
	fluids = vector<Fluid>(fl_ct);

	//in_file.seekg(0, std::ios::beg);
	in_file.close();
	in_file.open(file_name);
	
	if(in_file) {
		nd_ct = max_int;
		el_ct = max_int;
		ns_ct = max_int;
		es_ct = max_int;
		sec_ct = max_int;
		mat_ct = max_int;
		fl_ct = max_int;
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[0] == "nodes" && data_len == 4) {
				nd_ct = stoi(data[0]);
				nodes[nd_ct].label = nd_ct;
				doub_inp[0] = stod(data[1]);
				doub_inp[1] = stod(data[2]);
				doub_inp[2] = stod(data[3]);
				nodes[nd_ct].set_crd(doub_inp);
			} else if(headings[0] == "elements") {
				if(headings[1] == "type" && data_len == 1) {
					if(data[0] == "tet4") {
						el_type = 4;
					} else if(data[0] == "wedge6") {
						el_type = 6;
					} else if(data[0] == "brick8") {
						el_type = 8;
					} else if(data[0] == "brickIM") {
						el_type = 81;
					}
					else if (data[0] == "tet10") {
						el_type = 10;
					}
					else if (data[0] == "shell3") {
						el_type = 3;
					} else if(data[0] == "shell4") {
						el_type = 41;
					} else if(data[0] == "beam2") {
						el_type = 2;
					}
					else if (data[0] == "frcFld") {
						el_type = 21;
					}
					else if (data[0] == "mass") {
						el_type = 1;
					}
					else if (data[0] == "fl4") {
						el_type = 400;
					}
					else if (data[0] == "fl6") {
						el_type = 600;
					}
					else if (data[0] == "fl8") {
						el_type = 800;
					}
					else if (data[0] == "fl10") {
						el_type = 1000;
					}
					else {
						string err_st = "Error: unrecognized Element type, " + data[0];
						throw invalid_argument(err_st);
					}
				} else if(headings[1] == "connectivity" && data_len > 1) {
					el_ct = stoi(data[0]);
					Element& new_el = elements[el_ct];
					new_el.initialize_type(el_type);
					new_el.label = el_ct;
					for (i1 = 1; i1 < data_len; i1++) {
						int_inp[i1-1] = stoi(data[i1]);
					}
					new_el.set_nodes(int_inp);
				} 
			} else if(headings[0] == "sets") {
				if(headings[1] == "node") {
					if(headings[2] == "name" && data_len == 1) {
						if (ns_ct == max_int) {
							ns_ct = 0;
						}
						else {
							ns_ct++;
						}
						Set& new_set = node_sets[ns_ct];
						new_set.name = data[0];
					} else if(headings[2] == "labels") {
						if(data_len == 1) {
							node_sets[ns_ct].labels.push_back(stoi(data[0]));
						} else if(data_len > 1) {
							i2 = stoi(data[0]);
							i3 = stoi(data[1]);
							if(data_len == 3) {
								i4 = stoi(data[2]);
							} else {
								i4 = 1;
							}
							for (i1 = i2; i1 < i3; i1+= i4) {
								node_sets[ns_ct].labels.push_back(i1);
							}
						}
					}
				} else if(headings[1] == "element") {
					if(headings[2] == "name" && data_len == 1) {
						if (es_ct == max_int) {
							es_ct = 0;
						}
						else {
							es_ct++;
						}
						Set& new_set = element_sets[es_ct];
						new_set.name = data[0];
					} else if(headings[2] == "labels") {
						if(data_len == 1) {
							element_sets[es_ct].labels.push_back(stoi(data[0]));
						} else if(data_len > 1) {
							i2 = stoi(data[0]);
							i3 = stoi(data[1]);
							if(data_len == 3) {
								i4 = stoi(data[2]);
							} else {
								i4 = 1;
							}
							for (i1 = i2; i1 < i3; i1+= i4) {
								element_sets[es_ct].labels.push_back(i1);
							}
						}
					}
				}
			} else if(headings[0] == "sections") {
				if(headings[1] == "type" && data_len == 1) {
					if (sec_ct == max_int) {
						sec_ct = 0;
					}
					else {
						sec_ct++;
					}
					sections[sec_ct].type = data[0];
				} else if(headings[1] == "material" && data_len == 1) {
					sections[sec_ct].mat_name = data[0];
				}
				else if (headings[1] == "fluid" && data_len == 1) {
					sections[sec_ct].fl_name = data[0];
				}
				else if (headings[1] == "orientation" && data_len == 6) {
					for (i1 = 0; i1 < 6; i1++) {
						doub_inp[i1] = stod(data[i1]);
					}
					sections[sec_ct].set_orientation(doub_inp);
				} else if(headings[1] == "layup") {
					if(headings[2] == "zOffset" && data_len == 1) {
						sections[sec_ct].z_offset = stod(data[0]);
					} else if(headings[2] == "layers") {
						if(headings[3] == "material" && data_len == 1) {
							Layer new_lay;
							new_lay.mat_name = data[0];
							sections[sec_ct].layers.push_back(new_lay);
						} else if(headings[3] == "thickness" && data_len == 1) {
							auto lay_pt = sections[sec_ct].layers.end();
							--lay_pt;
							lay_pt->thickness = stod(data[0]);
						} else if(headings[3] == "angle" && data_len == 1) {
							auto lay_pt = sections[sec_ct].layers.end();
							--lay_pt;
							lay_pt->angle = stod(data[0]);
						}
					}
				} else if(headings[1] == "beamProps") {
					if(headings[2] == "area" && data_len == 1) {
						sections[sec_ct].area = stod(data[0]);
					} else if(headings[2] == "I" && data_len == 5) {
						for(i1 = 0; i1 < 5; i1++) {
							doub_inp[i1] = stod(data[i1]);
						}
						sections[sec_ct].set_area_moment(doub_inp);
					} else if(headings[2] == "J" && data_len == 1) {
						sections[sec_ct].polar_moment = stod(data[1]);
					} else if(headings[2] == "stiffness" && data_len == 3) {
						int_inp[0] = stoi(data[0]) - 1;
						int_inp[1] = stoi(data[1]) - 1;
						doub_inp[0] = stod(data[2]);
						sections[sec_ct].set_stiffness(int_inp[0],int_inp[1],doub_inp[0]);
					} else if(headings[2] == "mass" && data_len == 3) {
						int_inp[0] = stoi(data[0]) - 1;
						int_inp[1] = stoi(data[1]) - 1;
						doub_inp[0] = stod(data[2]);
						sections[sec_ct].set_mass(int_inp[0],int_inp[1],doub_inp[0]);
					}
					else if (headings[2] == "damping" && data_len == 3) {
						int_inp[0] = stoi(data[0]) - 1;
						int_inp[1] = stoi(data[1]) - 1;
						doub_inp[0] = stod(data[2]);
						sections[sec_ct].set_damping(int_inp[0], int_inp[1], doub_inp[0]);
					}
					else if (headings[2] == "expLoadCoef" && data_len == 6) {
						for(i1 = 0; i1 < 6; i1++) {
							doub_inp[i1] = stod(data[i1]);
						}
						sections[sec_ct].set_exp_ld(doub_inp);
					} else if(headings[2] == "conductivity" && data_len == 1) {
						sections[sec_ct].conductivity = stod(data[0]);
					} else if(headings[2] == "specHeat" && data_len == 1) {
						sections[sec_ct].spec_heat = stod(data[0]);
					}
				}
				else if (headings[1] == "potField") {
					if (headings[2] == "coef" && data_len == 1) {
						sections[sec_ct].pot_coef = stod(data[0]);
					}
					else if (headings[2] == "exp" && data_len == 1) {
						sections[sec_ct].pot_exp = stod(data[0]);
					}
				}
				else if (headings[1] == "dampField") {
					if (headings[2] == "coef" && data_len == 1) {
						sections[sec_ct].damp_coef = stod(data[0]);
					}
					else if (headings[2] == "exp" && data_len == 1) {
						sections[sec_ct].damp_exp = stod(data[0]);
					}
				}
				else if (headings[1] == "thermField") {
					if (headings[2] == "condCoef" && data_len == 1) {
						sections[sec_ct].cond_coef = stod(data[0]);
					}
					else if (headings[2] == "radCoef" && data_len == 1) {
						sections[sec_ct].rad_coef = stod(data[0]);
					}
					else if (headings[2] == "refTemp" && data_len == 1) {
						sections[sec_ct].ref_temp = stod(data[0]);
					}
				}
				else if (headings[1] == "massPerEl" && data_len == 1) {
					sections[sec_ct].mass_per_el = stod(data[0]);
				}
				else if (headings[1] == "specHeat" && data_len == 1) {
					sections[sec_ct].spec_heat = stod(data[0]);
				}
				else if (headings[1] == "fluidParams") {
					if (headings[2] == "denVisCoef" && data_len == 1) {
						sections[sec_ct].den_vis_coef = stod(data[0]);
					}
					else if (headings[2] == "tempVisCoef" && data_len == 1) {
						sections[sec_ct].temp_vis_coef = stod(data[0]);
					}
					else if (headings[2] == "turbVisCoef" && data_len == 1) {
						sections[sec_ct].turb_vis_coef = stod(data[0]);
					}
					else if (headings[2] == "gradVTurbCoef" && data_len == 1) {
						sections[sec_ct].grad_vturb_coef = stod(data[0]);
					}
					else if (headings[2] == "dissTurbCoef" && data_len == 1) {
						sections[sec_ct].diss_turb_coef = stod(data[0]);
					}
					else if (headings[2] == "enthCoef" && data_len == 1) {
						sections[sec_ct].enth_coef = stod(data[0]);
					}
					else if (headings[2] == "enthExp" && data_len == 1) {
						sections[sec_ct].enth_exp = stod(data[0]);
					}
					else if (headings[2] == "presCoef" && data_len == 1) {
						sections[sec_ct].pres_coef = stod(data[0]);
					}
					else if (headings[2] == "presExp" && data_len == 1) {
						sections[sec_ct].pres_exp = stod(data[0]);
					}
					else if (headings[2] == "refTemp" && data_len == 1) {
						sections[sec_ct].ref_temp = stod(data[0]);
					}
					else if (headings[2] == "refDen" && data_len == 1) {
						sections[sec_ct].ref_den = stod(data[0]);
					}
					else if (headings[2] == "refTurbE" && data_len == 1) {
						sections[sec_ct].ref_turb_e = stod(data[0]);
					}
					else if (headings[2] == "refEnth" && data_len == 1) {
						sections[sec_ct].ref_enth = stod(data[0]);
					}
				}
				else if(headings[1] == "elementSet" && data_len == 1) {
					sections[sec_ct].el_set_name = data[0];
				}
			} else if(headings[0] == "materials") {
				if(headings[1] == "name" && data_len == 1) {
					if (mat_ct == max_int) {
						mat_ct = 0;
					}
					else {
						mat_ct++;
					}
					materials[mat_ct].name = data[0];
				} else if(headings[1] == "density" && data_len == 1) {
					materials[mat_ct].density = stod(data[0]);
				} else if(headings[1] == "elastic") {
					if(headings[2] == "E" && data_len > 0) {
						Material& this_mat = materials[mat_ct];
						this_mat.modulus[0] = stod(data[0]);
						if(data_len == 3) {
							this_mat.modulus[1] = stod(data[1]);
							this_mat.modulus[2] = stod(data[2]);
						} else {
							this_mat.modulus[1] = this_mat.modulus[0];
							this_mat.modulus[2] = this_mat.modulus[0];
						}
					} else if(headings[2] == "nu" && data_len > 0) {
						Material& this_mat = materials[mat_ct];
						this_mat.poisson_ratio[0] = stod(data[0]);
						if(data_len == 3) {
							this_mat.poisson_ratio[1] = stod(data[1]);
							this_mat.poisson_ratio[2] = stod(data[2]);
						} else {
							this_mat.poisson_ratio[1] = this_mat.poisson_ratio[0];
							this_mat.poisson_ratio[2] = this_mat.poisson_ratio[0];
						}
					} else if(headings[2] == "G" && data_len > 0) {
						Material& this_mat = materials[mat_ct];
						this_mat.shear_mod[0] = stod(data[0]);
						if(data_len == 3) {
							this_mat.shear_mod[1] = stod(data[1]);
							this_mat.shear_mod[2] = stod(data[2]);
						} else {
							this_mat.shear_mod[1] = this_mat.shear_mod[0];
							this_mat.shear_mod[2] = this_mat.shear_mod[0];
						}
					} else if(headings[2] == "stiffness" && data_len == 3) {
						int_inp[0] = stoi(data[0]) - 1;
						int_inp[1] = stoi(data[1]) - 1;
						doub_inp[0] = stod(data[2]);
						materials[mat_ct].set_stiffness(int_inp[0],int_inp[1],doub_inp[0]);
					} 
				} else if(headings[1] == "thermal") {
					if(headings[2] == "conductivity" && data_len == 6) {
						Material& this_mat = materials[mat_ct];
						for(i1 = 0; i1 < 6; i1++) {
							this_mat.conductivity[i1] = stod(data[i1]);
						}
					} else if(headings[2] == "expansion" && data_len == 6) {
						Material& this_mat = materials[mat_ct];
						for(i1 = 0; i1 < 6; i1++) {
							this_mat.expansion[i1] = stod(data[i1]);
						}
					} else if (headings[2] == "specHeat" && data_len == 1) {
						materials[mat_ct].spec_heat = stod(data[0]);
					}
				}
				else if (headings[1] == "damping" && data_len == 3) {
					int_inp[0] = stoi(data[0]) - 1;
					int_inp[1] = stoi(data[1]) - 1;
					doub_inp[0] = stod(data[2]);
					materials[mat_ct].set_damping(int_inp[0], int_inp[1], doub_inp[0]);
				}
				else if (headings[1] == "failure") {
					if(headings[2] == "maxStress") {
						if(headings[3] == "tensile" && data_len > 0) {
							Material& n_mat = materials[mat_ct];
							n_mat.max_ten_stress[0] = stod(data[0]);
							if(data_len == 3) {
								n_mat.max_ten_stress[1] = stod(data[1]);
								n_mat.max_ten_stress[2] = stod(data[2]);
							} else {
								n_mat.max_ten_stress[1] = n_mat.max_ten_stress[0];
								n_mat.max_ten_stress[2] = n_mat.max_ten_stress[0];
							}
						} else if(headings[3] == "compressive" && data_len > 0) {
							Material& n_mat = materials[mat_ct];
							n_mat.max_comp_stress[0] = stod(data[0]);
							if(data_len == 3) {
								n_mat.max_comp_stress[1] = stod(data[1]);
								n_mat.max_comp_stress[2] = stod(data[2]);
							} else {
								n_mat.max_comp_stress[1] = n_mat.max_comp_stress[0];
								n_mat.max_comp_stress[2] = n_mat.max_comp_stress[0];
							}
						} else if(headings[3] == "shear" && data_len > 0) {
							Material& n_mat = materials[mat_ct];
							n_mat.max_shear_stress[0] = stod(data[0]);
							if(data_len == 3) {
								n_mat.max_shear_stress[1] = stod(data[1]);
								n_mat.max_shear_stress[2] = stod(data[2]);
							} else {
								n_mat.max_shear_stress[1] = n_mat.max_shear_stress[0];
								n_mat.max_shear_stress[2] = n_mat.max_shear_stress[0];
							}
						}
					} else if(headings[2] == "maxStrain") {
						if(headings[3] == "tensile" && data_len > 0) {
							Material& n_mat = materials[mat_ct];
							n_mat.max_ten_strain[0] = stod(data[0]);
							if(data_len == 3) {
								n_mat.max_ten_strain[1] = stod(data[1]);
								n_mat.max_ten_strain[2] = stod(data[2]);
							} else {
								n_mat.max_ten_strain[1] = n_mat.max_ten_strain[0];
								n_mat.max_ten_strain[2] = n_mat.max_ten_strain[0];
							}
						} else if(headings[3] == "compressive" && data_len > 0) {
							Material& n_mat = materials[mat_ct];
							n_mat.max_comp_strain[0] = stod(data[0]);
							if(data_len == 3) {
								n_mat.max_comp_strain[1] = stod(data[1]);
								n_mat.max_comp_strain[2] = stod(data[2]);
							} else {
								n_mat.max_comp_strain[1] = n_mat.max_comp_strain[0];
								n_mat.max_comp_strain[2] = n_mat.max_comp_strain[0];
							}
						} else if(headings[3] == "shear" && data_len > 0) {
							Material& n_mat = materials[mat_ct];
							n_mat.max_shear_strain[0] = stod(data[0]);
							if(data_len == 3) {
								n_mat.max_shear_strain[1] = stod(data[1]);
								n_mat.max_shear_strain[2] = stod(data[2]);
							} else {
								n_mat.max_shear_strain[1] = n_mat.max_shear_strain[0];
								n_mat.max_shear_strain[2] = n_mat.max_shear_strain[0];
							}
						}
					} else if(headings[2] == "maxStrainEnergy" && data_len == 1) {
						materials[mat_ct].max_str_eng = stod(data[0]);
					} else if(headings[2] == "maxMises" && data_len == 1) {
						materials[mat_ct].max_mises = stod(data[0]);
					}
				}
			}
			else if (headings[0] == "fluids") {
				if (headings[1] == "name" && data_len == 1) {
					if (fl_ct == max_int) {
						fl_ct = 0;
					}
					else {
						fl_ct++;
					}
					fluids[fl_ct].name = data[0];
				}
				else if (headings[1] == "viscosity" && data_len == 1) {
					fluids[fl_ct].viscosity = stod(data[0]);
				}
				else if (headings[1] == "thermal") {
					if (headings[2] == "conductivity" && data_len == 1) {
						fluids[fl_ct].therm_cond = stod(data[0]);
					}
					else if (headings[2] == "specHeat" && data_len == 1) {
						fluids[fl_ct].spec_heat = stod(data[0]);
					}
				}
				else if (headings[1] == "idealGasConst" && data_len == 1) {
					fluids[fl_ct].ideal_gas = stod(data[0]);
				}
			}
		}
	} else {
		string err_st = "Error: could not open Model input file: " + file_name;
		throw invalid_argument(err_st);
	}

	// create all and individual nd/element sets

	i1 = nodes.size();
	for (i2 = 0; i2 < i1; i2++) {
		ns_ct++;
		node_sets[ns_ct].name = to_string(i2);
		node_sets[ns_ct].labels.push_back(i2);
	}

	ns_ct++;
	node_sets[ns_ct].name = "all";
	list<int>& labs = node_sets[ns_ct].labels;
	for (i2 = 0; i2 < i1; i2++) {
		labs.push_back(i2);
	}

	i1 = elements.size();
	for (i2 = 0; i2 < i1; i2++) {
		es_ct++;
		element_sets[es_ct].name = to_string(i2);
		element_sets[es_ct].labels.push_back(i2);
	}

	es_ct++;
	Set& new_set = element_sets[es_ct];
	new_set.name = "all";
	for (i2 = 0; i2 < i1; i2++) {
		new_set.labels.push_back(i2);
	}

	// populate Set arrays

	i1 = node_sets.size();

	i2 = 0;
	for (auto& ns : node_sets) {
		ns_map.insert(make_pair(ns.name, i2));
		i2++;
	}

	i1 = element_sets.size();

	i2 = 0;
	for (auto& es : element_sets) {
		es_map.insert(make_pair(es.name, i2));
		i2++;
	}
	
	in_file.close();
	
	return;
}

void Model::read_constraint_input(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	
	int i1;
	string all_types = "displacement temperature";
	string err_st;
	
	int ec_ct = 0;
	int tc_ct = 0;
	string this_type;
	int curr_tm;

	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "constraints") {
				if (headings[1] == "type" && data_len == 1) {
					if (data[0] == "displacement") {
						ec_ct++;
					}
					else if (data[0] == "temperature") {
						tc_ct++;
					}
				}
			}
		}
	}

	elastic_const.const_vec = vector<Constraint>(ec_ct);
	thermal_const.const_vec = vector<Constraint>(tc_ct);

	//in_file.seekg(0, std::ios::beg);
	in_file.close();
	in_file.open(file_name);

	Constraint* new_con = nullptr;
	if(in_file) {
		ec_ct = max_int;
		tc_ct = max_int;
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[0] == "constraints") {
				if (headings[1] == "type" && data_len == 1) {
					if (data[0] == "displacement") {
						if (ec_ct == max_int) {
							ec_ct = 0;
						}
						else {
							ec_ct++;
						}
						new_con = &elastic_const.const_vec[ec_ct];
						new_con->type = "displacement";
					}
					else if (data[0] == "temperature") {
						if (tc_ct == max_int) {
							tc_ct = 0;
						}
						else {
							tc_ct++;
						}
						new_con = &thermal_const.const_vec[tc_ct];
						new_con->type = "temperature";
					}
					else {
						err_st = "Error: " + data[0] + " is not a valid Constraint type.  Allowable values are " + all_types;
						throw invalid_argument(err_st);
					}
				}
				else if(headings[1] == "terms") {
					if(headings[2] == "nodeSet" && data_len == 1) {
						ConstraintTerm new_cn;
						new_cn.node_set = data[0];
						new_con->terms.push_back(new_cn);
					} else if(headings[2] == "dof" && data_len == 1) {
						auto tm_pt = new_con->terms.end();
						--tm_pt;
						tm_pt->dof = stoi(data[0]);
					} else if(headings[2] == "coef" && data_len == 1) {
						auto tm_pt = new_con->terms.end();
						--tm_pt;
						tm_pt->coef = stod(data[0]);
					}
				} else if(headings[1] == "rhs" && data_len == 1) {
					new_con->rhs = stod(data[0]);
				}
			}
		}
	} else {
		err_st = "Error: could not open Constraint input file: " + file_name;
		throw invalid_argument(err_st);
	}

	in_file.close();

	return;
}

void Model::read_load_input(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	
	int i1;
	double doub_inp[10] = { 0,0,0,0,0,0,0,0,0,0 };
	string elastic_list = "nodalForce bodyForce gravitational centrifugal surfacePressure surfaceTraction";
	string thermal_list = "bodyHeadGen surfaceFlux";
	Load *new_ld = nullptr;

	int e_ld_ct = 0;
	int t_ld_ct = 0;
	
	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "loads") {
				if (headings[1] == "type" && data_len == 1) {
					i1 = elastic_list.find(data[0]);
					if (i1 > -1) {
						e_ld_ct++;
					}
					i1 = thermal_list.find(data[0]);
					if (i1 > -1) {
						t_ld_ct++;
					}
				}
			}
		}
	}

	elastic_loads = vector<Load>(e_ld_ct);
	thermal_loads = vector<Load>(t_ld_ct);

	//in_file.seekg(0, std::ios::beg);
	in_file.close();
	in_file.open(file_name);
	if(in_file) {
		e_ld_ct = max_int;
		t_ld_ct = max_int;
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[0] == "loads") {
				if(headings[1] == "type" && data_len == 1) {
					i1 = elastic_list.find(data[0]);
					if (i1 > -1) {
						if (e_ld_ct == max_int) {
							e_ld_ct = 0;
						}
						else {
							e_ld_ct++;
						}
						elastic_loads[e_ld_ct].type = data[0];
						new_ld = &elastic_loads[e_ld_ct];
					}
					i1 = thermal_list.find(data[0]);
					if(i1 > -1) {
						if (t_ld_ct == max_int) {
							t_ld_ct = 0;
						}
						else {
							t_ld_ct++;
						}
						thermal_loads[t_ld_ct].type = data[0];
						new_ld = &thermal_loads[t_ld_ct];
					}
				} else if(headings[1] == "activeTime" && data_len > 0) {
					doub_inp[0] = stod(data[0]);
					if(data_len == 2) {
						doub_inp[1] = stod(data[1]);
					} else {
						doub_inp[1] = 1.0e+100;
					}
					new_ld->set_act_time(doub_inp);
				} else if(headings[1] == "load" && data_len > 0) {
					for (i1 = 0; i1 < 6; i1++) {
						if(i1 < data_len) {
							doub_inp[i1] = stod(data[i1]);
						} else {
							doub_inp[i1] = 0.0;
						}
					}
					new_ld->set_load(doub_inp);
				} else if(headings[1] == "nodeSet" && data_len == 1) {
					new_ld->node_set = data[0];
				} else if(headings[1] == "elementSet" && data_len == 1) {
					new_ld->element_set = data[0];
				} else if(headings[1] == "normDir" && data_len == 3) {
					doub_inp[0] = stod(data[0]);
					doub_inp[1] = stod(data[1]);
					doub_inp[2] = stod(data[2]);
					new_ld->set_norm_dir(doub_inp);
				} else if(headings[1] == "normTolerance" && data_len == 1) {
					new_ld->norm_tol = stod(data[0]);
				} else if(headings[1] == "center" && data_len == 3) {
					doub_inp[0] = stod(data[0]);
					doub_inp[1] = stod(data[1]);
					doub_inp[2] = stod(data[2]);
					new_ld->set_center(doub_inp);
				} else if(headings[1] == "axis" && data_len == 3) {
					doub_inp[0] = stod(data[0]);
					doub_inp[1] = stod(data[1]);
					doub_inp[2] = stod(data[2]);
					new_ld->set_axis(doub_inp);
				} else if(headings[1] == "angularVelocity") {
					new_ld->angular_vel = stod(data[0]);
				}
			}
		}
	} else {
		string err_st = "Error: could not open Load input file: " + file_name;
		throw invalid_argument(err_st);
	}
	in_file.close();
	
	return;
}

void Model::read_initial_state(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	
	int i1;
	int i2;
	int i3;
	int seti;
	double doub_inp[10] = { 0,0,0,0,0,0,0,0,0,0 };
	string disp_hdings = " displacement velocity acceleration";
	
	in_file.open(file_name);
	if(in_file) {
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[0] == "initialState") {
				i3 = disp_hdings.find(headings[1]);
				if(i3 > -1 && data_len > 3) {
					seti = ns_map.at(data[0]);
					for (auto& ndi : node_sets[seti].labels) {
						Node& this_nd = nodes[ndi];
						i2 = 1;
						for (i1 = 0; i1 < 6; i1++) {
							if (i2 < data_len) {
								doub_inp[i1] = stod(data[i2]);
							}
							else {
								doub_inp[i1] = 0.0;
							}
							i2++;
						}
						if (headings[1] == "displacement") {
							this_nd.set_initial_disp(doub_inp);
						}
						else if (headings[1] == "velocity") {
							this_nd.set_initial_vel(doub_inp);
						}
						else if (headings[1] == "acceleration") {
							this_nd.set_initial_acc(doub_inp);
						}
					}
				} else if(headings[1] == "temperature" && data_len == 2) {
					seti = ns_map.at(data[0]);
					for (auto& ndi : node_sets[seti].labels) {
						nodes[ndi].initial_temp = stod(data[1]);
					}
				} else if(headings[1] == "tdot" && data_len == 2) {
					seti = ns_map.at(data[0]);
					for (auto& ndi : node_sets[seti].labels) {
						nodes[ndi].initial_tdot = stod(data[1]);
					}
				}
			}
		}
	} else {
		string err_st = "Error: could not open initial state input file: " + file_name;
		throw invalid_argument(err_st);
	}

	in_file.close();
	
	return;
}

void Model::read_des_var_input(string file_name) {
    ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	
	int i1;
	double doub_inp[10] = {0,0,0,0,0,0,0,0,0,0};
	int int_inp[10] = {0,0,0,0,0,0,0,0,0,0};
	
	int dv_ct = 0;
	
	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "designVariables") {
				if (headings[1] == "category" && data_len == 1) {
					dv_ct++;
				}
			}
		}
	}

	design_vars = vector<DesignVariable>(dv_ct);

	//in_file.seekg(0, std::ios::beg);
	in_file.close();
	in_file.open(file_name);

	if(in_file) {
		dv_ct = max_int;
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[0] == "designVariables") {
				if(headings[1] == "category" && data_len == 1) {
					if (dv_ct == max_int) {
						dv_ct = 0;
					}
					else {
						dv_ct++;
					}
					design_vars[dv_ct].category = data[0];
				} else if(headings[1] == "elementSet" && data_len == 1) {
					design_vars[dv_ct].el_set_name = data[0];
				} else if(headings[1] == "nodeSet" && data_len == 1) {
					design_vars[dv_ct].nd_set_name = data[0];
				} else if(headings[1] == "activeTime" && data_len > 0) {
					doub_inp[0] = stod(data[0]);
					if(data_len == 2) {
						doub_inp[1] = stod(data[1]);
					} else {
						doub_inp[1] = 1.0e+100;
					}
					design_vars[dv_ct].set_active_time(doub_inp);
				} else if(headings[1] == "component") {
					if(data_len == 1) {
					    design_vars[dv_ct].component = stoi(data[0]);
					} else if(data_len == 2) {
						int_inp[0] = stoi(data[0]) - 1;
						int_inp[1] = stoi(data[1]) - 1;
						if(int_inp[0] >= int_inp[1]) {
							i1 = 6*int_inp[1] + int_inp[0];
						} else {
							i1 = 6*int_inp[0] + int_inp[1];
						}
						design_vars[dv_ct].component = i1;
					}
				} else if(headings[1] == "layer" && data_len == 1) {
					design_vars[dv_ct].layer = stoi(data[0]);
				} else if(headings[1] == "coefficients" && data_len == 1) {
					design_vars[dv_ct].coefs.push_back(stod(data[0]));
				}
			}
		}
	} else {
		cout << "entering !inFile block" << endl;
		string err_st = "Error: could not open design variable input file: " + file_name;
		throw invalid_argument(err_st);
	}
	in_file.close();
	
	return;
}

void Model::read_objective_input(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	
	double doub_inp[10] = { 0,0,0,0,0,0,0,0,0,0 };
	
	int ob_ct = 0;
	
	in_file.open(file_name);

	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "objectiveTerms") {
				if (headings[1] == "category" && data_len == 1) {
					ob_ct++;
				}
			}
		}
	}

	obj.terms = vector<ObjectiveTerm>(ob_ct);

	//in_file.seekg(0, std::ios::beg);
	in_file.close();
	in_file.open(file_name);

	if(in_file) {
		ob_ct = max_int;
		while(!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(headings[0] == "objectiveTerms") {
				if(headings[1] == "category" && data_len == 1) {
					if (ob_ct == max_int) {
						ob_ct = 0;
					}
					else {
						ob_ct++;
					}
					obj.terms[ob_ct].category = data[0];
				} else if(headings[1] == "operator" && data_len == 1) {
					obj.terms[ob_ct].optr = data[0];
				} else if(headings[1] == "activeTime" && data_len > 0) {
					doub_inp[0] = stod(data[0]);
					if(data_len == 2) {
						doub_inp[1] = stod(data[1]);
					} else {
						doub_inp[1] = 1.0e+100;
					}
					obj.terms[ob_ct].set_active_time(doub_inp);
				} else if(headings[1] == "component" && data_len == 1) {
					obj.terms[ob_ct].component = stoi(data[0]);
				} else if(headings[1] == "layer" && data_len == 1) {
					obj.terms[ob_ct].layer = stoi(data[0]);
				} else if(headings[1] == "coefficient" && data_len == 1) {
					obj.terms[ob_ct].coef = stod(data[0]);
				} else if(headings[1] == "exponent" && data_len == 1) {
					obj.terms[ob_ct].expnt = stod(data[0]);
				} else if(headings[1] == "elementSet" && data_len == 1) {
					obj.terms[ob_ct].el_set_name = data[0];
				} else if(headings[1] == "nodeSet" && data_len == 1) {
					obj.terms[ob_ct].nd_set_name = data[0];
				} else if(headings[1] == "targetValue" && data_len == 1) {
					if (is_doub(data[0])) {
						doub_inp[0] = stod(data[0]);
						obj.terms[ob_ct].tgt_vals.push_back(doub_inp[0]);
					} 
					else {
						obj.terms[ob_ct].tgt_tag = data[0];
					}
				}
			}
		}
	} else {
		string err_st = "Error: could not open Objective input file: " + file_name;
		throw invalid_argument(err_st);
	}

	in_file.close();
	
	return;
}

void Model::read_des_var_values(string file_name) {
	ifstream in_file;
	string file_line;
	string headings[4] = {"","","",""};
	int hd_ld_space[4] = {0,0,0,0};
	string data[11];
	int data_len;
	int label;
	double value;
	
	in_file.open(file_name);
	if(in_file) {
		while(!in_file.eof()) {
			getline(in_file, file_line);
		    read_input_line(file_line,headings,hd_ld_space,data,data_len);
			if(data_len == 2) {
				label = stoi(data[0]);
				value = stod(data[1]);
				design_vars[label].value.set_val(value);
			}
		}
	} else {
		string err_st = "Error: could not open design variable value input file: " + file_name;
		throw invalid_argument(err_st);		
	}

	in_file.close();
	
	return;
}

void Model::read_node_results(string file_name) {
	int i1;
	int i2;
	int nd;
	double nd_dat[6];
	Node* this_nd;
	ifstream in_file;
	string file_line;
	string headings[4] = { "","","","" };
	int hd_ld_space[4] = { 0,0,0,0 };
	string data[11];
	int data_len;
	string disp_fields = "displacement velocity acceleration";
	string thrm_fields = "temperature tdot";

	in_file.open(file_name);
	if (in_file) {
		while (!in_file.eof()) {
			getline(in_file, file_line);
			read_input_line(file_line, headings, hd_ld_space, data, data_len);
			if (headings[0] == "nodeResults" && data_len > 1) {
				nd = stoi(data[0]);
				Node& this_nd = nodes[nd];
				i2 = disp_fields.find(headings[1]);
				if (i2 > -1) {
					for (i1 = 0; i1 < this_nd.num_dof; i1++) {
						nd_dat[i1] = stod(data[i1 + 1]);
					}
					if (headings[1] == "displacement") {
						this_nd.set_displacement(nd_dat);
					}
					else if (headings[1] == "velocity") {
						this_nd.set_velocity(nd_dat);
					}
					else {
						this_nd.set_acceleration(nd_dat);
					}
				}
				i2 = thrm_fields.find(headings[1]);
				if (i2 > -1) {
					nd_dat[0] = stod(data[1]);
					if (headings[1] == "temperature") {
						this_nd.temperature = nd_dat[0];
					}
					else {
						this_nd.temp_change_rate = nd_dat[0];
					}
				}
			}
		}
	}
	else {
		string err_st = "Error: could not open Node result input file: " + file_name;
		throw invalid_argument(err_st);
	}

	in_file.close();
	
	return;
}

void Model::read_time_step_soln(int t_step) {
	int i1;
	int dof_per_nd;
	int num_int_dof;
	double* in_dat;
	char in_ln[72];
	string full_file = job[solve_cmd].file_name + "/solnTStep" + to_string(t_step) + ".out";
	ifstream in_file;
	in_file.open(full_file, std::ifstream::binary);

	if (in_file) {
		in_dat = reinterpret_cast<double*>(&in_ln[0]);
		for (auto& nd : nodes) {
			if (job[solve_cmd].thermal) {
				in_file.read(&in_ln[0], 8);
				nd.prev_temp = *in_dat;
				in_file.read(&in_ln[0], 8);
				nd.prev_tdot = *in_dat;
			}
			if (job[solve_cmd].elastic) {
				dof_per_nd = nd.num_dof;
				i1 = 8 * dof_per_nd;
				in_file.read(&in_ln[0], i1);
				nd.set_prev_disp(in_dat);
				in_file.read(&in_ln[0], i1);
				nd.set_prev_vel(in_dat);
				in_file.read(&in_ln[0], i1);
				nd.set_prev_acc(in_dat);
			}
		}

		if (job[solve_cmd].elastic) {
			for (auto& el : elements) {
				num_int_dof = el.num_int_dof;
				if (num_int_dof > 0) {
					i1 = 8 * num_int_dof;
					in_file.read(&in_ln[0], i1);
					el.set_int_prev_disp(in_dat);
				}
			}
		}

		in_file.close();
		return;
	}
	else {
		string err_st = "Error: could not open solutioin history file for time step: " + to_string(t_step);
		throw invalid_argument(err_st);
	}
}