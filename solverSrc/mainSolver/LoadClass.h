#ifndef loadclass
#define loadclass
#include <string>
#include <cmath>
#include "SetClass.h"

class Load {
	public:
	    std::string type;
		double active_time[2];
		std::string node_set;
		int nd_set_ptr;
		std::string element_set;
		int el_set_ptr;
		double load[6];
		double normal_dir[3];
		double norm_tol;
		double center[3];
		double axis[3];
		double angular_vel;
		
	    Load();
		
		void set_act_time(double new_at[]);
		
		void set_load(double new_ld[]);
		
		void set_norm_dir(double new_ndir[]);
		
		void set_center(double new_cent[]);
		
		void set_axis(double new_axis[]);
};

#endif