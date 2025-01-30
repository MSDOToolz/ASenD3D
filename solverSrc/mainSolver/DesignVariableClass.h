#ifndef dvarclass
#define dvarclass
#include <string>
#include <vector>
#include "ListEntClass.h"
#include "DiffDoubClass.h"
#include "SetClass.h"

class DesignVariable {
	public:
	    std::string category;
		int component;
		int layer;
		double active_time[2];
		std::string el_set_name;
		int el_set_ptr;
		std::string nd_set_name;
		int nd_set_ptr;
		std::list<double> coefs;
		DiffDoub0 value;
		DiffDoub1 diff_val;
		std::list<int> comp_el_list;
		
	    DesignVariable();
		
		void set_active_time(double new_at[]);
		
		void get_value_dfd0(DiffDoub0& inp);
		
		void get_value_dfd1(DiffDoub1& inp);

		void add_comp_el(int eli);

};

#endif