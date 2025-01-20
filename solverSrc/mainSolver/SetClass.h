#ifndef setclass
#define setclass
#include <string>
#include <vector>
#include "ListEntClass.h"

class Set {
	public:
	    std::string name;
		std::list<int> labels;
		
	    Set();

		bool add_if_absent(int new_i);
};

#endif