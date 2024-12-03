#ifndef SETCLASS
#define SETCLASS
#include <string>
#include <vector>
#include "ListEntClass.h"

class Set {
	public:
	    std::string name;
		std::list<int> labels;
		
	    Set();

		bool addIfAbsent(int newI);
};

#endif