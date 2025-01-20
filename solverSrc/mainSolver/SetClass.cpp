#include "SetClass.h"
#include <string>
#include <vector>
#include "ListEntClass.h"

using namespace std;

Set::Set() {
    name = "";
	labels.clear();
	return;
}

bool Set::add_if_absent(int new_i) {
	for (auto& i : labels) {
		if (i == new_i) {
			return false;
		}
	}
	labels.push_back(new_i);
	return true;
}

