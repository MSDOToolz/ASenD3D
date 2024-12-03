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

bool Set::addIfAbsent(int newI) {
	for (auto& i : labels) {
		if (i == newI) {
			return false;
		}
	}
	labels.push_back(newI);
	return true;
}

