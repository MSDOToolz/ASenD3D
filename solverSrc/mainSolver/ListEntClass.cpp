#include <cmath>
#include <cstddef>
#include <string>
#include <fstream>
#include "ListEntClass.h"

using namespace std;

IntListEnt::IntListEnt(int newVal) {
	value = newVal;
	next = nullptr;
}


IntList::IntList() {
	len = 0;
	first = nullptr;
	last = nullptr;
}

int IntList::getLength() {
	return len;
}

IntListEnt* IntList::getFirst() {
	return first;
}

IntListEnt* IntList::getLast() {
	return last;
}

void IntList::addEntry(int newInt) {
	IntListEnt *newEnt = new IntListEnt(newInt);
	if(!first) {
		first = newEnt;
		last = newEnt;
	} else {
		last->next = newEnt;
		last = newEnt;
	}
	len++;
	return;
}

void IntList::addIfAbsent(int newInt) {
	if(!first) {
		IntListEnt *newEnt = new IntListEnt(newInt);
		first = newEnt;
		last = newEnt;
		len++;
	} else {
		IntListEnt *thisEnt = first;
		bool inserted = false;
		while(thisEnt && !inserted) {
			if(thisEnt->value == newInt) {
				inserted = true;
			}
			thisEnt = thisEnt->next;
		}
		if(!inserted) {
			IntListEnt *newEnt = new IntListEnt(newInt);
			last->next = newEnt;
			last = newEnt;
			len++;
		}
	}
	return;
}

void IntList::destroy() {
	IntListEnt *thisEnt = first;
	IntListEnt *nextEnt;
	while(thisEnt) {
		nextEnt = thisEnt->next;
		delete thisEnt;
		thisEnt = nextEnt;
	}
	first = nullptr;
	last = nullptr;
	len = 0;
	return;
}


DoubListEnt::DoubListEnt(double newVal) {
	value = newVal;
	next = nullptr;
}

DoubList::DoubList() {
	first = nullptr;
	last = nullptr;
	len = 0;
}

int DoubList::getLength() {
	return len;
}

DoubListEnt* DoubList::getFirst() {
	return first;
}

void DoubList::addEntry(double newDoub) {
	DoubListEnt *newEnt = new DoubListEnt(newDoub);
	if(!first) {
		first = newEnt;
		last = newEnt;
	} else {
		last->next = newEnt;
		last = newEnt;
	}
	len++;
	return;
}

void DoubList::addIfAbsent(double newDoub) {
	if(!first) {
		DoubListEnt *newEnt = new DoubListEnt(newDoub);
		first = newEnt;
		last = newEnt;
		len++;
	} else {
		DoubListEnt *thisEnt = first;
		bool inserted = false;
		while(thisEnt && !inserted) {
			if(thisEnt->value == newDoub) {
				inserted = true;
			}
			thisEnt = thisEnt->next;
		}
		if(!inserted) {
			DoubListEnt *newEnt = new DoubListEnt(newDoub);
			last->next = newEnt;
			last = newEnt;
			len++;
		}
	}
	return;
}

void DoubList::destroy() {
	DoubListEnt *thisEnt = first;
	DoubListEnt *nextEnt;
	while(thisEnt) {
		nextEnt = thisEnt->next;
		delete thisEnt;
		thisEnt = nextEnt;
	}
	first = nullptr;
	last = nullptr;
	len = 0;
	return;
}

StringListEnt::StringListEnt(string newStr) {
	value = newStr;
	next = nullptr;
	return;
}

StringList::StringList() {
	len = 0;
	first = nullptr;
	last = nullptr;
	return;
}

int StringList::getLength() {
	return len;
}

StringListEnt* StringList::getFirst() {
	return first;
}

void StringList::addEntry(string newStr) {
	StringListEnt *newEnt = new StringListEnt(newStr);
	if(!first) {
		first = newEnt;
		last = newEnt;
	} else {
		last->next = newEnt;
		last = newEnt;
	}
	len++;
	return;
}

void StringList::destroy() {
	StringListEnt *thisEnt = first;
	StringListEnt *nextEnt;
	while(thisEnt) {
		nextEnt = thisEnt->next;
		delete thisEnt;
		thisEnt = nextEnt;
	}
	first = nullptr;
	last = nullptr;
	len = 0;
	return;
}

MatrixEnt::MatrixEnt(int newRow, int newCol, double newVal) {
	row = newRow;
	col = newCol;
	value = newVal;
	nextEnt = nullptr;
}


MEPtr::MEPtr () {
	ptr = nullptr;
}


SparseMat::SparseMat() {
	dim = 0;
	matrix = nullptr;
	return;
}

void SparseMat::setDim(int newDim) {
	dim = newDim;
	matrix = new MEPtr[newDim];
	return;
}

void SparseMat::zeroAll() {
	int i1;
	MatrixEnt *thisEnt;
	
	for (i1 = 0; i1 < dim; i1++) {
		thisEnt = matrix[i1].ptr;
		while(thisEnt) {
		    thisEnt->value = 0.0;
			thisEnt = thisEnt->nextEnt;
		}
	}
	
	return;
}

void SparseMat::addEntry(int row, int col, double val) {
	MatrixEnt *thisEnt = matrix[row].ptr;
	if(!thisEnt) {
		matrix[row].ptr = new MatrixEnt(row,col,val);
	} else {
		bool inserted = false;
		MatrixEnt *prevEnt = nullptr;
		while(thisEnt && !inserted) {
			if(thisEnt->col == col) {
				thisEnt->value+= val;
				inserted = true;
			}
			prevEnt = thisEnt;
			thisEnt = thisEnt->nextEnt;
		}
		if(!inserted) {
			MatrixEnt *newEnt = new MatrixEnt(row,col,val);
			prevEnt->nextEnt = newEnt;
		}
	}
	return;
}

int SparseMat::getDim() {
	return dim;
}

MatrixEnt* SparseMat::getFirstEnt(int row) {
	return matrix[row].ptr;
}

void SparseMat::vectorMultiply(double prod[], double inpVec[], bool transpose) {
	int i1;
	int col;
	MatrixEnt *thisEnt;
	if(transpose) {
		for (i1=0; i1<dim; i1++) {
			thisEnt = matrix[i1].ptr;
			while(thisEnt) {
				col = thisEnt->col;
				prod[col] += thisEnt->value*inpVec[i1];
				thisEnt = thisEnt->nextEnt;
			}
		}
	} else {
		for (i1=0; i1<dim; i1++) {
			thisEnt = matrix[i1].ptr;
			while(thisEnt) {
				col = thisEnt->col;
				prod[i1] += thisEnt->value*inpVec[col];
				thisEnt = thisEnt->nextEnt;
			}
		}
	}
	return;
}

double SparseMat::getMaxAbsVal() {
	int i1;
	MatrixEnt* thisEnt;
	double thisVal;
	double maxVal = 0.0;
	for (i1 = 0; i1 < dim; i1++) {
		thisEnt = matrix[i1].ptr;
		while (thisEnt) {
			thisVal = abs(thisEnt->value);
			if (thisVal > maxVal) {
				maxVal = thisVal;
			}
			thisEnt = thisEnt->nextEnt;
		}
	}
	return maxVal;
}

void SparseMat::writeToFile(ofstream& outFile) {
	int i1;
	MatrixEnt* thisEnt;
	for (i1 = 0; i1 < dim; i1++) {
		thisEnt = matrix[i1].ptr;
		while (thisEnt) {
			outFile << "    - [" << thisEnt->row << ", " << thisEnt->col << ", " << thisEnt->value << "]\n";
			thisEnt = thisEnt->nextEnt;
		}
	}
	return;
}

void SparseMat::destroy() {
	MatrixEnt *thisEnt;
	MatrixEnt *nextEnt;
	int i1;
	for (i1 = 0; i1 < dim; i1++){
		thisEnt = matrix[i1].ptr;
		while(thisEnt) {
			nextEnt = thisEnt->nextEnt;
			delete thisEnt;
			thisEnt = nextEnt;
		}
	}
	delete[] matrix;
	return;
}
