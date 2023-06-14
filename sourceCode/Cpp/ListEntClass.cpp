#include <cstddef>
#include <string>
#include "ListEntClass.h"

using namespace std;

IntListEnt::IntListEnt(int newVal) {
	value = newVal;
	next = NULL;
}


IntList::IntList() {
	len = 0;
	first = NULL;
	last = NULL;
}

int IntList::getLength() {
	return len;
}

IntListEnt* IntList::getFirst() {
	return first;
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
	first = NULL;
	last = NULL;
	len = 0;
	return;
}


DoubListEnt::DoubListEnt(double newVal) {
	value = newVal;
	next = NULL;
}

DoubList::DoubList() {
	first = NULL;
	last = NULL;
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
	first = NULL;
	last = NULL;
	len = 0;
	return;
}

StringListEnt::StringListEnt(string newStr) {
	value = newStr;
	next = NULL;
	return;
}

StringList::StringList() {
	len = 0;
	first = NULL;
	last = NULL;
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
	StringList *thisEnt = first;
	StringList *nextEnt;
	while(thisEnt) {
		nextEnt = thisEnt->next;
		delete thisEnt;
		thisEnt = nextEnt;
	}
	first = NULL;
	last = NULL;
	len = 0;
	return;
}

MatrixEnt::MatrixEnt(int newRow, int newCol, double newVal) {
	row = newRow;
	col = newCol;
	value = newVal;
	nextEnt = NULL;
}


MEPtr::MEPtr () {
	ptr = NULL;
}


SparseMat::SparseMat() {
	matrix = NULL;
	return;
}

void SparseMat::setDim(int newDim) {
	dim = newDim;
	matrix = new MEPtr[newDim];
	return;
}

void SparseMat::addEntry(int row, int col, double val) {
	MatrixEnt *thisEnt = matrix[row].ptr;
	if(!thisEnt) {
		matrix[row].ptr = new MatrixEnt(row,col,val);
	} else {
		bool inserted = false;
		MatrixEnt *prevEnt;
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

void SparseMat::vectorMultipy(double prod[], double inpVec[]) {
	int i1;
	int col;
	MatrixEnt *thisEnt;
	for (i1=0; i1<dim; i1++) {
		prod[i1] = 0.0;
		thisEnt = matrix[i1].ptr;
		while(thisEnt) {
			col = thisEnt->col;
			prod[i1] += thisEnt->value*inpVec[col];
			thisEnt = thisEnt->nextEnt;
		}
	}
	return;
}

void SparseMat::destroy() {
	MatrixEnt *thisEnt;
	MatrixEnt *nextEnt;
	int i1;
	for (i1 = 0; i1 < dim; i1++ ){
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
