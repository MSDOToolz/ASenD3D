#ifndef LISTENT
#define LISTENT
#include <string>

class IntListEnt {
    public:
	    int value;
		IntListEnt *next;
	
	    IntListEnt(int newVal);
};

class IntList {
	private:
	    int len;
	    IntListEnt *first;
		IntListEnt *last;
		
    public:
	    IntList();
		
		int getLength();
		
		IntListEnt* getFirst();
		
		IntListEnt* getLast();
		
		void addEntry(int newInt);
		
		void addIfAbsent(int newInt);
		
		void destroy();
};

class DoubListEnt {
    public:
	    double value;
		DoubListEnt *next;
		
		DoubListEnt(double newVal);
};

class DoubList {
	private:
	    int len;
	    DoubListEnt *first;
		DoubListEnt *last;
		
    public:
	    DoubList();
		
		int getLength();
		
		DoubListEnt* getFirst();
		
		void addEntry(double newDoub);
		
		void addIfAbsent(double newDoub);
		
		void destroy();
};

class StringListEnt {
	public:
	    std::string value;
		StringListEnt *next;
		
		StringListEnt(std::string newStr);
};

class StringList {
	private:
	    int len;
		StringListEnt *first;
		StringListEnt *last;
		
	public:
	    StringList();
		
		int getLength();
		
		StringListEnt* getFirst();
		
		void addEntry(std::string newStr);
		
		void destroy();
};

class MatrixEnt {
    public:
	    int row;
		int col;
		double value;
		MatrixEnt *nextEnt;
		
	    MatrixEnt(int newRow, int newCol, double newVal);
};

class MEPtr {
	public:
	    MatrixEnt *ptr;
		
		MEPtr();
};

class SparseMat {
	private:
	    int dim;
	    MEPtr *matrix;
		
	public:
	    SparseMat();
		
		void setDim(int newDim);
		
		void addEntry(int row, int col, double val);
		
		int getDim();
		
		MatrixEnt* getFirstEnt(int row);
		
		void vectorMultipy(double prod[], double inpVec[]);
		
		void destroy();
};
#endif