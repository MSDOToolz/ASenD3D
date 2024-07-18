#ifndef LISTENT
#define LISTENT
#include <string>
#include <fstream>

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
		
		~IntList();
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
		
		~DoubList();
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
		
		~StringList();
};

class MatrixEnt {
    public:
	    int row;
		int col;
		double value;
		MatrixEnt *nextEnt;
		
	    MatrixEnt(int newRow, int newCol, double newVal);
};

class SparseMat {
	private:
	    int dim;
	    MatrixEnt** matrix;
		
	public:
	    SparseMat();
		
		void setDim(int newDim);
		
		void zeroAll();
		
		void addEntry(int row, int col, double val);

		void addMatrix(SparseMat& inpMat);
		
		int getDim();
		
		MatrixEnt* getFirstEnt(int row);
		
		void vectorMultiply(double prod[], double inpVec[], bool transpose);

		double getMaxAbsVal();

		void writeToFile(std::ofstream& outFile);
		
		~SparseMat();
};
#endif