#include "MatrixEntClass.cpp"

class SparseMat {
	private:
	    int dim;
	    MEPtr *matrix;
		
	public:
	    SparseMat(int matDim) {
			dim = matDim;
			matrix = new MEPtr[matDim];
		}
		
		void addEntry(int row, int col, double val) {
			MatrixEnt *thisEnt = matrix[row].ptr;
			if(!thisEnt) {
				matrix[row].ptr = new MatrixEnt(row,col,val);
			} else {
				bool inserted = false;
				MatrixEnt *prevEnt;
				while(thisEnt && !inserted) {
					if(thisEnt->getCol() == col) {
						thisEnt->addToValue(val);
						inserted = true;
					}
					prevEnt = thisEnt;
					thisEnt = thisEnt->getNext();
				}
				if(!inserted) {
					MatrixEnt *newEnt = new MatrixEnt(row,col,val);
					prevEnt->setNext(newEnt);
				}
			}
			return;
		}
		
		void vectorMultipy(double *prod, double *inpVec) {
			int i1;
			int col;
			for (i1=0; i1<dim; i1++) {
				prod[i1] = 0.0;
				MatrixEnt *thisEnt = matrix[i1].ptr;
				while(thisEnt) {
					col = thisEnt->getCol();
					prod[i1] += thisEnt->getValue()*inpVec[col];
				}
			}
			return;
		}
};