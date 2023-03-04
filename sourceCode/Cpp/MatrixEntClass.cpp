class MatrixEnt {
	private:
	    int row;
		int col;
		double value;
		MatrixEnt *nextEnt;
		
    public:
	    MatrixEnt(int newRow, int newCol, double newVal) {
			row = newRow;
			col = newCol;
			value = newVal;
			nextEnt = NULL;
		}
		
		MatrixEnt* getNext() {
			return nextEnt;
		}
		
		int getRow() {
			return row;
		}
		
		int getCol() {
			return col;
		}
		
		double getValue() {
			return value;
		}
		
		void addToValue(double val) {
			value += val;
		}
		
		void setNext(MatrixEnt* newEnt) {
			nextEnt = newEnt;
		}
};

class MEPtr {
	public:
	    MatrixEnt *ptr;
		
		MEPtr () {
			ptr = NULL;
		}
};