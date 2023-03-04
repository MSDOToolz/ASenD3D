class DoubListEnt {
	private:
	    double value;
		DoubListEnt *next;
		
    public:
	    DoubListEnt(double newVal) {
			value = newVal;
			next = NULL;
		}
		
		double getValue() {
			return value;
		}
		
		DoubListEnt* getNext() {
			return next;
		}
		
		void setNext(DoubListEnt* newNext) {
			next = newNext;
			return;
		}
};