class IntListEnt {
	private:
	    int value;
		IntListEnt *next;
		
    public:
	    IntListEnt(int newVal) {
			value = newVal;
			next = NULL;
		}
		
		int getValue() {
			return value;
		}
		
		IntListEnt* getNext() {
			return next;
		}
		
		void setNext(IntListEnt* newNext) {
			next = newNext;
			return;
		}
};