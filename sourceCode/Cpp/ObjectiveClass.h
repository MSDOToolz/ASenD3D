#ifndef OBJECTIVE
#define OBJECTIVE

class ObjectiveTerm {
	private:
	    string category;
		string optr;
		double activeTime[2];
		int component;
		int layer;
		double coef;
		double expnt;
		string elSetName;
		string ndSetName;
		string tgtTag;
		DoubList tgtVals;
		
	public:
	    ObjectiveTerm(string newCat);
};

class Objective {
	private:
	    ObjectiveTerm *firstTerm;
		ObjectiveTerm *lastTerm;
		int length;
		
	public:
	    Objective();
		
		void addTerm(ObjectiveTerm *newTerm);
		
		int getLength();
		
		ObjectiveTerm* getFirst();
};

#endif
