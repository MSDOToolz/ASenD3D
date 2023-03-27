#ifndef OBJECTIVE
#define OBJECTIVE

class ObjectiveTerm {
	private:
	    std::string category;
		std::string optr;
		double activeTime[2];
		int component;
		int layer;
		double coef;
		double expnt;
		std::string elSetName;
		std::string ndSetName;
		std::string tgtTag;
		DoubList tgtVals;
		
	public:
	    ObjectiveTerm(std::string newCat);
		
		void setOperator(std::string newOp);
		
		void setActiveTime(double newAt[]);
		
		void setComponent(int newComp);
		
		void setLayer(int newLay);
		
		void setCoef(double newCoef);
		
		void setExponent(double newExp);
		
		void setElset(std::string newElset);
		
		void setNdset(std::string newNdset);
		
		void setTgtTag(std::string newTag);
		
		void addTargetValue(double newTgt);
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
