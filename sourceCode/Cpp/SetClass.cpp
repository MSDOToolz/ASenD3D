#include <string>
#include "IntListEntClass.cpp"

class Set {
	private:
	    string name;
		IntListEnt *firstEnt;
		IntListEnt *lastEnt;
		Set *nextSet;
		
    public:
	     Set(string newNm) {
			 name = newNm;
			 firstEnt = NULL;
			 lastEnt = NULL;
			 nextSet = NULL;
		 }
		 
		 void addEntry(int newLabel) {
			 IntListEnt *newEnt = new IntListEnt();
			 if(!firstEnt) {
				 firstEnt = newEnt;
				 lastEnt = newEnt;
			 } else {
				 lastEnt->setNext(newEnt);
				 lastEnt = newEnt;
			 }
		 }
		 
		 string getName() {
			 return name;
		 }
		 
		 IntListEnt* getFirstEntry() {
			 return firstEnt;
		 }
		 
		 Set* getNext() {
			 return nextSet;
		 }
		 
		 void setNext(Set* newNext) {
			 nextSet = newSet;
			 return;
		 }
};