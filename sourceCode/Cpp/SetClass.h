#ifndef SETCLASS
#define SETCLASS
#include <string>
#include "ListEntClass.h"

class Set {
	private:
	    std::string name;
		IntList labels;
		Set *nextSet;
		
    public:
	     Set(std::string newNm);
		 
		 void addEntry(int newLabel);
		 
		 std::string getName();
		 
		 int getLength();
		 
		 IntListEnt* getFirstEntry();
		 
		 Set* getNext();
		 
		 void setNext(Set* newNext);
};

class SetList {
	private:
	    Set *firstSet;
		Set *lastSet;
		int length;
		
	public:
	    SetList();
		
		void addSet(std::string newNm);
		
		int getLength();
		
		Set* getFirst();
};
#endif