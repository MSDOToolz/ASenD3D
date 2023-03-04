#include <string>
#include "ConstantVals.cpp"
#include "DoubListEntClass.cpp"
using namespace std;

class DesignVariable {
	private:
	    string category;
		int component;
		int layer;
		double activeTime[2];
		string elSetName;
		string ndlSetName;
		DoubListEnt *firstCoef;
		DoubListEnt *lastCoef;
		double r_value;
		DesignVariable *nextDV;
		
    public:
	    DesignVariable(string newCat,int newComp=1,int newLay=0,string newES="",string newNS="") {
			category = newCat;
			component = newComp;
			layer = newLay;
			elSetName = newES;
			ndlSetName = newNS;
			activeTime[0] = 0.0;
			activeTime[1] = 1.0e+10;
			r_value = 0.0;
			firstCoef = NULL;
			lastCoef = NULL;
			nextDV = NULL;
		}
		
		void setActiveTime(double st, double fn) {
			activeTime[0] = st;
			activeTime[1] = fn;
		}
		
		string getCategory() {
			return category;
		}
		
		double r_getValue() {
			return r_value;
		}
		
		int getComponent() {
			return component;
		}
		
		DesignVariable* getNext() {
			return nextDV;
		}
		
		void r_setValue(double r_newVal) {
			r_value = r_newVal;
		}
		
		void addCoefficient(double newCoef) {
			DoubListEnt *newPt = new DoubListEnt(newCoef);
			if(!firstCoef) {
				firstCoef = newPt;
				lastCoef = newPt;
			} else {
				lastCoef->setNext(newPt);
				lastCoef = newPt;
			}
		}
};

class DVPt {
	public:
	    DesignVariable *ptr;
		
		DVPt() {
			ptr = NULL;
		}
};