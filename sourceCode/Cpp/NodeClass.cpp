#include "IntListEntClass.cpp"
#include "DesignVariableClass.cpp"


class Node {
	private:
	    int label;
		int dofIndex[6];
	    double coord[3];
		double displacement[6];
		double velocity[6];
		double acceleration[6];
		double temperature;
		double tempChangeRate;
		IntListEnt *firstDVar;
		IntListEnt *lastDVar;
		DoubListEnt *firstCoef;
		DoubListEnt *lastCoef;
		Node *nextNd;
		
    public:
	    Node(int newLab, double newCrd[3]) {
			label = newLab;
			coord[0] = newCrd[0];
			coord[1] = newCrd[1];
			coord[2] = newCrd[2];
			int i1;
			for (i1=0; i1 < 6; i1++) {
				displacement[i1] = 0.0;
				velocity[i1] = 0.0;
				acceleration[i1] = 0.0;
				dofIndex[i1] = 0;
			}
			firstDVar = NULL;
			lastDVar = NULL;
			firstCoef = NULL;
			lastCoef = NULL;
			nextNd = NULL;
		}
		
		void addDesignVariable(int dIndex, double coef) {
			if(!firstDVar) {
				firstDVar = new IntListEnt(dIndex);
				lastDVar = firstDVar;
				firstCoef = new DoubListEnt(coef);
				lastCoef = firstCoef;
			} else {
				IntListEnt *newD = new IntListEnt(dIndex);
				lastDVar->setNext(newD);
				lastDVar = newD;
				DoubListEnt *newCoef = new DoubListEnt(coef);
				lastCoef->setNext(newCoef);
				lastCoef = newCoef;
			}
		}
		
		void r_getCrd(double r_crdOut[3], DVPt *dvAr) {
			r_crdOut[0] = r_1*coord[0];
			r_crdOut[1] = r_1*coord[1];
			r_crdOut[2] = r_1*coord[2];
			IntListEnt *currD = firstDVar;
			DoubListEnt *currCoef = firstCoef;
			int dIndex;
			DesignVariable *dPtr;
			double r_dVal;
			int comp;
			string cat;
			double coef;
			while(currD) {
				dIndex = currD->getValue();
				dPtr = dvAr[dIndex].ptr;
				r_dVal = dPtr->r_getValue();
				cat = dPtr->getCategory();
				comp = dPtr->getComponent() - 1;
				if(cat == "nodeCoord") {
					coef = currCoef->getValue();
					r_crdOut[comp] = r_crdOut[comp] + coef*r_dVal;
				}
				currD = currD->getNext();
				currCoef = currCoef->getNext();
			}
			return;
		}
		
		void destroy() {
			IntListEnt *currD = firstDVar;
			IntListEnt *cDNext;
			DoubListEnt *currCoef = firstCoef;
			DoubListEnt *cCNext;
			while(currD) {
				cDNext = currD->getNext();
				cCNext = currCoef->getNext();
				delete[] currD;
				delete[] currCoef;
				currD = cDNext;
				currCoef = cCNext;
			}
		}
};