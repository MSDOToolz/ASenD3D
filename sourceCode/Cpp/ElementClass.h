#ifndef ELEMENT
#define ELEMENT
#include "DiffDoubClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"

class Element {
	private:
	    int type;
		int label;
	    int *nodes;
		int numNds;
		int dofPerNd;
		int numIntDof;
		int nDim;
		int *dofTable;  //dofTable[i][0] = nd/basis function of dof i, dofTable[i][1] = component of dof i
		double *intPts;
		double *ipWt;
		int numIP;
		int *faceNds;
		int numFaces;
		double *internalDisp;
		IntList designVars;
		DoubList dvCoef;
		Section *sectPtr;
		Element *nextEl;
		
	public:
	    Element(int newType, int newLabel, int newNodes[]);
		
		void setNext(Element *newEl);
//dup1
        void getStiffMat(Doub& Cmat, DVPt& dvAr);

        void getNdCrds(Doub& xGlob, NdPt& ndAr, DVPt& dvAr);
		
		void getNdDisp(Doub& globDisp, NdPt& ndAr);

		void evalN(Doub& nVec, Doub& dNds, double spt[]);
		
		void getIpData(Doub& nVec, Doub& dNdx, Doub& detJ, Doub& locNds, double spt[]);
		
		void getInstOri(Doub& instOriMat, Doub& locOri, Doub& globDisp, bool stat);
		
		void getInstDisp(Doub& instDisp, Doub& globDisp, Doub& instOriMat, Doub& locOri, Doub& xGlob, int dv1, int dv2);
		
		void getSectionDef(Doub& secDef, Doub& globDisp,  Doub& instOriMat, Doub& locOri, Doub& xGlob, Doub& dNdx, Doub& nVec, int dv1, int dv2);
		
		void getSolidStrain(Doub& strain, Doub& ux, Doub& dNdx, Doub& locOri, int dv1, int dv2, bool nLGeom);
//end dup
};

class ElPt {
	public:
	    Element *ptr;
		
		ElPt();
};

class ElementList {
	private:
	    Element *firstEl;
		Element *lastEl;
		int length;
		
	public:
	    ElementList();
		
		void addElement(int newType, int newLabel, int newNodes[]);
		
		int getLength();
		
		Element* getFirst();
};

#endif