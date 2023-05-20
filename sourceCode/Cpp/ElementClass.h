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
		int defDim;
		int *dofTable;  //dofTable[i][0] = nd/basis function of dof i, dofTable[i][1] = component of dof i
		double *intPts;
		double *ipWt;
		int numIP;
		int *faceNds;
		int numFaces;
		double *internalDisp;
		Doub *internalRu;
		double *internalMat;
		IntList designVars;
		DoubList dvCoef;
		Section *sectPtr;
		Element *nextEl;
		
	public:
	    Element(int newType);
		
		void setLabel(int newLab);

        void setNodes(int newNds[]);	
		
		void setNext(Element *newEl);
		
		Element* getNext();
//dup1

// Properties
        void getThickOffset(Doub& thickness, Doub& zOffset, DVPt dvAr[]);
		
        void getStiffMat(Doub Cmat[], DVPt dvAr[]);

        void getNdCrds(Doub xGlob[], NdPt ndAr[], DVPt dvAr[]);
		
		void getLocOri(Doub locOri[], DVPt dvAr[]);
		
		void correctOrient(Doub locOri[], Doub xGlob[]);

// Solution Fields		
		void getNdDisp(Doub globDisp[], NdPt ndAr[]);

		void evalN(Doub nVec[], Doub dNds[], double spt[]);
		
		void getIpData(Doub nVec[], Doub dNdx[], Doub& detJ, Doub locNds[], double spt[]);
		
		void getInstOri(Doub instOriMat[], Doub locOri[], Doub globDisp[], bool stat);
		
		void getInstDisp(Doub instDisp[], Doub globDisp[], Doub instOriMat[], Doub locOri[], Doub xGlob[], int dv1, int dv2);
		
		void getSectionDef(Doub secDef[], Doub globDisp[],  Doub instOriMat[], Doub locOri[], Doub xGlob[], Doub dNdx[], Doub nVec[], int dv1, int dv2);
		
		void getSolidStrain(Doub strain[], Doub ux[], Doub dNdx[], Doub locOri[], int dv1, int dv2, bool nLGeom);
		
// Equations

        void getRuk(Doub Rvec[], double dRdu[], bool getMatrix, bool nLGeom, NdPt ndAr[], DVPt dvAr[]);
		
		void condenseMat(double mat[]);
		
		void getRu(Doub globR[], SparseMat& globdRdu, bool getMatrix, bool dyn, bool nLGeom, NdPt ndAr[], DVPt dvAr[]);
		
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
		
		void addElement(Element *newEl);
		
		int getLength();
		
		Element* getFirst();
};

#endif