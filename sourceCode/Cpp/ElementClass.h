#ifndef ELEMENT
#define ELEMENT
#include "DiffDoubClass.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"
#include "FaceClass.h"

//dup1
class DoubStressPrereq {
public:
	Doub globNds[30];
	Doub locNds[30];
	Doub locOri[9];
	Doub instOri[720];
	Doub globDisp[60];
	Doub Cmat[81];
	Doub* layerZ;
	Doub* layerThk;
	Doub* layerAng;
	Doub* layerQ;
	int currentLayLen;

	DoubStressPrereq();

	void destroy();
};
//end dup
//skip 

//DiffDoub versions: 
class DiffDoubStressPrereq {
public:
	DiffDoub globNds[30];
	DiffDoub locNds[30];
	DiffDoub locOri[9];
	DiffDoub instOri[720];
	DiffDoub globDisp[60];
	DiffDoub Cmat[81];
	DiffDoub* layerZ;
	DiffDoub* layerThk;
	DiffDoub* layerAng;
	DiffDoub* layerQ;
	int currentLayLen;

	DiffDoubStressPrereq();

	void destroy();
};

//end skip 

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
		int *dofTable;  //dofTable[2*i] = nd/basis function of dof i, dofTable[2*i+1] = component of dof i
		int intDofIndex;
		double* intPts;
		double* ipWt;
		int numIP;
		FaceList faces;
		double *internalDisp;
		double *internaldLdu;
		double *internalAdj;
		DiffDoub *internalRu;
		double *internalMat;
		IntList designVars;
		DoubList dvCoef;
		IntList compDVars;
		Section *sectPtr;
		Element *nextEl;
		
	public:
	    Element(int newType);
		
		void setLabel(int newLab);

        void setNodes(int newNds[]);

        void setSectPtr(Section *newSec);

        void initializeFaces();	

		void setIntDofIndex(int newInd);

		void setIntdLdU(double globdLdU[]);
		
		void setNext(Element *newEl);

		int getType();
		
		int getLabel();
		
		int getNumNds();
		
		int getDofPerNd();
		
		int getNumIntDof();

		int getNumIP();

		double* getIP();

		int getNumLayers();
		
		int* getNodes();
		
		IntList* getDesignVars();

		IntListEnt* getFirstCompDV();
		
		Face* getFirstFace();
		
		Element* getNext();
		
		void addDesignVariable(int dIndex, double coef);

		void addCompDVar(int dIndex);

		void destroy();
//dup1

// Properties
		void getLayerThkZ(Doub layThk[], Doub layZ[], Doub& zOffset, DVPt dvAr[]);

		void getLayerQ(Doub layQ[], DVPt dvAr[]);

		void getLayerAngle(Doub layAng[], DVPt dvAr[]);

		void transformStrain(Doub stnNew[], Doub stnOrig[], Doub& angle);

		void transformQ(Doub qNew[], Doub qOrig[], Doub& angle);

		void getSolidStiff(Doub Cmat[], DVPt dvAr[]);

		void getABD(Doub Cmat[], Doub layThk[], Doub layZ[], Doub layQ[], Doub layAng[]);

		void getBeamStiff(Doub Cmat[], DVPt dvAr[]);

		void getDensity(Doub& den, int layer, DVPt dvAr[]);

        void getNdCrds(Doub xGlob[], NdPt ndAr[], DVPt dvAr[]);
		
		void getLocOri(Doub locOri[], DVPt dvAr[]);
		
		void correctOrient(Doub locOri[], Doub xGlob[]);

// Solution Fields
		void getNdDisp(Doub globDisp[], NdPt ndAr[]);

		//double getDispSqr(Doub globDisp[]);

		void evalN(Doub nVec[], Doub dNds[], double spt[]);
		
		void getIpData(Doub nVec[], Doub dNdx[], Doub& detJ, Doub locNds[], double spt[]);
		
		void getInstOri(Doub instOriMat[], Doub locOri[], Doub globDisp[], bool stat);
		
		void getInstDisp(Doub instDisp[], Doub globDisp[], Doub instOriMat[], Doub locOri[], Doub xGlob[], int dv1, int dv2);

		void getStressPrereq(DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getVolume(Doub& vol, DoubStressPrereq& pre, int layer);
		
		void getSectionDef(Doub secDef[], Doub globDisp[],  Doub instOriMat[], Doub locOri[], Doub xGlob[], Doub dNdx[], Doub nVec[], int dv1, int dv2);
		
		void getSolidStrain(Doub strain[], Doub ux[], Doub dNdx[], Doub locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(Doub stress[], Doub strain[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre);

		void dStressStraindU(Doub dsdU[], Doub dedU[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, Doub elQVec[], int matRow, NdPt ndAr[]);
		
//end dup
//skip 
 
//DiffDoub versions: 

		void getLayerThkZ(DiffDoub layThk[], DiffDoub layZ[], DiffDoub& zOffset, DVPt dvAr[]);

		void getLayerQ(DiffDoub layQ[], DVPt dvAr[]);

		void getLayerAngle(DiffDoub layAng[], DVPt dvAr[]);

		void transformStrain(DiffDoub stnNew[], DiffDoub stnOrig[], DiffDoub& angle);

		void transformQ(DiffDoub qNew[], DiffDoub qOrig[], DiffDoub& angle);

		void getSolidStiff(DiffDoub Cmat[], DVPt dvAr[]);

		void getABD(DiffDoub Cmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layQ[], DiffDoub layAng[]);

		void getBeamStiff(DiffDoub Cmat[], DVPt dvAr[]);

		void getDensity(DiffDoub& den, int layer, DVPt dvAr[]);

        void getNdCrds(DiffDoub xGlob[], NdPt ndAr[], DVPt dvAr[]);
		
		void getLocOri(DiffDoub locOri[], DVPt dvAr[]);
		
		void correctOrient(DiffDoub locOri[], DiffDoub xGlob[]);

		void getNdDisp(DiffDoub globDisp[], NdPt ndAr[]);


		void evalN(DiffDoub nVec[], DiffDoub dNds[], double spt[]);
		
		void getIpData(DiffDoub nVec[], DiffDoub dNdx[], DiffDoub& detJ, DiffDoub locNds[], double spt[]);
		
		void getInstOri(DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub globDisp[], bool stat);
		
		void getInstDisp(DiffDoub instDisp[], DiffDoub globDisp[], DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], int dv1, int dv2);

		void getStressPrereq(DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getVolume(DiffDoub& vol, DiffDoubStressPrereq& pre, int layer);
		
		void getSectionDef(DiffDoub secDef[], DiffDoub globDisp[],  DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], DiffDoub dNdx[], DiffDoub nVec[], int dv1, int dv2);
		
		void getSolidStrain(DiffDoub strain[], DiffDoub ux[], DiffDoub dNdx[], DiffDoub locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub stress[], DiffDoub strain[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre);

		void dStressStraindU(DiffDoub dsdU[], DiffDoub dedU[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, DiffDoub elQVec[], int matRow, NdPt ndAr[]);
		
 
//end skip 
		
// Equations
		void condenseMat(double mat[]);
		
		void updateExternal(double extVec[], int forSoln, NdPt ndAr[]);
		
		void updateInternal(double extVec[], int forSoln, NdPt ndAr[]);

		double getIntAdjdRdD();

//dup1
        void getRuk(Doub Rvec[], double dRdu[], bool getMatrix, bool nLGeom, NdPt ndAr[], DVPt dvAr[]);
		
		void getRu(Doub globR[], SparseMat& globdRdu, bool getMatrix, bool dyn, bool nLGeom, NdPt ndAr[], DVPt dvAr[]);
		
//end dup
//skip 
 
//DiffDoub versions: 
        void getRuk(DiffDoub Rvec[], double dRdu[], bool getMatrix, bool nLGeom, NdPt ndAr[], DVPt dvAr[]);
		
		void getRu(DiffDoub globR[], SparseMat& globdRdu, bool getMatrix, bool dyn, bool nLGeom, NdPt ndAr[], DVPt dvAr[]);
		
 
//end skip 
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

		void destroy();
};



#endif