#ifndef ELEMENT
#define ELEMENT
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"
#include "FaceClass.h"
#include "LoadClass.h"
#include "JobClass.h"

//dup1
class DoubStressPrereq {
public:
	Doub* globNds;  //globNds[30];
	Doub* locNds;// [30] ;
	Doub* locOri;// [9] ;
	Doub* instOri;// [720] ;
	Doub* globDisp;// [60] ;
	Doub* globVel;
	Doub* globAcc;// [30];
	Doub* globTemp;
	Doub* globTdot;
	Doub* Cmat;// [81] ;
	Doub* Mmat;// [36];
	Doub* Dmat;
	Doub* thermExp;
	Doub* Einit;
	Doub* TCmat;
	Doub SpecHeat;
	Doub* BMat;
	Doub* CBMat;
	Doub* layerZ;
	Doub* layerThk;
	Doub* layerAng;
	Doub* layerQ;
	Doub* layerD;
	Doub* layerTE;
	Doub* layerE0;
	Doub* layerDen;
	Doub* layerTC;
	Doub* layerSH;
	Doub* frcFldCoef;
	Doub* frcFldExp;
	Doub massPerEl;

	int currentLayLen;

	DoubStressPrereq();

	void destroy();
};
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1
class DiffDoubStressPrereq {
public:
	DiffDoub* globNds;  //globNds[30];
	DiffDoub* locNds;// [30] ;
	DiffDoub* locOri;// [9] ;
	DiffDoub* instOri;// [720] ;
	DiffDoub* globDisp;// [60] ;
	DiffDoub* globVel;
	DiffDoub* globAcc;// [30];
	DiffDoub* globTemp;
	DiffDoub* globTdot;
	DiffDoub* Cmat;// [81] ;
	DiffDoub* Mmat;// [36];
	DiffDoub* Dmat;
	DiffDoub* thermExp;
	DiffDoub* Einit;
	DiffDoub* TCmat;
	DiffDoub SpecHeat;
	DiffDoub* BMat;
	DiffDoub* CBMat;
	DiffDoub* layerZ;
	DiffDoub* layerThk;
	DiffDoub* layerAng;
	DiffDoub* layerQ;
	DiffDoub* layerD;
	DiffDoub* layerTE;
	DiffDoub* layerE0;
	DiffDoub* layerDen;
	DiffDoub* layerTC;
	DiffDoub* layerSH;
	DiffDoub* frcFldCoef;
	DiffDoub* frcFldExp;
	DiffDoub massPerEl;

	int currentLayLen;

	DiffDoubStressPrereq();

	void destroy();
};
//end dup
 
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
		double *intPrevDisp;
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

		void setIntDisp(double newDisp[]);

		void setIntPrevDisp(double newDisp[]);

		void advanceIntDisp();

		void backstepIntDisp();

		void setIntdLdU(double globdLdU[]);
		
		void setNext(Element *newEl);

		int getType();
		
		int getLabel();
		
		int getNumNds();
		
		int getDofPerNd();
		
		int getNumIntDof();

		void getIntDisp(double dispOut[]);

		void getIntPrevDisp(double dispOut[]);

		int getNumIP();

		double* getIP();

		double* getIntAdj();

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

		void getLayerD(Doub layD[], DVPt dvAr[]);

		void getLayerAngle(Doub layAng[], DVPt dvAr[]);

		void getLayerThExp(Doub layThExp[], DVPt dvAr[]);

		void getLayerEinit(Doub layEinit[], DVPt dvAr[]);

		void getLayerDen(Doub layerDen[], DVPt dvAr[]);

		void getLayerCond(Doub layCond[], DVPt dvAr[]);

		void getLayerSpecHeat(Doub laySH[], DVPt dvAr[]);

		void transformStrain(Doub stnNew[], Doub stnOrig[], Doub& angle);

		void transformQ(Doub qNew[], Doub qOrig[], Doub& angle);

		void getSolidStiff(Doub Cmat[], DVPt dvAr[]);

		void getABD(Doub Cmat[], Doub layThk[], Doub layZ[], Doub layQ[], Doub layAng[]);

		void getBeamStiff(Doub Cmat[], DVPt dvAr[]);

		void getThermalExp(Doub thExp[], Doub Einit[], DVPt dvAr[]);

		void getShellExpLoad(Doub expLd[], Doub E0Ld[], Doub layThk[], Doub layZ[], Doub layQ[], Doub layThExp[], Doub layEInit[], Doub layAng[]);

		void getBeamExpLoad(Doub expLd[], Doub E0Ld[], DVPt dvAr[]);

		void getDensity(Doub& den, int layer, DVPt dvAr[]);

		void getShellMass(Doub Mmat[], Doub layThk[], Doub layZ[], Doub layDen[], DVPt dvAr[]);

		void getBeamMass(Doub Mmat[], DVPt dvAr[]);

		void getSolidDamp(Doub Dmat[], DVPt dvAr[]);

		void getShellDamp(Doub Dmat[], Doub layThk[], Doub layZ[], Doub layD[], Doub layAng[]);

		void getBeamDamp(Doub Dmat[], DVPt dvAr[]);

		void getConductivity(Doub tCond[], DVPt dvAr[]);

		void getShellCond(Doub tCond[], Doub layThk[], Doub layAng[], Doub layCond[], DVPt dvAr[]);

		void getBeamCond(Doub tCond[], DVPt dvAr[]);

		void getSpecificHeat(Doub& specHeat, DVPt dvAr[]);

		void getShellSpecHeat(Doub& specHeat, Doub layThk[], Doub laySH[], Doub layDen[]);

		void getBeamSpecHeat(Doub& specHeat, DVPt dvAr[]);

        void getNdCrds(Doub xGlob[], NdPt ndAr[], DVPt dvAr[]);
		
		void getLocOri(Doub locOri[], DVPt dvAr[]);
		
		void correctOrient(Doub locOri[], Doub xGlob[]);

		void getFrcFldConst(Doub coef[], Doub exp[], DVPt dvAr[]);

		void getMassPerEl(Doub& massPerEl, DVPt dvAr[]);

// Solution Fields
		void getNdDisp(Doub globDisp[], NdPt ndAr[]);

		void getNdVel(Doub globVel[], NdPt ndAr[]);

		void getNdAcc(Doub globAcc[], NdPt ndAr[]);

		void getNdTemp(Doub globTemp[], NdPt ndAr[]);

		void getNdTdot(Doub globTdot[], NdPt ndAr[]);

		void evalN(Doub nVec[], Doub dNds[], double spt[]);
		
		void getIpData(Doub nVec[], Doub dNdx[], Doub& detJ, Doub locNds[], double spt[]);
		
		void getInstOri(Doub instOriMat[], Doub locOri[], Doub globDisp[], int stat);
		
		void getInstDisp(Doub instDisp[], Doub globDisp[], Doub instOriMat[], Doub locOri[], Doub xGlob[], bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getVolume(Doub& vol, DoubStressPrereq& pre, int layer);
		
		void getSectionDef(Doub secDef[], Doub globDisp[],  Doub instOriMat[], Doub locOri[], Doub xGlob[], Doub dNdx[], Doub nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(Doub strain[], Doub ux[], Doub dNdx[], Doub locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(Doub stress[], Doub strain[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre);

		void dStressStraindU(Doub dsdU[], Doub dedU[], Doub dsdT[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre);

		void getDefFrcMom(Doub def[], Doub frcMom[], double spt[], bool nLGeom, DoubStressPrereq& pre);

		void dDefFrcMomdU(Doub dDefdU[], Doub dFrcMomdU[], Doub dFrcMomdT[], double spt[], bool nLGeom, DoubStressPrereq& pre);

		void getFluxTGrad(Doub flux[], Doub tGrad[], double spt[], int layer, DoubStressPrereq& pre);

		void dFluxTGraddT(Doub dFdT[], Doub dTG[], double spt[], int layer, DoubStressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, Doub elQVec[], bool forTherm, int matRow, NdPt ndAr[]);
		
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1

// Properties
		void getLayerThkZ(DiffDoub layThk[], DiffDoub layZ[], DiffDoub& zOffset, DVPt dvAr[]);

		void getLayerQ(DiffDoub layQ[], DVPt dvAr[]);

		void getLayerD(DiffDoub layD[], DVPt dvAr[]);

		void getLayerAngle(DiffDoub layAng[], DVPt dvAr[]);

		void getLayerThExp(DiffDoub layThExp[], DVPt dvAr[]);

		void getLayerEinit(DiffDoub layEinit[], DVPt dvAr[]);

		void getLayerDen(DiffDoub layerDen[], DVPt dvAr[]);

		void getLayerCond(DiffDoub layCond[], DVPt dvAr[]);

		void getLayerSpecHeat(DiffDoub laySH[], DVPt dvAr[]);

		void transformStrain(DiffDoub stnNew[], DiffDoub stnOrig[], DiffDoub& angle);

		void transformQ(DiffDoub qNew[], DiffDoub qOrig[], DiffDoub& angle);

		void getSolidStiff(DiffDoub Cmat[], DVPt dvAr[]);

		void getABD(DiffDoub Cmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layQ[], DiffDoub layAng[]);

		void getBeamStiff(DiffDoub Cmat[], DVPt dvAr[]);

		void getThermalExp(DiffDoub thExp[], DiffDoub Einit[], DVPt dvAr[]);

		void getShellExpLoad(DiffDoub expLd[], DiffDoub E0Ld[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layQ[], DiffDoub layThExp[], DiffDoub layEInit[], DiffDoub layAng[]);

		void getBeamExpLoad(DiffDoub expLd[], DiffDoub E0Ld[], DVPt dvAr[]);

		void getDensity(DiffDoub& den, int layer, DVPt dvAr[]);

		void getShellMass(DiffDoub Mmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layDen[], DVPt dvAr[]);

		void getBeamMass(DiffDoub Mmat[], DVPt dvAr[]);

		void getSolidDamp(DiffDoub Dmat[], DVPt dvAr[]);

		void getShellDamp(DiffDoub Dmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layD[], DiffDoub layAng[]);

		void getBeamDamp(DiffDoub Dmat[], DVPt dvAr[]);

		void getConductivity(DiffDoub tCond[], DVPt dvAr[]);

		void getShellCond(DiffDoub tCond[], DiffDoub layThk[], DiffDoub layAng[], DiffDoub layCond[], DVPt dvAr[]);

		void getBeamCond(DiffDoub tCond[], DVPt dvAr[]);

		void getSpecificHeat(DiffDoub& specHeat, DVPt dvAr[]);

		void getShellSpecHeat(DiffDoub& specHeat, DiffDoub layThk[], DiffDoub laySH[], DiffDoub layDen[]);

		void getBeamSpecHeat(DiffDoub& specHeat, DVPt dvAr[]);

        void getNdCrds(DiffDoub xGlob[], NdPt ndAr[], DVPt dvAr[]);
		
		void getLocOri(DiffDoub locOri[], DVPt dvAr[]);
		
		void correctOrient(DiffDoub locOri[], DiffDoub xGlob[]);

		void getFrcFldConst(DiffDoub coef[], DiffDoub exp[], DVPt dvAr[]);

		void getMassPerEl(DiffDoub& massPerEl, DVPt dvAr[]);

// Solution Fields
		void getNdDisp(DiffDoub globDisp[], NdPt ndAr[]);

		void getNdVel(DiffDoub globVel[], NdPt ndAr[]);

		void getNdAcc(DiffDoub globAcc[], NdPt ndAr[]);

		void getNdTemp(DiffDoub globTemp[], NdPt ndAr[]);

		void getNdTdot(DiffDoub globTdot[], NdPt ndAr[]);

		void evalN(DiffDoub nVec[], DiffDoub dNds[], double spt[]);
		
		void getIpData(DiffDoub nVec[], DiffDoub dNdx[], DiffDoub& detJ, DiffDoub locNds[], double spt[]);
		
		void getInstOri(DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub globDisp[], int stat);
		
		void getInstDisp(DiffDoub instDisp[], DiffDoub globDisp[], DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getVolume(DiffDoub& vol, DiffDoubStressPrereq& pre, int layer);
		
		void getSectionDef(DiffDoub secDef[], DiffDoub globDisp[],  DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], DiffDoub dNdx[], DiffDoub nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(DiffDoub strain[], DiffDoub ux[], DiffDoub dNdx[], DiffDoub locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub stress[], DiffDoub strain[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre);

		void dStressStraindU(DiffDoub dsdU[], DiffDoub dedU[], DiffDoub dsdT[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre);

		void getDefFrcMom(DiffDoub def[], DiffDoub frcMom[], double spt[], bool nLGeom, DiffDoubStressPrereq& pre);

		void dDefFrcMomdU(DiffDoub dDefdU[], DiffDoub dFrcMomdU[], DiffDoub dFrcMomdT[], double spt[], bool nLGeom, DiffDoubStressPrereq& pre);

		void getFluxTGrad(DiffDoub flux[], DiffDoub tGrad[], double spt[], int layer, DiffDoubStressPrereq& pre);

		void dFluxTGraddT(DiffDoub dFdT[], DiffDoub dTG[], double spt[], int layer, DiffDoubStressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, DiffDoub elQVec[], bool forTherm, int matRow, NdPt ndAr[]);
		
//end dup
 
//end skip 
 
 
		void getElVec(double elVec[], double globVec[], bool forTherm, bool intnl, NdPt ndAr[]);

		void addToGlobVec(double elVec[], double globVec[], bool forTherm, bool intnl, NdPt ndAr[]);
 
 
// Equations
		void condenseMat(double mat[]);
		
		void updateExternal(double extVec[], int forSoln, NdPt ndAr[]);
		
		void updateInternal(double extVec[], int forSoln, NdPt ndAr[]);

		double getIntAdjdRdD();

//dup1
        void getRuk(Doub Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getRum(Doub Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);
		
		void getRud(Doub Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getRu(Doub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);
		
		void getRtk(Doub Rvec[], double dRdT[], bool getMatrix, DoubStressPrereq& pre);

		void getRtm(Doub Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DoubStressPrereq& pre);

		void getRt(Doub globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, NdPt ndAr[]);

		void getRuFrcFld(Doub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, NdPt ndAr[]);

		void getAppLoad(Doub AppLd[], Load* ldPt, bool nLGeom, DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getAppThermLoad(Doub AppLd[], Load* ldPt, DoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1
        void getRuk(DiffDoub Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getRum(DiffDoub Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);
		
		void getRud(DiffDoub Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getRu(DiffDoub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);
		
		void getRtk(DiffDoub Rvec[], double dRdT[], bool getMatrix, DiffDoubStressPrereq& pre);

		void getRtm(DiffDoub Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DiffDoubStressPrereq& pre);

		void getRt(DiffDoub globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, NdPt ndAr[]);

		void getRuFrcFld(DiffDoub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, NdPt ndAr[]);

		void getAppLoad(DiffDoub AppLd[], Load* ldPt, bool nLGeom, DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

		void getAppThermLoad(DiffDoub AppLd[], Load* ldPt, DiffDoubStressPrereq& pre, NdPt ndAr[], DVPt dvAr[]);

//end dup
 
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