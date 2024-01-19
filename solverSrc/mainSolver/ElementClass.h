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
	double* scratch;

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
	double* scratch;

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
		FaceList* faces;
		double *internalDisp;
		double *intPrevDisp;
		double *internaldLdu;
		double *internalAdj;
		DiffDoub *internalRu;
		double *internalMat;
		IntList* designVars;
		DoubList* dvCoef;
		IntList* compDVars;
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
		void getLayerThkZ(Doub layThk[], Doub layZ[], Doub& zOffset, DesignVariable* dvAr[]);

		void getLayerQ(Doub layQ[], DesignVariable* dvAr[]);

		void getLayerD(Doub layD[], DesignVariable* dvAr[]);

		void getLayerAngle(Doub layAng[], DesignVariable* dvAr[]);

		void getLayerThExp(Doub layThExp[], DesignVariable* dvAr[]);

		void getLayerEinit(Doub layEinit[], DesignVariable* dvAr[]);

		void getLayerDen(Doub layerDen[], DesignVariable* dvAr[]);

		void getLayerCond(Doub layCond[], DesignVariable* dvAr[]);

		void getLayerSpecHeat(Doub laySH[], DesignVariable* dvAr[]);

		void transformStrain(Doub stnNew[], Doub stnOrig[], Doub& angle);

		void transformQ(Doub qNew[], Doub qOrig[], Doub& angle);

		void getSolidStiff(Doub Cmat[], DesignVariable* dvAr[]);

		void getABD(Doub Cmat[], Doub layThk[], Doub layZ[], Doub layQ[], Doub layAng[]);

		void getBeamStiff(Doub Cmat[], DesignVariable* dvAr[]);

		void getThermalExp(Doub thExp[], Doub Einit[], DesignVariable* dvAr[]);

		void getShellExpLoad(Doub expLd[], Doub E0Ld[], Doub layThk[], Doub layZ[], Doub layQ[], Doub layThExp[], Doub layEInit[], Doub layAng[]);

		void getBeamExpLoad(Doub expLd[], Doub E0Ld[], DesignVariable* dvAr[]);

		void getDensity(Doub& den, int layer, DesignVariable* dvAr[]);

		void getShellMass(Doub Mmat[], Doub layThk[], Doub layZ[], Doub layDen[], DesignVariable* dvAr[]);

		void getBeamMass(Doub Mmat[], DesignVariable* dvAr[]);

		void getSolidDamp(Doub Dmat[], DesignVariable* dvAr[]);

		void getShellDamp(Doub Dmat[], Doub layThk[], Doub layZ[], Doub layD[], Doub layAng[]);

		void getBeamDamp(Doub Dmat[], DesignVariable* dvAr[]);

		void getConductivity(Doub tCond[], DesignVariable* dvAr[]);

		void getShellCond(Doub tCond[], Doub layThk[], Doub layAng[], Doub layCond[], DesignVariable* dvAr[]);

		void getBeamCond(Doub tCond[], DesignVariable* dvAr[]);

		void getSpecificHeat(Doub& specHeat, DesignVariable* dvAr[]);

		void getShellSpecHeat(Doub& specHeat, Doub layThk[], Doub laySH[], Doub layDen[]);

		void getBeamSpecHeat(Doub& specHeat, DesignVariable* dvAr[]);

        void getNdCrds(Doub xGlob[], Node* ndAr[], DesignVariable* dvAr[]);
		
		void getLocOri(Doub locOri[], DesignVariable* dvAr[]);
		
		void correctOrient(Doub locOri[], Doub xGlob[]);

		void getFrcFldConst(Doub coef[], Doub exp[], DesignVariable* dvAr[]);

		void getMassPerEl(Doub& massPerEl, DesignVariable* dvAr[]);

// Solution Fields
		void getNdDisp(Doub globDisp[], Node* ndAr[]);

		void getNdVel(Doub globVel[], Node* ndAr[]);

		void getNdAcc(Doub globAcc[], Node* ndAr[]);

		void getNdTemp(Doub globTemp[], Node* ndAr[]);

		void getNdTdot(Doub globTdot[], Node* ndAr[]);

		void evalN(Doub nVec[], Doub dNds[], double spt[]);
		
		void getIpData(Doub nVec[], Doub dNdx[], Doub& detJ, Doub locNds[], double spt[]);
		
		void getInstOri(Doub instOriMat[], Doub locOri[], Doub globDisp[], int stat);
		
		void getInstDisp(Doub instDisp[], Doub globDisp[], Doub instOriMat[], Doub locOri[], Doub xGlob[], bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getVolume(Doub& vol, DoubStressPrereq& pre, int layer);
		
		void getSectionDef(Doub secDef[], Doub globDisp[],  Doub instOriMat[], Doub locOri[], Doub xGlob[], Doub dNdx[], Doub nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(Doub strain[], Doub ux[], Doub dNdx[], Doub locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(Doub stress[], Doub strain[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre);

		void dStressStraindU(Doub dsdU[], Doub dedU[], Doub dsdT[], double spt[], int layer, bool nLGeom, DoubStressPrereq& pre);

		void getDefFrcMom(Doub def[], Doub frcMom[], double spt[], bool nLGeom, DoubStressPrereq& pre);

		void dDefFrcMomdU(Doub dDefdU[], Doub dFrcMomdU[], Doub dFrcMomdT[], double spt[], bool nLGeom, DoubStressPrereq& pre);

		void getFluxTGrad(Doub flux[], Doub tGrad[], double spt[], int layer, DoubStressPrereq& pre);

		void dFluxTGraddT(Doub dFdT[], Doub dTG[], double spt[], int layer, DoubStressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, Doub elQVec[], bool forTherm, int matRow, Node* ndAr[]);
		
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1

// Properties
		void getLayerThkZ(DiffDoub layThk[], DiffDoub layZ[], DiffDoub& zOffset, DesignVariable* dvAr[]);

		void getLayerQ(DiffDoub layQ[], DesignVariable* dvAr[]);

		void getLayerD(DiffDoub layD[], DesignVariable* dvAr[]);

		void getLayerAngle(DiffDoub layAng[], DesignVariable* dvAr[]);

		void getLayerThExp(DiffDoub layThExp[], DesignVariable* dvAr[]);

		void getLayerEinit(DiffDoub layEinit[], DesignVariable* dvAr[]);

		void getLayerDen(DiffDoub layerDen[], DesignVariable* dvAr[]);

		void getLayerCond(DiffDoub layCond[], DesignVariable* dvAr[]);

		void getLayerSpecHeat(DiffDoub laySH[], DesignVariable* dvAr[]);

		void transformStrain(DiffDoub stnNew[], DiffDoub stnOrig[], DiffDoub& angle);

		void transformQ(DiffDoub qNew[], DiffDoub qOrig[], DiffDoub& angle);

		void getSolidStiff(DiffDoub Cmat[], DesignVariable* dvAr[]);

		void getABD(DiffDoub Cmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layQ[], DiffDoub layAng[]);

		void getBeamStiff(DiffDoub Cmat[], DesignVariable* dvAr[]);

		void getThermalExp(DiffDoub thExp[], DiffDoub Einit[], DesignVariable* dvAr[]);

		void getShellExpLoad(DiffDoub expLd[], DiffDoub E0Ld[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layQ[], DiffDoub layThExp[], DiffDoub layEInit[], DiffDoub layAng[]);

		void getBeamExpLoad(DiffDoub expLd[], DiffDoub E0Ld[], DesignVariable* dvAr[]);

		void getDensity(DiffDoub& den, int layer, DesignVariable* dvAr[]);

		void getShellMass(DiffDoub Mmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layDen[], DesignVariable* dvAr[]);

		void getBeamMass(DiffDoub Mmat[], DesignVariable* dvAr[]);

		void getSolidDamp(DiffDoub Dmat[], DesignVariable* dvAr[]);

		void getShellDamp(DiffDoub Dmat[], DiffDoub layThk[], DiffDoub layZ[], DiffDoub layD[], DiffDoub layAng[]);

		void getBeamDamp(DiffDoub Dmat[], DesignVariable* dvAr[]);

		void getConductivity(DiffDoub tCond[], DesignVariable* dvAr[]);

		void getShellCond(DiffDoub tCond[], DiffDoub layThk[], DiffDoub layAng[], DiffDoub layCond[], DesignVariable* dvAr[]);

		void getBeamCond(DiffDoub tCond[], DesignVariable* dvAr[]);

		void getSpecificHeat(DiffDoub& specHeat, DesignVariable* dvAr[]);

		void getShellSpecHeat(DiffDoub& specHeat, DiffDoub layThk[], DiffDoub laySH[], DiffDoub layDen[]);

		void getBeamSpecHeat(DiffDoub& specHeat, DesignVariable* dvAr[]);

        void getNdCrds(DiffDoub xGlob[], Node* ndAr[], DesignVariable* dvAr[]);
		
		void getLocOri(DiffDoub locOri[], DesignVariable* dvAr[]);
		
		void correctOrient(DiffDoub locOri[], DiffDoub xGlob[]);

		void getFrcFldConst(DiffDoub coef[], DiffDoub exp[], DesignVariable* dvAr[]);

		void getMassPerEl(DiffDoub& massPerEl, DesignVariable* dvAr[]);

// Solution Fields
		void getNdDisp(DiffDoub globDisp[], Node* ndAr[]);

		void getNdVel(DiffDoub globVel[], Node* ndAr[]);

		void getNdAcc(DiffDoub globAcc[], Node* ndAr[]);

		void getNdTemp(DiffDoub globTemp[], Node* ndAr[]);

		void getNdTdot(DiffDoub globTdot[], Node* ndAr[]);

		void evalN(DiffDoub nVec[], DiffDoub dNds[], double spt[]);
		
		void getIpData(DiffDoub nVec[], DiffDoub dNdx[], DiffDoub& detJ, DiffDoub locNds[], double spt[]);
		
		void getInstOri(DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub globDisp[], int stat);
		
		void getInstDisp(DiffDoub instDisp[], DiffDoub globDisp[], DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getVolume(DiffDoub& vol, DiffDoubStressPrereq& pre, int layer);
		
		void getSectionDef(DiffDoub secDef[], DiffDoub globDisp[],  DiffDoub instOriMat[], DiffDoub locOri[], DiffDoub xGlob[], DiffDoub dNdx[], DiffDoub nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(DiffDoub strain[], DiffDoub ux[], DiffDoub dNdx[], DiffDoub locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub stress[], DiffDoub strain[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre);

		void dStressStraindU(DiffDoub dsdU[], DiffDoub dedU[], DiffDoub dsdT[], double spt[], int layer, bool nLGeom, DiffDoubStressPrereq& pre);

		void getDefFrcMom(DiffDoub def[], DiffDoub frcMom[], double spt[], bool nLGeom, DiffDoubStressPrereq& pre);

		void dDefFrcMomdU(DiffDoub dDefdU[], DiffDoub dFrcMomdU[], DiffDoub dFrcMomdT[], double spt[], bool nLGeom, DiffDoubStressPrereq& pre);

		void getFluxTGrad(DiffDoub flux[], DiffDoub tGrad[], double spt[], int layer, DiffDoubStressPrereq& pre);

		void dFluxTGraddT(DiffDoub dFdT[], DiffDoub dTG[], double spt[], int layer, DiffDoubStressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, DiffDoub elQVec[], bool forTherm, int matRow, Node* ndAr[]);
		
//end dup
 
//end skip 
 
 
 
 
 
		void getElVec(double elVec[], double globVec[], bool forTherm, bool intnl, Node* ndAr[]);

		void addToGlobVec(double elVec[], double globVec[], bool forTherm, bool intnl, Node* ndAr[]);
 
 
// Equations
		void condenseMat(double mat[]);
		
		void updateExternal(double extVec[], int forSoln, Node* ndAr[]);
		
		void updateInternal(double extVec[], int forSoln, Node* ndAr[]);

		double getIntAdjdRdD();

//dup1
        void getRuk(Doub Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRum(Doub Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRud(Doub Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRu(Doub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRtk(Doub Rvec[], double dRdT[], bool getMatrix, DoubStressPrereq& pre);

		void getRtm(Doub Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DoubStressPrereq& pre);

		void getRt(Doub globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, Node* ndAr[]);

		void getRuFrcFld(Doub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DoubStressPrereq& pre, Node* ndAr[]);

		void getAppLoad(Doub AppLd[], Load* ldPt, bool nLGeom, DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getAppThermLoad(Doub AppLd[], Load* ldPt, DoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1
        void getRuk(DiffDoub Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRum(DiffDoub Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRud(DiffDoub Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRu(DiffDoub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRtk(DiffDoub Rvec[], double dRdT[], bool getMatrix, DiffDoubStressPrereq& pre);

		void getRtm(DiffDoub Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DiffDoubStressPrereq& pre);

		void getRt(DiffDoub globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, Node* ndAr[]);

		void getRuFrcFld(DiffDoub globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoubStressPrereq& pre, Node* ndAr[]);

		void getAppLoad(DiffDoub AppLd[], Load* ldPt, bool nLGeom, DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getAppThermLoad(DiffDoub AppLd[], Load* ldPt, DiffDoubStressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

//end dup
 
//end skip 
 
 
 
 
 
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