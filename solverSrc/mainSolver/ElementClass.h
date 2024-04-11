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
class DiffDoub0StressPrereq {
public:
	DiffDoub0* globNds;  //globNds[30];
	DiffDoub0* locNds;// [30] ;
	DiffDoub0* locOri;// [9] ;
	DiffDoub0* instOri;// [720] ;
	DiffDoub0* globDisp;// [60] ;
	DiffDoub0* globVel;
	DiffDoub0* globAcc;// [30];
	DiffDoub0* globTemp;
	DiffDoub0* globTdot;
	DiffDoub0* Cmat;// [81] ;
	DiffDoub0* Mmat;// [36];
	DiffDoub0* Dmat;
	DiffDoub0* thermExp;
	DiffDoub0* Einit;
	DiffDoub0* TCmat;
	DiffDoub0 SpecHeat;
	DiffDoub0* BMat;
	DiffDoub0* CBMat;
	DiffDoub0* layerZ;
	DiffDoub0* layerThk;
	DiffDoub0* layerAng;
	DiffDoub0* layerQ;
	DiffDoub0* layerD;
	DiffDoub0* layerTE;
	DiffDoub0* layerE0;
	DiffDoub0* layerDen;
	DiffDoub0* layerTC;
	DiffDoub0* layerSH;
	DiffDoub0* frcFldCoef;
	DiffDoub0* frcFldExp;
	DiffDoub0 massPerEl;
	double* scratch;

	int currentLayLen;

	DiffDoub0StressPrereq();

	void allocateLayers(int numLayers);

	~DiffDoub0StressPrereq();
};
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
class DiffDoub1StressPrereq {
public:
	DiffDoub1* globNds;  //globNds[30];
	DiffDoub1* locNds;// [30] ;
	DiffDoub1* locOri;// [9] ;
	DiffDoub1* instOri;// [720] ;
	DiffDoub1* globDisp;// [60] ;
	DiffDoub1* globVel;
	DiffDoub1* globAcc;// [30];
	DiffDoub1* globTemp;
	DiffDoub1* globTdot;
	DiffDoub1* Cmat;// [81] ;
	DiffDoub1* Mmat;// [36];
	DiffDoub1* Dmat;
	DiffDoub1* thermExp;
	DiffDoub1* Einit;
	DiffDoub1* TCmat;
	DiffDoub1 SpecHeat;
	DiffDoub1* BMat;
	DiffDoub1* CBMat;
	DiffDoub1* layerZ;
	DiffDoub1* layerThk;
	DiffDoub1* layerAng;
	DiffDoub1* layerQ;
	DiffDoub1* layerD;
	DiffDoub1* layerTE;
	DiffDoub1* layerE0;
	DiffDoub1* layerDen;
	DiffDoub1* layerTC;
	DiffDoub1* layerSH;
	DiffDoub1* frcFldCoef;
	DiffDoub1* frcFldExp;
	DiffDoub1 massPerEl;
	double* scratch;

	int currentLayLen;

	DiffDoub1StressPrereq();

	void allocateLayers(int numLayers);

	~DiffDoub1StressPrereq();
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
		double sCent[3];
		int numIP;
		FaceList faces;
		double *internalDisp;
		double *intPrevDisp;
		double *internaldLdu;
		double *internalAdj;
		DiffDoub1* internalRu;
		double* internalMat;
		IntList designVars;
		DoubList dvCoef;
		IntList compDVars;
		Section* sectPtr;
		Element* nextEl;
		
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

		double* getSCent();

		double* getIntAdj();

		int getNumLayers();
		
		int* getNodes();
		
		IntList* getDesignVars();

		IntListEnt* getFirstCompDV();
		
		Face* getFirstFace();
		
		Element* getNext();
		
		void addDesignVariable(int dIndex, double coef);

		void addCompDVar(int dIndex);

		~Element();
//dup1

// Properties
		void getLayerThkZ(DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0& zOffset, DesignVariable* dvAr[]);

		void getLayerQ(DiffDoub0 layQ[], DesignVariable* dvAr[]);

		void getLayerD(DiffDoub0 layD[], DesignVariable* dvAr[]);

		void getLayerAngle(DiffDoub0 layAng[], DesignVariable* dvAr[]);

		void getLayerThExp(DiffDoub0 layThExp[], DesignVariable* dvAr[]);

		void getLayerEinit(DiffDoub0 layEinit[], DesignVariable* dvAr[]);

		void getLayerDen(DiffDoub0 layerDen[], DesignVariable* dvAr[]);

		void getLayerCond(DiffDoub0 layCond[], DesignVariable* dvAr[]);

		void getLayerSpecHeat(DiffDoub0 laySH[], DesignVariable* dvAr[]);

		void transformStrain(DiffDoub0 stnNew[], DiffDoub0 stnOrig[], DiffDoub0& angle);

		void transformQ(DiffDoub0 qNew[], DiffDoub0 qOrig[], DiffDoub0& angle);

		void getSolidStiff(DiffDoub0 Cmat[], DesignVariable* dvAr[]);

		void getABD(DiffDoub0 Cmat[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layQ[], DiffDoub0 layAng[]);

		void getBeamStiff(DiffDoub0 Cmat[], DesignVariable* dvAr[]);

		void getThermalExp(DiffDoub0 thExp[], DiffDoub0 Einit[], DesignVariable* dvAr[]);

		void getShellExpLoad(DiffDoub0 expLd[], DiffDoub0 E0Ld[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layQ[], DiffDoub0 layThExp[], DiffDoub0 layEInit[], DiffDoub0 layAng[]);

		void getBeamExpLoad(DiffDoub0 expLd[], DiffDoub0 E0Ld[], DesignVariable* dvAr[]);

		void getDensity(DiffDoub0& den, int layer, DesignVariable* dvAr[]);

		void getShellMass(DiffDoub0 Mmat[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layDen[], DesignVariable* dvAr[]);

		void getBeamMass(DiffDoub0 Mmat[], DesignVariable* dvAr[]);

		void getSolidDamp(DiffDoub0 Dmat[], DesignVariable* dvAr[]);

		void getShellDamp(DiffDoub0 Dmat[], DiffDoub0 layThk[], DiffDoub0 layZ[], DiffDoub0 layD[], DiffDoub0 layAng[]);

		void getBeamDamp(DiffDoub0 Dmat[], DesignVariable* dvAr[]);

		void getConductivity(DiffDoub0 tCond[], DesignVariable* dvAr[]);

		void getShellCond(DiffDoub0 tCond[], DiffDoub0 layThk[], DiffDoub0 layAng[], DiffDoub0 layCond[], DesignVariable* dvAr[]);

		void getBeamCond(DiffDoub0 tCond[], DesignVariable* dvAr[]);

		void getSpecificHeat(DiffDoub0& specHeat, DesignVariable* dvAr[]);

		void getShellSpecHeat(DiffDoub0& specHeat, DiffDoub0 layThk[], DiffDoub0 laySH[], DiffDoub0 layDen[]);

		void getBeamSpecHeat(DiffDoub0& specHeat, DesignVariable* dvAr[]);

        void getNdCrds(DiffDoub0 xGlob[], Node* ndAr[], DesignVariable* dvAr[]);
		
		void getLocOri(DiffDoub0 locOri[], DesignVariable* dvAr[]);
		
		void correctOrient(DiffDoub0 locOri[], DiffDoub0 xGlob[]);

		void getFrcFldConst(DiffDoub0 coef[], DiffDoub0 exp[], DesignVariable* dvAr[]);

		void getMassPerEl(DiffDoub0& massPerEl, DesignVariable* dvAr[]);

// Solution Fields
		void getNdDisp(DiffDoub0 globDisp[], Node* ndAr[]);

		void getNdVel(DiffDoub0 globVel[], Node* ndAr[]);

		void getNdAcc(DiffDoub0 globAcc[], Node* ndAr[]);

		void getNdTemp(DiffDoub0 globTemp[], Node* ndAr[]);

		void getNdTdot(DiffDoub0 globTdot[], Node* ndAr[]);

		void evalN(DiffDoub0 nVec[], DiffDoub0 dNds[], double spt[]);
		
		void getIpData(DiffDoub0 nVec[], DiffDoub0 dNdx[], DiffDoub0& detJ, DiffDoub0 locNds[], double spt[]);
		
		void getInstOri(DiffDoub0 instOriMat[], DiffDoub0 locOri[], DiffDoub0 globDisp[], int stat);
		
		void getInstDisp(DiffDoub0 instDisp[], DiffDoub0 globDisp[], DiffDoub0 instOriMat[], DiffDoub0 locOri[], DiffDoub0 xGlob[], bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getVolume(DiffDoub0& vol, DiffDoub0StressPrereq& pre, int layer);
		
		void getSectionDef(DiffDoub0 secDef[], DiffDoub0 globDisp[],  DiffDoub0 instOriMat[], DiffDoub0 locOri[], DiffDoub0 xGlob[], DiffDoub0 dNdx[], DiffDoub0 nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(DiffDoub0 strain[], DiffDoub0 ux[], DiffDoub0 dNdx[], DiffDoub0 locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub0 stress[], DiffDoub0 strain[], double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre);

		void dStressStraindU(DiffDoub0 dsdU[], DiffDoub0 dedU[], DiffDoub0 dsdT[], double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre);

		void getDefFrcMom(DiffDoub0 def[], DiffDoub0 frcMom[], double spt[], bool nLGeom, DiffDoub0StressPrereq& pre);

		void dDefFrcMomdU(DiffDoub0 dDefdU[], DiffDoub0 dFrcMomdU[], DiffDoub0 dFrcMomdT[], double spt[], bool nLGeom, DiffDoub0StressPrereq& pre);

		void getFluxTGrad(DiffDoub0 flux[], DiffDoub0 tGrad[], double spt[], int layer, DiffDoub0StressPrereq& pre);

		void dFluxTGraddT(DiffDoub0 dFdT[], DiffDoub0 dTG[], double spt[], int layer, DiffDoub0StressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, DiffDoub0 elQVec[], bool forTherm, int matRow, Node* ndAr[]);
		
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

// Properties
		void getLayerThkZ(DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1& zOffset, DesignVariable* dvAr[]);

		void getLayerQ(DiffDoub1 layQ[], DesignVariable* dvAr[]);

		void getLayerD(DiffDoub1 layD[], DesignVariable* dvAr[]);

		void getLayerAngle(DiffDoub1 layAng[], DesignVariable* dvAr[]);

		void getLayerThExp(DiffDoub1 layThExp[], DesignVariable* dvAr[]);

		void getLayerEinit(DiffDoub1 layEinit[], DesignVariable* dvAr[]);

		void getLayerDen(DiffDoub1 layerDen[], DesignVariable* dvAr[]);

		void getLayerCond(DiffDoub1 layCond[], DesignVariable* dvAr[]);

		void getLayerSpecHeat(DiffDoub1 laySH[], DesignVariable* dvAr[]);

		void transformStrain(DiffDoub1 stnNew[], DiffDoub1 stnOrig[], DiffDoub1& angle);

		void transformQ(DiffDoub1 qNew[], DiffDoub1 qOrig[], DiffDoub1& angle);

		void getSolidStiff(DiffDoub1 Cmat[], DesignVariable* dvAr[]);

		void getABD(DiffDoub1 Cmat[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layQ[], DiffDoub1 layAng[]);

		void getBeamStiff(DiffDoub1 Cmat[], DesignVariable* dvAr[]);

		void getThermalExp(DiffDoub1 thExp[], DiffDoub1 Einit[], DesignVariable* dvAr[]);

		void getShellExpLoad(DiffDoub1 expLd[], DiffDoub1 E0Ld[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layQ[], DiffDoub1 layThExp[], DiffDoub1 layEInit[], DiffDoub1 layAng[]);

		void getBeamExpLoad(DiffDoub1 expLd[], DiffDoub1 E0Ld[], DesignVariable* dvAr[]);

		void getDensity(DiffDoub1& den, int layer, DesignVariable* dvAr[]);

		void getShellMass(DiffDoub1 Mmat[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layDen[], DesignVariable* dvAr[]);

		void getBeamMass(DiffDoub1 Mmat[], DesignVariable* dvAr[]);

		void getSolidDamp(DiffDoub1 Dmat[], DesignVariable* dvAr[]);

		void getShellDamp(DiffDoub1 Dmat[], DiffDoub1 layThk[], DiffDoub1 layZ[], DiffDoub1 layD[], DiffDoub1 layAng[]);

		void getBeamDamp(DiffDoub1 Dmat[], DesignVariable* dvAr[]);

		void getConductivity(DiffDoub1 tCond[], DesignVariable* dvAr[]);

		void getShellCond(DiffDoub1 tCond[], DiffDoub1 layThk[], DiffDoub1 layAng[], DiffDoub1 layCond[], DesignVariable* dvAr[]);

		void getBeamCond(DiffDoub1 tCond[], DesignVariable* dvAr[]);

		void getSpecificHeat(DiffDoub1& specHeat, DesignVariable* dvAr[]);

		void getShellSpecHeat(DiffDoub1& specHeat, DiffDoub1 layThk[], DiffDoub1 laySH[], DiffDoub1 layDen[]);

		void getBeamSpecHeat(DiffDoub1& specHeat, DesignVariable* dvAr[]);

        void getNdCrds(DiffDoub1 xGlob[], Node* ndAr[], DesignVariable* dvAr[]);
		
		void getLocOri(DiffDoub1 locOri[], DesignVariable* dvAr[]);
		
		void correctOrient(DiffDoub1 locOri[], DiffDoub1 xGlob[]);

		void getFrcFldConst(DiffDoub1 coef[], DiffDoub1 exp[], DesignVariable* dvAr[]);

		void getMassPerEl(DiffDoub1& massPerEl, DesignVariable* dvAr[]);

// Solution Fields
		void getNdDisp(DiffDoub1 globDisp[], Node* ndAr[]);

		void getNdVel(DiffDoub1 globVel[], Node* ndAr[]);

		void getNdAcc(DiffDoub1 globAcc[], Node* ndAr[]);

		void getNdTemp(DiffDoub1 globTemp[], Node* ndAr[]);

		void getNdTdot(DiffDoub1 globTdot[], Node* ndAr[]);

		void evalN(DiffDoub1 nVec[], DiffDoub1 dNds[], double spt[]);
		
		void getIpData(DiffDoub1 nVec[], DiffDoub1 dNdx[], DiffDoub1& detJ, DiffDoub1 locNds[], double spt[]);
		
		void getInstOri(DiffDoub1 instOriMat[], DiffDoub1 locOri[], DiffDoub1 globDisp[], int stat);
		
		void getInstDisp(DiffDoub1 instDisp[], DiffDoub1 globDisp[], DiffDoub1 instOriMat[], DiffDoub1 locOri[], DiffDoub1 xGlob[], bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getVolume(DiffDoub1& vol, DiffDoub1StressPrereq& pre, int layer);
		
		void getSectionDef(DiffDoub1 secDef[], DiffDoub1 globDisp[],  DiffDoub1 instOriMat[], DiffDoub1 locOri[], DiffDoub1 xGlob[], DiffDoub1 dNdx[], DiffDoub1 nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(DiffDoub1 strain[], DiffDoub1 ux[], DiffDoub1 dNdx[], DiffDoub1 locOri[], int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub1 stress[], DiffDoub1 strain[], double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre);

		void dStressStraindU(DiffDoub1 dsdU[], DiffDoub1 dedU[], DiffDoub1 dsdT[], double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre);

		void getDefFrcMom(DiffDoub1 def[], DiffDoub1 frcMom[], double spt[], bool nLGeom, DiffDoub1StressPrereq& pre);

		void dDefFrcMomdU(DiffDoub1 dDefdU[], DiffDoub1 dFrcMomdU[], DiffDoub1 dFrcMomdT[], double spt[], bool nLGeom, DiffDoub1StressPrereq& pre);

		void getFluxTGrad(DiffDoub1 flux[], DiffDoub1 tGrad[], double spt[], int layer, DiffDoub1StressPrereq& pre);

		void dFluxTGraddT(DiffDoub1 dFdT[], DiffDoub1 dTG[], double spt[], int layer, DiffDoub1StressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, DiffDoub1 elQVec[], bool forTherm, int matRow, Node* ndAr[]);
		
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
        void getRuk(DiffDoub0 Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRum(DiffDoub0 Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRud(DiffDoub0 Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRu(DiffDoub0 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRtk(DiffDoub0 Rvec[], double dRdT[], bool getMatrix, DiffDoub0StressPrereq& pre);

		void getRtm(DiffDoub0 Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DiffDoub0StressPrereq& pre);

		void getRt(DiffDoub0 globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[]);

		void getRuFrcFld(DiffDoub0 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub0StressPrereq& pre, Node* ndAr[]);

		void getAppLoad(DiffDoub0 AppLd[], Load* ldPt, bool nLGeom, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getAppThermLoad(DiffDoub0 AppLd[], Load* ldPt, DiffDoub0StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
        void getRuk(DiffDoub1 Rvec[], double dRdu[], double dRdT[], bool getMatrix, bool nLGeom, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRum(DiffDoub1 Rvec[], double dRdA[], bool getMatrix, bool actualProps, bool nLGeom, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRud(DiffDoub1 Rvec[], double dRdV[], bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getRu(DiffDoub1 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);
		
		void getRtk(DiffDoub1 Rvec[], double dRdT[], bool getMatrix, DiffDoub1StressPrereq& pre);

		void getRtm(DiffDoub1 Rvec[], double dRdTdot[], bool getMatrix, bool actualProps, DiffDoub1StressPrereq& pre);

		void getRt(DiffDoub1 globR[], SparseMat& globdRdT, bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[]);

		void getRuFrcFld(DiffDoub1 globR[], SparseMat& globdRdu, bool getMatrix, JobCommand* cmd, DiffDoub1StressPrereq& pre, Node* ndAr[]);

		void getAppLoad(DiffDoub1 AppLd[], Load* ldPt, bool nLGeom, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

		void getAppThermLoad(DiffDoub1 AppLd[], Load* ldPt, DiffDoub1StressPrereq& pre, Node* ndAr[], DesignVariable* dvAr[]);

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

		~ElementList();
};


#endif