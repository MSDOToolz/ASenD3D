#ifndef ELEMENT
#define ELEMENT
#include <vector>
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
	std::vector<DiffDoub0> globNds;  //globNds[30];
	std::vector<DiffDoub0> locNds;// [30] ;
	std::vector<DiffDoub0> locOri;// [9] ;
	std::vector<DiffDoub0> instOri;// [720] ;
	std::vector<DiffDoub0> globDisp;// [60] ;
	std::vector<DiffDoub0> globVel;
	std::vector<DiffDoub0> globAcc;// [30];
	std::vector<DiffDoub0> globTemp;
	std::vector<DiffDoub0> globTdot;
	std::vector<DiffDoub0> Cmat;// [81] ;
	std::vector<DiffDoub0> Mmat;// [36];
	std::vector<DiffDoub0> Dmat;
	std::vector<DiffDoub0> thermExp;
	std::vector<DiffDoub0> Einit;
	std::vector<DiffDoub0> TCmat;
	DiffDoub0 SpecHeat;
	std::vector<DiffDoub0> BMat;
	std::vector<DiffDoub0> CBMat;
	std::vector<DiffDoub0> layerZ;
	std::vector<DiffDoub0> layerThk;
	std::vector<DiffDoub0> layerAng;
	std::vector<DiffDoub0> layerQ;
	std::vector<DiffDoub0> layerD;
	std::vector<DiffDoub0> layerTE;
	std::vector<DiffDoub0> layerE0;
	std::vector<DiffDoub0> layerDen;
	std::vector<DiffDoub0> layerTC;
	std::vector<DiffDoub0> layerSH;
	std::vector<DiffDoub0> frcFldCoef;
	std::vector<DiffDoub0> frcFldExp;
	std::vector<DiffDoub0> thrmFldCoef;
	DiffDoub0 refTemp;
	DiffDoub0 massPerEl;
	std::vector<double> scrMat1;
	std::vector<double> scrMat2;
	std::vector<double> scrMat3;
	std::vector<double> scrMat4;
	std::vector<double> scrMat5;
	std::vector<DiffDoub0> scrVec1;
	std::vector<DiffDoub0> scrVec2;
	std::vector<DiffDoub0> scrVec3;
	std::vector<DiffDoub0> scrVec4;
	std::vector<DiffDoub0> scrVec5;

	int currentLayLen;

	DiffDoub0StressPrereq();

	void allocateLayers(int numLayers);

};

class DiffDoub0FlPrereq {
public:
	std::vector<DiffDoub0> globNds;
	std::vector<DiffDoub0> globDisp;
	std::vector<DiffDoub0> globVel;
	std::vector<DiffDoub0> flDen;
	std::vector<DiffDoub0> flVel;
	std::vector<DiffDoub0> flTemp;
	std::vector<DiffDoub0> flTurbE;
	std::vector<DiffDoub0> flDenDot;
	std::vector<DiffDoub0> flVelDot;
	std::vector<DiffDoub0> flTDot;
	std::vector<DiffDoub0> flTurbEDot;
	DiffDoub0 refTurbE;
	DiffDoub0 gradVTurbCoef;
	DiffDoub0 dissTurbCoef;
	DiffDoub0 refVisc;
	DiffDoub0 denVisCoef;
	DiffDoub0 tempVisCoef;
	DiffDoub0 turbVisCoef;
	DiffDoub0 refEnth;
	DiffDoub0 denEnthCoef;
	DiffDoub0 denEnthExp;
	DiffDoub0 denPresCoef;
	DiffDoub0 denPresExp;
	DiffDoub0 refDen;
	DiffDoub0 refTemp;
	DiffDoub0 thermCond;
	DiffDoub0 specHeat;
	DiffDoub0 iGConst;
	std::vector<double> scratch;

	DiffDoub0FlPrereq();

};

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
class DiffDoub1StressPrereq {
public:
	std::vector<DiffDoub1> globNds;  //globNds[30];
	std::vector<DiffDoub1> locNds;// [30] ;
	std::vector<DiffDoub1> locOri;// [9] ;
	std::vector<DiffDoub1> instOri;// [720] ;
	std::vector<DiffDoub1> globDisp;// [60] ;
	std::vector<DiffDoub1> globVel;
	std::vector<DiffDoub1> globAcc;// [30];
	std::vector<DiffDoub1> globTemp;
	std::vector<DiffDoub1> globTdot;
	std::vector<DiffDoub1> Cmat;// [81] ;
	std::vector<DiffDoub1> Mmat;// [36];
	std::vector<DiffDoub1> Dmat;
	std::vector<DiffDoub1> thermExp;
	std::vector<DiffDoub1> Einit;
	std::vector<DiffDoub1> TCmat;
	DiffDoub1 SpecHeat;
	std::vector<DiffDoub1> BMat;
	std::vector<DiffDoub1> CBMat;
	std::vector<DiffDoub1> layerZ;
	std::vector<DiffDoub1> layerThk;
	std::vector<DiffDoub1> layerAng;
	std::vector<DiffDoub1> layerQ;
	std::vector<DiffDoub1> layerD;
	std::vector<DiffDoub1> layerTE;
	std::vector<DiffDoub1> layerE0;
	std::vector<DiffDoub1> layerDen;
	std::vector<DiffDoub1> layerTC;
	std::vector<DiffDoub1> layerSH;
	std::vector<DiffDoub1> frcFldCoef;
	std::vector<DiffDoub1> frcFldExp;
	std::vector<DiffDoub1> thrmFldCoef;
	DiffDoub1 refTemp;
	DiffDoub1 massPerEl;
	std::vector<double> scrMat1;
	std::vector<double> scrMat2;
	std::vector<double> scrMat3;
	std::vector<double> scrMat4;
	std::vector<double> scrMat5;
	std::vector<DiffDoub1> scrVec1;
	std::vector<DiffDoub1> scrVec2;
	std::vector<DiffDoub1> scrVec3;
	std::vector<DiffDoub1> scrVec4;
	std::vector<DiffDoub1> scrVec5;

	int currentLayLen;

	DiffDoub1StressPrereq();

	void allocateLayers(int numLayers);

};

class DiffDoub1FlPrereq {
public:
	std::vector<DiffDoub1> globNds;
	std::vector<DiffDoub1> globDisp;
	std::vector<DiffDoub1> globVel;
	std::vector<DiffDoub1> flDen;
	std::vector<DiffDoub1> flVel;
	std::vector<DiffDoub1> flTemp;
	std::vector<DiffDoub1> flTurbE;
	std::vector<DiffDoub1> flDenDot;
	std::vector<DiffDoub1> flVelDot;
	std::vector<DiffDoub1> flTDot;
	std::vector<DiffDoub1> flTurbEDot;
	DiffDoub1 refTurbE;
	DiffDoub1 gradVTurbCoef;
	DiffDoub1 dissTurbCoef;
	DiffDoub1 refVisc;
	DiffDoub1 denVisCoef;
	DiffDoub1 tempVisCoef;
	DiffDoub1 turbVisCoef;
	DiffDoub1 refEnth;
	DiffDoub1 denEnthCoef;
	DiffDoub1 denEnthExp;
	DiffDoub1 denPresCoef;
	DiffDoub1 denPresExp;
	DiffDoub1 refDen;
	DiffDoub1 refTemp;
	DiffDoub1 thermCond;
	DiffDoub1 specHeat;
	DiffDoub1 iGConst;
	std::vector<double> scratch;

	DiffDoub1FlPrereq();

};

//end dup
 
//end skip 
 
 
 
 
 
 
class Element {
	public:
	    int type;
		int label;
	    std::vector<int> nodes;
		int numNds;
		int dofPerNd;
		int numIntDof;
		int nDim;
		int defDim;
		std::vector<int> dofTable;  //dofTable[2*i] = nd/basis function of dof i, dofTable[2*i+1] = component of dof i
		int intDofIndex;
		std::vector<double> intPts;
		std::vector<double> ipWt;
		std::vector<double> ndSpts;
		double sCent[3];
		int numIP;
		int numFaces;
		std::list<int> faces;
		std::vector<double> internalDisp;
		std::vector<double> intPrevDisp;
		std::vector<double> internaldLdu;
		std::vector<double> internalAdj;
		std::vector<DiffDoub1> internalRu;
		std::vector<double> internalMat;
		std::list<IDCapsule> designVars;
		std::list<int> compDVars;
		int sectPtr;
		
		Element();

	    void initializeType(int newType);

        void setNodes(int newNds[]);

        void initializeFaces(std::vector<Face>& globFcLst, int& fi);	

		void setIntDisp(double newDisp[]);

		void setIntPrevDisp(double newDisp[]);

		void advanceIntDisp();

		void backstepIntDisp();

		void setIntdLdU(std::vector<double>& globdLdU);

		int getNumLayers(std::vector<Section>& secLst);
		
		void addDesignVariable(int dIndex, double coef);

		void addCompDVar(int dIndex);

//dup1

// Properties
		void getGenProp(DiffDoub0& prop, std::string propKey, std::vector<DesignVariable>& dvAr);

		void getLayerThkZ(std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& layZ, DiffDoub0& zOffset, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getLayerQ(std::vector<DiffDoub0>& layQ, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerD(std::vector<DiffDoub0>& layD, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerAngle(std::vector<DiffDoub0>& layAng, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getLayerThExp(std::vector<DiffDoub0>& layThExp, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerEinit(std::vector<DiffDoub0>& layEinit, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getLayerDen(std::vector<DiffDoub0>& layerDen, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerCond(std::vector<DiffDoub0>& layCond, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerSpecHeat(std::vector<DiffDoub0>& laySH, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void transformStrain(DiffDoub0 stnNew[], DiffDoub0 stnOrig[], DiffDoub0& angle);

		void transformQ(DiffDoub0 qNew[], DiffDoub0 qOrig[], DiffDoub0& angle);

		void getSolidStiff(std::vector<DiffDoub0>& Cmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getABD(std::vector<DiffDoub0>& Cmat, std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& layZ, std::vector<DiffDoub0>& layQ, std::vector<DiffDoub0>& layAng, std::vector<Section>& secAr);

		void getBeamStiff(std::vector<DiffDoub0>& Cmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getThermalExp(std::vector<DiffDoub0>& thExp, std::vector<DiffDoub0>& Einit, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellExpLoad(std::vector<DiffDoub0>& expLd, std::vector<DiffDoub0>& E0Ld, std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& layZ, std::vector<DiffDoub0>& layQ, std::vector<DiffDoub0>& layThExp, std::vector<DiffDoub0>& layEInit, std::vector<DiffDoub0>& layAng, std::vector<Section>& secAr);

		void getBeamExpLoad(std::vector<DiffDoub0>& expLd, std::vector<DiffDoub0>& E0Ld, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getDensity(DiffDoub0& den, int layer, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellMass(std::vector<DiffDoub0>& Mmat, std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& layZ, std::vector<DiffDoub0>& layDen, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getBeamMass(std::vector<DiffDoub0>& Mmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getSolidDamp(std::vector<DiffDoub0>& Dmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellDamp(std::vector<DiffDoub0>& Dmat, std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& layZ, std::vector<DiffDoub0>& layD, std::vector<DiffDoub0>& layAng, std::vector<Section>& secAr);

		void getBeamDamp(std::vector<DiffDoub0>& Dmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getConductivity(std::vector<DiffDoub0>& tCond, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellCond(std::vector<DiffDoub0>& tCond, std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& layAng, std::vector<DiffDoub0>& layCond, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getBeamCond(std::vector<DiffDoub0>& tCond, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getSpecificHeat(DiffDoub0& specHeat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellSpecHeat(DiffDoub0& specHeat, std::vector<DiffDoub0>& layThk, std::vector<DiffDoub0>& laySH, std::vector<DiffDoub0>& layDen, std::vector<Section>& secAr);

		void getBeamSpecHeat(DiffDoub0& specHeat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

        void getNdCrds(std::vector<DiffDoub0>& xGlob, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);
		
		void getLocOri(std::vector<DiffDoub0>& locOri, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);
		
		void correctOrient(std::vector<DiffDoub0>& locOri, std::vector<DiffDoub0>& xGlob);

		void getFrcFldConst(std::vector<DiffDoub0>& coef, std::vector<DiffDoub0>& exp, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getThrmFldConst(std::vector<DiffDoub0>& coef, DiffDoub0& refT, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getMassPerEl(DiffDoub0& massPerEl, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

// Solution Fields
		void getNdDisp(std::vector<DiffDoub0>& globDisp, std::vector<Node>& ndAr);

		void getNdVel(std::vector<DiffDoub0>& globVel, std::vector<Node>& ndAr);

		void getNdAcc(std::vector<DiffDoub0>& globAcc, std::vector<Node>& ndAr);

		void getNdFlVel(std::vector<DiffDoub0>& flVel, std::vector<Node>& ndAr);

		void getNdFlVDot(std::vector<DiffDoub0>& flVDot, std::vector<Node>& ndAr);

		void getNdTemp(std::vector<DiffDoub0>& globTemp, std::vector<Node>& ndAr);

		void getNdTdot(std::vector<DiffDoub0>& globTdot, std::vector<Node>& ndAr);

		void getNdFlDen(std::vector<DiffDoub0>& flDen, std::vector<Node>& ndAr);

		void getNdFlDenDot(std::vector<DiffDoub0>& flDenDot, std::vector<Node>& ndAr);

		void getNdTurbE(std::vector<DiffDoub0>& turbE, std::vector<Node>& ndAr);

		void getNdTurbEDot(std::vector<DiffDoub0>& turbEDot, std::vector<Node>& ndAr);

		void evalN(DiffDoub0 nVec[], DiffDoub0 dNds[], double spt[]);
		
		void getIpData(DiffDoub0 nVec[], DiffDoub0 dNdx[], DiffDoub0& detJ, std::vector<DiffDoub0>& locNds, double spt[]);
		
		void getInstOri(std::vector<DiffDoub0>& instOriMat, std::vector<DiffDoub0>& locOri, std::vector<DiffDoub0>& globDisp, int stat);
		
		void getInstDisp(DiffDoub0 instDisp[], std::vector<DiffDoub0>& globDisp, std::vector<DiffDoub0>& instOriMat, std::vector<DiffDoub0>& locOri, std::vector<DiffDoub0>& xGlob, bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DiffDoub0StressPrereq& pre, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getFluidPrereq(DiffDoub0FlPrereq& pre, std::vector<Section>& secAr, std::vector<Fluid>& flAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getVolume(DiffDoub0& vol, DiffDoub0StressPrereq& pre, int layer, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);
		
		void getSectionDef(DiffDoub0 secDef[], std::vector<DiffDoub0>& globDisp, std::vector<DiffDoub0>& instOriMat, std::vector<DiffDoub0>& locOri, std::vector<DiffDoub0>& xGlob, DiffDoub0 dNdx[], DiffDoub0 nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(DiffDoub0 strain[], DiffDoub0 ux[], DiffDoub0 dNdx[], std::vector<DiffDoub0>& locOri, int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub0 stress[], DiffDoub0 strain[], double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre);

		void dStressStraindU(std::vector<DiffDoub0>& dsdU, std::vector<DiffDoub0>& dedU, std::vector<DiffDoub0>& dsdT, double spt[], int layer, bool nLGeom, DiffDoub0StressPrereq& pre);

		void getDefFrcMom(DiffDoub0 def[], DiffDoub0 frcMom[], double spt[], bool nLGeom, DiffDoub0StressPrereq& pre);

		void dDefFrcMomdU(std::vector<DiffDoub0>& dDefdU, std::vector<DiffDoub0>& dFrcMomdU, std::vector<DiffDoub0>& dFrcMomdT, double spt[], bool nLGeom, DiffDoub0StressPrereq& pre);

		void getFluxTGrad(DiffDoub0 flux[], DiffDoub0 tGrad[], double spt[], int layer, DiffDoub0StressPrereq& pre);

		void dFluxTGraddT(std::vector<DiffDoub0>& dFdT, std::vector<DiffDoub0>& dTG, double spt[], int layer, DiffDoub0StressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, std::vector<DiffDoub0>& elQVec, bool forTherm, int matRow, std::vector<Node>& ndAr);
		
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

// Properties
		void getGenProp(DiffDoub1& prop, std::string propKey, std::vector<DesignVariable>& dvAr);

		void getLayerThkZ(std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& layZ, DiffDoub1& zOffset, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getLayerQ(std::vector<DiffDoub1>& layQ, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerD(std::vector<DiffDoub1>& layD, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerAngle(std::vector<DiffDoub1>& layAng, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getLayerThExp(std::vector<DiffDoub1>& layThExp, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerEinit(std::vector<DiffDoub1>& layEinit, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getLayerDen(std::vector<DiffDoub1>& layerDen, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerCond(std::vector<DiffDoub1>& layCond, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getLayerSpecHeat(std::vector<DiffDoub1>& laySH, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void transformStrain(DiffDoub1 stnNew[], DiffDoub1 stnOrig[], DiffDoub1& angle);

		void transformQ(DiffDoub1 qNew[], DiffDoub1 qOrig[], DiffDoub1& angle);

		void getSolidStiff(std::vector<DiffDoub1>& Cmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getABD(std::vector<DiffDoub1>& Cmat, std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& layZ, std::vector<DiffDoub1>& layQ, std::vector<DiffDoub1>& layAng, std::vector<Section>& secAr);

		void getBeamStiff(std::vector<DiffDoub1>& Cmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getThermalExp(std::vector<DiffDoub1>& thExp, std::vector<DiffDoub1>& Einit, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellExpLoad(std::vector<DiffDoub1>& expLd, std::vector<DiffDoub1>& E0Ld, std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& layZ, std::vector<DiffDoub1>& layQ, std::vector<DiffDoub1>& layThExp, std::vector<DiffDoub1>& layEInit, std::vector<DiffDoub1>& layAng, std::vector<Section>& secAr);

		void getBeamExpLoad(std::vector<DiffDoub1>& expLd, std::vector<DiffDoub1>& E0Ld, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getDensity(DiffDoub1& den, int layer, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellMass(std::vector<DiffDoub1>& Mmat, std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& layZ, std::vector<DiffDoub1>& layDen, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getBeamMass(std::vector<DiffDoub1>& Mmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getSolidDamp(std::vector<DiffDoub1>& Dmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellDamp(std::vector<DiffDoub1>& Dmat, std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& layZ, std::vector<DiffDoub1>& layD, std::vector<DiffDoub1>& layAng, std::vector<Section>& secAr);

		void getBeamDamp(std::vector<DiffDoub1>& Dmat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getConductivity(std::vector<DiffDoub1>& tCond, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellCond(std::vector<DiffDoub1>& tCond, std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& layAng, std::vector<DiffDoub1>& layCond, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getBeamCond(std::vector<DiffDoub1>& tCond, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getSpecificHeat(DiffDoub1& specHeat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

		void getShellSpecHeat(DiffDoub1& specHeat, std::vector<DiffDoub1>& layThk, std::vector<DiffDoub1>& laySH, std::vector<DiffDoub1>& layDen, std::vector<Section>& secAr);

		void getBeamSpecHeat(DiffDoub1& specHeat, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<DesignVariable>& dvAr);

        void getNdCrds(std::vector<DiffDoub1>& xGlob, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);
		
		void getLocOri(std::vector<DiffDoub1>& locOri, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);
		
		void correctOrient(std::vector<DiffDoub1>& locOri, std::vector<DiffDoub1>& xGlob);

		void getFrcFldConst(std::vector<DiffDoub1>& coef, std::vector<DiffDoub1>& exp, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getThrmFldConst(std::vector<DiffDoub1>& coef, DiffDoub1& refT, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

		void getMassPerEl(DiffDoub1& massPerEl, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);

// Solution Fields
		void getNdDisp(std::vector<DiffDoub1>& globDisp, std::vector<Node>& ndAr);

		void getNdVel(std::vector<DiffDoub1>& globVel, std::vector<Node>& ndAr);

		void getNdAcc(std::vector<DiffDoub1>& globAcc, std::vector<Node>& ndAr);

		void getNdFlVel(std::vector<DiffDoub1>& flVel, std::vector<Node>& ndAr);

		void getNdFlVDot(std::vector<DiffDoub1>& flVDot, std::vector<Node>& ndAr);

		void getNdTemp(std::vector<DiffDoub1>& globTemp, std::vector<Node>& ndAr);

		void getNdTdot(std::vector<DiffDoub1>& globTdot, std::vector<Node>& ndAr);

		void getNdFlDen(std::vector<DiffDoub1>& flDen, std::vector<Node>& ndAr);

		void getNdFlDenDot(std::vector<DiffDoub1>& flDenDot, std::vector<Node>& ndAr);

		void getNdTurbE(std::vector<DiffDoub1>& turbE, std::vector<Node>& ndAr);

		void getNdTurbEDot(std::vector<DiffDoub1>& turbEDot, std::vector<Node>& ndAr);

		void evalN(DiffDoub1 nVec[], DiffDoub1 dNds[], double spt[]);
		
		void getIpData(DiffDoub1 nVec[], DiffDoub1 dNdx[], DiffDoub1& detJ, std::vector<DiffDoub1>& locNds, double spt[]);
		
		void getInstOri(std::vector<DiffDoub1>& instOriMat, std::vector<DiffDoub1>& locOri, std::vector<DiffDoub1>& globDisp, int stat);
		
		void getInstDisp(DiffDoub1 instDisp[], std::vector<DiffDoub1>& globDisp, std::vector<DiffDoub1>& instOriMat, std::vector<DiffDoub1>& locOri, std::vector<DiffDoub1>& xGlob, bool nLGeom, int dv1, int dv2);

		void getStressPrereq(DiffDoub1StressPrereq& pre, std::vector<Section>& secAr, std::vector<Material>& matAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getFluidPrereq(DiffDoub1FlPrereq& pre, std::vector<Section>& secAr, std::vector<Fluid>& flAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getVolume(DiffDoub1& vol, DiffDoub1StressPrereq& pre, int layer, std::vector<Section>& secAr, std::vector<DesignVariable>& dvAr);
		
		void getSectionDef(DiffDoub1 secDef[], std::vector<DiffDoub1>& globDisp, std::vector<DiffDoub1>& instOriMat, std::vector<DiffDoub1>& locOri, std::vector<DiffDoub1>& xGlob, DiffDoub1 dNdx[], DiffDoub1 nVec[], bool nLGeom, int dv1, int dv2);
		
		void getSolidStrain(DiffDoub1 strain[], DiffDoub1 ux[], DiffDoub1 dNdx[], std::vector<DiffDoub1>& locOri, int dv1, int dv2, bool nLGeom);

		void getStressStrain(DiffDoub1 stress[], DiffDoub1 strain[], double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre);

		void dStressStraindU(std::vector<DiffDoub1>& dsdU, std::vector<DiffDoub1>& dedU, std::vector<DiffDoub1>& dsdT, double spt[], int layer, bool nLGeom, DiffDoub1StressPrereq& pre);

		void getDefFrcMom(DiffDoub1 def[], DiffDoub1 frcMom[], double spt[], bool nLGeom, DiffDoub1StressPrereq& pre);

		void dDefFrcMomdU(std::vector<DiffDoub1>& dDefdU, std::vector<DiffDoub1>& dFrcMomdU, std::vector<DiffDoub1>& dFrcMomdT, double spt[], bool nLGeom, DiffDoub1StressPrereq& pre);

		void getFluxTGrad(DiffDoub1 flux[], DiffDoub1 tGrad[], double spt[], int layer, DiffDoub1StressPrereq& pre);

		void dFluxTGraddT(std::vector<DiffDoub1>& dFdT, std::vector<DiffDoub1>& dTG, double spt[], int layer, DiffDoub1StressPrereq& pre);

		void putVecToGlobMat(SparseMat& qMat, std::vector<DiffDoub1>& elQVec, bool forTherm, int matRow, std::vector<Node>& ndAr);
		
//end dup
 
//end skip 
 
 
 
 
 
 
		void getElVec(std::vector<double>& elVec, std::vector<double>& globVec, bool forTherm, bool intnl, std::vector<Node>& ndAr);

		void addToGlobVec(std::vector<double>& elVec, std::vector<double>& globVec, bool forTherm, bool intnl, std::vector<Node>& ndAr);
 
 
// Equations
		void condenseMat(std::vector<double>& mat, std::vector<double>& scr1, std::vector<double>& scr2);
		
		void updateExternal(std::vector<double>& extVec, int forSoln, std::vector<Node>& ndAr, std::vector<double>& scr1, std::vector<double>& scr2);
		
		void updateInternal(std::vector<double>& extVec, int forSoln, std::vector<Node>& ndAr, std::vector<double>& scr1, std::vector<double>& scr2);

		double getIntAdjdRdD();

//dup1
        void getRuk(std::vector<DiffDoub0>& Rvec, std::vector<double>& dRdu, std::vector<double>& dRdT, bool getMatrix, bool nLGeom, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getRum(std::vector<DiffDoub0>& Rvec, std::vector<double>& dRdA, bool getMatrix, bool actualProps, bool nLGeom, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);
		
		void getRud(std::vector<DiffDoub0>& Rvec, std::vector<double>& dRdV, bool getMatrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getRu(std::vector<DiffDoub0>& globR, SparseMat& globdRdu, bool getMatrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);
		
		void getRtk(std::vector<DiffDoub0>& Rvec, std::vector<double>& dRdT, bool getMatrix, DiffDoub0StressPrereq& pre);

		void getRtm(std::vector<DiffDoub0>& Rvec, std::vector<double>& dRdTdot, bool getMatrix, bool actualProps, DiffDoub0StressPrereq& pre);

		void getRt(std::vector<DiffDoub0>& globR, SparseMat& globdRdT, bool getMatrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr);

		void getRuFrcFld(std::vector<DiffDoub0>& globR, SparseMat& globdRdu, bool getMatrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr);

		void getRtFrcFld(std::vector<DiffDoub0>& globR, SparseMat& globdRdT, std::vector<double>& dRdU, std::vector<double>& dRdV, bool getMatrix, JobCommand& cmd, DiffDoub0StressPrereq& pre, std::vector<Node>& ndAr);

		void getAppLoad(std::vector<DiffDoub0>& AppLd, Load& ldPt, bool nLGeom, DiffDoub0StressPrereq& pre, std::vector<Section>& secAr, std::vector<Face>& fcAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getAppThermLoad(std::vector<DiffDoub0>& AppLd, Load& ldPt, DiffDoub0StressPrereq& pre, std::vector<Section>& secAr, std::vector<Face>& fcAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getRfUnst(std::vector<DiffDoub0>& Rvec, DiffDoub0FlPrereq& pre);

		void getRfConv(std::vector<DiffDoub0>& Rvec, DiffDoub0FlPrereq& pre);

		void getRfVisc(std::vector<DiffDoub0>& Rvec, DiffDoub0FlPrereq& pre);

		void getRfPres(std::vector<DiffDoub0>& Rvec, DiffDoub0FlPrereq& pre);

		void getRfLoad(std::vector<DiffDoub0>& Rvec, DiffDoub0FlPrereq& pre, std::vector<DiffDoub0>& bodyF, DiffDoub0& bodyQ);

		void getRf(std::vector<DiffDoub0>& globR, SparseMat& globdRdV, bool getMatrix, JobCommand& cmd, DiffDoub0FlPrereq& pre);

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
        void getRuk(std::vector<DiffDoub1>& Rvec, std::vector<double>& dRdu, std::vector<double>& dRdT, bool getMatrix, bool nLGeom, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getRum(std::vector<DiffDoub1>& Rvec, std::vector<double>& dRdA, bool getMatrix, bool actualProps, bool nLGeom, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);
		
		void getRud(std::vector<DiffDoub1>& Rvec, std::vector<double>& dRdV, bool getMatrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getRu(std::vector<DiffDoub1>& globR, SparseMat& globdRdu, bool getMatrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);
		
		void getRtk(std::vector<DiffDoub1>& Rvec, std::vector<double>& dRdT, bool getMatrix, DiffDoub1StressPrereq& pre);

		void getRtm(std::vector<DiffDoub1>& Rvec, std::vector<double>& dRdTdot, bool getMatrix, bool actualProps, DiffDoub1StressPrereq& pre);

		void getRt(std::vector<DiffDoub1>& globR, SparseMat& globdRdT, bool getMatrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr);

		void getRuFrcFld(std::vector<DiffDoub1>& globR, SparseMat& globdRdu, bool getMatrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr);

		void getRtFrcFld(std::vector<DiffDoub1>& globR, SparseMat& globdRdT, std::vector<double>& dRdU, std::vector<double>& dRdV, bool getMatrix, JobCommand& cmd, DiffDoub1StressPrereq& pre, std::vector<Node>& ndAr);

		void getAppLoad(std::vector<DiffDoub1>& AppLd, Load& ldPt, bool nLGeom, DiffDoub1StressPrereq& pre, std::vector<Section>& secAr, std::vector<Face>& fcAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getAppThermLoad(std::vector<DiffDoub1>& AppLd, Load& ldPt, DiffDoub1StressPrereq& pre, std::vector<Section>& secAr, std::vector<Face>& fcAr, std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		void getRfUnst(std::vector<DiffDoub1>& Rvec, DiffDoub1FlPrereq& pre);

		void getRfConv(std::vector<DiffDoub1>& Rvec, DiffDoub1FlPrereq& pre);

		void getRfVisc(std::vector<DiffDoub1>& Rvec, DiffDoub1FlPrereq& pre);

		void getRfPres(std::vector<DiffDoub1>& Rvec, DiffDoub1FlPrereq& pre);

		void getRfLoad(std::vector<DiffDoub1>& Rvec, DiffDoub1FlPrereq& pre, std::vector<DiffDoub1>& bodyF, DiffDoub1& bodyQ);

		void getRf(std::vector<DiffDoub1>& globR, SparseMat& globdRdV, bool getMatrix, JobCommand& cmd, DiffDoub1FlPrereq& pre);

//end dup
 
//end skip 
 
 
 
 
 
 
};


#endif