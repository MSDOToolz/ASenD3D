#ifndef MODEL
#define MODEL
#include <fstream>
#include <iostream>
#include <map>
#include "JobClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "FaceClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "LoadClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"
#include "DiffDoubClass.h"

class Model {
	private:
	    NodeList nodes;
		Node** nodeArray;
		ElementList elements;
		Element** elementArray;
		SetList nodeSets;
		Set** nsArray;
		std::map<std::string, int> nsMap;
		SetList elementSets;
		Set** esArray;
		std::map<std::string, int> esMap;
		SectionList sections;
		MaterialList materials;
		ConstraintList elasticConst;
		ConstraintList thermalConst;
		LoadList elasticLoads;
		LoadList thermalLoads;
		DesVarList designVars;
		DesignVariable** dVarArray;
		Objective objective;
		Job job;
		
		int elMatDim;
		int totGlobDof;
		bool anPrepRun;
		int timeStepsSaved;
		JobCommand* solveCmd;
		JobCommand* modalCmd;
		DiffDoub0StressPrereq d0Pre;
		DiffDoub1StressPrereq d1Pre;
		
		SparseMat elasticMat;
		LowerTriMat elasticLT;
		bool elasticScaled;
		SparseMat thermMat;
		LowerTriMat thermLT;
		bool thermScaled;

		double* eigVecs;
		double* eigVals;
		double* diagMass;
		double* loadFact;
		
		double* tempV1;
		double* tempV2;
		double* tempV3;
		double* tempV4;
		double* tempV5;
		double* tempV6;

		DiffDoub0 *tempD1;

		double* dLdU;
		double* dLdV;
		double* dLdA;
		double* dLdT;
		double* dLdTdot;
		double* uAdj;
		double* vAdj;
		double* aAdj;
		double* tAdj;
		double* tdotAdj;

		DiffDoub1* dRudD;
		DiffDoub1* dRtdD;

		int* elInD;

		double* dLdD;
	
	public:
	    Model();
		
		NodeList* getNodes();
		
		ElementList* getElements();
		
		SetList* getNodeSets();
		
		SetList* getElementSets();
		
		SectionList* getSections();
		
		MaterialList* getMaterials();
		
		ConstraintList* getElasticConstraints();
		
		DesVarList* getDesignVars();
		
		void executeJob();

        ~Model();	// Finish	
		
		// Input
		
		void readInputLine(std::ifstream& inFile,std::string& fileLine,std::string headings[],int hdLdSpace[], std::string data[], int& dataLen);
		
		void readJob(std::string fileName);
		
		void readModelInput(std::string fileName);
		
		void readConstraintInput(std::string fileName);
		
		void readLoadInput(std::string fileName);
		
		void readInitialState(std::string fileName);
		
		void readDesVarInput(std::string fileName);
		
		void readObjectiveInput(std::string fileName);
		
		void readDesVarValues(std::string fileName);

		void readNodeResults(std::string fileName);

		void readTimeStepSoln(int tStep);
		
		// Analysis
		
		void reorderNodes(int blockDim);

		void buildConstraintMats();
		
		void updateReference();
		
		void findSurfaceFaces();

		void prepMatrixFactorizations();
		
		void analysisPrep();
		
		void buildElasticAppLoad(double appLd[], double time);

		void buildThermalAppLoad(double appLd[], double time);
		
		void buildElasticSolnLoad(double solnLd[], bool buildMat);

		void buildThermalSolnLoad(double solnLd[], bool buildMat);

		void scaleElasticConst();

		void scaleThermalConst();

		void buildElasticConstLoad(double constLd[]);

		void buildThermalConstLoad(double constLd[]);
		
		void solveStep(JobCommand *cmd, double time, double appLdFact);
		
		void solve(JobCommand *cmd);

		void zeroSolution(StringList& fields);
		
		void eigenSolve(JobCommand* cmd);

	/*	void backupElastic();

		void restoreElastic();*/

		void setSolnToMode(std::string field, int mode, double maxVal);

		void augmentdLdU();

		void solveForAdjoint(double time);

		void dRthermaldD(int dVarNum);

		void dRelasticdD(int dVarNum);

		void getObjective();
		
		void getObjGradient();
		
		// Output

		void writeTimeStepSoln(int tStep);
		
		void writeNodeResults(std::string fileName, std::string nodeSet, StringList& fields, int timeStep);

		void writeElementResults(std::string fileName, std::string elSet, StringList& fields, std::string position, int timeStep);

		void writeModalResults(std::string fileName, bool writeModes);

		void writeObjective(std::string fileName, StringList& includeFields, bool writeGrad);
		
		//
};

#endif