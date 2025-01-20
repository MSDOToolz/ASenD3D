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
	public:
	    std::vector<Node> nodes;
		std::vector<Element> elements;
		std::vector<Face> faces;
		std::vector<Set> nodeSets;
		std::map<std::string, int> nsMap;
		std::vector<Set> elementSets;
		std::map<std::string, int> esMap;
		std::vector<Section> sections;
		std::vector<Material> materials;
		std::vector<Fluid> fluids;
		ConstraintList elasticConst;
		ConstraintList thermalConst;
		std::vector<Load> elasticLoads;
		std::vector<Load> thermalLoads;
		std::vector<DesignVariable> designVars;
		Objective objective;
		std::vector<JobCommand> job;
		
		int elMatDim;
		int totGlobDof;
		bool anPrepRun;
		int timeStepsSaved;
		int solveCmd;
		int modalCmd;
		DiffDoub0StressPrereq d0Pre;
		DiffDoub1StressPrereq d1Pre;
		
		SparseMat elasticMat;
		LowerTriMat elasticLT;
		SparseMat nonFrcElMat;
		bool elasticScaled;
		SparseMat thermMat;
		LowerTriMat thermLT;
		bool thermScaled;

		std::vector<double> eigVecs;
		std::vector<double> eigVals;
		std::vector<double> diagMass;
		std::vector<double> loadFact;
		
		std::vector<double> tempV1;
		std::vector<double> tempV2;
		std::vector<double> tempV3;
		std::vector<double> tempV4;
		std::vector<double> tempV5;
		std::vector<double> tempV6;

		std::vector<DiffDoub0> tempD1;

		std::vector<double> dLdU;
		std::vector<double> dLdV;
		std::vector<double> dLdA;
		std::vector<double> dLdT;
		std::vector<double> dLdTdot;
		std::vector<double> uAdj;
		std::vector<double> vAdj;
		std::vector<double> aAdj;
		std::vector<double> tAdj;
		std::vector<double> tdotAdj;

		std::vector<DiffDoub1> dRudD;
		std::vector<DiffDoub1> dRtdD;

		std::vector<int> elInD;

		std::vector<double> dLdD;
	
	    Model();

		bool key_in_map(std::map<std::string, int>& inMap, std::string& key);

		bool is_int(std::string& inStr);

		bool is_doub(std::string& inStr);
		
		void executeJob();
		
		// Input
		
		void readInputLine(std::string& fileLine,std::string headings[],int hdLdSpace[], std::string data[], int& dataLen);
		
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
		
		void buildElasticAppLoad(std::vector<double>& appLd, double time);

		void buildThermalAppLoad(std::vector<double>& appLd, double time);
		
		void buildElasticSolnLoad(std::vector<double>& solnLd, bool buildMat, bool fullRef);

		void buildThermalSolnLoad(std::vector<double>& solnLd, bool buildMat);

		void scaleElasticConst();

		void scaleThermalConst();

		void buildElasticConstLoad(std::vector<double>& constLd);

		void buildThermalConstLoad(std::vector<double>& constLd);
		
		void solveStep(double time, double appLdFact, bool fullRef);
		
		void solve();

		void zeroSolution(std::list<std::string>& fields);
		
		void eigenSolve();

	/*	void backupElastic();

		void restoreElastic();*/

		void setSolnToMode(std::string field, int mode, double maxVal);

		void augmentdLdU();

		void solveForAdjoint(double time, bool fullRef);

		void dRthermaldD(int dVarNum);

		void dRelasticdD(int dVarNum);

		void getObjective();
		
		void getObjGradient();
		
		// Output

		void writeTimeStepSoln(int tStep);
		
		void writeNodeResults(std::string fileName, std::string nodeSet, std::list<std::string>& fields, int timeStep);

		void writeElementResults(std::string fileName, std::string elSet, std::list<std::string>& fields, std::string position, int timeStep);

		void writeModalResults(std::string fileName, bool writeModes);

		void writeObjective(std::string fileName, std::list<std::string>& includeFields, bool writeGrad);
		
		//
};

#endif