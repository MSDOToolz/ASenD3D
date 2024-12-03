#include <cmath>
#include <iostream>
#include <iomanip>
#include "ModelClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "LowerTriMatClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"
#include "FaceClass.h"
#include "DiffDoubClass.h"
#include "LoadClass.h"
#include "matrixFunctions.h"

using namespace std;

void Model::reorderNodes(int blockDim) {
	int i1;
	int i2;
	int i3;
	int nd1;
	int nd2;
	int minCt;
	int minNd;
	double minDist;
	int sinceRestart;
	double dist;
	int elNumNds;
	int elDofPerNd;
	int numNodes = nodes.size();
	
	vector<Set> nodalConn(numNodes);
	vector<int> nodeInserted(numNodes);
	
	// Build nodal connectivity
	for (auto& el : elements) {
		elNumNds = el.numNds;
		elDofPerNd = el.dofPerNd;
		for (i1 = 0; i1 < elNumNds; i1++) {
			nd1 = el.nodes[i1];
			if (elDofPerNd > 3) {
				nodes[nd1].numDof = elDofPerNd;
			}
			for (i2 = i1+1; i2 < elNumNds; i2++) {
				nd2 = el.nodes[i2];
				nodalConn[nd1].addIfAbsent(nd2);
				nodalConn[nd2].addIfAbsent(nd1);
			}
		}
	}

	for (auto& thisConst : elasticConst.constVec) {
		for (auto& term1 : thisConst.terms) {
			list<int>& t1Labs = nodeSets[term1.nsPtr].labels;
			i1 = t1Labs.size();
			for(auto& term2 : thisConst.terms) {
				list<int>& t2Labs = nodeSets[term2.nsPtr].labels;
				i2 = t2Labs.size();
				if (i1 > 1 && i2 > 1) {
					auto iter1 = t1Labs.begin();
					auto end1 = t1Labs.end();
					auto iter2 = t2Labs.begin();
					auto end2 = t2Labs.end();
					while (iter1 != end1 && iter2 != end2) {
						nodalConn[*iter1].addIfAbsent(*iter2);
						nodalConn[*iter2].addIfAbsent(*iter1);
						++iter1;
						++iter2;
					}
				}
				else {
					for (auto& nd1 : t1Labs) {
						for (auto& nd2 : t2Labs) {
							nodalConn[nd1].addIfAbsent(nd2);
							nodalConn[nd2].addIfAbsent(nd1);
						}
					}
				}
			}
		}
	}

	for (auto& thisConst : thermalConst.constVec) {
		for (auto& term1 : thisConst.terms) {
			list<int>& t1Labs = nodeSets[term1.nsPtr].labels;
			for (auto& nd1 : t1Labs) {
				for (auto& term2 : thisConst.terms) {
					list<int>& t2Labs = nodeSets[term2.nsPtr].labels;
					for (auto& nd2 : t2Labs) {
						nodalConn[nd1].addIfAbsent(nd2);
						nodalConn[nd2].addIfAbsent(nd1);
					}
				}
			}
		}
	}
	
	// Find node with least connectivity
	minCt = numNodes;
	minNd = 0;
	for (i1 = 0; i1 < numNodes; i1++) {
		i2 = nodalConn[i1].labels.size();
		if(i2 < minCt) {
			minNd = i1;
			minCt = i2;
		}
		nodeInserted[i1] = 0;
	}
	
	// Put nodes into an integer list in level order.
	list<int> orderedNds;
	
	orderedNds.push_back(minNd);
	auto thisNd = orderedNds.begin();
	nodeInserted[minNd] = 1;
	sinceRestart = 0;
	
	while(orderedNds.size() < numNodes) {
		nd1 = *thisNd;
		auto neighborNd = nodalConn[nd1].labels.begin();
		while(neighborNd != nodalConn[nd1].labels.end() && sinceRestart < blockDim) {
			nd2 = *neighborNd;
			if(nodeInserted[nd2] == 0) {
				orderedNds.push_back(nd2);
				nodeInserted[nd2] = 1;
				sinceRestart++;
			}
			++neighborNd;
		}
		if(sinceRestart >= blockDim || thisNd == orderedNds.end()) {
			thisNd = orderedNds.end();
			--thisNd;
			nd1 = *thisNd;
			minDist = 1.0e+100;
			for (i1 = 0; i1 < numNodes; i1++) {
				if(nodeInserted[i1] == 0) {
					dist = getDist(nodes[nd1].coord, nodes[i1].coord);
					if(dist < minDist) {
						minDist = dist;
						minNd = i1;
					}
				}
			}
			orderedNds.push_back(minNd);
			nodeInserted[minNd] = 1;
			sinceRestart++;
			if(sinceRestart >= blockDim) {
				sinceRestart = 0;
			}
		}
		++thisNd;
	}
	
	// Update the global degree of freedom indexes for the nodes
	
    i2 = 0; // index in elastic matrix
	i3 = 0; // sorted rank for solid nodes
	for (auto& ndi : orderedNds) {
		Node& thisNd = nodes[ndi];
		if (!thisNd.fluid) {
			thisNd.sortedRank = i3;
			i3++;
			thisNd.dofIndex[0] = i2;
			i2++;
			thisNd.dofIndex[1] = i2;
			i2++;
			thisNd.dofIndex[2] = i2;
			i2++;
			if (thisNd.numDof == 6) {
				thisNd.dofIndex[3] = i2;
				i2++;
				thisNd.dofIndex[4] = i2;
				i2++;
				thisNd.dofIndex[5] = i2;
				i2++;
			}
		}
	}
	elMatDim = i2;
	elasticMat.setDim(elMatDim);
	nonFrcElMat.setDim(elMatDim);
	thermMat.setDim(nodes.size());

	JobCommand& scmd = job[solveCmd];
	if (scmd.solverMethod == "iterative" && scmd.maxIt == 0) {
		scmd.maxIt = elMatDim;
	}

	for (auto& thisEl : elements) {
		i3 = thisEl.numIntDof;
		if (i3 > 0) {
			thisEl.intDofIndex = i2;
			i2 += i3;
		}
	}
	totGlobDof = i2;
	
	tempV1 = vector<double>(elMatDim);
	tempV2 = vector<double>(elMatDim);
	tempV3 = vector<double>(elMatDim);
	tempV4 = vector<double>(elMatDim);
	tempV5 = vector<double>(elMatDim);
	tempV6 = vector<double>(elMatDim);
	tempD1 = vector<DiffDoub0>(elMatDim);

	i3 = nodes.size();
	dLdU = vector<double>(totGlobDof);
	dLdV = vector<double>(elMatDim);
	dLdA = vector<double>(elMatDim);
	dLdT = vector<double>(i3);
	dLdTdot = vector<double>(i3);
	uAdj = vector<double>(elMatDim);
	vAdj = vector<double>(elMatDim);
	aAdj = vector<double>(elMatDim);
	tAdj = vector<double>(i3);
	tdotAdj = vector<double>(i3);

	dRudD = vector<DiffDoub1>(elMatDim);
	dRtdD = vector<DiffDoub1>(i3);

	elInD = vector<int>(elements.size());

	i3 = designVars.size();
	if (i3 > 0) {
		dLdD = vector<double>(i3);
	}

	return;
}

void Model::buildConstraintMats() {
	for (auto& thisConst : elasticConst.constVec) {
		thisConst.buildMat(nodes,nodeSets);
	}
	for (auto& thisConst : thermalConst.constVec) {
		thisConst.buildMat(nodes,nodeSets);
	}
	return;
}

void Model::updateReference() {
	int i1;
	int i2;
	int i3;
	// Set the section pointers for all elements and material pointers for all sections
	string elSet;
	Set *esPtr;
	string matName;
	i2 = 0;
	for (auto& thisSec : sections) {
		elSet = thisSec.elSetName;
		i1 = esMap.at(elSet);
		for (auto& eli : elementSets[i1].labels) {
			elements[eli].sectPtr = i2;
		}
		i3 = 0;
		for (auto& thisMat : materials) {
			matName = thisMat.name;
			if(matName == thisSec.matName) {
				thisSec.matPtr = i3;
			}
			for (auto& thisLay : thisSec.layers) {
				if(matName == thisLay.matName) {
					thisLay.matPtr = i3;
				}
			}
			i3++;
		}
		i2++;
	}
	
	// Set node & element set pointers in loads, constraints, design variables and objectives
	string ndSet;
	Set *nsPtr;
	for (auto& thisLoad : elasticLoads) {
		try {
			ndSet = thisLoad.nodeSet;
			i1 = nsMap.at(ndSet);
			thisLoad.ndSetPtr = i1;
		}
		catch (...) {
			elSet = thisLoad.elementSet;
			i1 = esMap.at(elSet);
			thisLoad.elSetPtr = i1;
		}
	}
	for (auto& thisLoad : thermalLoads) {
		try {
			ndSet = thisLoad.nodeSet;
			i1 = nsMap.at(ndSet);
			thisLoad.ndSetPtr = i1;
		}
		catch (...) {
			elSet = thisLoad.elementSet;
			i1 = esMap.at(elSet);
			thisLoad.elSetPtr = i1;
		}
	}

	for (auto& thisConst : elasticConst.constVec) {
		for (auto& thisCTerm : thisConst.terms) {
			ndSet = thisCTerm.nodeSet;
			i1 = nsMap.at(ndSet);
			thisCTerm.nsPtr = i1;
		}
	}
	for (auto& thisConst : thermalConst.constVec) {
		for (auto& thisCTerm : thisConst.terms) {
			ndSet = thisCTerm.nodeSet;
			i1 = nsMap.at(ndSet);
			thisCTerm.nsPtr = i1;
		}
	}

	for (auto& thisDV : designVars) {
		try {
			ndSet = thisDV.ndSetName;
			i1 = nsMap.at(ndSet);
			thisDV.ndSetPtr = i1;
		}
		catch (...) {
			elSet = thisDV.elSetName;
			i1 = esMap.at(elSet);
			thisDV.elSetPtr = i1;
		}
	}

	for (auto& thisTerm : objective.terms) {
		try {
			elSet = thisTerm.elSetName;
			i1 = esMap.at(elSet);
			thisTerm.elSetPtr = i1;
		}
		catch (...) {
			ndSet = thisTerm.ndSetName;
			i1 = nsMap.at(ndSet);
			thisTerm.ndSetPtr = i1;
		}
	}
	
	// Build DV reference list for nodes and elements
	int coefLen;
	double constCoef;
	int DVi = 0;
	for (auto& thisDV : designVars) {
		coefLen = thisDV.coefs.size();
		if(thisDV.elSetPtr > -1) {
			if(coefLen < 2) {
				if(coefLen == 0) {
					constCoef = 1.0;
				} else {
					constCoef = *thisDV.coefs.begin();
				}
				for (auto& eli : elementSets[thisDV.elSetPtr].labels) {
					elements[eli].addDesignVariable(DVi,constCoef);
				}
			}
			else {
				list<int>& setLabs = elementSets[thisDV.elSetPtr].labels;
				auto setIter = setLabs.begin();
				auto setEnd = setLabs.end();
				auto coefIter = thisDV.coefs.begin();
				auto coefEnd = thisDV.coefs.end();
				while (setIter != setEnd && coefIter != coefEnd) {
					elements[*setIter].addDesignVariable(DVi, *coefIter);
					++setIter;
					++coefIter;
				}
		    }
		}
		
		if (thisDV.ndSetPtr > -1) {
			if (coefLen < 2) {
				if (coefLen == 0) {
					constCoef = 1.0;
				}
				else {
					constCoef = *thisDV.coefs.begin();
				}
				for (auto& ndi : nodeSets[thisDV.ndSetPtr].labels) {
					nodes[ndi].addDesignVariable(DVi, constCoef);
				}
			}
			else {
				list<int>& setLabs = nodeSets[thisDV.ndSetPtr].labels;
				auto setIter = setLabs.begin();
				auto setEnd = setLabs.end();
				auto coefIter = thisDV.coefs.begin();
				auto coefEnd = thisDV.coefs.end();
				while (setIter != setEnd && coefIter != coefEnd) {
					nodes[*setIter].addDesignVariable(DVi, *coefIter);
					++setIter;
					++coefIter;
				}
			}
		}
		DVi++;
	}
	
	// Build comprehensive element list for each design variable
	
	int elLabel;
	int elNumNds;
	for (auto& thisEl : elements) {
		elLabel = thisEl.label;
		for (auto& dv : thisEl.designVars) {
			designVars[dv.intDat].addCompEl(elLabel);
			thisEl.addCompDVar(dv.intDat);
		}
		elNumNds = thisEl.numNds;
		for (i1 = 0; i1 < elNumNds; i1++) {
			Node& thisNd = nodes[thisEl.nodes[i1]];
			for (auto& dv : thisNd.dVarLst) {
				DesignVariable& thisDV = designVars[dv.intDat];
				if (thisDV.category == "nodeCoord") {
					thisDV.addCompEl(elLabel);
					thisEl.addCompDVar(dv.intDat);
				}
			}
		}
	}
	
	return;
}

void Model::findSurfaceFaces() {
	int i1;
	int lowNd;
	bool added;

	int fcCt = 0;
	for (auto& el : elements) {
		fcCt += el.numFaces;
	}

	faces = vector<Face>(fcCt);

	i1 = 0;
	for (auto& el : elements) {
		el.initializeFaces(faces,i1);
	}
	
	int numNodes = nodes.size();
	vector<FacePtList> fLArray(numNodes);
	for (auto& thisEl : elements) {
		if (thisEl.dofPerNd == 3) {
			for (auto& thisFc : thisEl.faces) {
				lowNd = faces[thisFc].getLowNd();
				added = fLArray[lowNd].addIfAbsent(thisFc,faces);
			}
		}
	}
	
	return;
}

void Model::prepMatrixFactorizations() {
	double zeroAr[9] = { 0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };

	JobCommand& scmd = job[solveCmd];
	if (scmd.thermal) {
		for (auto& thisNd : nodes) {
			thisNd.initializeTemp();
			if (scmd.dynamic) {
				thisNd.updateTdot(scmd.newmarkGamma, scmd.timeStep);
			}
		}
		if (!thermLT.isAllocated()) {
			buildThermalSolnLoad(tempV1, true);
			thermLT.allocateFromSparseMat(thermMat, thermalConst, scmd.solverBlockDim);
		}
		if (!thermScaled) {
			scaleThermalConst();
		}
	}

	if (scmd.elastic) {
		for (auto& thisNd : nodes) {
			thisNd.initializeDisp();
			if (scmd.dynamic) {
				thisNd.updateVelAcc(scmd.newmarkBeta, scmd.newmarkGamma, scmd.timeStep);
			}
		}
		for (auto& thisEl : elements) {
			thisEl.setIntDisp(zeroAr);
			thisEl.setIntPrevDisp(zeroAr);
		}
		if (!elasticLT.isAllocated()) {
			buildElasticSolnLoad(tempV1, true, true);
			scaleElasticConst();
			elasticLT.allocateFromSparseMat(elasticMat, elasticConst, 6 * scmd.solverBlockDim);
		}
	}

	return;
}

void Model::analysisPrep() {
	int i1;
	int i2;
	int numNds;
	int blockDim;
	double time;
	double tInc;

	if (solveCmd > -1) {
		JobCommand& scmd = job[solveCmd];
		if (scmd.solverMethod == "direct") {
			scmd.solverBlockDim = 2000000000;
		}
		else {
			i1 = 6;
			numNds = nodes.size();
			while ((i1 * i1) < numNds) {
				i1 += 6;
			}
			if (scmd.solverBlockDim == 2000000000) {
				scmd.solverBlockDim = i1;
			}
		}
		blockDim = scmd.solverBlockDim;

		if (scmd.staticLoadTime.size() == 0) {
			scmd.staticLoadTime.push_back(0.0);
		}
	}
	else {
		blockDim = 2000000000;
	}

	i1 = 0;
	for (auto& sec : sections) {
		i2 = sec.layers.size();
		if (i2 > i1) {
			i1 = i2;
		}
	}
	if (i1 > 0) {
		d0Pre.allocateLayers(i1);
		d1Pre.allocateLayers(i1);
	}

	updateReference();
	reorderNodes(blockDim);
	buildConstraintMats();
	prepMatrixFactorizations();
	findSurfaceFaces();

	anPrepRun = true;
	return;
}

void Model::buildElasticAppLoad(vector<double>& appLd, double time) {
	// Construct loads from input file
	int i1;
	int numDof;
	int dofInd;
	string ldType;
	DiffDoub0 ndDVLd[6];

	for (i1 = 0; i1 < elMatDim; i1++) {
		tempD1[i1].setVal(0.0);
	}

	//Loads from model input file
	for (auto& thisLoad : elasticLoads) {
		ldType = thisLoad.type;
		if(time >= thisLoad.activeTime[0] && time <= thisLoad.activeTime[1]) {
			if(ldType == "nodalForce") {
				for (auto& ndi : nodeSets[thisLoad.ndSetPtr].labels) {
					Node& thisNd = nodes[ndi];
					numDof = thisNd.numDof;
					for (i1 = 0; i1 < numDof; i1++) {
						dofInd = thisNd.dofIndex[i1];
						appLd[dofInd] += thisLoad.load[i1];
					}
				}
			}
			else {
				for (auto& eli : elementSets[thisLoad.elSetPtr].labels) {
					Element& thisEl = elements[eli];
					thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
					JobCommand& scmd = job[solveCmd];
					thisEl.getAppLoad(tempD1, thisLoad, scmd.nonlinearGeom, d0Pre, sections, faces, nodes, designVars);
				}
			}
		}
	}

	for (i1 = 0; i1 < elMatDim; i1++) {
		appLd[i1] += tempD1[i1].val;
	}
	
	// Design variable dependent loads.
	for (auto& thisNd : nodes) {
		thisNd.getElasticDVLoad(ndDVLd,designVars);
		numDof = thisNd.numDof;
		for (i1 = 0; i1 < numDof; i1++) {
			dofInd = thisNd.dofIndex[i1];
			appLd[dofInd] += ndDVLd[i1].val;
		}
	}
	
	return;
}

void Model::buildThermalAppLoad(vector<double>& appLd, double time) {
	// Construct loads from input file
	int i1;
	int totNodes;
	int numDof;
	int dofInd;
	Node* thisNd;
	Element* thisEl;
	string ldType;
	double actTime[2];
	double ndLoad[6];
	DiffDoub0 ndDVLd;
	Set* thisSet;

	totNodes = nodes.size();
	for (i1 = 0; i1 < totNodes; i1++) {
		tempD1[i1].setVal(0.0);
	}

	//Loads from model input file
	for (auto& thisLoad : thermalLoads) {
		ldType = thisLoad.type;
		if (time >= thisLoad.activeTime[0] && time <= thisLoad.activeTime[1]) {
			if (ldType == "nodalHeatGen") {
				for (auto& ndi : nodeSets[thisLoad.ndSetPtr].labels) {
					Node& thisNd = nodes[ndi];
					dofInd = thisNd.sortedRank;
					appLd[dofInd] += ndLoad[0];
				}
			}
			else {
				for (auto& eli : elementSets[thisLoad.elSetPtr].labels) {
					Element& thisEl = elements[eli];
					thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
					thisEl.getAppThermLoad(tempD1, thisLoad, d0Pre, sections, faces, nodes, designVars);
				}
			}
		}
	}

	for (i1 = 0; i1 < totNodes; i1++) {
		appLd[i1] += tempD1[i1].val;
	}

	// Design variable dependent loads.
	for (auto& thisNd : nodes) {
		thisNd.getThermalDVLoad(ndDVLd, designVars);
		dofInd = thisNd.sortedRank;
		appLd[dofInd] += ndDVLd.val;
	}

	return;
}

void Model::buildElasticSolnLoad(vector<double>& solnLd, bool buildMat, bool fullRef) {
	int i1;
	
	for (i1 = 0; i1 < elMatDim; i1++) {
		tempD1[i1].setVal(0.0);
	}
	if(buildMat) {
		elasticMat.zeroAll();
	}
	if (fullRef) {
		nonFrcElMat.zeroAll();
	}
	
	JobCommand& scmd = job[solveCmd];
	for (auto& thisEl : elements) {
		thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
		if (thisEl.type == 21) {
			thisEl.getRu(tempD1, elasticMat, buildMat, scmd, d0Pre, nodes, designVars);
		}
		else {
			if (fullRef) {
				thisEl.getRu(tempD1, nonFrcElMat, true, scmd, d0Pre, nodes, designVars);
			}
			else {
				thisEl.getRu(tempD1, nonFrcElMat, false, scmd, d0Pre, nodes, designVars);
			}
		}
	}
	
	if (buildMat) {
		elasticMat.addMatrix(nonFrcElMat);
	}
	
	for (i1 = 0; i1 < elMatDim; i1++) {
		solnLd[i1]-= tempD1[i1].val;
	}
	
	return;
}

void Model::buildThermalSolnLoad(vector<double>& solnLd, bool buildMat) {
	int i1;
	int numNodes;

	numNodes = nodes.size();
	for (i1 = 0; i1 < numNodes; i1++) {
		tempD1[i1].setVal(0.0);
	}

	if (buildMat) {
		thermMat.zeroAll();
	}

	JobCommand& scmd = job[solveCmd];
	for (auto& thisEl : elements) {
		thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
		thisEl.getRt(tempD1, thermMat, buildMat, scmd, d0Pre, nodes);
	}

	for (i1 = 0; i1 < numNodes; i1++) {
		solnLd[i1] -= tempD1[i1].val;
	}

	return;
}

void Model::scaleElasticConst() {
	double scaleFact = 100000.0*elasticMat.getMaxAbsVal();
	elasticConst.setScaleFact(scaleFact);
	elasticScaled = true;
	return;
}

void Model::scaleThermalConst() {
	double scaleFact = 100000.0 * thermMat.getMaxAbsVal();
	thermalConst.setScaleFact(scaleFact);
	thermScaled = true;
	return;
}

void Model::buildElasticConstLoad(vector<double>& constLd) {
	int i1;
	int ndDof;
	int globInd;
	for (auto& thisNd : nodes) {
		ndDof = thisNd.numDof;
		for (i1 = 0; i1 < ndDof; i1++) {
			globInd = thisNd.dofIndex[i1];
			tempV2[globInd] = thisNd.displacement[i1];
		}
	}
	elasticConst.getTotalLoad(constLd, tempV2, tempV3, elMatDim);
	return;
}

void Model::buildThermalConstLoad(vector<double>& constLd) {
	int i1;
	int numNodes;
	double ndTemp;
	int globInd;
	for (auto& thisNd : nodes) {
		ndTemp = thisNd.temperature;
		globInd = thisNd.sortedRank;
		tempV2[globInd] = ndTemp;
	}
	numNodes = nodes.size();
	thermalConst.getTotalLoad(constLd, tempV2, tempV3, numNodes);
	return;
}

void Model::solveStep(double time, double appLdFact, bool fullRef) {
	int i1;
	int i2;
	int numNodes;
	int ideb;
	int maxNLit;
	double dUnorm;
	double dUtol;
	double absdU;
	int ndof;
	int dofInd;
	double ndDelDisp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	Node *thisNd;
	Element *thisEl;

	JobCommand& cmd = job[solveCmd];
	if(cmd.thermal) {
		numNodes = nodes.size();
		for (i1 = 0; i1 < numNodes; i1++) {
			tempV1[i1] = 0.0;
			tempV2[11] = 0.0;
		}
		buildThermalAppLoad(tempV1, time);
		buildThermalSolnLoad(tempV1, false);
		if (!thermScaled) {
			scaleThermalConst();
		}
		buildThermalConstLoad(tempV1);
		if (cmd.solverMethod == "direct") {
			thermLT.ldlSolve(tempV2, tempV1);
		}
		else {
			conjGradSparse(tempV2, thermMat, thermalConst, thermLT, tempV1, cmd.convTol, cmd.maxIt);
			//gMResSparse(tempV2, *thermMat, *thermalConst, *thermLT, tempV1, cmd->convTol, cmd->maxIt, cmd->solverBlockDim);
		}
		for (auto& thisNd : nodes) {
			dofInd = thisNd.sortedRank;
			thisNd.temperature += tempV2[dofInd];
			if (cmd.dynamic) {
				thisNd.updateTdot(cmd.newmarkGamma, cmd.timeStep);
			}
		}
	}
	
	if(cmd.elastic) {
		if(cmd.nonlinearGeom) {
			maxNLit = 50;
		} else {
			maxNLit = 1;
		}
		dUtol = 1.0e-12;
		dUnorm = 1.0;
		i2 = 0;
		while(i2 < maxNLit && dUnorm > dUtol) {
			for (i1 = 0; i1 < elMatDim; i1++) {
				tempV1[i1] = 0.0;
				tempV2[i1] = 0.0;
			}
			buildElasticAppLoad(tempV1,time);
			for (i1 = 0; i1 < elMatDim; i1++) {
				tempV1[i1] *= appLdFact;
			}
			buildElasticSolnLoad(tempV1,cmd.nonlinearGeom,fullRef);
			if (!elasticScaled) {
				scaleElasticConst();
			}
			buildElasticConstLoad(tempV1);
			for (auto& thisEl : elements) {
				if(thisEl.numIntDof > 0) {
				    thisEl.updateExternal(tempV1,1,nodes,d0Pre.scrMat1,d0Pre.scrMat2);
				}
			}
			if (cmd.nonlinearGeom) {
				elasticLT.populateFromSparseMat(elasticMat, elasticConst);
				elasticLT.ldlFactor();
			}
			if(cmd.solverMethod == "direct") {
				elasticLT.ldlSolve(tempV2,tempV1);
			}
			else {
				conjGradSparse(tempV2, elasticMat, elasticConst, elasticLT, tempV1, cmd.convTol, cmd.maxIt);
				//gMResSparse(tempV2, *elasticMat, *elasticConst, *elasticLT, tempV1, cmd->convTol, cmd->maxIt, 6*cmd->solverBlockDim);
			}
			for (auto& thisNd : nodes) {
				ndof = thisNd.numDof;
				for (i1 = 0; i1 < ndof; i1++) {
					dofInd = thisNd.dofIndex[i1];
					ndDelDisp[i1] = tempV2[dofInd];
				}
				thisNd.addToDisplacement(ndDelDisp);
				if(cmd.dynamic) {
					thisNd.updateVelAcc(cmd.newmarkBeta,cmd.newmarkGamma,cmd.timeStep);
				}
			}
			for (auto& thisEl : elements) {
				if(thisEl.numIntDof > 0) {
				    thisEl.updateInternal(tempV2,1,nodes,d0Pre.scrMat1,d0Pre.scrMat2);
				}
			}
			dUnorm = 0.0;
			for (i1 = 0; i1 < elMatDim; i1++) {
				absdU = abs(tempV2[i1]);
				if(absdU > dUnorm) {
					dUnorm = absdU;
				}
			}
			if (cmd.nonlinearGeom) {
				cout << "Nonlinear iteration: " << i2 << ", max solution step: " << dUnorm << endl;
			}
			i2++;
		}
	}

	
	return;
}

void Model::solve() {
	JobCommand& cmd = job[solveCmd];
	int i1;
	int i2;
	int i3;
	int sinceRef;
	double appLdFact;
	double time;
	int ldSteps = cmd.loadRampSteps;
	double c1 = 1.0/(ldSteps*ldSteps);
	double zeroAr[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	if(!anPrepRun) {
		analysisPrep();
	}

	if(cmd.thermal) {
		thermLT.populateFromSparseMat(thermMat, thermalConst);
		thermLT.ldlFactor();
	}
	
	if(cmd.elastic && !cmd.nonlinearGeom) {
		elasticLT.populateFromSparseMat(elasticMat, elasticConst);
		cout << "factoring stiffness matrix" << endl;
		elasticLT.ldlFactor();
		cout << "finished factoring" << endl;
	}
	
	if(cmd.dynamic) {
		if (cmd.saveSolnHist) {
			writeTimeStepSoln(0);
		}
		time = 0.0;
		i1 = 1;
		sinceRef = 0;
		while (time < cmd.simPeriod) {
			sinceRef++;
			if (sinceRef == cmd.fullReform) {
				solveStep(time, 1.0, true);
				sinceRef = 0;
			}
			else {
				solveStep(time, 1.0, false);
			}
			for (auto& thisNd : nodes) {
				if (cmd.thermal) {
					thisNd.advanceTemp();
					thisNd.updateTdot(cmd.newmarkGamma,cmd.timeStep);
				}
				if (cmd.elastic) {
					thisNd.advanceDisp();
					thisNd.updateVelAcc(cmd.newmarkBeta,cmd.newmarkGamma,cmd.timeStep);
				}
			}
			if (cmd.elastic) {
				for (auto& thisEl : elements) {
					thisEl.advanceIntDisp();
				}
			}
			if (cmd.saveSolnHist) {
				writeTimeStepSoln(i1);
				timeStepsSaved = i1;
				i1++;
			}
			time += cmd.timeStep;
		}
	} else {
		i3 = 0;
		for (auto& thisLd : cmd.staticLoadTime) {
			if (cmd.nonlinearGeom) {
				for (i1 = 0; i1 < ldSteps; i1++) {
					i2 = ldSteps - i1 - 1;
					appLdFact = 1.0 - c1 * i2 * i2;
					solveStep(thisLd, appLdFact, true);
				}
			}
			else {
				appLdFact = 1.0;
				solveStep(thisLd, appLdFact, true);
			}
			if (cmd.saveSolnHist) {
				for (auto& thisNd : nodes) {
					if (cmd.thermal) {
						thisNd.advanceTemp();
					}
					if (cmd.elastic) {
						thisNd.advanceDisp();
					}
				}
				if (cmd.elastic) {
					for (auto& thisEl : elements) {
						thisEl.advanceIntDisp();
					}
				}
				writeTimeStepSoln(i3);
				timeStepsSaved = i3 + 1;
			}
			i3++;
		}
	}
	
	return;
}

void Model::zeroSolution(list<string>& fields) {
	double zeroAr[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	int nIntDof;

	bool dispInc = false;
	for (auto& thisField : fields) {
		if (thisField == "displacement") {
			dispInc = true;
		}
	}

	for (auto& thisNd : nodes) {
		for (auto& thisField : fields) {
			if (thisField == "temperature") {
				thisNd.temperature = 0.0;
			}
			else if (thisField == "tdot") {
				thisNd.tempChangeRate = 0.0;
			}
			else if (thisField == "displacement") {
				thisNd.setDisplacement(zeroAr);
			}
			else if (thisField == "velocity") {
				thisNd.setVelocity(zeroAr);
			}
			else if (thisField == "acceleration") {
				thisNd.setAcceleration(zeroAr);
			}
		}
	}

	if (dispInc) {
		for (auto& thisEl : elements) {
			nIntDof = thisEl.numIntDof;
			if (nIntDof > 0) {
				thisEl.setIntDisp(zeroAr);
			}
		}
	}

	return;
}

void Model::eigenSolve() {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int nDof;
	int globInd;
	int stepInc;
	double ndDisp[6];
	vector<DiffDoub0> Rvec(33);
	vector<double> Rv2(33);
	vector<double> dRdA(1089);
	double zeros[9] = {0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};
	double shft;
	double vKv;
	double dp;
	double mAvg;
	double mMin;

	JobCommand& cmd = job[modalCmd];
	JobCommand& scmd = job[solveCmd];

	if (!anPrepRun) {
		analysisPrep();
	}

	if (eigVecs.size() == 0) {
		i1 = cmd.numModes;
		eigVals = vector<double>(i1);
		loadFact = vector<double>(i1);
		i1 *= elMatDim;
		eigVecs = vector<double>(i1);
		diagMass = vector<double>(elMatDim);
	}

	if (cmd.type == "buckling") {
		for (i1 = 0; i1 < elMatDim; i1++) {
			diagMass[i1] = 1.0;
		}
	}
	else {
		for (i1 = 0; i1 < elMatDim; i1++) {
			diagMass[i1] = 0.0;
		}
		for (auto& thisEl : elements) {
			thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
			i2 = thisEl.numNds * thisEl.dofPerNd;
			for (i1 = 0; i1 < i2; i1++) {
				d0Pre.globAcc[i1].setVal(1.0);
			}
			thisEl.getRum(Rvec, dRdA, false, true, scmd.nonlinearGeom, d0Pre, nodes, designVars);
			for (i1 = 0; i1 < i2; i1++) {
				Rv2[i1] = Rvec[i1].val;
			}
			thisEl.addToGlobVec(Rv2, diagMass, false, false, nodes);
		}
	}
	
	mAvg = 0.0;
	for (i1 = 0; i1 < elMatDim; i1++) {
		mAvg += abs(diagMass[i1]);
	}
	mAvg /= elMatDim;
	mMin = 1.0e-6 * mAvg;
	for (i1 = 0; i1 < elMatDim; i1++) {
		if (diagMass[i1] < mMin) {
			cout << "Warning: small or negative value in diagonal mass matrix: " << diagMass[i1] << endl;
			cout << "Adjusting to minimum value" << mMin << endl;
			diagMass[i1] = mMin;
		}
	}

	cmd.nonlinearGeom = scmd.nonlinearGeom;
	cmd.dynamic = scmd.dynamic;
	scmd.nonlinearGeom = true;
	scmd.dynamic = false;
	cout << "building stiffness matrix" << endl;
	buildElasticSolnLoad(tempV1, true, true);
	cout << "finished building matrix" << endl;
	if (cmd.type == "buckling" || cmd.type == "frequency") {
		for (i1 = 0; i1 < elMatDim; i1++) {
			shft = -cmd.tgtEval * diagMass[i1];
			elasticMat.addEntry(i1, i1, shft);
		}
		if (cmd.solverMethod == "direct") {
			elasticLT.populateFromSparseMat(elasticMat, elasticConst);
			cout << "factoring stiffness matrix" << endl;
			elasticLT.ldlFactor();
			cout << "finished factoring.  beginning eigensolve" << endl;
			eigenSparseDirect(eigVals, eigVecs, cmd.numModes, elasticLT, diagMass, elMatDim);
			cout << "finished eigensolve" << endl;
		}
		else {

		}
		for (i1 = 0; i1 < cmd.numModes; i1++) {
			eigVals[i1] += cmd.tgtEval;
		}
	}
	else {
		//stepInc = timeStepsSaved / cmd.numModes;
		//if (stepInc == 0) {
		//	string erStr = "Error: not enough solution time steps have been saved to find the requested number of active eigenmodes.\n";
		//	erStr += "Be sure to set the saveSolnHist option to true in the dynamic solve command.\n";
		//	throw invalid_argument(erStr);
		//}
		//for (i1 = 0; i1 < cmd.numModes; i1++) {
		//	i2 = (i1 + 1) * stepInc;
		//	i3 = i1 * elMatDim;
		//	readTimeStepSoln(i2);
		//	thisNd = nodes.getFirst();
		//	while (thisNd) {
		//		thisNd->getPrevDisp(ndDisp);
		//		nDof = thisNd->getNumDof();
		//		for (i4 = 0; i4 < nDof; i4++) {
		//			globInd = thisNd->getDofIndex(i4);
		//			eigVecs[i3 + globInd] = ndDisp[i4];
		//		}
		//		thisNd = thisNd->getNext();
		//	}
		//}
		////getNearestEvecRQ(*elasticMat, *elasticConst, diagMass, eigVecs, eigVals, cmd->numModes, cmd->maxIt);
		//getNearestEvecSubspace(elasticMat, elasticConst, diagMass, eigVecs, eigVals, cmd->numModes);
	}

	if (cmd.type == "buckling") {
		writeTimeStepSoln(-1);
		for (auto& thisNd : nodes) {
			thisNd.advanceDisp();
			thisNd.setDisplacement(zeros);
		}
		for (auto& thisEl : elements) {
			thisEl.advanceIntDisp();
			thisEl.setIntDisp(zeros);
		}
		buildElasticSolnLoad(tempV1, true, true);
		i2 = 0;
		for (i1 = 0; i1 < cmd.numModes; i1++) {
			for (i3 = 0; i3 < elMatDim; i3++) {
				tempV1[i3] = 0.0;
			}
			subVec(tempV2, eigVecs, i2, i2 + elMatDim);
			elasticMat.vectorMultiply(tempV1, tempV2, false);
			vKv = 0.0;
			for (i3 = 0; i3 < elMatDim; i3++) {
				vKv += eigVecs[i2 + i3] * tempV1[i3];
			}
			try {
				loadFact[i1] = vKv / (vKv - eigVals[i1]);
			}
			catch (...) {
				loadFact[i1] = 1.0e+100;
			}
			i2 += elMatDim;
		}
		for (auto& thisNd : nodes) {
			thisNd.backstepDisp();
		}
		for (auto& thisEl : elements) {
			thisEl.backstepIntDisp();
		}
		readTimeStepSoln(-1);
	}
	else {
		for (i1 = 0; i1 < cmd.numModes; i1++) {
			try {
				loadFact[i1] = 0.159154943091895335 * sqrt(eigVals[i1]);
			}
			catch (...) {
				loadFact[i1] = -1.0;
			}
		}
	}

	scmd.nonlinearGeom = cmd.nonlinearGeom;
	scmd.dynamic = cmd.dynamic;

	return;
}


void Model::setSolnToMode(string field, int mode, double maxVal) {
	int i1;
	int i2;
	int modeSt;
	int nDof;
	int globInd;
	string dispFields = "displacement velocity acceleration";
	string erStr;
	double ndDat[6];
	double scaleFact;
	double abVal;
	
	modeSt = mode * elMatDim;
	i2 = modeSt;
	scaleFact = 0.0;
	for (i1 = 0; i1 < elMatDim; i1++) {
		abVal = abs(eigVecs[i2]);
		if (abVal > scaleFact) {
			scaleFact = abVal;
		}
		i2++;
	}
	scaleFact = maxVal / scaleFact;

	i1 = dispFields.find(field);
	if (i1 > -1) {
	    for (auto& thisNd : nodes) {
			nDof = thisNd.numDof;
			for (i2 = 0; i2 < nDof; i2++) {
				globInd = thisNd.dofIndex[i2];
				ndDat[i2] = scaleFact * eigVecs[modeSt + globInd];
			}
			if (field == "displacement") {
				thisNd.setDisplacement(ndDat);
			}
			else if (field == "velocity") {
				thisNd.setVelocity(ndDat);
			}
			else {
				thisNd.setAcceleration(ndDat);
			}
		}
	}
	else {
		erStr = "Warning: " + field + " is not a valid field to set to a mode. Aborting setSolnToMode()";
		cout << erStr << endl;
	}

	return;
}

void Model::augmentdLdU() {
	int i1;
	int i2;
	int i3;
	int numNds;
	int ndDof;
	JobCommand& scmd = job[solveCmd];
	double dt = scmd.timeStep;
	double bet = scmd.newmarkBeta;
	double gam = scmd.newmarkGamma;
	double dTdotdTp = -1.0 / (gam * dt);
	double dTdotdTdotp = -(1.0 - gam) / gam;
	double dAdUp = 1.0 / (dt * dt * (bet - gam));
	double dAdVp = dAdUp * dt;
	double dAdAp = dAdUp * dt * dt * (0.5 + bet - gam);
	double dVdUp = dt * gam * dAdUp;
	double dVdVp = 1.0 + dt * gam * dAdVp;
	double dVdAp = dt * (1.0 - gam) + dt * gam * dAdAp;

	double c1 = dAdUp;
	double c2 = dAdAp;
	double c3 = dt * (1.0 - gam);
	double c4 = 1.0 / (dt*(bet - gam));

	vector<DiffDoub0> Rvec(30);
	vector<double> Mmat(900);
	vector<double> Dmat(900);
	vector<double> elAdj(30);
	vector<double> eldLdT(10);
	vector<double> eldLdTdot(10);
	vector<double> eldLdU(30);
	vector<double> eldLdV(30);
	vector<double> eldLdA(30);

	Element* thisEl;

	if (scmd.thermal) {
		for (auto& thisEl : elements) {
			thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
			thisEl.getRtm(Rvec, Mmat, true, true, d0Pre);
			thisEl.getElVec(elAdj, tAdj, true, false, nodes);
			numNds = thisEl.numNds;
			i3 = 0;
			for (i1 = 0; i1 < numNds; i1++) {
				eldLdT[i1] = 0.0;
				eldLdTdot[i1] = 0.0;
				for (i2 = 0; i2 < numNds; i2++) {
					eldLdT[i1] -= dTdotdTp * Mmat[i3] * elAdj[i2];
					eldLdTdot[i1] -= dTdotdTdotp * Mmat[i3] * elAdj[i2];
					i3++;
				}
			}
			thisEl.addToGlobVec(eldLdT, dLdT, true, false, nodes);
			thisEl.addToGlobVec(eldLdTdot, dLdTdot, true, false, nodes);
		}
		i2 = nodes.size();
		for (i1 = 0; i1 < i2; i1++) {
			dLdT[i1] += tdotAdj[i1];
			dLdTdot[i1] += c3 * tdotAdj[i1];
		}
	}

	if (scmd.elastic) {
		for (auto& thisEl : elements) {
			thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
			thisEl.getRum(Rvec, Mmat, true, true, scmd.nonlinearGeom, d0Pre, nodes, designVars);
			thisEl.getRud(Rvec, Dmat, true, scmd, d0Pre, nodes, designVars);
			thisEl.getElVec(elAdj, uAdj, false, false, nodes);
			ndDof = thisEl.numNds * thisEl.dofPerNd;
			i3 = 0;
			for (i1 = 0; i1 < ndDof; i1++) {
				eldLdU[i1] = 0.0;
				eldLdV[i1] = 0.0;
				eldLdA[i1] = 0.0;
				for (i2 = 0; i2 < ndDof; i2++) {
					eldLdU[i1] -= (dAdUp * Mmat[i3] + dVdUp * Dmat[i3]) * elAdj[i2];
					eldLdA[i1] -= (dAdAp * Mmat[i3] + dVdAp * Dmat[i3]) * elAdj[i2];
					eldLdV[i1] -= (dAdVp * Mmat[i3] + dVdVp * Dmat[i3]) * elAdj[i2];
					i3++;
				}
			}
			thisEl.addToGlobVec(eldLdU, dLdU, false, false, nodes);
			thisEl.addToGlobVec(eldLdA, dLdA, false, false, nodes);
			thisEl.addToGlobVec(eldLdV, dLdV, false, false, nodes);
		}
		for (i1 = 0; i1 < elMatDim; i1++) {
			dLdU[i1] += c1 * aAdj[i1];
			dLdA[i1] += (c2 * aAdj[i1] + c3 * vAdj[i1]);
			dLdV[i1] += c4 * aAdj[i1] + vAdj[i1];
		}
	}

	return;
}

void Model::solveForAdjoint(double time, bool fullRef) {
	int i1;
	int i2;
	int i3;
	int totNodes;
	int elNumNds;
	int elDofPerNd;
	int elNdDof;
	int elIntDof;
	double c1;
	double c2;
	double c3;
	Element* thisEl;
	vector<DiffDoub0> Rvec(33);
	vector<double> dRdU(1089);
	vector<double> dRdT(330);
	vector<double> eldLdT(10);
	vector<double> elAdj(33);
	double* intAdj;
	vector<double> scr1(33);
	vector<double> scr2(33);
	JobCommand& scmd = job[solveCmd];

	objective.calculatedLdU(dLdU, dLdV, dLdA, dLdT, dLdTdot, time, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d0Pre);
    
	if (scmd.elastic) {
		if (scmd.dynamic) {
			c1 = scmd.timeStep * scmd.newmarkGamma;
			c2 = scmd.timeStep;
			c2 = 1.0 / (c2 * c2 * (scmd.newmarkBeta - scmd.newmarkGamma));
			for (i1 = 0; i1 < elMatDim; i1++) {
				vAdj[i1] = dLdV[i1];
				aAdj[i1] = (dLdA[i1] + c1 * vAdj[i1]);
				dLdU[i1] -= c2 * aAdj[i1];
			}
		}
		if (scmd.nonlinearGeom) {
			buildElasticSolnLoad(tempV1, true, fullRef);
			elasticLT.populateFromSparseMat(elasticMat, elasticConst);
			elasticLT.ldlFactor();
		}
		for (auto& thisEl : elements) {
			thisEl.setIntdLdU(dLdU);
			thisEl.updateExternal(dLdU, 0, nodes, scr1, scr2);
		}
		if (scmd.solverMethod == "direct") {
			elasticLT.ldlSolve(uAdj, dLdU);
		}
		else {
			conjGradSparse(uAdj, elasticMat, elasticConst, elasticLT, dLdU, scmd.convTol, scmd.maxIt);
			//gMResSparse(uAdj, *elasticMat, *elasticConst, *elasticLT, dLdU, solveCmd->convTol, solveCmd->maxIt, 6*solveCmd->solverBlockDim);
		}
		for (auto& thisEl : elements) {
			thisEl.updateInternal(uAdj, 0, nodes, scr1, scr2);
		}
	}

	if (scmd.thermal) {
		if (scmd.elastic) {
			for (auto& thisEl : elements) {
				thisEl.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
				thisEl.getRuk(Rvec, dRdU, dRdT, true, scmd.nonlinearGeom, d0Pre, nodes, designVars);
				thisEl.getElVec(elAdj, uAdj, false, false, nodes);
				elNumNds = thisEl.numNds;
				elDofPerNd = thisEl.dofPerNd;
				elIntDof = thisEl.numIntDof;
				elNdDof = elNumNds * elDofPerNd;
				for (i1 = 0; i1 < elNumNds; i1++) {
					eldLdT[i1] = 0.0;
					i3 = i1;
					for (i2 = 0; i2 < elNdDof; i2++) {
						//i3 = i2 * elNumNds + i1;
						eldLdT[i1] -= dRdT[i3] * elAdj[i2];
						i3 += elNumNds;
					}
					for (i2 = 0; i2 < elIntDof; i2++) {
						eldLdT[i1] -= dRdT[i3] * elAdj[i2];
						i3 += elNumNds;
					}
				}
				thisEl.addToGlobVec(eldLdT, dLdT, true, false, nodes);
			}
		}
		totNodes = nodes.size();
		if (scmd.dynamic) {
			c3 = -1.0 / (scmd.timeStep * scmd.newmarkGamma);
			for (i1 = 0; i1 < totNodes; i1++) {
				tdotAdj[i1] = c3 * dLdTdot[i1];
				dLdT[i1] -= tdotAdj[i1];
			}
		}
		if (scmd.solverMethod == "direct") {
			thermLT.ldlSolve(tAdj, dLdT);
		}
		else {
			conjGradSparse(tAdj, thermMat, thermalConst, thermLT, dLdT, scmd.convTol, scmd.maxIt);
			//gMResSparse(tAdj, *thermMat, *thermalConst, *thermLT, dLdT, solveCmd->convTol, solveCmd->maxIt, solveCmd->solverBlockDim);
		}
	}
	
	return;
}

void Model::dRthermaldD(int dVarNum) {
	int i1;
	int totNodes;
	int numDof;
	int globInd;
	DiffDoub0 dvVal;

	JobCommand& scmd = job[solveCmd];
	DesignVariable& thisDV = designVars[dVarNum];
	thisDV.getValue(dvVal);
	thisDV.diffVal.setVal(dvVal.val,1.0);

	totNodes = nodes.size();
	for (i1 = 0; i1 < totNodes; i1++) {
		dRtdD[i1].setVal(0.0);
	}

	for (i1 = 0; i1 < elements.size(); i1++) {
		elInD[i1] = 0;
	}

	for (auto& eli : thisDV.compElList) {
		elInD[eli] = 1;
	}
	
	// Applied contribution
	for (auto& thisLd : thermalLoads) {
		if (thisLd.type != "nodalHeatGen") {
			for (auto& eli : elementSets[thisLd.elSetPtr].labels) {
				if (elInD[eli]) {
					Element& thisEl = elements[eli];
					thisEl.getStressPrereq(d1Pre, sections, materials, nodes, designVars);
					thisEl.getAppThermLoad(dRtdD, thisLd, d1Pre, sections, faces, nodes, designVars);
				}
			}
		}
	}

	for (i1 = 0; i1 < totNodes; i1++) {
		dRtdD[i1].neg();
	}

	// Solution-dependent contribution of load
	for (auto& eli : thisDV.compElList) {
		Element& thisEl = elements[eli];
		thisEl.getStressPrereq(d1Pre, sections, materials, nodes, designVars);
		thisEl.getRt(dRtdD,thermMat,false,scmd,d1Pre,nodes);
	}


	// Design variable dependent contribution
	DiffDoub1 ndLd;
	for (auto& ndi : nodeSets[thisDV.ndSetPtr].labels) {
		Node& thisNd = nodes[ndi];
		thisNd.getThermalDVLoad(ndLd, designVars);
		globInd = thisNd.sortedRank;
		dRtdD[globInd].sub(ndLd);
	}

	thisDV.diffVal.setVal(dvVal.val, 0.0);
	
	return;
}

void Model::dRelasticdD(int dVarNum) {
	int i1;
	int numDof;
	int globInd;
	DiffDoub0 dvVal;

	JobCommand& scmd = job[solveCmd];
	DesignVariable& thisDV = designVars[dVarNum];
	thisDV.getValue(dvVal);
	thisDV.diffVal.setVal(dvVal.val, 1.0);

	for (i1 = 0; i1 < elMatDim; i1++) {
		dRudD[i1].setVal(0.0);
	}

	for (i1 = 0; i1 < elements.size(); i1++) {
		elInD[i1] = 0;
	}

	for (auto& eli : thisDV.compElList) {
		elInD[eli] = 1;
	}
	
	// Applied contribution

    for (auto& thisLd : elasticLoads) {
		if (thisLd.type != "nodalForce") {
			for (auto eli : elementSets[thisLd.elSetPtr].labels) {
				if (elInD[eli] > 0) {
					Element& thisEl = elements[eli];
					thisEl.getStressPrereq(d1Pre, sections, materials, nodes, designVars);
					thisEl.getAppLoad(dRudD, thisLd, scmd.nonlinearGeom, d1Pre, sections, faces, nodes, designVars);
				}
			}
		}
	}

	for (i1 = 0; i1 < elMatDim; i1++) {
		dRudD[i1].neg();
	}

	// Solution-dependent contribution of load
	for (auto& eli : thisDV.compElList) {
		Element& thisEl = elements[eli];
		thisEl.getStressPrereq(d1Pre, sections, materials, nodes, designVars);
		thisEl.getRu(dRudD, elasticMat, false, scmd, d1Pre, nodes, designVars);
	}


	// Design variable dependent contribution
	DiffDoub1 ndLd[6];
	for (auto& ndi : nodeSets[thisDV.ndSetPtr].labels) {
		Node& thisNd = nodes[ndi];
		thisNd.getElasticDVLoad(ndLd, designVars);
		numDof = thisNd.numDof;
		for (i1 = 0; i1 < numDof; i1++) {
			globInd = thisNd.dofIndex[i1];
			ndLd[i1].neg();
			dRudD[globInd].add(ndLd[i1]);
		}
	}

	thisDV.diffVal.setVal(dvVal.val, 0.0);

	return;
}

void Model::getObjective() {
	int i1;
	double time;
	Node* thisNd;
	Element* thisEl;
	string erStr;

	JobCommand& scmd = job[solveCmd];
	if (!scmd.saveSolnHist) {
		erStr = "Error: The solution history must be saved to calculate the objective gradient.\n";
		erStr += "Save history by setting saveSolnHist=True in the solve command.\n";
		throw invalid_argument(erStr);
	}

	if (scmd.dynamic) {
		objective.clearValues();
		readTimeStepSoln(timeStepsSaved);
		i1 = timeStepsSaved - 1;
		time = scmd.timeStep * timeStepsSaved;
		while (i1 >= 0) {
			for (auto& thisNd : nodes) {
				if (scmd.elastic) {
					thisNd.backstepDisp();
				}
				if (scmd.thermal) {
					thisNd.backstepTemp();
				}
			}
			if (scmd.elastic) {
				for (auto& thisEl : elements) {
					thisEl.backstepIntDisp();
				}
			}
			readTimeStepSoln(i1);
			objective.calculateTerms(time, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d0Pre);
			i1--;
			time -= scmd.timeStep;
		}
	}
	else {
		objective.clearValues();
		i1 = 0;
		for (auto& thisLdTm : scmd.staticLoadTime) {
			readTimeStepSoln(i1);
			for (auto& thisNd : nodes) {
				if (scmd.elastic) {
					thisNd.backstepDisp();
				}
				if (scmd.thermal) {
					thisNd.backstepTemp();
				}
			}
			if (scmd.elastic) {
				for (auto& thisEl : elements) {
					thisEl.backstepIntDisp();
				}
			}
			objective.calculateTerms(thisLdTm, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d0Pre);
			i1++;
		}
	}
	return;
}

void Model::getObjGradient() {
	int i1;
	int i2;
	int i3;
	int i4;
	double time;
	int sinceRef;
	int numDV = designVars.size();
	Element* thisEl;
	Node* thisNd;
	string erStr;

	JobCommand& scmd = job[solveCmd];
	if (!scmd.saveSolnHist) {
		erStr = "Error: The solution history must be saved to calculate the objective gradient.\n";
		erStr += "Save history by setting saveSolnHist=True in the solve command.\n";
		throw invalid_argument(erStr);
	}

	objective.clearValues();
	for (i1 = 0; i1 < numDV; i1++) {
		dLdD[i1] = 0.0;
	}

	if (scmd.dynamic) {
		time = scmd.timeStep * timeStepsSaved;
		i1 = timeStepsSaved;
		sinceRef = 0;
		readTimeStepSoln(i1);
		while (i1 > 0) {
			for (i2 = 0; i2 < nodes.size(); i2++) {
				dLdT[i2] = 0.0;
				dLdTdot[i2] = 0.0;
			}
			for (i2 = 0; i2 < elMatDim; i2++) {
				dLdA[i2] = 0.0;
				dLdV[i2] = 0.0;
			}
			for (i2 = 0; i2 < totGlobDof; i2++) {
				dLdU[i2] = 0.0;
			}
			if (i1 < timeStepsSaved) {
				augmentdLdU();
			}
			for (auto& thisNd : nodes) {
				if (scmd.elastic) {
					thisNd.backstepDisp();
				}
				if (scmd.thermal) {
					thisNd.backstepTemp();
				}
			}
			for (auto& thisEl : elements) {
				thisEl.backstepIntDisp();
			}
			readTimeStepSoln(i1 - 1);
			objective.calculateTerms(time, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d0Pre);
			sinceRef++;
			if (sinceRef == scmd.fullReform) {
				solveForAdjoint(time, true);
				sinceRef = 0;
			}
			else {
				solveForAdjoint(time, false);
			}
			objective.calculatedLdD(dLdD, time, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d1Pre);
			for (i2 = 0; i2 < numDV; i2++) {
				if (scmd.thermal) {
					dRthermaldD(i2);
					for (i3 = 0; i3 < nodes.size(); i3++) {
						dLdD[i2] -= tAdj[i3] * dRtdD[i3].dval;
					}
				}
				if (scmd.elastic) {
					dRelasticdD(i2);
					for (i3 = 0; i3 < elMatDim; i3++) {
						dLdD[i2] -= uAdj[i3] * dRudD[i3].dval;
					}
					for (auto& thisEl : elements) {
						dLdD[i2] -= thisEl.getIntAdjdRdD();
					}
				}
			}
			i1--;
		}
	}
	else {
		i4 = 0;
		for (auto& thisLd : scmd.staticLoadTime) {
			readTimeStepSoln(i4);
			for (auto& thisNd : nodes) {
				if (scmd.elastic) {
					thisNd.backstepDisp();
				}
				if (scmd.thermal) {
					thisNd.backstepTemp();
				}
			}
			for (auto& thisEl : elements) {
				thisEl.backstepIntDisp();
			}
			objective.calculateTerms(thisLd, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d0Pre);
			for (i2 = 0; i2 < nodes.size(); i2++) {
				dLdT[i2] = 0.0;
				dLdTdot[i2] = 0.0;
			}
			for (i2 = 0; i2 < elMatDim; i2++) {
				dLdA[i2] = 0.0;
				dLdV[i2] = 0.0;
			}
			for (i2 = 0; i2 < totGlobDof; i2++) {
				dLdU[i2] = 0.0;
			}
			solveForAdjoint(thisLd,true);
			objective.calculatedLdD(dLdD, thisLd, scmd.nonlinearGeom, nodes, elements, nodeSets, elementSets, sections, materials, designVars, d1Pre);
			for (i1 = 0; i1 < numDV; i1++) {
				if (scmd.thermal) {
					dRthermaldD(i1);
					for (i3 = 0; i3 < nodes.size(); i3++) {
						dLdD[i1] -= tAdj[i3] * dRtdD[i3].dval;
					}
				}
				if (scmd.elastic) {
					dRelasticdD(i1);
					for (i2 = 0; i2 < elMatDim; i2++) {
						dLdD[i1] -= uAdj[i2] * dRudD[i2].dval;
					}
					for (auto& thisEl : elements) {
						dLdD[i1] -= thisEl.getIntAdjdRdD();
					}
				}
			}
			i4++;
		}
	}
	return;
}