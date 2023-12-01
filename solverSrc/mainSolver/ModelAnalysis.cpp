#include <cmath>
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
	double p1[3];
	double p2[3];
	double dist;
	Element *thisEl;
	int *elNodes;
	int elNumNds;
	int elDofPerNd;
	int numNodes = nodes.getLength();
	
	IntList *nodalConn = new IntList[numNodes];
	int *nodeInserted = new int[numNodes];
	
	// Build nodal connectivity
	thisEl = elements.getFirst();
	while(thisEl) {
		elNumNds = thisEl->getNumNds();
		elDofPerNd = thisEl->getDofPerNd();
		elNodes = thisEl->getNodes();
		for (i1 = 0; i1 < elNumNds; i1++) {
			nd1 = elNodes[i1];
			if (elDofPerNd > 3) {
				nodeArray[nd1].ptr->setNumDof(elDofPerNd);
			}
			for (i2 = i1+1; i2 < elNumNds; i2++) {
				nd2 = elNodes[i2];
				nodalConn[nd1].addIfAbsent(nd2);
				nodalConn[nd2].addIfAbsent(nd1);
			}
		}
		thisEl = thisEl->getNext();
	}
	
	// Find node with least connectivity
	minCt = numNodes;
	minNd = 0;
	for (i1 = 0; i1 < numNodes; i1++) {
		i2 = nodalConn[i1].getLength();
		if(i2 < minCt) {
			minNd = i1;
			minCt = i2;
		}
		nodeInserted[i1] = 0;
	}
	
	// Put nodes into an integer list in level order.
	IntList orderedNds;
	IntListEnt *thisNd;
	IntListEnt *neighborNd;
	
	orderedNds.addEntry(minNd);
	thisNd = orderedNds.getFirst();
	nodeInserted[minNd] = 1;
	sinceRestart = 0;
	
	while(orderedNds.getLength() < numNodes) {
		nd1 = thisNd->value;
		neighborNd = nodalConn[nd1].getFirst();
		while(neighborNd && sinceRestart < blockDim) {
			nd2 = neighborNd->value;
			if(nodeInserted[nd2] == 0) {
				orderedNds.addEntry(nd2);
				nodeInserted[nd2] = 1;
				sinceRestart++;
			}
			neighborNd = neighborNd->next;
		}
		if(sinceRestart >= blockDim || !thisNd->next) {
			thisNd = orderedNds.getLast();
			nd1 = thisNd->value;
			nodeArray[nd1].ptr->getCrd(p1);
			minDist = 1.0e+100;
			for (i1 = 0; i1 < numNodes; i1++) {
				if(nodeInserted[i1] == 0) {
					nodeArray[i1].ptr->getCrd(p2);
					dist = getDist(p1, p2);
					if(dist < minDist) {
						minDist = dist;
						minNd = i1;
					}
				}
			}
			orderedNds.addEntry(minNd);
			if(sinceRestart >= blockDim) {
				sinceRestart = 0;
			}
		}
		thisNd = thisNd->next;
	}
	
	// Update the global degree of freedom indexes for the nodes
	
    i2 = 0;
	i3 = 0;
	thisNd = orderedNds.getFirst();
	Node *thisPt;
	while(thisNd) {
		nd1 = thisNd->value;
		thisPt = nodeArray[nd1].ptr;
		thisPt->setSortedRank(i3);
		i3++;
		thisPt->setDofIndex(0,i2);
		i2++;
		thisPt->setDofIndex(1,i2);
		i2++;
		thisPt->setDofIndex(2,i2);
		i2++;
		if(thisPt->getNumDof() == 6) {
			thisPt->setDofIndex(3,i2);
			i2++;
			thisPt->setDofIndex(4,i2);
			i2++;
			thisPt->setDofIndex(5,i2);
			i2++;			
		}
		thisNd = thisNd->next;
	}
	elMatDim = i2;
	elasticMat.setDim(elMatDim);
	thermMat.setDim(nodes.getLength());

	thisEl = elements.getFirst();
	while (thisEl) {
		i3 = thisEl->getNumIntDof();
		if (i3 > 0) {
			thisEl->setIntDofIndex(i2);
			i2 += i3;
		}
		thisEl = thisEl->getNext();
	}
	totGlobDof = i2;
	
	tempV1 = new double[elMatDim];
	tempV2 = new double[elMatDim];
	tempV3 = new double[elMatDim];
	tempV4 = new double[elMatDim];
	tempV5 = new double[elMatDim];
	tempV6 = new double[elMatDim];
	tempD1 = new Doub[elMatDim];

	i3 = nodes.getLength();
	dLdU = new double[totGlobDof];
	dLdV = new double[elMatDim];
	dLdA = new double[elMatDim];
	dLdT = new double[i3];
	dLdTdot = new double[i3];
	uAdj = new double[elMatDim];
	vAdj = new double[elMatDim];
	aAdj = new double[elMatDim];
	tAdj = new double[i3];
	tdotAdj = new double[i3];

	dRudD = new DiffDoub[elMatDim];
	dRtdD = new DiffDoub[i3];

	elInD = new int[elements.getLength()];

	i3 = designVars.getLength();
	if (i3 > 0) {
		dLdD = new double[i3];
	}
	else {
		dLdD = nullptr;
	}

    for (i1 = 0; i1 < numNodes; i1++) {
		nodalConn[i1].destroy();
	}

    orderedNds.destroy();	
	
	delete[] nodalConn;
	delete[] nodeInserted;
	return;
}

void Model::buildConstraintMats() {
	Constraint* thisConst = elasticConst.getFirst();
	while (thisConst) {
		thisConst->buildMat(nodeArray);
		thisConst = thisConst->getNext();
	}
	thisConst = thermalConst.getFirst();
	while (thisConst) {
		thisConst->buildMat(nodeArray);
		thisConst = thisConst->getNext();
	}
	return;
}

void Model::updateReference() {
	int i1;
	// Set the section pointers for all elements and material pointers for all sections
	string elSet;
	Set *esPtr;
	Section *thisSec = sections.getFirst();
	string matName;
	Material *thisMat;
	Layer *thisLay;
	IntListEnt *thisInt;
	while(thisSec) {
		elSet = thisSec->getElset();
		i1 = esMap.at(elSet);
		esPtr = esArray[i1].ptr;
		thisInt = esPtr->getFirstEntry();
		while(thisInt) {
			elementArray[thisInt->value].ptr->setSectPtr(thisSec);
			thisInt = thisInt->next;
		}
		thisMat = materials.getFirst();
		while(thisMat) {
			matName = thisMat->getName();
			if(matName == thisSec->getMaterial()) {
				thisSec->setMatPtr(thisMat);
			}
			thisLay = thisSec->getFirstLayer();
			while(thisLay) {
				if(matName == thisLay->getMatName()) {
					thisLay->setMatPtr(thisMat);
				}
				thisLay = thisLay->getNext();
			}
			thisMat = thisMat->getNext();
		}
		thisSec = thisSec->getNext();
	}
	
	// Set node & element set pointers in loads, constraints, design variables and objectives
	string ndSet;
	Set *nsPtr;
	Load *thisLoad = elasticLoads.getFirst();
	while(thisLoad) {
		try {
			ndSet = thisLoad->getNodeSet();
			i1 = nsMap.at(ndSet);
			thisLoad->setNdSetPtr(nsArray[i1].ptr);
		}
		catch (...) {
			elSet = thisLoad->getElSet();
			i1 = esMap.at(elSet);
			thisLoad->setElSetPtr(esArray[i1].ptr);
		}
		thisLoad = thisLoad->getNext();
	}
	thisLoad = thermalLoads.getFirst();
	while (thisLoad) {
		try {
			ndSet = thisLoad->getNodeSet();
			i1 = nsMap.at(ndSet);
			thisLoad->setNdSetPtr(nsArray[i1].ptr);
		}
		catch (...) {
			elSet = thisLoad->getElSet();
			i1 = esMap.at(elSet);
			thisLoad->setElSetPtr(esArray[i1].ptr);
		}
		thisLoad = thisLoad->getNext();
	}

	Constraint* thisConst = elasticConst.getFirst();
	ConstraintTerm* thisCTerm;
	while (thisConst) {
		thisCTerm = thisConst->getFirst();
		while (thisCTerm) {
			ndSet = thisCTerm->getSetName();
			i1 = nsMap.at(ndSet);
			thisCTerm->setNsPtr(nsArray[i1].ptr);
			thisCTerm = thisCTerm->getNext();
		}
		thisConst = thisConst->getNext();
	}
	thisConst = thermalConst.getFirst();
	while (thisConst) {
		thisCTerm = thisConst->getFirst();
		while (thisCTerm) {
			ndSet = thisCTerm->getSetName();
			i1 = nsMap.at(ndSet);
			thisCTerm->setNsPtr(nsArray[i1].ptr);
			thisCTerm = thisCTerm->getNext();
		}
		thisConst = thisConst->getNext();
	}

	DesignVariable* thisDV = designVars.getFirst();
	while (thisDV) {
		try {
			ndSet = thisDV->getNdSet();
			i1 = nsMap.at(ndSet);
			thisDV->setNdsetPtr(nsArray[i1].ptr);
		}
		catch (...) {
			elSet = thisDV->getElSet();
			i1 = esMap.at(elSet);
			thisDV->setElsetPtr(esArray[i1].ptr);
		}
		thisDV = thisDV->getNext();
	}

	ObjectiveTerm* thisTerm = objective.getFirst();
	while (thisTerm) {
		try {
			elSet = thisTerm->getElsetName();
			i1 = esMap.at(elSet);
			thisTerm->setElsetPtr(esArray[i1].ptr);
		}
		catch (...) {
			ndSet = thisTerm->getNdsetName();
			i1 = nsMap.at(ndSet);
			thisTerm->setNdsetPtr(nsArray[i1].ptr);
		}
		thisTerm = thisTerm->getNext();
	}
	
	// Build DV reference list for nodes and elements
	DoubList *coefs;
	DoubListEnt *thisDoub;
	int coefLen;
	double constCoef;
	int DVi = 0;
	thisDV = designVars.getFirst();
	while(thisDV) {
		coefs = thisDV->getCoefs();
		coefLen = coefs->getLength();
		esPtr = thisDV->getElsetPtr();
		if(esPtr) {
			if(coefLen < 2) {
				if(coefLen == 0) {
					constCoef = 1.0;
				} else {
					thisDoub = coefs->getFirst();
					constCoef = thisDoub->value;
				}
				thisInt = esPtr->getFirstEntry();
				while(thisInt) {
					elementArray[thisInt->value].ptr->addDesignVariable(DVi,constCoef);
					thisInt = thisInt->next;
				}
			}
			else {
				thisInt = esPtr->getFirstEntry();
				thisDoub = coefs->getFirst();
				while (thisInt && thisDoub) {
					elementArray[thisInt->value].ptr->addDesignVariable(DVi, thisDoub->value);
					thisInt = thisInt->next;
					thisDoub = thisDoub->next;
				}
		    }
		}
		
		nsPtr = thisDV->getNdsetPtr();
		if (nsPtr) {
			if (coefLen < 2) {
				if (coefLen == 0) {
					constCoef = 1.0;
				}
				else {
					thisDoub = coefs->getFirst();
					constCoef = thisDoub->value;
				}
				thisInt = nsPtr->getFirstEntry();
				while (thisInt) {
					nodeArray[thisInt->value].ptr->addDesignVariable(DVi, constCoef);
					thisInt = thisInt->next;
				}
			}
			else {
				thisInt = nsPtr->getFirstEntry();
				thisDoub = coefs->getFirst();
				while (thisInt && thisDoub) {
					nodeArray[thisInt->value].ptr->addDesignVariable(DVi, thisDoub->value);
					thisInt = thisInt->next;
					thisDoub = thisDoub->next;
				}
			}
		}
		thisDV = thisDV->getNext();
		DVi++;
	}
	
	// Build comprehensive element list for each design variable
	
	int elLabel;
	int elNumNds;
	int *elNodes;
	Element *thisEl = elements.getFirst();
	IntList *desVars;
	while(thisEl) {
		elLabel = thisEl->getLabel();
		desVars = thisEl->getDesignVars();
		thisInt = desVars->getFirst();
		while(thisInt) {
			dVarArray[thisInt->value].ptr->addCompEl(elLabel);
			thisEl->addCompDVar(thisInt->value);
			thisInt = thisInt->next;
		}
		elNumNds = thisEl->getNumNds();
		elNodes = thisEl->getNodes();
		for (i1 = 0; i1 < elNumNds; i1++) {
			desVars = nodeArray[elNodes[i1]].ptr->getDesignVars();
			thisInt = desVars->getFirst();
			while(thisInt) {
				thisDV = dVarArray[thisInt->value].ptr;
				if (thisDV->getCategory() == "nodeCoord") {
					dVarArray[thisInt->value].ptr->addCompEl(elLabel);
					thisEl->addCompDVar(thisInt->value);
				}
				thisInt = thisInt->next;
			}
		}
		thisEl = thisEl->getNext();
	}
	
	return;
}

void Model::findSurfaceFaces() {
	int i1;
	int lowNd;
	bool added;
	
	Element *thisEl = elements.getFirst();
	while(thisEl) {
		thisEl->initializeFaces();
		thisEl = thisEl->getNext();
	}
	
	int numNodes = nodes.getLength();
	FacePtList* fLArray = new FacePtList[numNodes];
	thisEl = elements.getFirst();
	Face* thisFc;
	while(thisEl) {
		thisFc = thisEl->getFirstFace();
		while(thisFc) {
			lowNd = thisFc->getLowNd();
			added = fLArray[lowNd].addIfAbsent(thisFc);
			thisFc = thisFc->getNext();
		}
		thisEl = thisEl->getNext();
	}
	
	for (i1 = 0; i1 < numNodes; i1++) {
		fLArray[i1].destroy();
	}
	delete[] fLArray;
	
	return;
}

void Model::analysisPrep(int blockDim) {
	updateReference();
	reorderNodes(blockDim);
	buildConstraintMats();
	findSurfaceFaces();
	
	anPrepRun = true;
	return;
}

void Model::buildElasticAppLoad(double appLd[], double time) {
	// Construct loads from input file
	int i1;
	int numDof;
	int dofInd;
	Node *thisNd;
	Element* thisEl;
	Load *thisLoad = elasticLoads.getFirst();
	string ldType;
	double actTime[2];
	double ndLoad[6];
	Doub ndDVLd[6];
	Set *thisSet;
	IntListEnt *thisEnt;
	DoubStressPrereq pre;

	for (i1 = 0; i1 < elMatDim; i1++) {
		tempD1[i1].setVal(0.0);
	}

	//Loads from model input file
	while(thisLoad) {
		ldType = thisLoad->getType();
		thisLoad->getActTime(actTime);
		if(time >= actTime[0] && time <= actTime[1]) {
			if(ldType == "nodalForce") {
				thisLoad->getLoad(ndLoad);
				thisSet = thisLoad->getNdSetPtr();
				thisEnt = thisSet->getFirstEntry();
				while(thisEnt) {
					thisNd = nodeArray[thisEnt->value].ptr;
					numDof = thisNd->getNumDof();
					for (i1 = 0; i1 < numDof; i1++) {
						dofInd = thisNd->getDofIndex(i1);
						appLd[dofInd]+= ndLoad[i1];
					}
					thisEnt = thisEnt->next;
				}
			}
			else {
				thisSet = thisLoad->getElSetPtr();
				thisEnt = thisSet->getFirstEntry();
				while (thisEnt) {
					thisEl = elementArray[thisEnt->value].ptr;
					thisEl->getStressPrereq(pre, nodeArray, dVarArray);
					thisEl->getAppLoad(tempD1, thisLoad, solveCmd->nonlinearGeom, pre, nodeArray, dVarArray);
					thisEnt = thisEnt->next;
				}
			}
		}
		thisLoad = thisLoad->getNext();
	}

	for (i1 = 0; i1 < elMatDim; i1++) {
		appLd[0] += tempD1[i1].val;
	}
	
	// Design variable dependent loads.
    thisNd = nodes.getFirst();
	while(thisNd) {
		thisNd->getElasticDVLoad(ndDVLd,dVarArray);
		numDof = thisNd->getNumDof();
		for (i1 = 0; i1 < numDof; i1++) {
			dofInd = thisNd->getDofIndex(i1);
			appLd[dofInd]+= ndDVLd[i1].val;
		}
		thisNd = thisNd->getNext();
	}

	pre.destroy();
	
	return;
}

void Model::buildThermalAppLoad(double appLd[], double time) {
	// Construct loads from input file
	int i1;
	int totNodes;
	int numDof;
	int dofInd;
	Node* thisNd;
	Element* thisEl;
	Load* thisLoad = thermalLoads.getFirst();
	string ldType;
	double actTime[2];
	double ndLoad[6];
	Doub ndDVLd;
	Set* thisSet;
	IntListEnt* thisEnt;
	DoubStressPrereq pre;

	totNodes = nodes.getLength();
	for (i1 = 0; i1 < totNodes; i1++) {
		tempD1[i1].setVal(0.0);
	}

	//Loads from model input file
	while (thisLoad) {
		ldType = thisLoad->getType();
		thisLoad->getActTime(actTime);
		if (time >= actTime[0] && time <= actTime[1]) {
			if (ldType == "nodalHeatGen") {
				thisLoad->getLoad(ndLoad);
				thisSet = thisLoad->getNdSetPtr();
				thisEnt = thisSet->getFirstEntry();
				while (thisEnt) {
					thisNd = nodeArray[thisEnt->value].ptr;
					dofInd = thisNd->getSortedRank();
					appLd[dofInd] += ndLoad[0];
					thisEnt = thisEnt->next;
				}
			}
			else {
				thisSet = thisLoad->getElSetPtr();
				thisEnt = thisSet->getFirstEntry();
				while (thisEnt) {
					thisEl = elementArray[thisEnt->value].ptr;
					thisEl->getStressPrereq(pre, nodeArray, dVarArray);
					thisEl->getAppThermLoad(tempD1, thisLoad, pre, nodeArray, dVarArray);
					thisEnt = thisEnt->next;
				}
			}
		}
		thisLoad = thisLoad->getNext();
	}

	for (i1 = 0; i1 < totNodes; i1++) {
		appLd[i1] += tempD1[i1].val;
	}

	// Design variable dependent loads.
	thisNd = nodes.getFirst();
	while (thisNd) {
		thisNd->getThermalDVLoad(ndDVLd, dVarArray);
		dofInd = thisNd->getSortedRank();
		appLd[dofInd] += ndDVLd.val;
		thisNd = thisNd->getNext();
	}

	pre.destroy();

	return;
}

void Model::buildElasticSolnLoad(double solnLd[], bool buildMat) {
	int i1;
	Element *thisEl;
	DoubStressPrereq pre;
	
	for (i1 = 0; i1 < elMatDim; i1++) {
		tempD1[i1].setVal(0.0);
	}
	
	if(buildMat) {
		elasticMat.zeroAll();
	}
	
	thisEl = elements.getFirst();
	while(thisEl) {
		thisEl->getStressPrereq(pre, nodeArray, dVarArray);
		thisEl->getRu(tempD1, elasticMat, buildMat, solveCmd, pre, nodeArray, dVarArray);
		thisEl = thisEl->getNext();
	}
	
	for (i1 = 0; i1 < elMatDim; i1++) {
		solnLd[i1]-= tempD1[i1].val;
	}

	pre.destroy();
	
	return;
}

void Model::buildThermalSolnLoad(double solnLd[], bool buildMat) {
	int i1;
	int numNodes;
	Element* thisEl;
	DoubStressPrereq pre;

	numNodes = nodes.getLength();
	for (i1 = 0; i1 < numNodes; i1++) {
		tempD1[i1].setVal(0.0);
	}

	if (buildMat) {
		thermMat.zeroAll();
	}

	thisEl = elements.getFirst();
	while (thisEl) {
		thisEl->getStressPrereq(pre, nodeArray, dVarArray);
		thisEl->getRt(tempD1, thermMat, buildMat, solveCmd, pre, nodeArray);
		thisEl = thisEl->getNext();
	}

	for (i1 = 0; i1 < numNodes; i1++) {
		solnLd[i1] -= tempD1[i1].val;
	}

	pre.destroy();

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

void Model::buildElasticConstLoad(double constLd[]) {
	int i1;
	double ndDisp[6];
	int ndDof;
	int globInd;
	Node* thisNd = nodes.getFirst();
	while (thisNd) {
		thisNd->getDisp(ndDisp);
		ndDof = thisNd->getNumDof();
		for (i1 = 0; i1 < ndDof; i1++) {
			globInd = thisNd->getDofIndex(i1);
			tempV2[globInd] = ndDisp[i1];
		}
		thisNd = thisNd->getNext();
	}
	elasticConst.getTotalLoad(constLd, tempV2, tempV3, elMatDim);
	return;
}

void Model::buildThermalConstLoad(double constLd[]) {
	int i1;
	int numNodes;
	double ndTemp;
	int globInd;
	Node* thisNd = nodes.getFirst();
	while (thisNd) {
		ndTemp = thisNd->getTemperature();
		globInd = thisNd->getSortedRank();
		tempV2[globInd] = ndTemp;
		thisNd = thisNd->getNext();
	}
	numNodes = nodes.getLength();
	thermalConst.getTotalLoad(constLd, tempV2, tempV3, numNodes);
	return;
}

void Model::solveStep(JobCommand *cmd, double time, double appLdFact) {
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
	
	if(cmd->thermal) {
		numNodes = nodes.getLength();
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
		if (cmd->solverMethod == "direct") {
			thermLT.ldlSolve(tempV2, tempV1);
		}
		thisNd = nodes.getFirst();
		while (thisNd) {
			dofInd = thisNd->getSortedRank();
			thisNd->addToTemperature(tempV2[dofInd]);
			if (cmd->dynamic) {
				thisNd->updateTdot(cmd->newmarkGamma, cmd->timeStep);
			}
			thisNd = thisNd->getNext();
		}
	}
		
	if(cmd->elastic) {
		if(cmd->nonlinearGeom) {
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
			buildElasticSolnLoad(tempV1,cmd->nonlinearGeom);
			if (!elasticScaled) {
				scaleElasticConst();
			}
			buildElasticConstLoad(tempV1);
			thisEl = elements.getFirst();
			while(thisEl) {
				if(thisEl->getNumIntDof() > 0) {
				    thisEl->updateExternal(tempV1,1,nodeArray);
				}
				thisEl = thisEl->getNext();
			}
			if(cmd->solverMethod == "direct") {
				if(cmd->nonlinearGeom) {
					elasticLT.populateFromSparseMat(elasticMat,elasticConst);
					elasticLT.ldlFactor();
				}
				elasticLT.ldlSolve(tempV2,tempV1);
			}
			thisNd = nodes.getFirst();
			while(thisNd) {
				ndof = thisNd->getNumDof();
				for (i1 = 0; i1 < ndof; i1++) {
					dofInd = thisNd->getDofIndex(i1);
					ndDelDisp[i1] = tempV2[dofInd];
				}
				thisNd->addToDisplacement(ndDelDisp);
				if(cmd->dynamic) {
					thisNd->updateVelAcc(cmd->newmarkBeta,cmd->newmarkGamma,cmd->timeStep);
				}
				thisNd = thisNd->getNext();
			}
			thisEl = elements.getFirst();
			while(thisEl) {
				if(thisEl->getNumIntDof() > 0) {
				    thisEl->updateInternal(tempV2,1,nodeArray);
				}
				thisEl = thisEl->getNext();
			}
			dUnorm = 0.0;
			for (i1 = 0; i1 < elMatDim; i1++) {
				absdU = abs(tempV2[i1]);
				if(absdU > dUnorm) {
					dUnorm = absdU;
				}
			}
			i2++;
		}
	}

	
	return;
}

void Model::solve(JobCommand *cmd) {
	int i1;
	int i2;
	double appLdFact;
	double time;
	int ldSteps = cmd->loadRampSteps;
	double c1 = 1.0/(ldSteps*ldSteps);
	Node *thisNd;
	Element* thisEl;
	double zeroAr[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	
	if(!anPrepRun) {
		analysisPrep(cmd->solverBlockDim);
	}
	
	if(cmd->thermal) {
		thisNd = nodes.getFirst();
		while (thisNd) {
			thisNd->initializeTemp();
			if (cmd->dynamic) {
				thisNd->updateTdot(cmd->newmarkGamma, cmd->timeStep);
			}
			thisNd = thisNd->getNext();
		}
		if (!thermLT.isAllocated()) {
			buildThermalSolnLoad(tempV1, true);
			thermLT.allocateFromSparseMat(thermMat, thermalConst, cmd->solverBandwidth);
		}
		if (!thermScaled) {
			scaleThermalConst();
		}
		thermLT.populateFromSparseMat(thermMat, thermalConst);
		thermLT.ldlFactor();
	}
	
	if(cmd->elastic) {
		thisNd = nodes.getFirst();
		while(thisNd) {
			thisNd->initializeDisp();
			if (cmd->dynamic) {
				thisNd->updateVelAcc(cmd->newmarkBeta, cmd->newmarkGamma, cmd->timeStep);
			}
			thisNd = thisNd->getNext();
		}
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->setIntDisp(zeroAr);
			thisEl->setIntPrevDisp(zeroAr);
			thisEl = thisEl->getNext();
		}
		if(!elasticLT.isAllocated()) {
			buildElasticSolnLoad(tempV1,true);
			elasticLT.allocateFromSparseMat(elasticMat,elasticConst,cmd->solverBandwidth);
		}
		if(!cmd->nonlinearGeom) {
			if (!elasticScaled) {
				scaleElasticConst();
			}
			elasticLT.populateFromSparseMat(elasticMat, elasticConst);
			elasticLT.ldlFactor();
		}
	}
	
	if(cmd->dynamic) {
		if (cmd->saveSolnHist) {
			writeTimeStepSoln(0);
		}
		time = 0.0;
		i1 = 1;
		while (time < cmd->simPeriod) {
			solveStep(cmd, time, 1.0);
			thisNd = nodes.getFirst();
			while (thisNd) {
				if (cmd->thermal) {
					thisNd->advanceTemp();
					thisNd->updateTdot(cmd->newmarkGamma,cmd->timeStep);
				}
				if (cmd->elastic) {
					thisNd->advanceDisp();
					thisNd->updateVelAcc(cmd->newmarkBeta,cmd->newmarkGamma,cmd->timeStep);
				}
				thisNd = thisNd->getNext();
			}
			if (cmd->elastic) {
				thisEl = elements.getFirst();
				while (thisEl) {
					thisEl->advanceIntDisp();
					thisEl = thisEl->getNext();
				}
			}
			if (cmd->saveSolnHist) {
				writeTimeStepSoln(i1);
				timeStepsSaved = i1;
				i1++;
			}
			time += cmd->timeStep;
		}
	} else {
		if(cmd->nonlinearGeom) {
			for (i1 = 0; i1 < ldSteps; i1++) {
				i2 = ldSteps - i1 - 1;
				appLdFact = 1.0 - c1*i2*i2;
				solveStep(cmd,cmd->staticLoadTime,appLdFact);
			}
		} else {
			appLdFact = 1.0;
			solveStep(cmd,cmd->staticLoadTime,appLdFact);
		}
	}
	
	return;
}

void Model::eigenSolve(JobCommand* cmd) {
	int i1;
	int i2;
	int i3;
	Element* thisEl;
	Node* thisNd;
	Doub Rvec[33];
	double Rv2[33];
	double dRdA[2];
	double zeros[9] = {0.0, 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0};
	double vKv;
	DoubStressPrereq pre;

	if (!anPrepRun) {
		analysisPrep(cmd->solverBlockDim);
	}
	if (!eigVecs) {
		i1 = cmd->numModes;
		eigVals = new double[i1];
		loadFact = new double[i1];
		i1 *= elMatDim;
		eigVecs = new double[i1];
		diagMass = new double[elMatDim];
	}
	if (cmd->type == "buckling") {
		for (i1 = 0; i1 < elMatDim; i1++) {
			diagMass[i1] = 1.0;
		}
	}
	else {
		for (i1 = 0; i1 < elMatDim; i1++) {
			diagMass[i1] = 0.0;
		}
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->getStressPrereq(pre, nodeArray, dVarArray);
			i2 = thisEl->getNumNds() * thisEl->getDofPerNd();
			for (i1 = 0; i1 < i2; i1++) {
				pre.globAcc[i1].setVal(1.0);
			}
			thisEl->getRum(Rvec, dRdA, false, true, solveCmd->nonlinearGeom, pre, nodeArray, dVarArray);
			for (i1 = 0; i1 < i2; i1++) {
				Rv2[i1] = Rvec[i1].val;
			}
			thisEl->addToGlobVec(Rv2, diagMass, false, false, nodeArray);
			thisEl = thisEl->getNext();
		}
	}
	cmd->nonlinearGeom = solveCmd->nonlinearGeom;
	cmd->dynamic = solveCmd->dynamic;
	solveCmd->nonlinearGeom = true;
	solveCmd->dynamic = false;
	buildElasticSolnLoad(tempV1, true);
	if (cmd->solverMethod == "direct") {
		elasticLT.populateFromSparseMat(elasticMat, elasticConst);
		elasticLT.ldlFactor();
		eigenSparseDirect(eigVals, eigVecs, cmd->numModes, elasticLT, diagMass, elMatDim);
	}
	else {

	}

	if (cmd->type == "buckling") {
		writeTimeStepSoln(-1);
		thisNd = nodes.getFirst();
		while (thisNd) {
			thisNd->advanceDisp();
			thisNd->setDisplacement(zeros);
			thisNd = thisNd->getNext();
		}
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->advanceIntDisp();
			thisEl->setIntDisp(zeros);
			thisEl = thisEl->getNext();
		}
		buildElasticSolnLoad(tempV1, true);
		i2 = 0;
		for (i1 = 0; i1 < cmd->numModes; i1++) {
			elasticMat.vectorMultiply(tempV1, &eigVecs[i2],false);
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
		thisNd = nodes.getFirst();
		while (thisNd) {
			thisNd->backstepDisp();
			thisNd = thisNd->getNext();
		}
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->backstepIntDisp();
			thisEl = thisEl->getNext();
		}
		readTimeStepSoln(-1);
	}
	else {
		for (i1 = 0; cmd->numModes; i1++) {
			try {
				loadFact[i1] = 0.159154943091895335 * sqrt(eigVals[i1]);
			}
			catch (...) {
				loadFact[i1] = -1.0;
			}
		}
	}

	solveCmd->nonlinearGeom = cmd->nonlinearGeom;
	solveCmd->dynamic = cmd->dynamic;

	pre.destroy();

	return;
}

//void Model::backupElastic() {
//	int i1;
//	int i2;
//	int numIntDof;
//	Node* thisNd = nodes.getFirst();
//	i1 = 0;
//	i2 = 0;
//	double zeros[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//	while (thisNd) {
//		thisNd->getPrevDisp(&tempV1[i1]);
//		thisNd->getPrevVel(&tempV2[i1]);
//		thisNd->getPrevAcc(&tempV3[i1]);
//		tempV4[i2] = thisNd->getTemperature();
//		thisNd->getDisp(&tempV5[i1]);
//		thisNd->setPrevDisp(zeros);
//		thisNd->setPrevVel(zeros);
//		thisNd->setPrevAcc(zeros);
//		thisNd->setTemperature(zeros[0]);
//		thisNd->setDisplacement(zeros);
//		thisNd = thisNd->getNext();
//		i1 += thisNd->getNumDof();
//		i2++;
//	}
//
//	Element* thisEl = elements.getFirst();
//	i1 = 0;
//	while (thisEl) {
//		numIntDof = thisEl->getNumIntDof();
//		if (numIntDof > 0) {
//			thisEl->getIntDisp(&tempV6[i1]);
//			i1 += numIntDof;
//		}
//		thisEl = thisEl->getNext();
//	}
//	return;
//}
//
//void Model::restoreElastic() {
//	int i1;
//	int i2;
//	Node* thisNd = nodes.getFirst();
//	i1 = 0;
//	i2 = 0;
//	double zeros[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//	while (thisNd) {
//		thisNd->setPrevDisp(&tempV1[i1]);
//		thisNd->setPrevVel(&tempV2[i1]);
//		thisNd->setPrevAcc(&tempV3[i1]);
//		thisNd->setTemperature(tempV4[i2]);
//		thisNd->setDisplacement(&tempV5[i1]);
//		thisNd = thisNd->getNext();
//		i1 += thisNd->getNumDof();
//		i2++;
//	}
//	return;
//}

void Model::setSolnToMode(string field, int mode, double maxVal) {
	int i1;
	int i2;
	int modeSt;
	int nDof;
	int globInd;
	string dispFields = "displacement velocity acceleration";
	string erStr;
	Node* thisNd;
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
		thisNd = nodes.getFirst();
		while (thisNd) {
			nDof = thisNd->getNumDof();
			for (i2 = 0; i2 < nDof; i2++) {
				globInd = thisNd->getDofIndex(i2);
				ndDat[i2] = scaleFact * eigVecs[modeSt + globInd];
			}
			if (field == "displacement") {
				thisNd->setDisplacement(ndDat);
			}
			else if (field == "velocity") {
				thisNd->setVelocity(ndDat);
			}
			else {
				thisNd->setAcceleration(ndDat);
			}
			thisNd = thisNd->getNext();
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
	double dt = solveCmd->timeStep;
	double bet = solveCmd->newmarkBeta;
	double gam = solveCmd->newmarkGamma;
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

	Doub Rvec[30];
	double Mmat[900];
	double Dmat[900];
	double elAdj[30];
	double eldLdT[10];
	double eldLdTdot[10];
	double eldLdU[30];
	double eldLdV[30];
	double eldLdA[30];
	DoubStressPrereq pre;

	Element* thisEl;

	if (solveCmd->thermal) {
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->getStressPrereq(pre, nodeArray, dVarArray);
			thisEl->getRtm(Rvec, Mmat, true, true, pre);
			thisEl->getElVec(elAdj, tAdj, true, false, nodeArray);
			numNds = thisEl->getNumNds();
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
			thisEl->addToGlobVec(eldLdT, dLdT, true, false, nodeArray);
			thisEl->addToGlobVec(eldLdTdot, dLdTdot, true, false, nodeArray);
			thisEl = thisEl->getNext();
		}
		i2 = nodes.getLength();
		for (i1 = 0; i1 < i2; i1++) {
			dLdT[i1] += tdotAdj[i1];
			dLdTdot[i1] += c3 * tdotAdj[i1];
		}
	}

	if (solveCmd->elastic) {
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->getStressPrereq(pre, nodeArray, dVarArray);
			thisEl->getRum(Rvec, Mmat, true, true, solveCmd->nonlinearGeom, pre, nodeArray, dVarArray);
			thisEl->getRud(Rvec, Dmat, true, solveCmd, pre, nodeArray, dVarArray);
			thisEl->getElVec(elAdj, uAdj, false, false, nodeArray);
			ndDof = thisEl->getNumNds() * thisEl->getDofPerNd();
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
			thisEl->addToGlobVec(eldLdU, dLdU, false, false, nodeArray);
			thisEl->addToGlobVec(eldLdA, dLdA, false, false, nodeArray);
			thisEl->addToGlobVec(eldLdV, dLdV, false, false, nodeArray);
			thisEl = thisEl->getNext();
		}
		for (i1 = 0; i1 < elMatDim; i1++) {
			dLdU[i1] += c1 * aAdj[i1];
			dLdA[i1] += (c2 * aAdj[i1] + c3 * vAdj[i1]);
			dLdV[i1] += c4 * aAdj[i1] + vAdj[i1];
		}
	}
	
	pre.destroy();

	return;
}

void Model::solveForAdjoint() {
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
	DoubStressPrereq pre;
	Doub Rvec[33];
	double dRdU[1089];
	double dRdT[330];
	double eldLdT[10];
	double elAdj[33];
	double* intAdj;

	objective.calculatedLdU(dLdU, dLdV, dLdA, dLdT, dLdTdot, solveCmd->staticLoadTime, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
    //
	cout << "dLdU:" << endl;
	for (i1 = 0; i1 < totGlobDof; i1++) {
		cout << dLdU[i1] << endl;
	}
	//
	if (solveCmd->elastic) {
		if (solveCmd->dynamic) {
			c1 = solveCmd->timeStep * solveCmd->newmarkGamma;
			c2 = solveCmd->timeStep;
			c2 = 1.0 / (c2 * c2 * (solveCmd->newmarkBeta - solveCmd->newmarkGamma));
			for (i1 = 0; i1 < elMatDim; i1++) {
				vAdj[i1] = dLdV[i1];
				aAdj[i1] = (dLdA[i1] + c1 * vAdj[i1]);
				dLdU[i1] -= c2 * aAdj[i1];
			}
		}
		if (solveCmd->nonlinearGeom) {
			buildElasticSolnLoad(tempV1, true);
			elasticLT.populateFromSparseMat(elasticMat, elasticConst);
			elasticLT.ldlFactor();
		}
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->setIntdLdU(dLdU);
			thisEl->updateExternal(dLdU, 0, nodeArray);
			thisEl = thisEl->getNext();
		}
		if (solveCmd->solverMethod == "direct") {
			elasticLT.ldlSolve(uAdj, dLdU);
		}
		thisEl = elements.getFirst();
		while (thisEl) {
			thisEl->updateInternal(uAdj, 0, nodeArray);
			thisEl = thisEl->getNext();
		}
	}

	if (solveCmd->thermal) {
		if (solveCmd->elastic) {
			thisEl = elements.getFirst();
			while (thisEl) {
				thisEl->getStressPrereq(pre, nodeArray, dVarArray);
				thisEl->getRuk(Rvec, dRdU, dRdT, true, solveCmd->nonlinearGeom, pre, nodeArray, dVarArray);
				thisEl->getElVec(elAdj, uAdj, false, false, nodeArray);
				elNumNds = thisEl->getNumNds();
				elDofPerNd = thisEl->getDofPerNd();
				elIntDof = thisEl->getNumIntDof();
				elNdDof = elNumNds * elDofPerNd;
				intAdj = thisEl->getIntAdj();
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
				thisEl->addToGlobVec(eldLdT, dLdT, true, false, nodeArray);
				thisEl = thisEl->getNext();
			}
		}
		totNodes = nodes.getLength();
		if (solveCmd->dynamic) {
			c3 = -1.0 / (solveCmd->timeStep * solveCmd->newmarkGamma);
			for (i1 = 0; i1 < totNodes; i1++) {
				tdotAdj[i1] = c3 * dLdTdot[i1];
				dLdT[i1] -= tdotAdj[i1];
			}
		}
		if (solveCmd->solverMethod == "direct") {
			thermLT.ldlSolve(tAdj, dLdT);
		}
	}
	
	pre.destroy();
	
	return;
}

void Model::dRthermaldD(int dVarNum) {
	int i1;
	int totNodes;
	int numDof;
	int globInd;
	Doub dvVal;
	DiffDoubStressPrereq pre;
	DesignVariable* thisDV;
	IntListEnt* thisEnt;
	Element* thisElPt;
	Load* thisLd;

	thisDV = dVarArray[dVarNum].ptr;
	thisDV->getValue(dvVal);
	thisDV->setDiffVal(dvVal.val, 1.0);

	totNodes = nodes.getLength();
	for (i1 = 0; i1 < totNodes; i1++) {
		dRtdD[i1].setVal(0.0);
	}

	for (i1 = 0; i1 < elements.getLength(); i1++) {
		elInD[i1] = 0;
	}

	thisEnt = thisDV->getFirstEl();
	while (thisEnt) {
		elInD[thisEnt->value] = 1;
		thisEnt = thisEnt->next;
	}
	
	// Applied contribution
	thisLd = thermalLoads.getFirst();
	while (thisLd) {
		if (thisLd->getType() != "nodalHeatGen") {
			thisEnt = thisLd->getElSetPtr()->getFirstEntry();
			while (thisEnt) {
				i1 = thisEnt->value;
				if (elInD[i1]) {
					thisElPt = elementArray[i1].ptr;
					thisElPt->getStressPrereq(pre, nodeArray, dVarArray);
					thisElPt->getAppThermLoad(dRtdD, thisLd, pre, nodeArray, dVarArray);
				}
				thisEnt = thisEnt->next;
			}
		}
		thisLd = thisLd->getNext();
	}

	for (i1 = 0; i1 < totNodes; i1++) {
		dRtdD[i1].neg();
	}

	// Solution-dependent contribution of load
	thisEnt = thisDV->getFirstEl();
	while (thisEnt) {
		thisElPt = elementArray[thisEnt->value].ptr;
		thisElPt->getStressPrereq(pre, nodeArray, dVarArray);
		thisElPt->getRt(dRtdD,thermMat,false,solveCmd,pre,nodeArray);
		thisEnt = thisEnt->next;
	}


	// Design variable dependent contribution
	thisEnt = thisDV->getFirstNd();
	Node* thisNdPt;
	DiffDoub ndLd;
	while (thisEnt) {
		thisNdPt = nodeArray[thisEnt->value].ptr;
		thisNdPt->getThermalDVLoad(ndLd, dVarArray);
		globInd = thisNdPt->getSortedRank();
		dRtdD[globInd].sub(ndLd);
		thisEnt = thisEnt->next;
	}

	thisDV->setDiffVal(dvVal.val, 0.0);
	
	pre.destroy();
	return;
}

void Model::dRelasticdD(int dVarNum) {
	int i1;
	int numDof;
	int globInd;
	Doub dvVal;
	DiffDoubStressPrereq pre;
	DesignVariable* thisDV;
	IntListEnt* thisEnt;
	Element* thisElPt;
	Load* thisLd;

	thisDV = dVarArray[dVarNum].ptr;
	thisDV->getValue(dvVal);
	thisDV->setDiffVal(dvVal.val, 1.0);

	for (i1 = 0; i1 < elMatDim; i1++) {
		dRudD[i1].setVal(0.0);
	}

	for (i1 = 0; i1 < elements.getLength(); i1++) {
		elInD[i1] = 0;
	}

	thisEnt = thisDV->getFirstEl();
	while (thisEnt) {
		elInD[thisEnt->value] = 1;
		thisEnt = thisEnt->next;
	}
	
	// Applied contribution

	thisLd = elasticLoads.getFirst();
	while (thisLd) {
		if (thisLd->getType() != "nodalForce") {
			thisEnt = thisLd->getElSetPtr()->getFirstEntry();
			while (thisEnt) {
				i1 = thisEnt->value;
				if (elInD[i1] > 0) {
					thisElPt = elementArray[i1].ptr;
					thisElPt->getStressPrereq(pre, nodeArray, dVarArray);
					thisElPt->getAppLoad(dRudD, thisLd, solveCmd->nonlinearGeom, pre, nodeArray, dVarArray);
				}
				thisEnt = thisEnt->next;
			}
		}
		thisLd = thisLd->getNext();
	}

	for (i1 = 0; i1 < elMatDim; i1++) {
		dRudD[i1].neg();
	}

	// Solution-dependent contribution of load
	thisEnt = thisDV->getFirstEl();
	while (thisEnt) {
		thisElPt = elementArray[thisEnt->value].ptr;
		thisElPt->getStressPrereq(pre, nodeArray, dVarArray);
		thisElPt->getRu(dRudD, elasticMat, false, solveCmd, pre, nodeArray, dVarArray);
		thisEnt = thisEnt->next;
	}


	// Design variable dependent contribution
	thisEnt = thisDV->getFirstNd();
	Node* thisNdPt;
	DiffDoub ndLd[6];
	while (thisEnt) {
		thisNdPt = nodeArray[thisEnt->value].ptr;
		thisNdPt->getElasticDVLoad(ndLd, dVarArray);
		numDof = thisNdPt->getNumDof();
		for (i1 = 0; i1 < numDof; i1++) {
			globInd = thisNdPt->getDofIndex(i1);
			ndLd[i1].neg();
			dRudD[globInd].add(ndLd[i1]);
		}
		thisEnt = thisEnt->next;
	}

	thisDV->setDiffVal(dvVal.val, 0.0);
	
	pre.destroy();

	return;
}

void Model::getObjGradient() {
	int i1;
	int i2;
	int i3;
	double time;
	int numDV = designVars.getLength();
	Element* thisEl;
	Node* thisNd;

	objective.clearValues();
	for (i1 = 0; i1 < numDV; i1++) {
		dLdD[i1] = 0.0;
	}

	if (solveCmd->dynamic) {
		time = solveCmd->timeStep * timeStepsSaved;
		i1 = timeStepsSaved;
		readTimeStepSoln(i1);
		while (i1 > 0) {
			for (i2 = 0; i2 < nodes.getLength(); i2++) {
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
			thisNd = nodes.getFirst();
			while (thisNd) {
				if (solveCmd->elastic) {
					thisNd->backstepDisp();
				}
				if (solveCmd->thermal) {
					thisNd->backstepTemp();
				}
				thisNd = thisNd->getNext();
			}
			thisEl = elements.getFirst();
			while (thisEl) {
				thisEl->backstepIntDisp();
				thisEl = thisEl->getNext();
			}
			readTimeStepSoln(i1 - 1);
			objective.calculateTerms(time, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
			solveForAdjoint();
			objective.calculatedLdD(dLdD, time, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
			for (i2 = 0; i2 < numDV; i2++) {
				if (solveCmd->thermal) {
					dRthermaldD(i2);
					for (i3 = 0; i3 < nodes.getLength(); i3++) {
						dLdD[i2] -= tAdj[i3] * dRtdD[i3].dval;
					}
				}
				if (solveCmd->elastic) {
					dRelasticdD(i2);
					for (i3 = 0; i3 < elMatDim; i3++) {
						dLdD[i2] -= uAdj[i3] * dRudD[i3].dval;
					}
					thisEl = elements.getFirst();
					while (thisEl) {
						dLdD[i2] -= thisEl->getIntAdjdRdD();
						thisEl = thisEl->getNext();
					}
				}
			}
			i1--;
		}
	}
	else {
		objective.calculateTerms(solveCmd->staticLoadTime, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
		for (i2 = 0; i2 < nodes.getLength(); i2++) {
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
		solveForAdjoint();
		//
		cout << "Adjoint:" << endl;
		for (i1 = 0; i1 < elMatDim; i1++) {
			cout << uAdj[i1] << endl;
		}
		//
		objective.calculatedLdD(dLdD, solveCmd->staticLoadTime, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
		//
		cout << "dLdD:" << endl;
		for (i1 = 0; i1 < numDV; i1++) {
			cout << dLdD[i1] << endl;
		}
		//
		for (i1 = 0; i1 < numDV; i1++) {
			if (solveCmd->thermal) {
				dRthermaldD(i1);
				for (i3 = 0; i3 < nodes.getLength(); i3++) {
					dLdD[i1] -= tAdj[i3] * dRtdD[i3].dval;
				}
			}
			if (solveCmd->elastic) {
				dRelasticdD(i1);
				cout << "dRdD:" << endl;
				for (i2 = 0; i2 < elMatDim; i2++) {
					dLdD[i1] -= uAdj[i2] * dRudD[i2].dval;
					cout << dRudD[i2].dval << endl;
				}
				thisEl = elements.getFirst();
				while (thisEl) {
					dLdD[i1] -= thisEl->getIntAdjdRdD();
					thisEl = thisEl->getNext();
				}
			}
		}
	}
	return;
}