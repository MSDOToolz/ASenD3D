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
			nodeArray[nd1].ptr->setNumDof(elDofPerNd);
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
	Set* firstNdSet = nodeSets.getFirst();
	Constraint* thisConst = constraints.getFirst();
	while (thisConst) {
		thisConst->buildMat(firstNdSet, nodeArray);
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
		esPtr = elementSets.getFirst();
		while(esPtr) {
			if(esPtr->getName() == elSet) {
				thisInt = esPtr->getFirstEntry();
				while(thisInt) {
					elementArray[thisInt->value].ptr->setSectPtr(thisSec);
					thisInt = thisInt->next;
				}
			}
			esPtr = esPtr->getNext();
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
	
	// Set node & element set pointers in loads
	string ndSet;
	Set *nsPtr;
	Load *thisLoad = loads.getFirst();
	while(thisLoad) {
		ndSet = thisLoad->getNodeSet();
		nsPtr = nodeSets.getFirst();
		while(nsPtr) {
			if(nsPtr->getName() == ndSet) {
				thisLoad->setNdSetPtr(nsPtr);
			}
			nsPtr = nsPtr->getNext();
		}
		elSet = thisLoad->getElSet();
		esPtr = elementSets.getFirst();
		while(esPtr) {
			if(esPtr->getName() == elSet) {
				thisLoad->setElSetPtr(esPtr);
			}
			esPtr = esPtr->getNext();
		}
		thisLoad = thisLoad->getNext();
	}
	
	// Build DV reference list for nodes and elements
	DesignVariable *thisDV = designVars.getFirst();
	DoubList *coefs;
	DoubListEnt *thisDoub;
	int coefLen;
	double constCoef;
	int DVi = 0;
	while(thisDV) {
		elSet = thisDV->getElSet();
		ndSet = thisDV->getNdSet();
		coefs = thisDV->getCoefs();
		coefLen = coefs->getLength();
		esPtr = elementSets.getFirst();
		while(esPtr) {
			if(esPtr->getName() == elSet) {
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
				} else {
					thisInt = esPtr->getFirstEntry();
					thisDoub = coefs->getFirst();
					while(thisInt && thisDoub) {
						elementArray[thisInt->value].ptr->addDesignVariable(DVi,thisDoub->value);
						thisInt = thisInt->next;
						thisDoub = thisDoub->next;
					}
				}
			}
			esPtr = esPtr->getNext();
		}
		
		nsPtr = nodeSets.getFirst();
		while(nsPtr) {
			if(nsPtr->getName() == ndSet) {
				if(coefLen < 2) {
					if(coefLen == 0) {
						constCoef = 1.0;
					} else {
						thisDoub = coefs->getFirst();
						constCoef = thisDoub->value;
					}
					thisInt = nsPtr->getFirstEntry();
					while(thisInt) {
						nodeArray[thisInt->value].ptr->addDesignVariable(DVi,constCoef);
						thisInt = thisInt->next;
					}
				} else {
					thisInt = nsPtr->getFirstEntry();
					thisDoub = coefs->getFirst();
					while(thisInt && thisDoub) {
						nodeArray[thisInt->value].ptr->addDesignVariable(DVi,thisDoub->value);
						thisInt = thisInt->next;
						thisDoub = thisDoub->next;
					}
				}
				thisDV->setNdsetPtr(nsPtr);
			}
			nsPtr = nsPtr->getNext();
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
				}
				thisInt = thisInt->next;
			}
		}
		thisEl = thisEl->getNext();
	}

	// Set the node and element set pointers in the objective terms

	ObjectiveTerm* thisTerm = objective.getFirst();
	while (thisTerm) {
		elSet = thisTerm->getElsetName();
		esPtr = elementSets.getFirst();
		while (esPtr) {
			if (esPtr->getName() == elSet) {
				thisTerm->setElsetPtr(esPtr);
			}
			esPtr = esPtr->getNext();
		}
		ndSet = thisTerm->getNdsetName();
		nsPtr = nodeSets.getFirst();
		while (nsPtr) {
			if (nsPtr->getName() == ndSet) {
				thisTerm->setNdsetPtr(nsPtr);
			}
			nsPtr = nsPtr->getNext();
		}
		thisTerm = thisTerm->getNext();
	}
	
	return;
}

void Model::findSurfaceFaces() {
	int i1;
	int lowNd;
	
	Element *thisEl = elements.getFirst();
	while(thisEl) {
		thisEl->initializeFaces();
		thisEl = thisEl->getNext();
	}
	
	int numNodes = nodes.getLength();
	FaceList *fLArray = new FaceList[numNodes];
	thisEl = elements.getFirst();
	Face *thisFc;
	while(thisEl) {
		thisFc = thisEl->getFirstFace();
		while(thisFc) {
			lowNd = thisFc->getLowNd();
			fLArray[lowNd].addIfAbsent(thisFc);
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
	Load *thisLoad = loads.getFirst();
	string ldType;
	double actTime[2];
	double ndLoad[6];
	Doub ndDVLd[6];
	Set *thisSet;
	IntListEnt *thisEnt;
	Node *thisNd;

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
		}
		thisLoad = thisLoad->getNext();
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
	
	return;
}

void Model::buildElasticSolnLoad(double solnLd[], bool buildMat, bool dyn, bool nLGeom) {
	int i1;
	Element *thisEl;
	
	for (i1 = 0; i1 < elMatDim; i1++) {
		tempD1[i1].setVal(0.0);
	}
	
	if(buildMat) {
		elasticMat.zeroAll();
	}
	
	thisEl = elements.getFirst();
	while(thisEl) {
		thisEl->getRu(tempD1,elasticMat,buildMat,dyn,nLGeom,nodeArray,dVarArray);
		thisEl = thisEl->getNext();
	}
	
	for (i1 = 0; i1 < elMatDim; i1++) {
		solnLd[i1]-= tempD1[i1].val;
	}
	
	return;
}

void Model::scaleElasticConst() {
	double scaleFact = 1000.0*elasticMat.getMaxAbsVal();
	constraints.scaleElastic(scaleFact);
	elasticScaled = true;
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
	constraints.getTotalLoad(constLd, tempV2, tempV3, elMatDim);
	return;
}

void Model::solveStep(JobCommand *cmd, double time, double appLdFact) {
	int i1;
	int i2;
	int maxNLit;
	double dUnorm;
	double dUtol;
	double absdU;
	int ndof;
	int dofInd;
	double ndDelDisp[6];
	Node *thisNd;
	Element *thisEl;
	
	if(cmd->thermal) {
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
			}
			buildElasticAppLoad(tempV1,time);
			for (i1 = 0; i1 < elMatDim; i1++) {
				tempV1[i1] *= appLdFact;
			}
			buildElasticSolnLoad(tempV1,cmd->nonlinearGeom,cmd->dynamic,cmd->nonlinearGeom);
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
					elasticLT.populateFromSparseMat(elasticMat,constraints);
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
	double ldSteps = cmd->loadRampSteps;
	double c1 = 1.0/(ldSteps*ldSteps);
	Node *thisNd;
	
	if(!anPrepRun) {
		analysisPrep(cmd->solverBlockDim);
	}
	
	if(cmd->thermal) {
	}
	
	if(cmd->elastic) {
		thisNd = nodes.getFirst();
		while(thisNd) {
			thisNd->initializeDisp();
			thisNd = thisNd->getNext();
		}
		if(!elasticLT.isAllocated()) {
			buildElasticSolnLoad(tempV1,true,cmd->dynamic,cmd->nonlinearGeom);
			elasticLT.allocateFromSparseMat(elasticMat,constraints,cmd->solverBandwidth);
		}
		if(!cmd->nonlinearGeom) {
			if (!elasticScaled) {
				scaleElasticConst();
			}
			elasticLT.populateFromSparseMat(elasticMat, constraints);
			elasticLT.ldlFactor();
		}
	}
	
	if(cmd->dynamic) {
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

void Model::solveForAdjoint() {
	int i1;
	Element* thisEl;

	if (solveCmd->elastic) {
		if (solveCmd->dynamic) {

		}
		else {
			for (i1 = 0; i1 < totGlobDof; i1++) {
				dLdU[i1] = 0.0;
			}
			objective.calculatedLdU(dLdU, dLdV, dLdA, dLdT, dLdTdot, solveCmd->staticLoadTime, solveCmd->nonlinearGeom, nodeArray,elementArray,dVarArray);
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
	}

	if (solveCmd->thermal) {

	}
	return;
}

void Model::dRelasticdD(int dVarNum) {
	int i1;
	int numDof;
	int globInd;
	Doub dvVal;
	DesignVariable* thisDV = dVarArray[dVarNum].ptr;
	thisDV->getValue(dvVal);
	thisDV->setDiffVal(dvVal.val, 1.0);

	for (i1 = 0; i1 < elMatDim; i1++) {
		dRudD[i1].setVal(0.0);
	}

	// Solution-dependent contribution of load
	IntListEnt* thisEnt = thisDV->getFirstEl();
	Element* thisElPt;
	while (thisEnt) {
		thisElPt = elementArray[thisEnt->value].ptr;
		thisElPt->getRu(dRudD, elasticMat, false, solveCmd->dynamic, solveCmd->nonlinearGeom, nodeArray, dVarArray);
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

	return;
}

void Model::getObjGradient() {
	int i1;
	int i2;
	int numDV = designVars.getLength();
	Element* thisEl;

	objective.clearValues();
	for (i1 = 0; i1 < numDV; i1++) {
		dLdD[i1] = 0.0;
	}

	if (solveCmd->dynamic) {

	}
	else {
		objective.calculateTerms(solveCmd->staticLoadTime, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
		solveForAdjoint();
		objective.calculatedLdD(dLdD, solveCmd->staticLoadTime, solveCmd->nonlinearGeom, nodeArray, elementArray, dVarArray);
		for (i1 = 0; i1 < numDV; i1++) {
			if (solveCmd->elastic) {
				dRelasticdD(i1);
				for (i2 = 0; i2 < elMatDim; i2++) {
					dLdD[i1] -= uAdj[i2] * dRudD[i2].dval;
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