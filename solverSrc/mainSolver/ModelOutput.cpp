#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "ModelClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "ListEntClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"
#include "JobClass.h"

using namespace std;

const int max_int = 2000000000;

void Model::writeTimeStepSoln(int tStep) {
	int i1;
	int dofPerNd;
	int numIntDof;
	double ndDat[6];
	double elDat[9];
	char* buf;
	JobCommand& scmd = job[solveCmd];
	string fullFile = scmd.fileName + "/solnTStep" + to_string(tStep) + ".out";
	ofstream outFile;
	outFile.open(fullFile, std::ofstream::binary);

	for (auto& thisNd : nodes) {
		if (scmd.thermal) {
			buf = reinterpret_cast<char*>(&thisNd.prevTemp);
			outFile.write(buf, 8);
			buf = reinterpret_cast<char*>(&thisNd.prevTdot);
			outFile.write(buf, 8);
		}
		if (scmd.elastic) {
			dofPerNd = thisNd.numDof;
			i1 = 8 * dofPerNd;
			buf = reinterpret_cast<char*>(&thisNd.prevDisp[0]);
			outFile.write(buf, i1);
			buf = reinterpret_cast<char*>(&thisNd.prevVel[0]);
			outFile.write(buf, i1);
			buf = reinterpret_cast<char*>(&thisNd.prevAcc[0]);
			outFile.write(buf, i1);
		}
	}

	if (scmd.elastic) {
		buf = reinterpret_cast<char*>(&elDat[0]);
		for (auto& thisEl : elements) {
			numIntDof = thisEl.numIntDof;
			if (numIntDof > 0) {
				i1 = 8 * numIntDof;
				vecToAr(elDat, thisEl.intPrevDisp, 0, numIntDof);
				outFile.write(buf, i1);
			}
		}
	}

	outFile.close();
	return;
}

void Model::writeNodeResults(string fileName, string nodeSet, list<string>& fields, int timeStep) {
	int i1;
	int globInd;
	int setPt;
	string errStr;
	string thisField;
	double ndDat[6];
	int ndof;
	double time;

	JobCommand& scmd = job[solveCmd];
	
	if(timeStep < max_int) {
		// Read the results from the time step file and store them in nodes
		readTimeStepSoln(timeStep);
		for (auto& ndPt : nodes) {
			if (scmd.thermal) {
				ndPt.backstepTemp();
			}
			if (scmd.elastic) {
				ndPt.backstepDisp();
			}
		}
		if (scmd.dynamic) {
			time = scmd.timeStep * timeStep;
		}
		else {
			auto thisLdTm = scmd.staticLoadTime.begin();
			i1 = 0;
			while (thisLdTm != scmd.staticLoadTime.end() && i1 < timeStep) {
				++thisLdTm;
				i1++;
			}
			time = *thisLdTm;
		}
	}
	else {
		time = -1.0;
	}
	
	if (key_in_map(nsMap,nodeSet)) {
		setPt = nsMap.at(nodeSet);
	}
	else {
		errStr = "Warning: there is no node set named " + nodeSet + ". Defaulting to all nodes in writeNodeResults";
		cout << errStr << endl;
		setPt = nsMap.at("all");
	}
	Set& thisSet = nodeSets[setPt];
	
	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);
	
	outFile << "nodeResults:\n";
	outFile << "    time: " << time << "\n";
	outFile << "    nodeSet: " << nodeSet << "\n";
	
	for (auto& thisField : fields) {
		outFile << "    " << thisField << ":\n";
		// -----------------
		// Calculate reaction force if necessary
		if (thisField == "reactionForce") {
			for (i1 = 0; i1 < elMatDim; i1++) {
				tempV1[i1] = 0.0;
			}
			buildElasticSolnLoad(tempV1, false, false);
			for (auto& ndLabel : thisSet.labels) {
				Node& ndPt = nodes[ndLabel];
				outFile << "        - [" << ndLabel;
				ndof = ndPt.numDof;
				for (i1 = 0; i1 < ndof; i1++) {
					globInd = ndPt.dofIndex[i1];
					outFile << ", " << -tempV1[globInd];
				}
				outFile << "]\n";
			}
		}
		else if (thisField == "reactionHeatGen") {
			for (i1 = 0; i1 < elMatDim; i1++) {
				tempV1[i1] = 0.0;
			}
			buildThermalSolnLoad(tempV1, false);
			for (auto& ndLabel : thisSet.labels) {
				Node& ndPt = nodes[ndLabel];
				outFile << "        - [" << ndLabel;
				globInd = ndPt.sortedRank;
				outFile << ", " << -tempV1[globInd] << "\n";
			}
		}
		else {
			for (auto& ndLabel : thisSet.labels) {
				Node& ndPt = nodes[ndLabel];
				outFile << "        - [" << ndLabel;
				if (thisField == "displacement") {
					ndof = ndPt.numDof;
					for (i1 = 0; i1 < ndof; i1++) {
						outFile << ", " << ndPt.displacement[i1];
					}
					outFile << "]\n";
				}
				else if (thisField == "velocity") {
					ndof = ndPt.numDof;
					for (i1 = 0; i1 < ndof; i1++) {
						outFile << ", " << ndPt.velocity[i1];
					}
					outFile << "]\n";
				}
				else if (thisField == "acceleration") {
					ndof = ndPt.numDof;
					for (i1 = 0; i1 < ndof; i1++) {
						outFile << ", " << ndPt.acceleration[i1];
					}
					outFile << "]\n";
				}
				else if (thisField == "temperature") {
					ndDat[0] = ndPt.temperature;
					outFile << ", " << ndDat[0] << "]\n";
				}
				else if (thisField == "tdot") {
					ndDat[0] = ndPt.tempChangeRate;
					outFile << ", " << ndDat[0] << "]\n";
				}
			}
		}
	}
	
	outFile.close();
	
	return;
}

void Model::writeElementResults(string fileName, string elSet, list<string>& fields, string position, int timeStep) {
	int i1;
	int i2;
	int i3;
	string errStr;
	string thisField;
	int setPt;
	int type;
	int numIP;
	int numLay;
	double intPts[24];
	double cent[3] = { 0.0,0.0,0.0 };
	DiffDoub0 strain[6];
	DiffDoub0 stress[6];
	double SEDen;
	DiffDoub0 def[9];
	DiffDoub0 frcMom[9];
	DiffDoub0 tGrad[3];
	DiffDoub0 flux[3];
	string fieldList;
	double time;

	JobCommand& scmd = job[solveCmd];

	if (timeStep < max_int) {
		// Read the results from the time step file and store them in nodes
		readTimeStepSoln(timeStep);
		for (auto& ndPt : nodes) {
			if (scmd.thermal) {
				ndPt.backstepTemp();
			}
			if (scmd.elastic) {
				ndPt.backstepDisp();
			}
		}
		if (scmd.elastic) {
			for (auto& elPt : elements) {
				elPt.backstepIntDisp();
			}
		}
		if (scmd.dynamic) {
			time = scmd.timeStep * timeStep;
		}
		else {
			auto thisLdTm = scmd.staticLoadTime.begin();
			i1 = 0;
			while (thisLdTm != scmd.staticLoadTime.end() && i1 < timeStep) {
				++thisLdTm;
				i1++;
			}
			time = *thisLdTm;
		}
	}
	else {
		time = -1.0;
	}

	if (key_in_map(esMap,elSet)) {
		setPt = esMap.at(elSet);
	}
	else {
		errStr = "Warning: there is no element set named " + elSet + ". Defaulting to all elements in writeNodeResults";
		cout << errStr << endl;
		setPt = esMap.at("all");
	}
	Set& thisSet = elementSets[setPt];

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "elementResults:\n";
	outFile << "    time: " << time << "\n";
	outFile << "    elSet: " << elSet << "\n";

	for (auto& thisField : fields) {
		outFile << "    " << thisField << ":\n";
		fieldList = "stress strain";
		i2 = fieldList.find(thisField);
		if (i2 > -1) {
			outFile << "    ##  - [element label, integration pt, layer, S11, S22, S33, S12, S13, S23]\n";
		}
		fieldList = "strainEnergyDen";
		i2 = fieldList.find(thisField);
		if (i2 > -1) {
			outFile << "    ##  - [element label, integration pt, layer, strain energy]\n";
		}
		fieldList = "sectionDef sectionFrcMom";
		i2 = fieldList.find(thisField);
		if (i2 > -1) {
			outFile << "    ## for shells:  - [element label, integration pt, S11, S22, S12, K11, K22, K12]\n";
			outFile << "    ## for beams:  - [element label, integration pt, S11, S12, S13, K11, K12, K13]\n";
		}
		if (thisField == "heatFlux") {
			outFile << "    ##  - [element label, integration pt, layer, f1, f2, f3]\n";
		}
		if (thisField == "tempGradient") {
			outFile << "    ##  - [element label, integration pt, layer, dT/dx, dT/dy, dT/dz]\n";
		}

		for (auto& elLabel : thisSet.labels) {
			Element& elPt = elements[elLabel];
			type = elPt.type;
			
			fieldList = "stress strain strainEnergyDen";
			i2 = fieldList.find(thisField);
			if (i2 > -1) {
				elPt.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
				if (position == "intPts") {
					numIP = elPt.numIP;
					vecToAr(intPts, elPt.intPts, 0, 3 * numIP);
				}
				else {
					numIP = 1;
					intPts[0] = elPt.sCent[0];
					intPts[1] = elPt.sCent[1];
					intPts[2] = elPt.sCent[2];
				}
				for (i1 = 0; i1 < numIP; i1++) {
					if (type == 3 || type == 41) {
						numLay = sections[elPt.sectPtr].layers.size();
					} else {
						numLay = 1;
					}
					for (i2 = 0; i2 < numLay; i2++) {
						elPt.getStressStrain(stress, strain, &intPts[3 * i1], i2, scmd.nonlinearGeom, d0Pre);
						outFile << "        - [" << elLabel << ", ";
						if (thisField == "strain") {
							outFile << i1 << ", " << i2 << ", " << strain[0].val;
							for (i3 = 1; i3 < 6; i3++) {
								outFile << ", " << strain[i3].val;
							}
							outFile << "]\n";
						} else if(thisField == "stress") {
							outFile << i1 << ", " << i2 << ", " << stress[0].val;
							for (i3 = 1; i3 < 6; i3++) {
								outFile << ", " << stress[i3].val;
							}
							outFile << "]\n";
						} else {
							SEDen = 0.0;
							for (i3 = 0; i3 < 6; i3++) {
								SEDen += stress[i3].val * strain[i3].val;
							}
							SEDen *= 0.5;
							outFile << i1 << ", " << i2 << ", " << SEDen << "]\n";
						}
					}
				}
			}

			fieldList = "sectionDef sectionFrcMom";
			i2 = fieldList.find(thisField);
			if (i2 > -1 && elPt.dofPerNd == 6) {
				elPt.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
				if (position == "intPts") {
					numIP = elPt.numIP;
					vecToAr(intPts, elPt.intPts, 0, 3 * numIP);
				}
				else {
					numIP = 1;
					intPts[0] = elPt.sCent[0];
					intPts[1] = elPt.sCent[1];
					intPts[2] = elPt.sCent[2];
				}
				for (i1 = 0; i1 < numIP; i1++) {
					type = elPt.type;
					//elPt->getStressStrain(stress, strain, &intPts[3 * i1], i2, solveCmd->nonlinearGeom, stPre);
					elPt.getDefFrcMom(def, frcMom, &intPts[3 * i1], scmd.nonlinearGeom, d0Pre);
					outFile << "        - [" << elLabel << ", ";
					if (thisField == "sectionDef") {
						outFile << i1 << ", " << def[0].val;
						for (i3 = 1; i3 < 6; i3++) {
							outFile << ", " << def[i3].val;
						}
						outFile << "]\n";
					}
					else if (thisField == "sectionFrcMom") {
						outFile << i1 << ", " << frcMom[0].val;
						for (i3 = 1; i3 < 6; i3++) {
							outFile << ", " << frcMom[i3].val;
						}
						outFile << "]\n";
					}
				}
			}

			fieldList = "tempGradient heatFlux";
			i2 = fieldList.find(thisField);
			if (i2 > -1) {
				elPt.getStressPrereq(d0Pre, sections, materials, nodes, designVars);
				if (position == "intPts") {
					numIP = elPt.numIP;
					vecToAr(intPts, elPt.intPts, 0, 3 * numIP);
				}
				else {
					numIP = 1;
					intPts[0] = elPt.sCent[0];
					intPts[1] = elPt.sCent[1];
					intPts[2] = elPt.sCent[2];
				}
				for (i1 = 0; i1 < numIP; i1++) {
					type = elPt.type;
					if (type == 3 || type == 41) {
						numLay = sections[elPt.sectPtr].layers.size();
					}
					else {
						numLay = 1;
					}
					for (i2 = 0; i2 < numLay; i2++) {
						//elPt->getStressStrain(stress, strain, &intPts[3 * i1], i2, solveCmd->nonlinearGeom, stPre);
						elPt.getFluxTGrad(flux, tGrad, &intPts[3 * i1], i2, d0Pre);
						outFile << "        - [" << elLabel << ", " << i1 << ", " << i2 << ", ";
						if (thisField == "tempGradient") {
							outFile << tGrad[0].val << ", " << tGrad[1].val << ", " << tGrad[2].val << "]\n";
						}
						else if (thisField == "heatFlux") {
							outFile << flux[0].val << ", " << flux[1].val << ", " << flux[2].val << "]\n";
						}
					}
				}
			}

		}
	}

	outFile.close();

	return;
}

void Model::writeModalResults(string fileName, bool writeModes) {
	int i1;
	int i2;
	int i3;
	int nd;
	int dofPerNd;
	int globInd;
	JobCommand& mcmd = job[modalCmd];
	int nMds = mcmd.numModes;
	Node* thisNd;

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "modalResults:\n";
	outFile << "    eigenValues:\n";
	for (i1 = 0; i1 < nMds; i1++) {
		outFile << "      - " << eigVals[i1] << "\n";
	}
	if (mcmd.type == "buckling") {
		outFile << "    loadFactors:\n";
	}
	else {
		outFile << "    frequencies:\n";
	}
	for (i1 = 0; i1 < nMds; i1++) {
		outFile << "      - " << loadFact[i1] << "\n";
	}

	if (mcmd.writeModes) {
		outFile << "    modes:\n";
		for (i1 = 0; i1 < nMds; i1++) {
			outFile << "      - mode: " << i1 << "\n";
			outFile << "        displacement:\n";
			for (auto& thisNd : nodes) {
				nd = thisNd.label;
				outFile << "          - [" << nd;
				dofPerNd = thisNd.numDof;
				for (i2 = 0; i2 < dofPerNd; i2++) {
					globInd = thisNd.dofIndex[i2];
					i3 = i1 * elMatDim + globInd;
					outFile << ", " << eigVecs[i3];
				}
				outFile << "]\n";
			}
		}
	}

	outFile.close();

	return;
}

void Model::writeObjective(string fileName, list<string>& includeFields, bool writeGrad) {
	int i1;
	int numDV;
	double totObj;
	double thisVal;
	string fldStr;

	ofstream outFile;
	outFile.open(fileName);
	outFile << setprecision(12);

	outFile << "objective:\n";
	outFile << "    terms:\n";
	totObj = 0.0;
	for (auto& thisTerm : objective.terms) {
		thisVal = thisTerm.value;
		outFile << "        - value: " << thisVal << "\n";
		for (auto& fldStr : includeFields) {
			if (fldStr == "category") {
				outFile << "          category: " << thisTerm.category << "\n";
			}
			else if (fldStr == "operator") {
				outFile << "          operator: " << thisTerm.optr << "\n";
			}
			else if (fldStr == "component") {
				outFile << "          component: " << thisTerm.component << "\n";
			}
			else if (fldStr == "layer") {
				outFile << "          layer: " << thisTerm.layer << "\n";
			}
			else if (fldStr == "coefficient") {
				outFile << "          coefficient: " << thisTerm.coef << "\n";
			}
			else if (fldStr == "exponent") {
				outFile << "          exponent: " << thisTerm.expnt << "\n";
			}
			else if (fldStr == "elementSet") {
				outFile << "          elementSet: " << thisTerm.elSetName << "\n";
			}
			else if (fldStr == "nodeSet") {
				outFile << "          nodeSet: " << thisTerm.ndSetName << "\n";
			}
			else if (fldStr == "activeTime") {
				outFile << "          activeTime: [" << thisTerm.activeTime[0] << ", " << thisTerm.activeTime[1] << "]\n";
			}
		}
		totObj += thisVal;
	}
	outFile << "    totalValue: " << totObj << "\n";
	if (writeGrad) {
		outFile << "objectiveGradient:\n";
		numDV = designVars.size();
		for (i1 = 0; i1 < numDV; i1++) {
			outFile << "    - [" << i1 << ", " << dLdD[i1] << "]\n";
		}
	}

	outFile.close();

	return;
}