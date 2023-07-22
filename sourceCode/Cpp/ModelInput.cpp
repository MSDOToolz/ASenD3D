#include <fstream>
#include <iostream>
#include "ModelClass.h"
#include "JobClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "SetClass.h"
#include "SectionClass.h"
#include "ConstraintClass.h"
#include "DesignVariableClass.h"
#include "ObjectiveClass.h"

using namespace std;

// Input

void Model::readInputLine(ifstream& inFile, string& fileLine, string headings[], int hdLdSpace[], string data[], int& dataLen) {
	int i1;
	int i2;
	int i3;
	int i4;
	int lnLen;
	int wrdLen;
	bool brk;
	getline(inFile,fileLine);
	i1 = fileLine.find("#");
	if(i1 > -1) {
		fileLine = fileLine.substr(0,i1);
	}
	lnLen = fileLine.length();
	i1 = fileLine.find(":");
	dataLen = 0;
	if(i1 > -1) {
		i2 = fileLine.find_first_not_of(" -\n\t");
		wrdLen = i1 - i2;
		if(headings[0] == "" || hdLdSpace[0] == i2) {
			headings[0] = fileLine.substr(i2,wrdLen);
			hdLdSpace[0] = i2;
			headings[1] = "";
			hdLdSpace[1] = 0;
			headings[2] = "";
			hdLdSpace[2] = 0;
			headings[3] = "";
			hdLdSpace[3] = 0;
		} else if(headings[1] == "" || hdLdSpace[1] == i2) {
			headings[1] = fileLine.substr(i2,wrdLen);
			hdLdSpace[1] = i2;
			headings[2] = "";
			hdLdSpace[2] = 0;
			headings[3] = "";
			hdLdSpace[3] = 0;
		} else if(headings[2] == "" || hdLdSpace[2] == i2) {
			headings[2] = fileLine.substr(i2,wrdLen);
			hdLdSpace[2] = i2;
			headings[3] = "";
			hdLdSpace[3] = 0;
		} else {
			headings[3] = fileLine.substr(i2,wrdLen);
			hdLdSpace[3] = i2;
		}
		i1++;
		while(i1 < lnLen) {
			fileLine = fileLine.substr(i1);
			i1 = fileLine.find_first_not_of(" ,[]\t\n");
			if(i1 > -1) {
				fileLine = fileLine.substr(i1);
				lnLen = fileLine.length();
				i1 = fileLine.find_first_of(" ,[]\t\n");
				if(i1 > -1) {
				    data[dataLen] = fileLine.substr(0,i1);
				    dataLen++;
				} else {
					i1 = lnLen;
				}
			} else {
				i1 = lnLen;
			}
		}
	} else {
		i1 = fileLine.find("- ");
		if(i1 > -1) {
			i1++;
			while(i1 < lnLen) {
				fileLine = fileLine.substr(i1);
				i1 = fileLine.find_first_not_of(" ,[]\t\n");
				if(i1 > -1) {
					fileLine = fileLine.substr(i1);
					lnLen = fileLine.length();
					i1 = fileLine.find_first_of(" ,[]\t\n");
					if(i1 > -1) {
						data[dataLen] = fileLine.substr(0,i1);
						dataLen++;
					} else {
						i1 = lnLen;
					}
				} else {
					i1 = lnLen;
				}
			}
		}
	}
	
	return;
}

void Model::readJob(string fileName) {
	int i1;
	int i2;
	int i3;
	int i4;
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	string prevLdHd = "";
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	string errSt;
	
	JobCommand *newCmd;
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] != prevLdHd) {
				newCmd = new JobCommand();
				newCmd->cmdString = headings[0];
				job.addCommand(newCmd);
			}
			if(headings[1] == "fileName" && dataLen == 1) {
				newCmd->fileName = data[0];
			} else if(headings[1] == "dynamic" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->dynamic = true;
				} else {
					newCmd->dynamic = false;
				}
			} else if(headings[1] == "elastic" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->elastic = true;
				} else {
					newCmd->elastic = false;
				}
			} else if(headings[1] == "loadRampSteps" && dataLen == 1) {
				newCmd->loadRampSteps = stoi(data[0]);
			} else if(headings[1] == "newmarkBeta" && dataLen == 1) {
				newCmd->newmarkBeta = stod(data[0]);
			} else if(headings[1] == "newmarkGamma" && dataLen == 1) {
				newCmd->newmarkGamma = stod(data[0]);
			} else if(headings[1] == "nonlinearGeom" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->nonlinearGeom = true;
				} else {
					newCmd->nonlinearGeom = false;
				}
			} else if(headings[1] == "saveSolnHist" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->saveSolnHist = true;
				} else {
					newCmd->saveSolnHist = false;
				}
			} else if(headings[1] == "simPeriod" && dataLen == 1) {
				newCmd->simPeriod = stod(data[0]);
			} else if(headings[1] == "solverBandwidth" && dataLen == 1) {
				newCmd->solverBandwidth = stoi(data[0]);
			} else if(headings[1] == "solverBlockDim" && dataLen == 1) {
				newCmd->solverBlockDim = stoi(data[0]);
			} else if(headings[1] == "solverMethod" && dataLen == 1) {
				newCmd->solverMethod = data[0];
			} else if(headings[1] == "staticLoadTime" && dataLen == 1) {
				newCmd->staticLoadTime = stod(data[0]);
			} else if(headings[1] == "thermal" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->thermal = true;
				} else {
					newCmd->thermal = false;
				}
			} else if(headings[1] == "timeStep" && dataLen == 1) {
				newCmd->timeStep = stod(data[0]);
			} else if(headings[1] == "type" && dataLen == 1) {
				newCmd->type = data[0];
			} else if(headings[1] == "numModes" && dataLen == 1) {
				newCmd->numModes = stoi(data[0]);
			} else if(headings[1] == "solnField" && dataLen == 1) {
				newCmd->solnField = data[0];
			} else if(headings[1] == "mode" && dataLen == 1) {
				newCmd->mode = stoi(data[0]);
			} else if(headings[1] == "maxAmplitude" && dataLen == 1) {
				newCmd->maxAmplitude = stod(data[0]);
			} else if(headings[1] == "nodeSet" && dataLen == 1) {
				newCmd->nodeSet = data[0];
			} else if(headings[1] == "fields" && dataLen == 1) {
				newCmd->fields.addEntry(data[0]);
			} else if(headings[1] == "timeSteps") {
				if(dataLen == 1) {
					if(data[0] == "all") {
						newCmd->timeStepTag = "all";
					} else if(data[0] == "last") {
						newCmd->timeStepTag = "last";
					} else {
						try {
							i1 = stoi(data[0]);
							newCmd->timeSteps.addEntry(stoi(data[0]));
						} catch(...) {
							errSt = "Warning: possible invalid entry for " + headings[0] + " timeSteps in job file " + fileName;
							cout << errSt << endl;
						}
					}
				} else if(dataLen > 1) {
					i2 = stoi(data[0]);
					i3 = stoi(data[1]);
					if(dataLen == 3) {
						i4 = stoi(data[2]);
					} else {
						i4 = 1;
					}
					for (i1 = i2; i1 < i3; i1+= i4) {
						newCmd->timeSteps.addEntry(i1);
					}
				}
			} else if(headings[1] == "elementSet" && dataLen == 1) {
				newCmd->elementSet = data[0];
			} else if(headings[1] == "writeModes" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->writeModes = true;
				} else {
					newCmd->writeModes = false;
				}
			} else if(headings[1] == "properties" && dataLen == 1) {
				newCmd->properties.addEntry(data[0]);
			} else if(headings[1] == "include" && dataLen == 1) {
				newCmd->objInclude.addEntry(data[0]);
			} else if(headings[1] == "writeGradient" && dataLen == 1) {
				if(data[0] == "yes") {
					newCmd->writeGradient = true;
				} else {
					newCmd->writeGradient = false;
				}
			}
			
			prevLdHd = headings[0];
		}
	} else {
		errSt = "Error: could not open job file: " + fileName;
		throw invalid_argument(errSt);
	}
	
	return;
}
		
void Model::readModelInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	int i1;
	int i2;
	double doubInp[10];
	int intInp[10];
	
	Node *newNd;
	Element *newEl;
	int elType;
	Set *newSet;
	Section *newSec;
	Layer *newLay;
	Material *newMat;
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "nodes" && dataLen == 4) {
				newNd = new Node(stoi(data[0]));
				doubInp[0] = stod(data[1]);
				doubInp[1] = stod(data[2]);
				doubInp[2] = stod(data[3]);
				newNd->setCrd(doubInp);
				nodes.addNode(newNd);
			} else if(headings[0] == "elements") {
				if(headings[1] == "type" && dataLen == 1) {
					if(data[0] == "tet4") {
						elType = 4;
					} else if(data[0] == "wedge6") {
						elType = 6;
					} else if(data[0] == "brick8") {
						elType = 8;
					} else if(data[0] == "brickIM") {
						elType = 81;
					} else if(data[0] == "shell3") {
						elType = 3;
					} else if(data[0] == "shell4") {
						elType = 41;
					} else if(data[0] == "beam2") {
						elType = 2;
					} else {
						string errSt = "Error: unrecognized element type, " + data[0];
						throw invalid_argument(errSt);
					}
				} else if(headings[1] == "connectivity" && dataLen > 1) {
					newEl = new Element(elType);
					newEl->setLabel(stoi(data[0]));
					for (i1 = 1; i1 < dataLen; i1++) {
						intInp[i1-1] = stoi(data[i1]);
					}
					newEl->setNodes(intInp);
					elements.addElement(newEl);
				} 
			} else if(headings[0] == "sets") {
				if(headings[1] == "node") {
					if(headings[2] == "name" && dataLen == 1) {
						newSet = new Set(data[0]);
						nodeSets.addSet(newSet);
					} else if(headings[2] == "labels") {
						if(dataLen == 1) {
							if(data[0] == "all") {
								i2 = nodes.getLength();
								for (i1 = 0; i1 < i2; i1++) {
									newSet->addEntry(i1);
								}
							} else {
								newSet->addEntry(stoi(data[0]));
							}
						} else if(dataLen > 1) {
							i2 = stoi(data[0]);
							i3 = stoi(data[1]);
							if(dataLen == 3) {
								i4 = stoi(data[2]);
							} else {
								i4 = 1;
							}
							for (i1 = i2; i1 < i3; i1+= i4) {
								newSet->addEntry(i1);
							}
						}
					}
				} else if(headings[1] == "element") {
					if(headings[2] == "name" && dataLen == 1) {
						newSet = new Set(data[0]);
						elementSets.addSet(newSet);
					} else if(headings[2] == "labels") {
						if(dataLen == 1) {
							if(data[0] == "all") {
								i2 = elements.getLength();
								for (i1 = 0; i1 < i2; i1++) {
									newSet->addEntry(i1);
								}
							} else {
								newSet->addEntry(stoi(data[0]));
							}
						} else if(dataLen > 1) {
							i2 = stoi(data[0]);
							i3 = stoi(data[1]);
							if(dataLen == 3) {
								i4 = stoi(data[2]);
							} else {
								i4 = 1;
							}
							for (i1 = i2; i1 < i3; i1+= i4) {
								newSet->addEntry(i1);
							}
						}
					}
				}
			} else if(headings[0] == "sections") {
				if(headings[1] == "type" && dataLen == 1) {
					newSec = new Section(data[0]);
					sections.addSection(newSec);
				} else if(headings[1] == "material" && dataLen == 1) {
					newSec->setMaterial(data[0]);
				} else if(headings[1] == "orientation" && dataLen == 6) {
					for (i1 = 0; i1 < 6; i1++) {
						doubInp[i1] = stod(data[i1]);
					}
					newSet->setOrientation(doubInp);
				} else if(headings[1] == "layup") {
					if(headings[2] == "zOffset" && dataLen == 1) {
						newSec->setZOffset(stod(data[0]));
					} else if(headings[2] == "layers") {
						if(headings[3] == "material" && dataLen == 1) {
							newLay = new Layer(data[0]);
							newSec->addLayer(newLay);
						} else if(headings[3] == "thickness" && dataLen == 1) {
							newLay->setThickness(stod(data[0]));
						} else if(headings[3] == "angle" && dataLen == 1) {
							newLay->setAngle(stod(data[0]));
						}
					}
				} else if(headings[1] == "beamProps") {
					if(headings[2] == "area" && dataLen == 1) {
						newSec->setArea(stod(data[0]));
					} else if(headings[2] == "I" && dataLen == 5) {
						for(i1 = 0; i1 < 5; i1++) {
							doubInp[i1] = stod(data[i1]);
						}
						newSec->setAreaMoment(doubInp);
					} else if(headings[2] == "J" && dataLen == 1) {
						newSec->setPolarMoment(stod(data[1]);
					} else if(headings[2] == "stiffness" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						newSec->setStiffness(intInp[0],intInp[1],doubInp[0]);
					} else if(headings[2] == "mass" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						newSec->setMass(intInp[0],intInp[1],doubInp[0]);
					} else if(headings[2] == "expLoadCoef" && dataLen == 6) {
						for(i1 = 0; i1 < 6; i1++) {
							doubInp[i1] = stod(data[i1]);
						}
						newSec->setExpLd(doubInp);
					} else if(headings[2] == "conductivity" && dataLen == 1) {
						newSec->setConductivity(stod(data[0]));
					} else if(headings[2] == "specHeat" && dataLen == 1) {
						newSec->setSpecHeat(stod(data[0]));
					}
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					newSec->setElset(data[0]);
				}
			} else if(headings[0] == "materials") {
				if(headings[1] == "name" && dataLen == 1) {
					newMat = new Material(data[0]);
					materials.addMaterial(newMat);
				} else if(headings[1] == "density" && dataLen == 1) {
					newMat->setDensity(stod(data[0]));
				} else if(headings[1] == "elastic") {
					if(headings[2] == "E" && dataLen > 0) {
						doubInp[0] = stod(data[0]);
						if(dataLen == 3) {
							doubInp[1] = stod(data[1]);
							doubInp[2] = stod(data[2]);
						} else {
							doubInp[1] = doubInp[0];
							doubInp[2] = doubInp[0];
						}
						newMat->setModulus(doubInp);
					} else if(headings[2] == "nu" && dataLen > 0) {
						doubInp[0] = stod(data[0]);
						if(dataLen == 3) {
							doubInp[1] = stod(data[1]);
							doubInp[2] = stod(data[2]);
						} else {
							doubInp[1] = doubInp[0];
							doubInp[2] = doubInp[0];
						}
						newMat->setPoissonRatio(doubInp);
					} else if(headings[2] == "G" && dataLen > 0) {
						doubInp[0] = stod(data[0]);
						if(dataLen == 3) {
							doubInp[1] = stod(data[1]);
							doubInp[2] = stod(data[2]);
						} else {
							doubInp[1] = doubInp[0];
							doubInp[2] = doubInp[0];
						}
						newMat->setShearMod(doubInp);
					} else if(headings[2] == "stiffness" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						newMat->setStiffness(intInp[0],intInp[1],doubInp[0]);
					} 
				} else if(headings[1] == "thermal") {
					if(headings[2] == "conductivity" && dataLen == 6) {
						for(i1 = 0; i1 < 6; i1++) {
							doubInp[i1] = stod(data[i1]);
						}
						newMat->setConductivity(doubInp);
					} else if(headings[2] == "expansion" && dataLen == 6) {
						for(i1 = 0; i1 < 6; i1++) {
							doubInp[i1] = stod(data[i1]);
						}
						newMat->setExpansion(doubInp);
					} else if(headings[2] == "specHeat" && dataLen == 1) {
						newMat->setSpecHeat(stod(data[0]));
				} else if(headings[1] == "failure") {
					if(headings[2] == "maxStress") {
						if(headings[3] == "tensile" && dataLen > 0) {
							doubInp[0] = stod(data[0]);
							if(dataLen == 3) {
								doubInp[1] = stod(data[1]);
								doubInp[2] = stod(data[2]);
							} else {
								doubInp[1] = doubInp[0];
								doubInp[2] = doubInp[0]:
							}
							newMat->setMaxTenStress(doubInp);
						} else if(headings[3] == "compressive" && dataLen > 0) {
							doubInp[0] = stod(data[0]);
							if(dataLen == 3) {
								doubInp[1] = stod(data[1]);
								doubInp[2] = stod(data[2]);
							} else {
								doubInp[1] = doubInp[0];
								doubInp[2] = doubInp[0]:
							}
							newMat->setMaxCompStress(doubInp);
						} else if(headings[3] == "shear" && dataLen > 0) {
							doubInp[0] = stod(data[0]);
							if(dataLen == 3) {
								doubInp[1] = stod(data[1]);
								doubInp[2] = stod(data[2]);
							} else {
								doubInp[1] = doubInp[0];
								doubInp[2] = doubInp[0]:
							}
							newMat->setMaxShearStress(doubInp);
						}
					} else if(headings[2] == "maxStrain") {
						if(headings[3] == "tensile" && dataLen > 0) {
							doubInp[0] = stod(data[0]);
							if(dataLen == 3) {
								doubInp[1] = stod(data[1]);
								doubInp[2] = stod(data[2]);
							} else {
								doubInp[1] = doubInp[0];
								doubInp[2] = doubInp[0]:
							}
							newMat->setMaxTenStrain(doubInp);
						} else if(headings[3] == "compressive" && dataLen > 0) {
							doubInp[0] = stod(data[0]);
							if(dataLen == 3) {
								doubInp[1] = stod(data[1]);
								doubInp[2] = stod(data[2]);
							} else {
								doubInp[1] = doubInp[0];
								doubInp[2] = doubInp[0]:
							}
							newMat->setMaxCompStrain(doubInp);
						} else if(headings[3] == "shear" && dataLen > 0) {
							doubInp[0] = stod(data[0]);
							if(dataLen == 3) {
								doubInp[1] = stod(data[1]);
								doubInp[2] = stod(data[2]);
							} else {
								doubInp[1] = doubInp[0];
								doubInp[2] = doubInp[0]:
							}
							newMat->setMaxShearStrain(doubInp);
						}
					} else if(headings[2] == "maxStrainEnergy" && dataLen == 1) {
						newMat->setMaxStrEng(stod(data[0]));
					} else if(headings[2] == "maxMises" && dataLen == 1) {
						newMat->setMaxMises(stod(data[0]));
					}
				}
			} 
		}
	} else {
		string errSt = "Error: could not open model input file: " + fileName;
		throw invalid_argument(errSt);
	}
	
	// Populate node array
	i1 = nodes.getLength();
	nodeArray = new NdPt[i1];
	newNd = nodes.getFirst();
	while(newNd) {
		i1 = newNd->getLabel();
		nodeArray[i1].ptr = newNd;
		newNd = newNd->getNext();
	}
	
	// Populate element array
	i1 = elements.getLength();
	elementArray = new ElPt[i1];
	newEl = elements.getFirst();
	while(newEl) {
		i1 = newEl->getLabel();
		elArray[i1].ptr = newEl;
		newEl = newEl->getNext();
	}
	
	inFile.close();
	
	return;
}

void Model::readConstraintInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	int i1;
	int i2;
	double doubInp[10];
	int intInp[10];
	string prevHd2 = "";
	
    Constraint *newConst = NULL;
	ConstraintTerm *newTerm;
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "constraints") {
				if(headings[1] == "terms") {
					if(headings[2] == "nodeSet" && dataLen == 1) {
						if(prevHd2 == "") {
							newConst = new Constraint;
							constraints.addConstraint(newConst);
						}
						newTerm = new ConstraintTerm(data[0]);
						newConst->addTerm(newTerm);
					} else if(headings[2] == "dof" && dataLen == 1) {
						newTerm->setDof(stoi(data[0]));
					} else if(headings[2] == "coef" && dataLen == 1) {
						newTerm->setCoef(stod(data[0]);
					}
				} else if(headings[1] == "rhs" && dataLen == 1) {
					newConst->setRhs(stod(data[0]));
				}
			}
			prevHd2 = headings[2];
		}
	} else {
		string errSt = "Error: could not open constraint input file: " + fileName;
		throw invalid_argument(errSt);
	}
	return;
}

void Model::readLoadInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	int i1;
	int i2;
	double doubInp[10];
	int intInp[10];
	
	Load *newLd;
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "loads") {
				if(headings[1] == "type" && dataLen == 1) {
					newLd = new Load(data[0]);
					loads.addLoad(newLd);
				} else if(headings[1] == "activeTime" && dataLen > 0) {
					doubInp[0] = stod(data[0]);
					if(dataLen == 2) {
						doubInp[1] = stod(data[1]);
					} else {
						doubInp[1] = 1.0e+100;
					}
					newLd->setActTime(doubInp);
				} else if(headings[1] == "load" && dataLen > 0) {
					for (i1 = 0; i1 < 6; i1++) {
						if(i1 < dataLen) {
							doubInp[i1] = stod(data[i1]);
						} else {
							doubInp[i1] = 0.0;
						}
					}
					newLd->setLoad(doubInp);
				} else if(headings[1] == "nodeSet" && dataLen == 1) {
					newLd->setNodeSet(data[0]);
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					newLd->setElSet(data[0]);
				} else if(headings[1] == "normDir" && dataLen == 3) {
					doubInp[0] = stod(data[0]);
					doubInp[1] = stod(data[1]);
					doubInp{2] = stod(data[2]);
					newLd->setNormDir(doubInp);
				} else if(headings[1] == "normTolerance" && dataLen == 1) {
					newLd->setNormTol(stod(data[0]));
				} else if(headings[1] == "center" && dataLen == 3) {
					doubInp[0] = stod(data[0]);
					doubInp[1] = stod(data[1]);
					doubInp{2] = stod(data[2]);
					newLd->setCenter(doubInp);
				} else if(headings[1] == "axis" && dataLen == 3) {
					doubInp[0] = stod(data[0]);
					doubInp[1] = stod(data[1]);
					doubInp{2] = stod(data[2]);
					newLd->setAxis(doubInp);
				} else if(headings[1] == "angularVelocity") {
					newLd->setAngVel(stod(data[0]));
				}
			}
		}
	} else {
		string errSt = "Error: could not open load input file: " + fileName;
		throw invalid_argument(errSt);
	}
	
	return;
}

void Model::readInitialState(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	int i1;
	int i2;
	int i3;
	int ndi;
	double doubInp[10];
	int intInp[10];
	string dispHdings = " displacement velocity acceleration";
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "initialState") {
				i3 = dispHdings.find(headings[1]);
				if(i3 > -1 && dataLen > 3) {
					ndi = stoi(data[0]);
					i2 = 1;
					for (i1 = 0; i1 < 6; i1++) {
						if(i2 < dataLen) {
							doubInp[i1] = stod(data[i2])
						} else {
							doubInp[i1] = 0.0;
						}
						i2++;
					}
					if(headings[1] == "displacement") {
						nodeArray[ndi].ptr->setInitialDisp(doubInp);
					} else if(headings[1] == "velocity") {
						nodeArray[ndi].ptr->setInitialVel(doubInp);
					} else if(headings[1] == "acceleration") {
						nodeArray[ndi].ptr->setInitialAcc(doubInp);
					}
				} else if(headings[1] == "temperature" && dataLen == 2) {
					ndi = stoi(data[0]);
					nodeArray[ndi].ptr->setInitialTemp(stod(data[1]));
				} else if(headings[1] == "tdot" && dataLen == 2) {
					ndi = stoi(data[0]);
					nodeArray[ndi].ptr->setInitialTdot(stod(data[1]));
				}
			}
		}
	} else {
		string errSt = "Error: could not open initial state input file: " + fileName;
		throw invalid_argument(errSt);
	}
	
	return;
}

void Model::readDesVarInput(string fileName) {
    ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	int i1;
	int i2;
	double doubInp[10];
	int intInp[10];
	
	DesignVariable *newDVar;
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "designVariables") {
				if(headings[1] == "category" && dataLen == 1) {
					newDVar = new DesignVariable(data[0]);
					designVars.addDVar(newDVar);
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					newDVar->setElset(data[0]);
				} else if(headings[1] == "nodeSet" && dataLen == 1) {
					newDVar->setNdset(data[0]);
				} else if(headings[1] == "activeTime" && dataLen > 0) {
					doubInp[0] = stod(data[0]);
					if(dataLen == 2) {
						doubInp[1] = stod(data[1]);
					} else {
						doubInp[1] = 1.0e+100;
					}
					newDVar->setActiveTime(doubInp);
				} else if(headings[1] == "component") {
					if(dataLen == 1) {
					    newDVar->setComponent(stoi(data[0]);
					} else if(dataLen == 2) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						if(intInp[0] >= intInp[1]) {
							i1 = 6*intInp[1] + intInp[0];
						} else {
							i1 = 6*intInp[0] + intInp[1];
						}
						newDVar->setComponent(i1);
					}
				} else if(headings[1] == "layer" && dataLen == 1) {
					newDVar->setLayer(stoi(data[0]));
				} else if(headings[1] == "coefficients" && dataLen == 1) {
					newDVar->addCoefficient(stod(data[0]));
				}
			}
		}
	} else {
		string errSt = "Error: could not open design variable input file: " + fileName;
		throw invalid_argument(errSt);
	}
	
	// Populate design variable array
	i1 = designVars.getLength();
	dVarArray = new DVPt[i1];
	newDVar = designVars.getFirst();
	i2 = 0;
	while(newDVar) {
		dVarArray[i2].ptr = newDVar;
		newDVar = newDVar->getNext();
		i2++;
	}
	
	return;
}

void Model::readObjectiveInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	int i1;
	int i2;
	double doubInp[10];
	int intInp[10];
	
	ObjectiveTerm *newTerm;
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "objectiveTerms") {
				if(headings[1] == "category" && dataLen == 1) {
					newTerm = new ObjectiveTerm(data[0]);
					objective.addTerm(newTerm);
				} else if(headings[1] == "operator" && dataLen == 1) {
					newTerm->setOperator(data[0]);
				} else if(headings[1] == "activeTime" && dataLen > 0) {
					doubInp[0] = stod(data[0]);
					if(dataLen == 2) {
						doubInp[1] = stod(data[1]);
					} else {
						doubInp[1] = 1.0e+100;
					}
					newTerm->setActiveTime(doubInp);
				} else if(headings[1] == "component" && dataLen == 1) {
					newTerm->setComponent(stoi(data[0]));
				} else if(headings[1] == "layer" && dataLen == 1) {
					newTerm->setLayer(stoi(data[0]));
				} else if(headings[1] == "coefficient" && dataLen == 1) {
					newTerm->setCoef(stod(data[0]));
				} else if(headings[1] == "exponent" && dataLen == 1) {
					newTerm->setExponent(stod(data[0]));
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					newTerm->setElset(data[0]);
				} else if(headings[1] == "nodeSet" && dataLen == 1) {
					newTerm->setNdset(data[0]);
				} else if(headings[1] == "targetValue" && dataLen == 1) {
					try {
						doubInp[0] = stod(data[0]);
						newTerm->addTargetValue(doubInp[0]);
					} catch(...) {
						newTerm->setTgtTag(data[0]);
					}
				}
			}
		}
	} else {
		string errSt = "Error: could not open objective input file: " + fileName;
		throw invalid_argument(errSt);
	}
	
	return;
}

void Model::readDesVarValues(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	int label;
	double value;
	
	infile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
		    readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(dataLen == 2) {
				label = stoi(data[0]);
				value = stod(data[1]);
				dVarArray[label].ptr->setValue(value);
			}
		}
	} else {
		string errSt = "Error: could not open design variable value input file: " + fileName;
		throw invalid_argument(errSt);		
	}
	
	return;
}