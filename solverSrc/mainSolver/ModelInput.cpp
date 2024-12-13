#include <fstream>
#include <iostream>
#include <string>
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
	int lnLen;
	int wrdLen;
	getline(inFile,fileLine);
	i1 = fileLine.find("#");
	if(i1 > -1) {
		fileLine = fileLine.substr(0,i1);
	}
	fileLine = fileLine + " ";
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
	int cmdCt;
	
	cmdCt = 0;
	inFile.open(fileName);
	while (!inFile.eof()) {
		readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
		if (headings[1] == "command" && dataLen == 1) {
			cmdCt++;
		}
	}

	job = vector<JobCommand>(cmdCt);
	//inFile.seekg(0, std::ios::beg);
	inFile.close();
	inFile.open(fileName);
	
	if(inFile) {
		cmdCt = -1;
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[1] == "command" && dataLen == 1) {
				cmdCt++;
				job[cmdCt].cmdString = data[0];
			} else if(headings[1] == "fileName" && dataLen == 1) {
				job[cmdCt].fileName = data[0];
			} else if(headings[1] == "dynamic" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].dynamic = true;
				} else {
					job[cmdCt].dynamic = false;
				}
			}
			else if (headings[1] == "elastic" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].elastic = true;
				} else {
					job[cmdCt].elastic = false;
				}
			} else if(headings[1] == "loadRampSteps" && dataLen == 1) {
				job[cmdCt].loadRampSteps = stoi(data[0]);
			} else if(headings[1] == "newmarkBeta" && dataLen == 1) {
				job[cmdCt].newmarkBeta = stod(data[0]);
			} else if(headings[1] == "newmarkGamma" && dataLen == 1) {
				job[cmdCt].newmarkGamma = stod(data[0]);
			} else if(headings[1] == "nonlinearGeom" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].nonlinearGeom = true;
				} else {
					job[cmdCt].nonlinearGeom = false;
				}
			} else if(headings[1] == "saveSolnHist" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].saveSolnHist = true;
				} else {
					job[cmdCt].saveSolnHist = false;
				}
			}
			else if (headings[1] == "solnHistDir" && dataLen == 1) {
				job[cmdCt].fileName = data[0];
			}
			else if (headings[1] == "lumpMass" && dataLen == 1) {
				if (data[0] == "yes") {
					job[cmdCt].lumpMass = true;
				}
				else {
					job[cmdCt].lumpMass = false;
				}
			}
			else if (headings[1] == "fullReformFreq") {
				job[cmdCt].fullReform = stoi(data[0]);
			}
			else if (headings[1] == "simPeriod" && dataLen == 1) {
				job[cmdCt].simPeriod = stod(data[0]);
			} else if(headings[1] == "solverBandwidth" && dataLen == 1) {
				job[cmdCt].solverBandwidth = stoi(data[0]);
			} else if(headings[1] == "solverBlockDim" && dataLen == 1) {
				job[cmdCt].solverBlockDim = stoi(data[0]);
			} else if(headings[1] == "solverMethod" && dataLen == 1) {
				job[cmdCt].solverMethod = data[0];
			}
			else if (headings[1] == "maxIterations" && dataLen == 1) {
				job[cmdCt].maxIt = stoi(data[0]);
			}
			else if (headings[1] == "convergenceTol" && dataLen == 1) {
				job[cmdCt].convTol = stod(data[0]);
			}
			else if (headings[1] == "staticLoadTime" && dataLen == 1) {
				job[cmdCt].staticLoadTime.push_back(stod(data[0]));
				//newCmd->staticLoadTime = stod(data[0]);
			} else if(headings[1] == "thermal" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].thermal = true;
				} else {
					job[cmdCt].thermal = false;
				}
			} else if(headings[1] == "timeStep" && dataLen == 1) {
				job[cmdCt].timeStep = stod(data[0]);
			} else if(headings[1] == "type" && dataLen == 1) {
				job[cmdCt].type = data[0];
			} else if(headings[1] == "numModes" && dataLen == 1) {
				job[cmdCt].numModes = stoi(data[0]);
			}
			else if (headings[1] == "targetEigenvalue" && dataLen == 1) {
				job[cmdCt].tgtEval = stod(data[0]);
			}
			else if (headings[1] == "solnField" && dataLen == 1) {
				job[cmdCt].solnField = data[0];
			} else if(headings[1] == "mode" && dataLen == 1) {
				job[cmdCt].mode = stoi(data[0]);
			} else if(headings[1] == "maxAmplitude" && dataLen == 1) {
				job[cmdCt].maxAmplitude = stod(data[0]);
			} else if(headings[1] == "nodeSet" && dataLen == 1) {
				job[cmdCt].nodeSet = data[0];
			} else if(headings[1] == "fields" && dataLen == 1) {
				job[cmdCt].fields.push_back(data[0]);
			} else if(headings[1] == "timeSteps") {
				if(dataLen == 1) {
					if(data[0] == "all") {
						job[cmdCt].timeStepTag = "all";
					} else if(data[0] == "last") {
						job[cmdCt].timeStepTag = "last";
					} else {
						try {
							i1 = stoi(data[0]);
							job[cmdCt].timeSteps.push_back(stoi(data[0]));
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
						job[cmdCt].timeSteps.push_back(i1);
					}
				}
			} else if(headings[1] == "elementSet" && dataLen == 1) {
				job[cmdCt].elementSet = data[0];
			}
			else if (headings[1] == "position" && dataLen == 1) {
				job[cmdCt].position = data[0];
			}
			else if (headings[1] == "writeModes" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].writeModes = true;
				} else {
					job[cmdCt].writeModes = false;
				}
			} else if(headings[1] == "properties" && dataLen == 1) {
				job[cmdCt].properties.push_back(data[0]);
			} else if(headings[1] == "include" && dataLen == 1) {
				job[cmdCt].objInclude.push_back(data[0]);
			} else if(headings[1] == "writeGradient" && dataLen == 1) {
				if(data[0] == "yes") {
					job[cmdCt].writeGradient = true;
				} else {
					job[cmdCt].writeGradient = false;
				}
			}
			
			prevLdHd = headings[0];
		}
	} else {
		cout << "entered !inFile block" << endl;
		errSt = "Error: could not open job file: " + fileName;
		throw invalid_argument(errSt);
	}

	inFile.close();
	
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
	int i3;
	int i4;
	double doubInp[10] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	int intInp[10] = {0,0,0,0,0,0,0,0,0,0};
	int elType;
	
	int ndCt = 0;
	int elCt = 0;
	int nsCt = 0;
	int esCt = 0;
	int secCt = 0;
	int matCt = 0;
	int flCt = 0;
	int layCt = 0;
	int setInd = 0;

	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "nodes" && dataLen == 4) {
				ndCt++;
			}
			else if (headings[0] == "elements") {
				if (headings[1] == "connectivity" && dataLen > 1) {
					elCt++;
				}
			}
			else if (headings[0] == "sets") {
				if (headings[1] == "node") {
					if (headings[2] == "name" && dataLen == 1) {
						nsCt++;
					}
				}
				else if (headings[1] == "element") {
					if (headings[2] == "name" && dataLen == 1) {
						esCt++;
					}
				}
			}
			else if (headings[0] == "sections") {
				if (headings[1] == "type" && dataLen == 1) {
					secCt++;
				}
			}
			else if (headings[0] == "materials") {
				if (headings[1] == "name" && dataLen == 1) {
					matCt++;
				}
			}
			else if (headings[0] == "fluids") {
				if (headings[1] == "name" && dataLen == 1) {
					flCt++;
				}
			}
		}
	}

	nodes = vector<Node>(ndCt);
	elements = vector<Element>(elCt);
	nodeSets = vector<Set>(nsCt + ndCt + 1);
	elementSets = vector<Set>(esCt + elCt + 1);
	sections = vector<Section>(secCt);
	materials = vector<Material>(matCt);
	fluids = vector<Fluid>(flCt);

	//inFile.seekg(0, std::ios::beg);
	inFile.close();
	inFile.open(fileName);
	
	if(inFile) {
		ndCt = -1;
		elCt = -1;
		nsCt = -1;
		esCt = -1;
		secCt = -1;
		matCt = -1;
		flCt = -1;
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "nodes" && dataLen == 4) {
				ndCt = stoi(data[0]);
				nodes[ndCt].label = ndCt;
				doubInp[0] = stod(data[1]);
				doubInp[1] = stod(data[2]);
				doubInp[2] = stod(data[3]);
				nodes[ndCt].setCrd(doubInp);
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
					}
					else if (data[0] == "tet10") {
						elType = 10;
					}
					else if (data[0] == "shell3") {
						elType = 3;
					} else if(data[0] == "shell4") {
						elType = 41;
					} else if(data[0] == "beam2") {
						elType = 2;
					}
					else if (data[0] == "frcFld") {
						elType = 21;
					}
					else if (data[0] == "mass") {
						elType = 1;
					}
					else if (data[0] == "fl4") {
						elType = 400;
					}
					else if (data[0] == "fl6") {
						elType = 600;
					}
					else if (data[0] == "fl8") {
						elType = 800;
					}
					else if (data[0] == "fl10") {
						elType = 1000;
					}
					else {
						string errSt = "Error: unrecognized element type, " + data[0];
						throw invalid_argument(errSt);
					}
				} else if(headings[1] == "connectivity" && dataLen > 1) {
					elCt = stoi(data[0]);
					Element& newEl = elements[elCt];
					newEl.initializeType(elType);
					newEl.label = elCt;
					for (i1 = 1; i1 < dataLen; i1++) {
						intInp[i1-1] = stoi(data[i1]);
					}
					newEl.setNodes(intInp);
				} 
			} else if(headings[0] == "sets") {
				if(headings[1] == "node") {
					if(headings[2] == "name" && dataLen == 1) {
						nsCt++;
						Set& newSet = nodeSets[nsCt];
						newSet.name = data[0];
					} else if(headings[2] == "labels") {
						if(dataLen == 1) {
							nodeSets[nsCt].labels.push_back(stoi(data[0]));
						} else if(dataLen > 1) {
							i2 = stoi(data[0]);
							i3 = stoi(data[1]);
							if(dataLen == 3) {
								i4 = stoi(data[2]);
							} else {
								i4 = 1;
							}
							for (i1 = i2; i1 < i3; i1+= i4) {
								nodeSets[nsCt].labels.push_back(i1);
							}
						}
					}
				} else if(headings[1] == "element") {
					if(headings[2] == "name" && dataLen == 1) {
						esCt++;
						Set& newSet = elementSets[esCt];
						newSet.name = data[0];
					} else if(headings[2] == "labels") {
						if(dataLen == 1) {
							elementSets[esCt].labels.push_back(stoi(data[0]));
						} else if(dataLen > 1) {
							i2 = stoi(data[0]);
							i3 = stoi(data[1]);
							if(dataLen == 3) {
								i4 = stoi(data[2]);
							} else {
								i4 = 1;
							}
							for (i1 = i2; i1 < i3; i1+= i4) {
								elementSets[esCt].labels.push_back(i1);
							}
						}
					}
				}
			} else if(headings[0] == "sections") {
				if(headings[1] == "type" && dataLen == 1) {
					secCt++;
					sections[secCt].type = data[0];
					layCt = -1;
				} else if(headings[1] == "material" && dataLen == 1) {
					sections[secCt].matName = data[0];
				}
				else if (headings[1] == "fluid" && dataLen == 1) {
					sections[secCt].flName = data[0];
				}
				else if (headings[1] == "orientation" && dataLen == 6) {
					for (i1 = 0; i1 < 6; i1++) {
						doubInp[i1] = stod(data[i1]);
					}
					sections[secCt].setOrientation(doubInp);
				} else if(headings[1] == "layup") {
					if(headings[2] == "zOffset" && dataLen == 1) {
						sections[secCt].zOffset = stod(data[0]);
					} else if(headings[2] == "layers") {
						if(headings[3] == "material" && dataLen == 1) {
							Layer newLay;
							newLay.matName = data[0];
							sections[secCt].layers.push_back(newLay);
						} else if(headings[3] == "thickness" && dataLen == 1) {
							auto layPt = sections[secCt].layers.end();
							--layPt;
							layPt->thickness = stod(data[0]);
						} else if(headings[3] == "angle" && dataLen == 1) {
							auto layPt = sections[secCt].layers.end();
							--layPt;
							layPt->angle = stod(data[0]);
						}
					}
				} else if(headings[1] == "beamProps") {
					if(headings[2] == "area" && dataLen == 1) {
						sections[secCt].area = stod(data[0]);
					} else if(headings[2] == "I" && dataLen == 5) {
						for(i1 = 0; i1 < 5; i1++) {
							doubInp[i1] = stod(data[i1]);
						}
						sections[secCt].setAreaMoment(doubInp);
					} else if(headings[2] == "J" && dataLen == 1) {
						sections[secCt].polarMoment = stod(data[1]);
					} else if(headings[2] == "stiffness" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						sections[secCt].setStiffness(intInp[0],intInp[1],doubInp[0]);
					} else if(headings[2] == "mass" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						sections[secCt].setMass(intInp[0],intInp[1],doubInp[0]);
					}
					else if (headings[2] == "damping" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						sections[secCt].setDamping(intInp[0], intInp[1], doubInp[0]);
					}
					else if (headings[2] == "expLoadCoef" && dataLen == 6) {
						for(i1 = 0; i1 < 6; i1++) {
							doubInp[i1] = stod(data[i1]);
						}
						sections[secCt].setExpLd(doubInp);
					} else if(headings[2] == "conductivity" && dataLen == 1) {
						sections[secCt].conductivity = stod(data[0]);
					} else if(headings[2] == "specHeat" && dataLen == 1) {
						sections[secCt].specHeat = stod(data[0]);
					}
				}
				else if (headings[1] == "potField") {
					if (headings[2] == "coef" && dataLen == 1) {
						sections[secCt].potCoef = stod(data[0]);
					}
					else if (headings[2] == "exp" && dataLen == 1) {
						sections[secCt].potExp = stod(data[0]);
					}
				}
				else if (headings[1] == "dampField") {
					if (headings[2] == "coef" && dataLen == 1) {
						sections[secCt].dampCoef = stod(data[0]);
					}
					else if (headings[2] == "exp" && dataLen == 1) {
						sections[secCt].dampExp = stod(data[0]);
					}
				}
				else if (headings[1] == "thermField") {
					if (headings[2] == "condCoef" && dataLen == 1) {
						sections[secCt].condCoef = stod(data[0]);
					}
					else if (headings[2] == "radCoef" && dataLen == 1) {
						sections[secCt].radCoef = stod(data[0]);
					}
					else if (headings[2] == "refTemp" && dataLen == 1) {
						sections[secCt].refTemp = stod(data[0]);
					}
				}
				else if (headings[1] == "massPerEl" && dataLen == 1) {
					sections[secCt].massPerEl = stod(data[0]);
				}
				else if (headings[1] == "specHeat" && dataLen == 1) {
					sections[secCt].specHeat = stod(data[0]);
				}
				else if (headings[1] == "fluidParams") {
					if (headings[2] == "denVisCoef" && dataLen == 1) {
						sections[secCt].denVisCoef = stod(data[0]);
					}
					else if (headings[2] == "tempVisCoef" && dataLen == 1) {
						sections[secCt].tempVisCoef = stod(data[0]);
					}
					else if (headings[2] == "turbVisCoef" && dataLen == 1) {
						sections[secCt].turbVisCoef = stod(data[0]);
					}
					else if (headings[2] == "gradVTurbCoef" && dataLen == 1) {
						sections[secCt].gradVTurbCoef = stod(data[0]);
					}
					else if (headings[2] == "dissTurbCoef" && dataLen == 1) {
						sections[secCt].dissTurbCoef = stod(data[0]);
					}
					else if (headings[2] == "enthCoef" && dataLen == 1) {
						sections[secCt].enthCoef = stod(data[0]);
					}
					else if (headings[2] == "enthExp" && dataLen == 1) {
						sections[secCt].enthExp = stod(data[0]);
					}
					else if (headings[2] == "presCoef" && dataLen == 1) {
						sections[secCt].presCoef = stod(data[0]);
					}
					else if (headings[2] == "presExp" && dataLen == 1) {
						sections[secCt].presExp = stod(data[0]);
					}
					else if (headings[2] == "refTemp" && dataLen == 1) {
						sections[secCt].refTemp = stod(data[0]);
					}
					else if (headings[2] == "refDen" && dataLen == 1) {
						sections[secCt].refDen = stod(data[0]);
					}
					else if (headings[2] == "refTurbE" && dataLen == 1) {
						sections[secCt].refTurbE = stod(data[0]);
					}
					else if (headings[2] == "refEnth" && dataLen == 1) {
						sections[secCt].refEnth = stod(data[0]);
					}
				}
				else if(headings[1] == "elementSet" && dataLen == 1) {
					sections[secCt].elSetName = data[0];
				}
			} else if(headings[0] == "materials") {
				if(headings[1] == "name" && dataLen == 1) {
					matCt++;
					materials[matCt].name = data[0];
				} else if(headings[1] == "density" && dataLen == 1) {
					materials[matCt].density = stod(data[0]);
				} else if(headings[1] == "elastic") {
					if(headings[2] == "E" && dataLen > 0) {
						Material& thisMat = materials[matCt];
						thisMat.modulus[0] = stod(data[0]);
						if(dataLen == 3) {
							thisMat.modulus[1] = stod(data[1]);
							thisMat.modulus[2] = stod(data[2]);
						} else {
							thisMat.modulus[1] = thisMat.modulus[0];
							thisMat.modulus[2] = thisMat.modulus[0];
						}
					} else if(headings[2] == "nu" && dataLen > 0) {
						Material& thisMat = materials[matCt];
						thisMat.poissonRatio[0] = stod(data[0]);
						if(dataLen == 3) {
							thisMat.poissonRatio[1] = stod(data[1]);
							thisMat.poissonRatio[2] = stod(data[2]);
						} else {
							thisMat.poissonRatio[1] = thisMat.poissonRatio[0];
							thisMat.poissonRatio[2] = thisMat.poissonRatio[0];
						}
					} else if(headings[2] == "G" && dataLen > 0) {
						Material& thisMat = materials[matCt];
						thisMat.shearMod[0] = stod(data[0]);
						if(dataLen == 3) {
							thisMat.shearMod[1] = stod(data[1]);
							thisMat.shearMod[2] = stod(data[2]);
						} else {
							thisMat.shearMod[1] = thisMat.shearMod[0];
							thisMat.shearMod[2] = thisMat.shearMod[0];
						}
					} else if(headings[2] == "stiffness" && dataLen == 3) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						doubInp[0] = stod(data[2]);
						materials[matCt].setStiffness(intInp[0],intInp[1],doubInp[0]);
					} 
				} else if(headings[1] == "thermal") {
					if(headings[2] == "conductivity" && dataLen == 6) {
						Material& thisMat = materials[matCt];
						for(i1 = 0; i1 < 6; i1++) {
							thisMat.conductivity[i1] = stod(data[i1]);
						}
					} else if(headings[2] == "expansion" && dataLen == 6) {
						Material& thisMat = materials[matCt];
						for(i1 = 0; i1 < 6; i1++) {
							thisMat.expansion[i1] = stod(data[i1]);
						}
					} else if (headings[2] == "specHeat" && dataLen == 1) {
						materials[matCt].specHeat = stod(data[0]);
					}
				}
				else if (headings[1] == "damping" && dataLen == 3) {
					intInp[0] = stoi(data[0]) - 1;
					intInp[1] = stoi(data[1]) - 1;
					doubInp[0] = stod(data[2]);
					materials[matCt].setDamping(intInp[0], intInp[1], doubInp[0]);
				}
				else if (headings[1] == "failure") {
					if(headings[2] == "maxStress") {
						if(headings[3] == "tensile" && dataLen > 0) {
							Material& nMat = materials[matCt];
							nMat.maxTenStress[0] = stod(data[0]);
							if(dataLen == 3) {
								nMat.maxTenStress[1] = stod(data[1]);
								nMat.maxTenStress[2] = stod(data[2]);
							} else {
								nMat.maxTenStress[1] = nMat.maxTenStress[0];
								nMat.maxTenStress[2] = nMat.maxTenStress[0];
							}
						} else if(headings[3] == "compressive" && dataLen > 0) {
							Material& nMat = materials[matCt];
							nMat.maxCompStress[0] = stod(data[0]);
							if(dataLen == 3) {
								nMat.maxCompStress[1] = stod(data[1]);
								nMat.maxCompStress[2] = stod(data[2]);
							} else {
								nMat.maxCompStress[1] = nMat.maxCompStress[0];
								nMat.maxCompStress[2] = nMat.maxCompStress[0];
							}
						} else if(headings[3] == "shear" && dataLen > 0) {
							Material& nMat = materials[matCt];
							nMat.maxShearStress[0] = stod(data[0]);
							if(dataLen == 3) {
								nMat.maxShearStress[1] = stod(data[1]);
								nMat.maxShearStress[2] = stod(data[2]);
							} else {
								nMat.maxShearStress[1] = nMat.maxShearStress[0];
								nMat.maxShearStress[2] = nMat.maxShearStress[0];
							}
						}
					} else if(headings[2] == "maxStrain") {
						if(headings[3] == "tensile" && dataLen > 0) {
							Material& nMat = materials[matCt];
							nMat.maxTenStrain[0] = stod(data[0]);
							if(dataLen == 3) {
								nMat.maxTenStrain[1] = stod(data[1]);
								nMat.maxTenStrain[2] = stod(data[2]);
							} else {
								nMat.maxTenStrain[1] = nMat.maxTenStrain[0];
								nMat.maxTenStrain[2] = nMat.maxTenStrain[0];
							}
						} else if(headings[3] == "compressive" && dataLen > 0) {
							Material& nMat = materials[matCt];
							nMat.maxCompStrain[0] = stod(data[0]);
							if(dataLen == 3) {
								nMat.maxCompStrain[1] = stod(data[1]);
								nMat.maxCompStrain[2] = stod(data[2]);
							} else {
								nMat.maxCompStrain[1] = nMat.maxCompStrain[0];
								nMat.maxCompStrain[2] = nMat.maxCompStrain[0];
							}
						} else if(headings[3] == "shear" && dataLen > 0) {
							Material& nMat = materials[matCt];
							nMat.maxShearStrain[0] = stod(data[0]);
							if(dataLen == 3) {
								nMat.maxShearStrain[1] = stod(data[1]);
								nMat.maxShearStrain[2] = stod(data[2]);
							} else {
								nMat.maxShearStrain[1] = nMat.maxShearStrain[0];
								nMat.maxShearStrain[2] = nMat.maxShearStrain[0];
							}
						}
					} else if(headings[2] == "maxStrainEnergy" && dataLen == 1) {
						materials[matCt].maxStrEng = stod(data[0]);
					} else if(headings[2] == "maxMises" && dataLen == 1) {
						materials[matCt].maxMises = stod(data[0]);
					}
				}
			}
			else if (headings[0] == "fluids") {
				if (headings[1] == "name" && dataLen == 1) {
					flCt++;
					fluids[flCt].name = data[0];
				}
				else if (headings[1] == "viscosity" && dataLen == 1) {
					fluids[flCt].viscosity = stod(data[0]);
				}
				else if (headings[1] == "thermal") {
					if (headings[2] == "conductivity" && dataLen == 1) {
						fluids[flCt].thermCond = stod(data[0]);
					}
					else if (headings[2] == "specHeat" && dataLen == 1) {
						fluids[flCt].specHeat = stod(data[0]);
					}
				}
				else if (headings[1] == "idealGasConst" && dataLen == 1) {
					fluids[flCt].idealGas = stod(data[0]);
				}
			}
		}
	} else {
		string errSt = "Error: could not open model input file: " + fileName;
		throw invalid_argument(errSt);
	}

	// Create all and individual nd/element sets

	i1 = nodes.size();
	for (i2 = 0; i2 < i1; i2++) {
		nsCt++;
		nodeSets[nsCt].name = to_string(i2);
		nodeSets[nsCt].labels.push_back(i2);
	}

	nsCt++;
	nodeSets[nsCt].name = "all";
	list<int>& labs = nodeSets[nsCt].labels;
	for (i2 = 0; i2 < i1; i2++) {
		labs.push_back(i2);
	}

	i1 = elements.size();
	for (i2 = 0; i2 < i1; i2++) {
		esCt++;
		elementSets[esCt].name = to_string(i2);
		elementSets[esCt].labels.push_back(i2);
	}

	esCt++;
	Set& newSet = elementSets[esCt];
	newSet.name = "all";
	for (i2 = 0; i2 < i1; i2++) {
		newSet.labels.push_back(i2);
	}

	// Populate set arrays

	i1 = nodeSets.size();

	i2 = 0;
	for (auto& ns : nodeSets) {
		nsMap.insert(make_pair(ns.name, i2));
		i2++;
	}

	i1 = elementSets.size();

	i2 = 0;
	for (auto& es : elementSets) {
		esMap.insert(make_pair(es.name, i2));
		i2++;
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
	string allTypes = "displacement temperature";
	string errSt;
	
	int ecCt = 0;
	int tcCt = 0;
	string thisType;
	int currTm;

	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "constraints") {
				if (headings[1] == "type" && dataLen == 1) {
					if (data[0] == "displacement") {
						ecCt++;
					}
					else if (data[0] == "temperature") {
						tcCt++;
					}
				}
			}
		}
	}

	elasticConst.constVec = vector<Constraint>(ecCt);
	thermalConst.constVec = vector<Constraint>(tcCt);

	//inFile.seekg(0, std::ios::beg);
	inFile.close();
	inFile.open(fileName);

	Constraint* newCon = nullptr;
	if(inFile) {
		ecCt = -1;
		tcCt = -1;
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "constraints") {
				if (headings[1] == "type" && dataLen == 1) {
					if (data[0] == "displacement") {
						ecCt++;
						newCon = &elasticConst.constVec[ecCt];
						newCon->type = "displacement";
					}
					else if (data[0] == "temperature") {
						tcCt++;
						newCon = &thermalConst.constVec[tcCt];
						newCon->type = "temperature";
					}
					else {
						errSt = "Error: " + data[0] + " is not a valid constraint type.  Allowable values are " + allTypes;
						throw invalid_argument(errSt);
					}
				}
				else if(headings[1] == "terms") {
					if(headings[2] == "nodeSet" && dataLen == 1) {
						ConstraintTerm newCn;
						newCn.nodeSet = data[0];
						newCon->terms.push_back(newCn);
					} else if(headings[2] == "dof" && dataLen == 1) {
						auto tmPt = newCon->terms.end();
						--tmPt;
						tmPt->dof = stoi(data[0]);
					} else if(headings[2] == "coef" && dataLen == 1) {
						auto tmPt = newCon->terms.end();
						--tmPt;
						tmPt->coef = stod(data[0]);
					}
				} else if(headings[1] == "rhs" && dataLen == 1) {
					newCon->rhs = stod(data[0]);
				}
			}
		}
	} else {
		errSt = "Error: could not open constraint input file: " + fileName;
		throw invalid_argument(errSt);
	}

	inFile.close();

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
	double doubInp[10] = { 0,0,0,0,0,0,0,0,0,0 };
	string elasticList = "nodalForce bodyForce gravitational centrifugal surfacePressure surfaceTraction";
	string thermalList = "bodyHeadGen surfaceFlux";
	Load *newLd = nullptr;

	int eLdCt = 0;
	int tLdCt = 0;
	
	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "loads") {
				if (headings[1] == "type" && dataLen == 1) {
					i1 = elasticList.find(data[0]);
					if (i1 > -1) {
						eLdCt++;
					}
					i1 = thermalList.find(data[0]);
					if (i1 > -1) {
						tLdCt++;
					}
				}
			}
		}
	}

	elasticLoads = vector<Load>(eLdCt);
	thermalLoads = vector<Load>(tLdCt);

	//inFile.seekg(0, std::ios::beg);
	inFile.close();
	inFile.open(fileName);
	if(inFile) {
		eLdCt = -1;
		tLdCt = -1;
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "loads") {
				if(headings[1] == "type" && dataLen == 1) {
					i1 = elasticList.find(data[0]);
					if (i1 > -1) {
						eLdCt++;
						elasticLoads[eLdCt].type = data[0];
						newLd = &elasticLoads[eLdCt];
					}
					i1 = thermalList.find(data[0]);
					if(i1 > -1) {
						tLdCt++;
						thermalLoads[tLdCt].type = data[0];
						newLd = &thermalLoads[tLdCt];
					}
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
					newLd->nodeSet = data[0];
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					newLd->elementSet = data[0];
				} else if(headings[1] == "normDir" && dataLen == 3) {
					doubInp[0] = stod(data[0]);
					doubInp[1] = stod(data[1]);
					doubInp[2] = stod(data[2]);
					newLd->setNormDir(doubInp);
				} else if(headings[1] == "normTolerance" && dataLen == 1) {
					newLd->normTol = stod(data[0]);
				} else if(headings[1] == "center" && dataLen == 3) {
					doubInp[0] = stod(data[0]);
					doubInp[1] = stod(data[1]);
					doubInp[2] = stod(data[2]);
					newLd->setCenter(doubInp);
				} else if(headings[1] == "axis" && dataLen == 3) {
					doubInp[0] = stod(data[0]);
					doubInp[1] = stod(data[1]);
					doubInp[2] = stod(data[2]);
					newLd->setAxis(doubInp);
				} else if(headings[1] == "angularVelocity") {
					newLd->angularVel = stod(data[0]);
				}
			}
		}
	} else {
		string errSt = "Error: could not open load input file: " + fileName;
		throw invalid_argument(errSt);
	}
	inFile.close();
	
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
	int seti;
	double doubInp[10] = { 0,0,0,0,0,0,0,0,0,0 };
	string dispHdings = " displacement velocity acceleration";
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "initialState") {
				i3 = dispHdings.find(headings[1]);
				if(i3 > -1 && dataLen > 3) {
					seti = nsMap.at(data[0]);
					for (auto& ndi : nodeSets[seti].labels) {
						Node& thisNd = nodes[ndi];
						i2 = 1;
						for (i1 = 0; i1 < 6; i1++) {
							if (i2 < dataLen) {
								doubInp[i1] = stod(data[i2]);
							}
							else {
								doubInp[i1] = 0.0;
							}
							i2++;
						}
						if (headings[1] == "displacement") {
							thisNd.setInitialDisp(doubInp);
						}
						else if (headings[1] == "velocity") {
							thisNd.setInitialVel(doubInp);
						}
						else if (headings[1] == "acceleration") {
							thisNd.setInitialAcc(doubInp);
						}
					}
				} else if(headings[1] == "temperature" && dataLen == 2) {
					seti = nsMap.at(data[0]);
					for (auto& ndi : nodeSets[seti].labels) {
						nodes[ndi].initialTemp = stod(data[1]);
					}
				} else if(headings[1] == "tdot" && dataLen == 2) {
					seti = nsMap.at(data[0]);
					for (auto& ndi : nodeSets[seti].labels) {
						nodes[ndi].initialTdot = stod(data[1]);
					}
				}
			}
		}
	} else {
		string errSt = "Error: could not open initial state input file: " + fileName;
		throw invalid_argument(errSt);
	}

	inFile.close();
	
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
	double doubInp[10] = {0,0,0,0,0,0,0,0,0,0};
	int intInp[10] = {0,0,0,0,0,0,0,0,0,0};
	
	int dvCt = 0;
	
	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "designVariables") {
				if (headings[1] == "category" && dataLen == 1) {
					dvCt++;
				}
			}
		}
	}

	designVars = vector<DesignVariable>(dvCt);

	//inFile.seekg(0, std::ios::beg);
	inFile.close();
	inFile.open(fileName);

	if(inFile) {
		dvCt = -1;
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "designVariables") {
				if(headings[1] == "category" && dataLen == 1) {
					dvCt++;
					designVars[dvCt].category = data[0];
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					designVars[dvCt].elSetName = data[0];
				} else if(headings[1] == "nodeSet" && dataLen == 1) {
					designVars[dvCt].ndSetName = data[0];
				} else if(headings[1] == "activeTime" && dataLen > 0) {
					doubInp[0] = stod(data[0]);
					if(dataLen == 2) {
						doubInp[1] = stod(data[1]);
					} else {
						doubInp[1] = 1.0e+100;
					}
					designVars[dvCt].setActiveTime(doubInp);
				} else if(headings[1] == "component") {
					if(dataLen == 1) {
					    designVars[dvCt].component = stoi(data[0]);
					} else if(dataLen == 2) {
						intInp[0] = stoi(data[0]) - 1;
						intInp[1] = stoi(data[1]) - 1;
						if(intInp[0] >= intInp[1]) {
							i1 = 6*intInp[1] + intInp[0];
						} else {
							i1 = 6*intInp[0] + intInp[1];
						}
						designVars[dvCt].component = i1;
					}
				} else if(headings[1] == "layer" && dataLen == 1) {
					designVars[dvCt].layer = stoi(data[0]);
				} else if(headings[1] == "coefficients" && dataLen == 1) {
					designVars[dvCt].coefs.push_back(stod(data[0]));
				}
			}
		}
	} else {
		cout << "entering !inFile block" << endl;
		string errSt = "Error: could not open design variable input file: " + fileName;
		throw invalid_argument(errSt);
	}
	inFile.close();
	
	return;
}

void Model::readObjectiveInput(string fileName) {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[11];
	int dataLen;
	
	double doubInp[10] = { 0,0,0,0,0,0,0,0,0,0 };
	
	int obCt = 0;
	
	inFile.open(fileName);

	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "objectiveTerms") {
				if (headings[1] == "category" && dataLen == 1) {
					obCt++;
				}
			}
		}
	}

	objective.terms = vector<ObjectiveTerm>(obCt);

	//inFile.seekg(0, std::ios::beg);
	inFile.close();
	inFile.open(fileName);

	if(inFile) {
		obCt = -1;
		while(!inFile.eof()) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(headings[0] == "objectiveTerms") {
				if(headings[1] == "category" && dataLen == 1) {
					obCt++;
					objective.terms[obCt].category = data[0];
				} else if(headings[1] == "operator" && dataLen == 1) {
					objective.terms[obCt].optr = data[0];
				} else if(headings[1] == "activeTime" && dataLen > 0) {
					doubInp[0] = stod(data[0]);
					if(dataLen == 2) {
						doubInp[1] = stod(data[1]);
					} else {
						doubInp[1] = 1.0e+100;
					}
					objective.terms[obCt].setActiveTime(doubInp);
				} else if(headings[1] == "component" && dataLen == 1) {
					objective.terms[obCt].component = stoi(data[0]);
				} else if(headings[1] == "layer" && dataLen == 1) {
					objective.terms[obCt].layer = stoi(data[0]);
				} else if(headings[1] == "coefficient" && dataLen == 1) {
					objective.terms[obCt].coef = stod(data[0]);
				} else if(headings[1] == "exponent" && dataLen == 1) {
					objective.terms[obCt].expnt = stod(data[0]);
				} else if(headings[1] == "elementSet" && dataLen == 1) {
					objective.terms[obCt].elSetName = data[0];
				} else if(headings[1] == "nodeSet" && dataLen == 1) {
					objective.terms[obCt].ndSetName = data[0];
				} else if(headings[1] == "targetValue" && dataLen == 1) {
					try {
						doubInp[0] = stod(data[0]);
						objective.terms[obCt].tgtVals.push_back(doubInp[0]);
					} catch(...) {
						objective.terms[obCt].tgtTag = data[0];
					}
				}
			}
		}
	} else {
		string errSt = "Error: could not open objective input file: " + fileName;
		throw invalid_argument(errSt);
	}

	inFile.close();
	
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
	
	inFile.open(fileName);
	if(inFile) {
		while(!inFile.eof()) {
		    readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			if(dataLen == 2) {
				label = stoi(data[0]);
				value = stod(data[1]);
				designVars[label].value.setVal(value);
			}
		}
	} else {
		string errSt = "Error: could not open design variable value input file: " + fileName;
		throw invalid_argument(errSt);		
	}

	inFile.close();
	
	return;
}

void Model::readNodeResults(string fileName) {
	int i1;
	int i2;
	int nd;
	double ndDat[6];
	Node* thisNd;
	ifstream inFile;
	string fileLine;
	string headings[4] = { "","","","" };
	int hdLdSpace[4] = { 0,0,0,0 };
	string data[11];
	int dataLen;
	string dispFields = "displacement velocity acceleration";
	string thrmFields = "temperature tdot";

	inFile.open(fileName);
	if (inFile) {
		while (!inFile.eof()) {
			readInputLine(inFile, fileLine, headings, hdLdSpace, data, dataLen);
			if (headings[0] == "nodeResults" && dataLen > 1) {
				nd = stoi(data[0]);
				Node& thisNd = nodes[nd];
				i2 = dispFields.find(headings[1]);
				if (i2 > -1) {
					for (i1 = 0; i1 < thisNd.numDof; i1++) {
						ndDat[i1] = stod(data[i1 + 1]);
					}
					if (headings[1] == "displacement") {
						thisNd.setDisplacement(ndDat);
					}
					else if (headings[1] == "velocity") {
						thisNd.setVelocity(ndDat);
					}
					else {
						thisNd.setAcceleration(ndDat);
					}
				}
				i2 = thrmFields.find(headings[1]);
				if (i2 > -1) {
					ndDat[0] = stod(data[1]);
					if (headings[1] == "temperature") {
						thisNd.temperature = ndDat[0];
					}
					else {
						thisNd.tempChangeRate = ndDat[0];
					}
				}
			}
		}
	}
	else {
		string errSt = "Error: could not open node result input file: " + fileName;
		throw invalid_argument(errSt);
	}

	inFile.close();
	
	return;
}

void Model::readTimeStepSoln(int tStep) {
	int i1;
	int dofPerNd;
	int numIntDof;
	double* inDat;
	char inLn[72];
	string fullFile = job[solveCmd].fileName + "/solnTStep" + to_string(tStep) + ".out";
	ifstream inFile;
	inFile.open(fullFile, std::ifstream::binary);

	inDat = reinterpret_cast<double*>(&inLn[0]);
	for (auto& nd : nodes) {
		if (job[solveCmd].thermal) {
			inFile.read(&inLn[0], 8);
			nd.prevTemp = *inDat;
			inFile.read(&inLn[0], 8);
			nd.prevTdot = *inDat;
		}
		if (job[solveCmd].elastic) {
			dofPerNd = nd.numDof;
			i1 = 8 * dofPerNd;
			inFile.read(&inLn[0], i1);
			nd.setPrevDisp(inDat);
			inFile.read(&inLn[0], i1);
			nd.setPrevVel(inDat);
			inFile.read(&inLn[0], i1);
			nd.setPrevAcc(inDat);
		}
	}

	if (job[solveCmd].elastic) {
		for (auto& el : elements) {
			numIntDof = el.numIntDof;
			if (numIntDof > 0) {
				i1 = 8 * numIntDof;
				inFile.read(&inLn[0], i1);
				el.setIntPrevDisp(inDat);
			}
		}
	}

	inFile.close();
	return;
}