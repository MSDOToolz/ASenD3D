#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void readInputLine(ifstream& inFile, string& fileLine, string headings[], int hdLdSpace[], string data[], int& dataLen);

int main() {
	ifstream inFile;
	string fileLine;
	string headings[4] = {"","","",""};
	int hdLdSpace[4] = {0,0,0,0};
	string data[10];
	int dataLen;
	int i1;
	int it;
	
	inFile.open("testYaml.yaml");
	if(inFile) {
		it = 0;
		while(!inFile.eof() && it < 100) {
			readInputLine(inFile,fileLine,headings,hdLdSpace,data,dataLen);
			cout << headings[0] << " " << headings[1] << " " << headings[2] << " " << headings[3] << " : ";
			for(i1 = 0; i1 < dataLen; i1++) {
				cout << data[i1] << ", ";
			}
			cout << endl;
			it++;
		}
	}
	inFile.close();
	return 0;
}

void readInputLine(ifstream& inFile, string& fileLine, string headings[], int hdLdSpace[], string data[], int& dataLen) {
	int i1;
	int i2;
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
		while(i1 < (lnLen)) {
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
			while(i1 < (lnLen)) {
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