#ifndef UTILITIES
#define UTILITIES
#include <string>
#include <fstream>

void crossProd(double prod[], double v1[], double v2[]);

void qRFactor(double mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(double xVec[], double mat[], double bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void readInputLine(std::ifstream& inFile, std::string& fileLine, std::string headings[], int hdLdSpace[], std::string data[], int& dataLen);

#endif
