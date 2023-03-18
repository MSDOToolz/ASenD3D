#ifndef MATRIXFUN
#define MATRIXFUN
#include "DiffDoubClass.h"
#include "ListEntClass.h"

void qRFactor(double& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(double& xVec, double& mat, double& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void symFactor(double& mat, double& qMat, matDim);

void getCharFun(DiffDoub& cFun, DiffDoub& mat, int matDim, DoubList& eVals; double lam, int triDiag);

void getEvals(DoubList& eVals, double& mat, int matDim, double lamInit, double convTol, int triDiag);

double getLowEval(double& mat, int matDim);

void eigenSolve(DoubList& eVals, double& eVecs, double& mat, int matDim, int triDiag);

//dup1
void qRFactor(Doub& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(Doub& xVec, Doub& mat, Doub& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(Doub& det, Doub& inv, Doub& mat, int colDim, int triDiag);
//end dup

//dup2
void matMul(Doub& prod, Doub& mat1, Doub& mat2, int m1Rows, int m1Cols, int m2Cols);

void transpose(Doub& matT, Doub& mat, int rowDim, int colDim);

void crossProd(Doub& prod, Doub& v1, Doub& v2);

void rotateOrient(Doub& instOri, Doub& locOri, Doub& rot);
//end dup

//dup1

void dOridThet(Doub& instOri, Doub& locOri, Doub& rot, int v1, int v2);

//end dup
#endif