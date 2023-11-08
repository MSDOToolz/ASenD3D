#ifndef MATRIXFUN
#define MATRIXFUN
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "LowerTriMatClass.h"

double getDist(double p1[], double p2[]);

void crossProd(double prod[], double v1[], double v2[]);

void qRFactor(double mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(double xVec[], double mat[], double bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void symFactor(double mat[], double qMat[], int matDim);

void getCharFun(DiffDoub& cFun, DiffDoub mat[], int matDim, DoubList& eVals, double lam, int triDiag);

void getEvals(double eVals[], double mat[], int matDim, double lamInit, double convTol, int triDiag);

double getLowEval(double mat[], int matDim);

void eigenSolve(double eVals[], double eVecs[], double mat[], int matDim, int triDiag);

void eigenSparseDirect(double eVals[], double eVecs[], int numPairs, LowerTriMat& mat, double massMat[], int matDim);

//dup1
void qRFactor(Doub mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(Doub xVec[], Doub mat[], Doub bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(Doub& det, Doub inv[], Doub mat[], int colDim, int triDiag, Doub xVec[], Doub bVec[]);
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1
void qRFactor(DiffDoub mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(DiffDoub xVec[], DiffDoub mat[], DiffDoub bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(DiffDoub& det, DiffDoub inv[], DiffDoub mat[], int colDim, int triDiag, DiffDoub xVec[], DiffDoub bVec[]);
//end dup
 
//end skip 
 
 
//dup2
void matMul(Doub prod[], Doub mat1[], Doub mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(Doub matT[], Doub mat[], int rowDim, int colDim);

void crossProd(Doub prod[], Doub v1[], Doub v2[]);

void rotateOrient(Doub instOri[], Doub locOri[], Doub rot[]);
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup2
void matMul(DiffDoub prod[], DiffDoub mat1[], DiffDoub mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(DiffDoub matT[], DiffDoub mat[], int rowDim, int colDim);

void crossProd(DiffDoub prod[], DiffDoub v1[], DiffDoub v2[]);

void rotateOrient(DiffDoub instOri[], DiffDoub locOri[], DiffDoub rot[]);
//end dup
 
//Diff2Doub versions: 
//dup2
void matMul(Diff2Doub prod[], Diff2Doub mat1[], Diff2Doub mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(Diff2Doub matT[], Diff2Doub mat[], int rowDim, int colDim);

void crossProd(Diff2Doub prod[], Diff2Doub v1[], Diff2Doub v2[]);

void rotateOrient(Diff2Doub instOri[], Diff2Doub locOri[], Diff2Doub rot[]);
//end dup
 
//end skip 
 
 
//dup1

void dOridThet(Doub instOri[], Doub locOri[], Doub rot[], int v1, int v2);

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1

void dOridThet(DiffDoub instOri[], DiffDoub locOri[], DiffDoub rot[], int v1, int v2);

//end dup
 
//end skip 
 
 
 
 
 
 
 
 
#endif