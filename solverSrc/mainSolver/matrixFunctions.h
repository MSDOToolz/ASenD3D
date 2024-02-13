#ifndef MATRIXFUN
#define MATRIXFUN
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "LowerTriMatClass.h"
#include "ConstraintClass.h"

double getDist(double p1[], double p2[]);

void crossProd(double prod[], double v1[], double v2[]);

void qRFactor(double mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(double xVec[], double mat[], double bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void gMResSparse(double soln[], SparseMat& mat, ConstraintList& cnst, LowerTriMat& pcMat, double rhs[], double convTol, int maxIt, int restart);

void symFactor(double mat[], double qMat[], int matDim);

void getCharFun(DiffDoub1& cFun, DiffDoub1 mat[], int matDim, DoubList& eVals, double lam, int triDiag);

void getEvals(double eVals[], double mat[], int matDim, double lamInit, double convTol, int triDiag);

double getLowEval(double mat[], int matDim);

void eigenSolve(double eVals[], double eVecs[], double mat[], int matDim, int triDiag);

void symEigenSolve(double eVals[], double eVecs[], double mat[], int matDim, int triDiag);

void eigenFull(double eVals[], double eVecs[], int numPairs, LowerTriMat& mat, double massMat[], int matDim);

void eigenSparseDirect(double eVals[], double eVecs[], int numPairs, LowerTriMat& mat, double massMat[], int matDim);

//dup1
void qRFactor(DiffDoub0 mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(DiffDoub0 xVec[], DiffDoub0 mat[], DiffDoub0 bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(DiffDoub0& det, DiffDoub0 inv[], DiffDoub0 mat[], int colDim, int triDiag, DiffDoub0 xVec[], DiffDoub0 bVec[]);
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1
void qRFactor(DiffDoub1 mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(DiffDoub1 xVec[], DiffDoub1 mat[], DiffDoub1 bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(DiffDoub1& det, DiffDoub1 inv[], DiffDoub1 mat[], int colDim, int triDiag, DiffDoub1 xVec[], DiffDoub1 bVec[]);
//end dup
 
//end skip 
 
 
 
 
 
 
//dup2
void matMul(DiffDoub0 prod[], DiffDoub0 mat1[], DiffDoub0 mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(DiffDoub0 matT[], DiffDoub0 mat[], int rowDim, int colDim);

void crossProd(DiffDoub0 prod[], DiffDoub0 v1[], DiffDoub0 v2[]);

void rotateOrient(DiffDoub0 instOri[], DiffDoub0 locOri[], DiffDoub0 rot[]);
//end dup
 
//skip 
 
//DiffDoub versions: 
//dup2
void matMul(DiffDoub1 prod[], DiffDoub1 mat1[], DiffDoub1 mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(DiffDoub1 matT[], DiffDoub1 mat[], int rowDim, int colDim);

void crossProd(DiffDoub1 prod[], DiffDoub1 v1[], DiffDoub1 v2[]);

void rotateOrient(DiffDoub1 instOri[], DiffDoub1 locOri[], DiffDoub1 rot[]);
//end dup
 
//Diff2Doub versions: 
//dup2
void matMul(DiffDoub2 prod[], DiffDoub2 mat1[], DiffDoub2 mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(DiffDoub2 matT[], DiffDoub2 mat[], int rowDim, int colDim);

void crossProd(DiffDoub2 prod[], DiffDoub2 v1[], DiffDoub2 v2[]);

void rotateOrient(DiffDoub2 instOri[], DiffDoub2 locOri[], DiffDoub2 rot[]);
//end dup
 
//end skip 
 
 
 
 
 
 
//dup1

void dOridThet(DiffDoub0 instOri[], DiffDoub0 locOri[], DiffDoub0 rot[], int v1, int v2);

//end dup
 
//skip 
 
//DiffDoub versions: 
//dup1

void dOridThet(DiffDoub1 instOri[], DiffDoub1 locOri[], DiffDoub1 rot[], int v1, int v2);

//end dup
 
//end skip 
 
 
 
 
 
 
#endif