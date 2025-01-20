#ifndef MATRIXFUN
#define MATRIXFUN
#include <vector>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "LUMatClass.h"
#include "LowerTriMatClass.h"
#include "ConstraintClass.h"

void subVec(std::vector<double>& subV, std::vector<double>& vIn, int st, int end);

void returnSV(std::vector<double>& subV, std::vector<double>& vIn, int st, int end);

void vecToAr(double ar[], std::vector<double>& vc, int st, int end);

double getDist(double p1[], double p2[]);

void crossProd(double prod[], double v1[], double v2[]);

void qRFactor(std::vector<double>& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(std::vector<double>& xVec, std::vector<double>& mat, std::vector<double>& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void conjGradSparse(std::vector<double>& soln, SparseMat& mat, ConstraintList& cnst, LowerTriMat& pcMat, std::vector<double>& rhs, double convTol, int maxIt);

void gMResSparse(std::vector<double>& soln, SparseMat& mat, ConstraintList& cnst, LUMat& pcMat, std::vector<double>& rhs, double convTol, int maxIt, int restart);

void symFactor(std::vector<double>& mat, std::vector<double>& qMat, int matDim);

void getCharFun(DiffDoub1& cFun, std::vector<DiffDoub1>& mat, int matDim, std::vector<double>& eVals, double lam, int triDiag);

void getEvals(std::vector<double>& eVals, std::vector<double>& mat, int matDim, double lamInit, double convTol, int triDiag);

double getLowEval(std::vector<double>& mat, int matDim);

void eigenSolve(std::vector<double>& eVals, std::vector<double>& eVecs, std::vector<double>& mat, int matDim, int triDiag);

void symEigenSolve(std::vector<double>& eVals, std::vector<double>& eVecs, std::vector<double>& mat, int matDim, int triDiag);

//void eigenFull(std::vector<double>& eVals, std::vector<double>& eVecs, int numPairs, LowerTriMat& mat, std::vector<double>& massMat, int matDim);

void eigenSparseDirect(std::vector<double>& eVals, std::vector<double>& eVecs, int numPairs, LowerTriMat& mat, std::vector<double>& massMat, int matDim);

double rayQuot(std::vector<double>& grad, std::vector<double>& Kv, std::vector<double>& Mv, SparseMat& mat, ConstraintList& cnst, std::vector<double>& massMat, std::vector<double>& inVec);

double unitizeVec(std::vector<double>& vec, int dim);

void getNearestEvecRQ(SparseMat& mat, ConstraintList& cnst, std::vector<double>& massMat, std::vector<double>& inVecs, std::vector<double>& eVals, int numVecs, int maxIt);

//void getNearestEvecSubspace(SparseMat& mat, ConstraintList& cnst, std::vector<double>& massMat, std::vector<double>& inVecs, std::vector<double>& eVals, int numVecs);

//dup1
void qRFactor(std::vector<DiffDoub0>& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void qRFactor(DiffDoub0 mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(std::vector<DiffDoub0>& xVec, std::vector<DiffDoub0>& mat, std::vector<DiffDoub0>& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(DiffDoub0 xVec[], DiffDoub0 mat[], DiffDoub0 bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(DiffDoub0& det, std::vector<DiffDoub0>& inv, std::vector<DiffDoub0>& mat, int colDim, int triDiag, std::vector<DiffDoub0>& xVec, std::vector<DiffDoub0>& bVec);

void getDetInv(DiffDoub0& det, DiffDoub0 inv[], DiffDoub0 mat[], int colDim, int triDiag, DiffDoub0 xVec[], DiffDoub0 bVec[]);
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1
void qRFactor(std::vector<DiffDoub1>& mat, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void qRFactor(DiffDoub1 mat[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(std::vector<DiffDoub1>& xVec, std::vector<DiffDoub1>& mat, std::vector<DiffDoub1>& bVec, int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void solveqRxEqb(DiffDoub1 xVec[], DiffDoub1 mat[], DiffDoub1 bVec[], int colDim, int stRow, int endRow, int stCol, int endCol, int triDiag);

void getDetInv(DiffDoub1& det, std::vector<DiffDoub1>& inv, std::vector<DiffDoub1>& mat, int colDim, int triDiag, std::vector<DiffDoub1>& xVec, std::vector<DiffDoub1>& bVec);

void getDetInv(DiffDoub1& det, DiffDoub1 inv[], DiffDoub1 mat[], int colDim, int triDiag, DiffDoub1 xVec[], DiffDoub1 bVec[]);
//end dup
 
//end skip 
 
 
//dup2
void subVec(std::vector<DiffDoub0>& subV, std::vector<DiffDoub0>& inV, int st, int end);

void returnSV(std::vector<DiffDoub0>& subV, std::vector<DiffDoub0>& inV, int st, int end);

void vecToAr(DiffDoub0 ar[], std::vector<DiffDoub0>& vc, int st, int end);

void arToVec(DiffDoub0 ar[], std::vector<DiffDoub0>& vc, int st, int end);

void matMul(std::vector<DiffDoub0>& prod, std::vector<DiffDoub0>& mat1, std::vector<DiffDoub0>& mat2, int m1Rows, int m1Cols, int m2Cols);

void matMul(DiffDoub0 prod[], DiffDoub0 mat1[], DiffDoub0 mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(std::vector<DiffDoub0>& matT, std::vector<DiffDoub0>& mat, int rowDim, int colDim);

void transpose(DiffDoub0 matT[], DiffDoub0 mat[], int rowDim, int colDim);

void crossProd(DiffDoub0 prod[], DiffDoub0 v1[], DiffDoub0 v2[]);

void rotateOrient(DiffDoub0 instOri[], DiffDoub0 locOri[], DiffDoub0 rot[]);
//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup2
void subVec(std::vector<DiffDoub1>& subV, std::vector<DiffDoub1>& inV, int st, int end);

void returnSV(std::vector<DiffDoub1>& subV, std::vector<DiffDoub1>& inV, int st, int end);

void vecToAr(DiffDoub1 ar[], std::vector<DiffDoub1>& vc, int st, int end);

void arToVec(DiffDoub1 ar[], std::vector<DiffDoub1>& vc, int st, int end);

void matMul(std::vector<DiffDoub1>& prod, std::vector<DiffDoub1>& mat1, std::vector<DiffDoub1>& mat2, int m1Rows, int m1Cols, int m2Cols);

void matMul(DiffDoub1 prod[], DiffDoub1 mat1[], DiffDoub1 mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(std::vector<DiffDoub1>& matT, std::vector<DiffDoub1>& mat, int rowDim, int colDim);

void transpose(DiffDoub1 matT[], DiffDoub1 mat[], int rowDim, int colDim);

void crossProd(DiffDoub1 prod[], DiffDoub1 v1[], DiffDoub1 v2[]);

void rotateOrient(DiffDoub1 instOri[], DiffDoub1 locOri[], DiffDoub1 rot[]);
//end dup
 
//DiffDoub2 versions: 
//dup2
void subVec(std::vector<DiffDoub2>& subV, std::vector<DiffDoub2>& inV, int st, int end);

void returnSV(std::vector<DiffDoub2>& subV, std::vector<DiffDoub2>& inV, int st, int end);

void vecToAr(DiffDoub2 ar[], std::vector<DiffDoub2>& vc, int st, int end);

void arToVec(DiffDoub2 ar[], std::vector<DiffDoub2>& vc, int st, int end);

void matMul(std::vector<DiffDoub2>& prod, std::vector<DiffDoub2>& mat1, std::vector<DiffDoub2>& mat2, int m1Rows, int m1Cols, int m2Cols);

void matMul(DiffDoub2 prod[], DiffDoub2 mat1[], DiffDoub2 mat2[], int m1Rows, int m1Cols, int m2Cols);

void transpose(std::vector<DiffDoub2>& matT, std::vector<DiffDoub2>& mat, int rowDim, int colDim);

void transpose(DiffDoub2 matT[], DiffDoub2 mat[], int rowDim, int colDim);

void crossProd(DiffDoub2 prod[], DiffDoub2 v1[], DiffDoub2 v2[]);

void rotateOrient(DiffDoub2 instOri[], DiffDoub2 locOri[], DiffDoub2 rot[]);
//end dup
 
//end skip 
 
 
//dup1

void dOridThet(DiffDoub0 instOri[], DiffDoub0 locOri[], DiffDoub0 rot[], int v1, int v2);

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

void dOridThet(DiffDoub1 instOri[], DiffDoub1 locOri[], DiffDoub1 rot[], int v1, int v2);

//end dup
 
//end skip 
 
 
#endif