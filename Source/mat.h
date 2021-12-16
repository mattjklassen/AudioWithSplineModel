/*
  ==============================================================================

    mat.h
    Created: 13 Nov 2021 1:48:59pm
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "../JuceLibraryCode/JuceHeader.h"

using namespace std;

void fillbmat(int k, int d, double *bmat, double *knots);
void fillmat(int k, int d, double *mat, double *inputs, double *knots);
void gaussElim(int n, Array<float>& A, Array<float>& B);
void gausselim(int n, double *mat, double *matinv);
int  getPivot(int n, int j, Array<float>& A);
int  getpivot(int n, int j, double *mat);
void getZero(int n, int i, int j, Array<float>& A, Array<float>& B);
void getzero(int n, int i, int j, double *mat, double *matinv);
void initbmat(int k, int d, double *bmat);
void matcopy(int n, double *mat, double *temp);
void printmat(int n, double *mat);
void printMatrix(int n, Array<float>& A);
void printmult(int n, double *mat, double *matinv);
void printMult(int n, Array<float>& A, Array<float>& B);
void rowMult(int n, int i, float C, Array<float>& A, Array<float>& B);
void rowmult(int n, int i, double C, double *mat, double *matinv);
void rowSwap(int n, int m, int r, Array<float>& A, Array<float>& B);
void rowswap(int n, int m, int r, double *mat, double *matinv);
void setidentity(int n, double *mat);
Array<float> multMatCol(int n, Array<float>& A, Array<float>& x);
Array<float> matMult(int n, Array<float>& A, Array<float>& B);
