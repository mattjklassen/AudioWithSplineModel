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
using namespace std;

void fillbmat(int k, int d, float *bmat, float *knots);
void fillmat(int k, int d, float *mat, float *inputs, float *knots);
void gausselim(int n, float *mat, float *matinv);
int  getpivot(int n, int j, float *mat);
void getzero(int n, int i, int j, float *mat, float *matinv);
void initbmat(int k, int d, float *bmat);
void matcopy(int n, float *mat, float *temp);
void printmat(int n, float *mat);
void printmult(int n, float *mat, float *matinv);
void rowmult(int n, int i, float C, float *mat, float *matinv);
void rowswap(int n, int m, int r, float *mat, float *matinv);
void setidentity(int n, float *mat);
