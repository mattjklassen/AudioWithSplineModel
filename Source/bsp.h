/*
  =============================================================================

    bsp.
    Created: 13 Nov 2021 1:39:23p
    Author:  Matt Klasse

  =============================================================================
*/

#pragma once

#include <fstream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "CycleSpline.h"

using namespace std;

void fillinputs(int d, int k, float *inputs);
void fillknots(int d, int k, float *knots);
void deboorfit(int k, int d, float u0, float SPP, float *bcoeffs, float *knots, short *samples3, float *DeBoor);
void truncfit(int k, int d, float u0, float SPP, float *trunccoeffs, float *knots, short *samples2 );
void getbcoeffs(int k, int d, float *bcoeffs, float *trunccoeffs, float *bmatinv);
void getwaveoutputs(int k, int d, float u0, float *inputs, float *outputs, float SPP, short *samples);
void gettrunccoeffs(int k, int d, float *trunccoeffs, float *matinv, float *outputs);
void bsplinefit(int k, int d, unsigned n, int NumZeros, short *samples, short *samples2, short *samples3, float *inputs, float *outputs, float *knots, float *SPP, float *trunccoeffs, float *matinv, float *bcoeffsteps, float *bcoeffs, float *bmatinv, float *DeBoor, short *BSPsamples);
unsigned GetWavDataSampleRate(fstream& in);
unsigned GetWavData(fstream& in, char *data);
unsigned GetWavDataSize(fstream& in);
void GetWavHeader(fstream& in, char *header);
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
float zerocross(int I, short *samples);//
float interp(float input, short *samples);//
int FindAllZeros(int Periods, float Freq, char *data, float *allzeros);//
int FindZerosClosestToPeriods(int Periods, float Freq, float *zeros, float *allzeros, float *SPP, float LastZero);
void truncfit2(float u0, float SPP, short *samples2, short *samples3);
void initBSPsamples(int k, int d, short *BSPsamples);
void writebcoeffs(int j, int k, int d, float *bcoeffs, short *BSPsamples);
void getbcoeffsteps(int k, int d, int NumZeros, float *bcoeffsteps, float *bcoeffs);
void getlinearbcoeffs(int k, int d, float *bcoeffsteps, float *bcoeffs);
void computeCycleBcoeffs(CycleSpline* cycle, short *samples);
void computeCycleSplineOutputs(CycleSpline* cycle, float* splineOutputs);

