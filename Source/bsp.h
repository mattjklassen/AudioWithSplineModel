/*
  =============================================================================

    bsp.h
    Created: 13 Nov 2021 1:39:23p
    Author:  Matt Klassen

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
#include "mat.h"

using namespace std;

void fillinputs(int d, int k, double *inputs);
void fillknots(int d, int k, double *knots);
void deboorfit(int k, int d, double u0, double SPP, double *bcoeffs, double *knots, short *samples3, double *DeBoor);
void truncfit(int k, int d, double u0, double SPP, double *trunccoeffs, double *knots, short *samples2 );
void getbcoeffs(int k, int d, double *bcoeffs, double *trunccoeffs, double *bmatinv);
void getwaveoutputs(int k, int d, double u0, double *inputs, double *outputs, double SPP, short *samples);
void gettrunccoeffs(int k, int d, double *trunccoeffs, double *matinv, double *outputs);
void bsplinefit(int k, int d, unsigned n, int NumZeros, short *samples, short *samples2, short *samples3, double *inputs, double *outputs, double *knots, double *SPP, double *trunccoeffs, double *matinv, double *bcoeffsteps, double *bcoeffs, double *bmatinv, double *DeBoor, short *BSPsamples);
unsigned GetWavDataSampleRate(fstream& in);
unsigned GetWavData(fstream& in, char *data);
unsigned GetWavDataSize(fstream& in);
void GetWavHeader(fstream& in, char *header);
void fillbmat(int k, int d, double *bmat, double *knots);
void fillmat(int k, int d, double *mat, double *inputs, double *knots);
void gausselim(int n, double *mat, double *matinv);
int  getpivot(int n, int j, double *mat);
void getzero(int n, int i, int j, double *mat, double *matinv);
void initbmat(int k, int d, double *bmat);
void matcopy(int n, double *mat, double *temp);
void printmat(int n, double *mat);
void printmult(int n, double *mat, double *matinv);
void rowmult(int n, int i, double C, double *mat, double *matinv);
void rowswap(int n, int m, int r, double *mat, double *matinv);
void setidentity(int n, double *mat);
//double zerocross(int I, short *samples);//
double interp(double input, short *samples);//
int FindAllZeros(int Periods, double Freq, char *data, double *allzeros);//
int FindZerosClosestToPeriods(int Periods, double Freq, double *zeros, double *allzeros, double *SPP, double LastZero);
void truncfit2(double u0, double SPP, short *samples2, short *samples3);
void initBSPsamples(int k, int d, short *BSPsamples);
void writebcoeffs(int j, int k, int d, double *bcoeffs, short *BSPsamples);
void getbcoeffsteps(int k, int d, int NumZeros, double *bcoeffsteps, double *bcoeffs);
void getlinearbcoeffs(int k, int d, double *bcoeffsteps, double *bcoeffs);
void computeCycleBcoeffs(CycleSpline& cycle, AudioBuffer<float>& samples);
void computeCycleBcoeffs2(CycleSpline& cycle, AudioBuffer<float>& samples);
void computeCycleSplineOutputs(CycleSpline& cycle);
float bSplineVal(int k, int j, float t);

