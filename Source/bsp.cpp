/*
  ==============================================================================

    bsp.cpp
    Created: 13 Nov 2021 1:35:49pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "bsp.h"
#include "zeros.h"
#include "mat.h"
#include "CycleSpline.h"
#include "MetaSpline.h"

// bsp.cpp
// 
// Interpolate a wav file with splines
// using varying period lengths
//
// Matt Klassen 3/15/15
//
// compile with g++ bsp.cpp mat.cpp wavio.cpp zeros.cpp -o bsp
// usage:
//   ./bsp <d> <k> <F> <P> <file>
// where:
//   <d> --  degree for polynomial spline interpolation
//   <k> --  number of subintervals for splines
//   <F> --  guess at Frequency
//   <P> --  number of periods to compute
//   <file> -- input file (WAVE, 16bit mono, 44100 kHz)
// output:
//   tspline.wav - wave data interpolated with splines,
//                 using truncated power functon basis.
//   DeBoor.wav - wave data interpolated with splines,
//                 using B-spline basis and DeBoor algorithm.

#include <fstream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

struct fileheader
{
	char riff_label[4]; // (00) = {'R','I','F','F'}
	unsigned riff_size; // (04) = 36 + data_size
	char file_tag[4]; // (08) = {'W','A','V','E'}
	char fmt_label[4]; // (12) = {'f','m','t',' '}
	unsigned fmt_size; // (16) = 16
	unsigned short audio_format; // (20) = 1
	unsigned short channel_count; // (22) = 1 or 2
	unsigned sampling_rate; // (24) = (anything)
	unsigned bytes_per_second; // (28) = (see above)
	unsigned short bytes_per_sample; // (32) = (see above)
	unsigned short bits_per_sample; // (34) = 8 or 16
	char data_label[4]; // (36) = {'d','a','t','a'}
	unsigned data_size; // (40) = # bytes of data
};

// this function not used in JUCE project, original for command line version
int old_main(int argc, char *argv[]) 
{
  // process command line arguments
//
//  if (argc != 6)
//    return 0;
//  fstream in(argv[5],ios_base::binary|ios_base::in);
//  if (!in)
//    return 0;
//
//  int d = atoi(argv[1]); // degree of splines
//  int k = atoi(argv[2]); // number of subintervals
//  int Freq = atoi(argv[3]);  // guess at frequency
//  int Periods = atoi(argv[4]);  // number of periods on which to do spline fit
//    int sampleRate = 44100;
//    cout << endl << "number of periods is: " << Periods << endl;
//  //  cout << "size of int is: " << sizeof(int) << endl;
//  //  cout << "size of short is: " << sizeof(short) << endl;
//  //  cout << "size of unsigned is: " << sizeof(unsigned) << endl;
//  //  cout << "size of float is: " << sizeof(float) << endl;
//  //  cout << "size of float is: " << sizeof(float) << endl;
//
//  // Globals
//  // int i=0;
//  int N=k+2*d;  // N is last index of knot sequence
//  int NumAllZeros=0, NumZeros=0;
//  double LastZero=0.0f;
//
//  // set up arrays
//  double *allzeros = new double[10*Periods];
//  // this allows for about 10 zeros per period on average
//  double *zeros = new double[2*Periods];
//  double *SPP = new double[2*Periods];
//  // input values for interpolation
//  double *inputs = new double[k+d];
//  // output values for interpolation
//  double *outputs = new double[k+d];
//  // knot values for spline functions
//  double *knots = new double[N+1];
//  // arrays for the B-spline coefficients, k per period
//  // as char
//  char *BSPcoeffs = new char[(k+d)*2*441];
//  // and also same array as shorts, there are (k+d)*P BSPsamples
//  // BSPsamples could be used to write B-spline coeffcients to a wav file,
//  // as a type of compression of audio data
//  short *BSPsamples = reinterpret_cast<short*>(BSPcoeffs);
//
//  // cout << "degree d is: " << d << endl;
//  // cout << "number of subintervals k is: " << k << endl;
//  // cout << "dimension of splines is: " << k+d << endl;
//
//  // Get wave file data and set up arrays for data
//  char header[44];
//  in.read(header,44);
//
//  fileheader* h = 0;   // buffer for header
//  h = (fileheader*)header;
//
//  unsigned size = *reinterpret_cast<unsigned*>(header+40),
//  // size is number of bytes of data
//           rate = *reinterpret_cast<unsigned*>(header+24),
//           count = size/2;
//  // int n = (int)count;
//  // count is number of shorts of data (samples)
//  // unsigned size = GetWavDataSize(in);
//  // unsigned rate = GetWavDataSampleRate(in);
//  // GetWavData(in, data);
//  cout << "sample rate from wav file is:  " << rate << endl;
//
//  // input data
//  char *data  = new char[size];
//  in.read(data,size);
//  // output data
//  char *data2 = new char[size];
//  char *data3 = new char[size];
//
//  // short is 2 bytes, 16 bit data
//  short *samples = reinterpret_cast<short*>(data);
//  short *samples2 = reinterpret_cast<short*>(data2);  // to build truncated power output
//  short *samples3 = reinterpret_cast<short*>(data3);  // to build bspline output
//
//  cout << "The data size n is:  " << count << " samples " << endl;
//
//  // coefficients for spline function
//  double *trunccoeffs = new double[k+d];
//  double *bcoeffs = new double[k+d];
//  double *DeBoor = new double[(k+d)*(k+d)];
//
//  // coefficients for spline function
//  // double *bspcoeffs = new double[k+d];
//
//  // matrix coefficients for linear system to solve for
//  // shifted power function
//  double *mat = new double[(k+d)*(k+d)];
//  // matrix coefficients for inverse
//  double *matinv = new double[(k+d)*(k+d)];
//
//  // matrix for change of basis from B-spline basis
//  // to truncated power function basis
//  double *bmat = new double[(k+d)*(k+d)];
//  // inverse matrix for change of basis from truncated
//  // power function basis to B-spline basis
//  double *bmatinv = new double[(k+d)*(k+d)];
//  double *temp = new double[(k+d)*(k+d)];
//
//  ////////////////////////////////////////////////////////////////////
//  // the rest of main is now calling functions to do the interpolation
//
//  // initialize BSPsamples to 0
//  initBSPsamples(k, d, BSPsamples);
//
//  // set up interpolation input values in [0,1]
//  fillinputs(d, k, inputs);
//
//  // set up knot values extending [0,1]
//  fillknots(d, k, knots);
//
//  // fill array allzeros[] with all zeros from data
////  NumAllZeros = FindAllZeros(sampleRate, Periods, Freq, data, allzeros);
//    NumAllZeros = 0;  //  this will not be called !!
//  // for (i=NumAllZeros-20; i<NumAllZeros+1; i++)
//  // {
//  //    cout << allzeros[i] << endl;
//  // }
//  LastZero = allzeros[NumAllZeros];
//  // cout << "last zero:  " <<  LastZero << endl;
//
//  // fill array zeros[] with zeros closest to periods
////  NumZeros = FindZerosClosestToPeriods(sampleRate, Periods, Freq, zeros, allzeros, SPP, LastZero);
//    NumZeros = 0;  //   THIS IS NOT BEING CALLED !!
//
//  int NumSteps = NumZeros-1-40;
//  double *bcoeffsteps = new double[NumSteps];
//
//  // Truncated power function basis:
//  // (t-knots[0])^d_+,...,(t-knots[N-d-1])^d_+,
//  // Interpolating function:
//  // f(t) = a0*(t-knots[0])^d_+ +...+ a(N-d-1)*(t-knots[N-d-1])^d_+
//
//  // linear system to solve for function f(t)
//  // one row for each input value inputs[0],...,inputs[k+d-1]
//
//  // matrix is k+d by k+d, say with entries aij
//  // set matrix entries to solve for trunccoeffs
//
//  fillmat(k, d, mat, inputs, knots);
//
//  // matrix bmat describes B-splines as column vectors or coordinate
//  // vectors with respect to the truncated power function basis.
//  // Each coordinate vector comes from a quotient of Vandermonde
//  // determinants which interprets the divided difference formulation
//  // of B-splines through Cramer's Rule.
//
//  // initialize bmat with zeros:
//  initbmat(k, d, bmat);
//
//  fillbmat(k, d, bmat, knots);
//
//  // set entries for inverse matrix of mat and bmat (as identity)
//
//  setidentity(k+d, matinv);
//  setidentity(k+d, bmatinv);
//
//  // To compute spline coefficients, need inverse matrix.
//  // Solve matrix system Ax=b with A=mat[] and b=outputs[].
//  // Since we will use this same matrix equation with different
//  // values for b on each cycle, we should find the inverse
//  // of A and solve for x on each cycle as x=A^{-1}b.
//  // After computing the truncated power function coefficients
//  // we can use the polar forms, or use the change of basis
//  // matrix coming from Cramer's rule applied to the B-splines,
//  // to compute the B-spline coefficients.
//
//  // example: knot sequence 0,1,2,3,4,5,6,7,8,  degree d=2,
//  // shifted power basis: (t-0)^2_+,...,(t-5)^2_+, dim = 6.
//  // basic interval [2,6], inputs: 2,2.5,3,4,5,6.
//  // Note: Schoenberg-Whitney Theorem guarantees a unique solution.
//  // matrix:
//  // inputs, coeffs
//  // 2:      (2-0)^2   (2-1)^2     0        0       0      0
//  // 2.5:    (2.5-0)^2 (2.5-1)^2 (2.5-2)^2  0       0      0
//  // 3:      (3-0)^2   (3-1)^2   (3-2)^2    0       0      0
//  // 4:      (4-0)^2   (4-1)^2   (4-2)^2  (4-3)^2   0      0
//  // 5:      (5-0)^2   (5-1)^2   (5-2)^2  (5-3)^2 (5-4)^2  0
//  // 6:      (6-0)^2   (6-1)^2   (6-2)^2  (6-3)^2 (6-4)^2 (6-5)^2
//  //
//  // This matrix would only require two row operations to reach
//  // lower triangular form.  This would be more efficient than
//  // finding the inverse.  It could also be used with back
//  // substitution to solve for the coefficients in each cycle.
//  // If we do find the inverse, then it could be computed most
//  // efficiently by working from bottom right up to top left.
//  // Or maybe it just doesn't matter much.
//
//  // basic row operations:  swap, multiply, getzero
//  // swap swaps two rows
//  // multiply multiplies one row by a nonzero constant
//  // getzero gets a zero in a row and column using a pivot
//
//  // find matinv[] = inverse of mat[] using gausselim
//
//  matcopy(d+k, mat, temp);
//  gausselim(d+k, temp, matinv);
//
//  // find bmatinv[] = inverse of bmat[] using gausselim
//
//  matcopy(d+k, bmat, temp);
//  gausselim(d+k, temp, bmatinv);
//
//  // printmat(d+k, mat);
//  // printmat(d+k, matinv);
//  // printmat(d+k, bmat);
//  // printmat(d+k, bmatinv);
//
//  // Test matrix inverse with multiplication mat*matinv
//  // cout << "Product matrix mat*matinv:" << endl << endl;
//  // printmult(k+d, mat, matinv);
//
//  // Test matrix inverse with multiplication bmat*bmatinv
//  // cout << "Product matrix bmat*bmatinv:" << endl << endl;
//  // printmult(k+d, bmat, bmatinv);
//
//  // Now loop through all Periods and do bspline fit
//
//  bsplinefit(k, d, count, NumZeros, samples, samples2, samples3, inputs, outputs, knots, SPP, trunccoeffs, matinv, bcoeffsteps, bcoeffs, bmatinv, DeBoor, BSPsamples);
//
//  // Using A450.wav, a guitar note played at approx 450 Hz
//  // but Audacity plot spectrum says it's about 452 Hz
//  // which is about 46 cents over 440.  The number of samples
//  // per period at 452 Hz would be 97.56, so we are rounding up
//  // to 98 samples per period, which has the nice division
//  // property:  44100/98 = 450.  So we'll just call it 450 Hz.
//  // Note: 44100 = 2^2*3^2*5^2*7^2.
//
//  // write to wave files
//
//  // for (i=0; i<44; i++)
//  // {
//  //     cout << i << "  " << header[i] << endl;
//  // }
//
//  fstream out("tspline.wav",ios_base::binary|ios_base::out);
//  out.write(header,44);
//  out.write(data2,size);
//
//  fstream out1("DeBoor.wav",ios_base::binary|ios_base::out);
//  out1.write(header,44);
//  out1.write(data3,size);
//
//
//  size = (unsigned) (k+d)*2*441;
//  h->data_size = size;
//  fstream out2("bcoeffs.wav",ios_base::binary|ios_base::out);
//  out2.write(header,44);
//  out2.write(BSPcoeffs,(k+d)*2*441);
//
//  // free
//
//  delete[] data;
//  delete[] data2;
//  delete[] data3;
//  delete[] mat;
//  delete[] matinv;
//  delete[] bmat;
//  delete[] bmatinv;
//  delete[] inputs;
//  delete[] outputs;
//  delete[] knots;
//  delete[] bcoeffs;
//  delete[] trunccoeffs;
//  delete[] DeBoor;
//  delete[] allzeros;
//  delete[] zeros;
//  delete[] SPP;
//  delete[] temp;
//  delete[] BSPcoeffs;
//  // delete[] header;
//
  return 0;

}  // end main

// new functions for JUCE version

// compute B-spline value at t in [0,1] for basis with k subintervals on [0,1], d=3

float bSplineVal(int k, int j, float t)
{
    int d = 3;
    int n = k + d;  // dimension of splines
    int N = n + d;  // last index of knot sequence
    juce::Array<float> controlCoeffs;
    for (int i=0; i<n*n; i++) {
        controlCoeffs.add(0);
    }
    controlCoeffs.set(n*j, 1);
    float incr = 1 / (float) k;
    juce::Array<float> knotVals;
    for (int i=0; i<N+1; i++) {
        knotVals.add((i-d) * incr);
    }
    int i, J;  // J index in DeBoor Algorithm
    int p;  // stage in DeBoor Algorithm
    float denom, fac1, fac2, fval;

    fval = 0.0f;
    for (i=1; i<N; i++)
    {
        if (t < knotVals[i]) {
            J = i-1;
            if (J > n-1) {
                J = n-1;
            }
//            DBG("DeBoor J: " << J);
            break;
        }
    }
    for (p=1; p<d+1; p++)
    {
        for (i=J-d+p; i<J+1; i++)
        {
           denom = (knotVals[i+d-(p-1)]-knotVals[i]);
           fac1 = (t-knotVals[i]) / denom;
           fac2 = (knotVals[i+d-(p-1)]-t) / denom;
           controlCoeffs.set((k+d)*i+p, fac1 * controlCoeffs[(k+d)*i+(p-1)]
              + fac2 * controlCoeffs[(k+d)*(i-1)+(p-1)]);
        }
    }
    fval = controlCoeffs[(k+d)*J+d];
    return fval;
}

// use piecewise linear graph between targets to compute outputs
// no bcoeffs are necessary here, just targets and outputs.size
// targets for metaspline i will be bcoeff i from each key cycle
// outputs for metaspline i will be bcoeff i for all cycles
void computeLinearMetaSplineOutputs(MetaSpline& spline)
{
    int P = spline.outputs.size();
    float incr = 1 / (float)P;
    spline.outputs.set(0, spline.targets[0]);
    float t = 0, y = 0;
    int j = 0;
    for (int i=1; i<spline.outputs.size(); i++) {
        t = incr * i;
        j = 0;
        while(t > spline.inputs[j]) {
            j += 1;
        } // linear interp between targets[j-1] and targets[j]:
        float t0 = spline.inputs[j-1], t1 = spline.inputs[j];
        float y0 = spline.targets[j-1], y1 = spline.targets[j];
        y = (t-t0)/(t1-t0) * y1 + (t1-t)/(t1-t0) * y0;
        spline.outputs.set(i, y);
    }
}

// compute the outputs of a meta-spline using DeBoor algorithm
// assuming meta-spline has bcoeffs and targets (which are key cycle bcoeffs)
void computeMetaSplineOutputs(MetaSpline& spline)
{
    int d = spline.d, n = spline.n,  N = n + d;  // k = n - d;
    juce::Array<float> controlCoeffs;
    for (int i=0; i<4*n; i++) {
        controlCoeffs.add(0);
    }
    int i, h, J;  // J index in DeBoor Algorithm
    int p;  // stage in DeBoor Algorithm
    float denom, fac1, fac2, t, fval;
    // numOutputs = number of outputs from meta-spline = numCycles

    spline.outputs.clear();
    int H = spline.numOutputs;
    float incr = 1/((float)(H-1));

    for (h=0; h<H; h++)  // loop on H uniformly distributed values in [0,1]
    {
        t = h * incr;  // t value in interval [0,1]
        fval = 0.0f;
        // set J value
        for (i=1; i<N; i++)
        {
            if (t < spline.knots[i])
            {
              J = i-1;
                if (J > n-1) {
                    J = n-1;
                }
              break;
            }
        }
        for (i=0; i<n; i++)
        {
            controlCoeffs.set(4*i+0, spline.bcoeffs[i]);
        }
        for (p=1; p<d+1; p++)
        {
            for (i=J-d+p; i<J+1; i++)
            {
              denom = (spline.knots[i+d-(p-1)]-spline.knots[i]);
              fac1 = (t-spline.knots[i]) / denom;
              fac2 = (spline.knots[i+d-(p-1)]-t) / denom;
              controlCoeffs.set(4*i+p, fac1 * controlCoeffs[4*i+(p-1)]
                  + fac2 * controlCoeffs[4*(i-1)+(p-1)]);
            }
        }
        fval = controlCoeffs[4*J+d];
        spline.outputs.set(h, fval);
    }
}

// compute the outputs of a spline on one cycle using DeBoor algorithm
void computeCycleSplineOutputs(CycleSpline& cycle)
{
    int k = cycle.k, d = cycle.d;
    int n = k + d;
    int N = n + d;
    juce::Array<float> controlCoeffs;
    for (int i=0; i<4*(k+d); i++) {
        controlCoeffs.add(0);
    }
    float incr = 1 / (float) k;
    juce::Array<float> knotVals;
    for (int i=0; i<N+1; i++) {
        knotVals.add((i-d) * incr);
    }
    float a = cycle.a, b = cycle.b;
    int i, h, J;  // J index in DeBoor Algorithm
    int p;  // stage in DeBoor Algorithm
    float length = b - a;
    float max = 0;
    int maxI = 0;
    int numSamples = (int) b - (int) a;
    float denom, fac1, fac2, t, fval;
    int Imin = (int)ceil(a);       // first sample index
    int Imax = Imin + numSamples - 1;  // last sample index
    int M = (int)(Imax - Imin) + 1;   // number of samples
    if (M<1) return;
    // we need outputs from spline to match number of samples in interval
    // initialize to zeros, then reset with values computed from bcoeffs
    cycle.outputs.clear();
    for (int i=0; i<Imax-Imin+1; i++) {
        cycle.outputs.add(0);
    }
    for (h=Imin; h<Imax+1; h++)  // loop on samples
    {
        t = ((float)h-a)/length;  // t value in interval [0,1]
        fval = 0.0f;
        // set K value
        for (i=1; i<N; i++)
        {
            if (t < knotVals[i])
            {
              J = i-1;
                if (J > n-1) {
                    J = n-1;
                }
              break;
            }
        }
        for (i=0; i<d+k; i++)
        {
            controlCoeffs.set(4*i+0, cycle.bcoeffs[i]);
        }
        for (p=1; p<d+1; p++)
        {
            for (i=J-d+p; i<J+1; i++)
            {
              denom = (knotVals[i+d-(p-1)]-knotVals[i]);
              fac1 = (t-knotVals[i]) / denom;
              fac2 = (knotVals[i+d-(p-1)]-t) / denom;
              controlCoeffs.set(4*i+p, fac1 * controlCoeffs[4*i+(p-1)]
                  + fac2 * controlCoeffs[4*(i-1)+(p-1)]);
            }
        }
        fval = controlCoeffs[4*J+d];
        if (abs(fval) > max) {
            max = abs(fval);
            maxI = h-Imin;
        }
        cycle.outputs.set(h-Imin, fval);  // outputs for int h = 0 to Imax-Imin
    }
    cycle.maxVal = max;
    cycle.maxIndex = maxI;
}

// computes the bcoeffs of meta-spline which will pass through the keycoeffs or targets
// associated to one index i. Assumes all data in MetaSpline is set, including targets,
// except for the bcoeffs (and of course the outputs)

void computeMetaSplineBcoeffs(MetaSpline& spline)
{
    int n = spline.n, d = spline.d, k = n - d;

    // Set up n by n linear system to solve for the meta-spline B-spline coefficients.
    // As of 1/13/22 we are changing meta-splines to be natural cubic splines.
    // Once we solve for the matrix inverse, then just multiply matrix * column vector
    // (vector of target values which are keycoeffs of key cycles, and also 2 zeros)
    // The 2 zeros come from the natural cubic spline conditions, which for B-spline coeffs
    // are given as: c(0)-2*c(1)+c(2)=0 and c(n-3)-2*c(n-2)+c(n-1)=0
    // Each row of matrix corredsponds to plugging in one input, say s_i, and setting
    // equal to target value, say y_i, or a second derivative at the ends = 0.
    // The function to solve for is a sum of Bsplines.
    
    float val = 0;
    juce::Array<float> A;
    juce::Array<float> B;
    juce::Array<float> temp;
    for (int i=0; i<n*n; i++) {
        A.add(0);
        B.add(0);
        temp.add(0);
    }
    // spline collocation matrix, up to last two rows
    for (int i=0; i<n-2; i++) {
        for (int j=0; j<n; j++) {
            val = bSplineVal(k, j, spline.inputs[i]);  // A[i,j]
            A.set(i*n+j, val);
            temp.set(i*n+j, val);
//            B.set(i*n+j, 0);
        }
    }
    // set last two rows according to second derivative = 0 at ends
    int i = n-2;
    int j = 0;
    A.set(i*n+j, 1);
    j = 1;
    A.set(i*n+j, -2);
    j = 2;
    A.set(i*n+j, 1);
    i = n-1;
    j = n-3;
    A.set(i*n+j, 1);
    j = n-2;
    A.set(i*n+j, -2);
    j = n-1;
    A.set(i*n+j, 1);
    // set matrix B as identity
    for (int i=0; i<n; i++) {
        B.set(i*n+i, 1.0);
    }
//    printMatrix(n, A);
//    printMatrix(n, B);
//    Now we need inverse of A, n x n, with n = k+d
    gaussElim(n, temp, B);
//    printMult(n, A, B);
    Array<float> x;
    Array<float> y;
    for (int i=0; i<n; i++) {
        x.add(0);
        y.add(0);
    }
    y = multMatCol(n, B, spline.targets);
    for (int i=0; i<n; i++) {
        spline.bcoeffs.set(i, y[i]);
    }
}

void computeCycleBcoeffs(CycleSpline& cycle, AudioBuffer<float>& samples)
{
    int k = cycle.k, d = cycle.d;
    int n = k + d;
    int N = n + d;  // N is last index of knot sequence, N - d = n
    juce::Array<float> inputVals;
    for (int i=0; i<n; i++) {
        inputVals.add(cycle.inputs[i]);
    }
    // outputVals are the target wav data values to match with splines
    juce::Array<float> outputVals;
    for (int i=0; i<n; i++) {
        outputVals.add(0);
        cycle.targets.add(0);
    }
    for (int i=1; i<n-1; i++) {
        float output = interpFloat(inputVals[i], samples);
        outputVals.set(i, output);
        cycle.targets.set(i, output);
    }
//    DBG("set outputVals");
    // now reset inputs to be on the interval [0,1] in the usual way
    // (ie. on the interval breakpoints 0,1/k,2/k,...,(k-1)/k,k and then
    // also the two other points 0.5/k and 1 - 0.5/k )
    float incr = 1 / (float) k;
    inputVals.set(0, 0.0);
    inputVals.set(1, 0.5 * incr);
    inputVals.set(2, incr);
    inputVals.set(n-2, 1 - 0.5 * incr);
    inputVals.set(n-1, 1.0);
    // Need n inputs, i=0...n-1, have 5, so another n-5 at indices 3 to n-3
    for (int i=3; i<n-2; i++) {
        inputVals.set(i, inputVals[i-1] + incr);  // +incr,...,(n-5)*incr
    }
//    DBG("inputVals: ");
//    for (int i=0; i<n; i++) {
//        DBG("inputVals[" << i << "] = " << inputVals[i]);
//    }
    juce::Array<float> knotVals;
    for (int i=0; i<N+1; i++) {
        knotVals.add((i-d) * incr);
    }
//    DBG("knotVals: ");
//    for (int i=0; i<N+1; i++) {
//        DBG("knotVals[" << i << "] = " << knotVals[i]);
//    }
    // Next, need to set up (k+d) by (k+d) linear system to solve for the B-spline coefficients.
    // Once we solve for the matrix inverse, then just multiply matrix * column vector
    // (vector of target values called outputs) to get coefficients.
    // Each row of matrix corredsponds to plugging in one input, say s_i, and setting
    // equal to target value, say y_i.  The function to solve for is a sum of Bsplines.
    float val = 0;
    juce::Array<float> A;
    juce::Array<float> B;
    juce::Array<float> temp;
    for (int i=0; i<n*n; i++) {
        A.add(0);
        B.add(0);
        temp.add(0);
    }
//    DBG("set matrices A and B");
    for (int i=0; i<n; i++) {
//        DBG("computing row i: " << i);
//        DBG("inputVals[" << i << "]: " << inputVals[i]);
        for (int j=0; j<n; j++) {
            val = bSplineVal(k, j, inputVals[i]);  // A[i,j]
//            DBG("computing val with k,j: " << k << ", " << j << " val: " << val);
            A.set(i*n+j, val);
            temp.set(i*n+j, val);
            B.set(i*n+j, 0);
        }
    }
//    DBG("reset matrices A and B");
    for (int i=0; i<n; i++) {
        B.set(i*n+i, 1.0);
    }
//    DBG("B set to identity");
//    DBG("matrix A:");
//    printMatrix(n, A);
//    DBG("matrix B:");
//    printMatrix(n, B);

    // Now we need inverse of A, n x n, with n = k+d
    gaussElim(n, temp, B);
//    printMult(n, A, B);
    Array<float> x;
    Array<float> y;
    for (int i=0; i<n; i++) {
        x.add(0);
        y.add(0);
    }
    y = multMatCol(n, B, outputVals);
    for (int i=0; i<n; i++) {
        cycle.bcoeffs.set(i, y[i]);
    }
}

void computeCycleBcoeffsOld(CycleSpline& cycle, AudioBuffer<float>& samples)
{
//    from input cycle we get: k, d, and other data:
//    int d;  // degree of B-splines, default d = 3
//    int k; // number of subintervals, can vary, use default 20
//    float a, b; // for interval [a,b] on time axis, a = left end point, b = right end point
//    Array<float> subintervals;  // default subinterval size = (b-a)/k, endpoints: u_0,...,u_k (k+1 values)
//    Array<float> inputs;  // default at a, b, and subinterval breakpoints, and midpoints of outer two subintervals
//    Array<float> knots;   // default uniform knot sequence includes subinterval endpoints, and d more on each end
//    Array<float> bcoeffs; // B-spline coefficients, d+k of them
    
    int k = cycle.k, d = cycle.d;
    juce::Array<float> inputVals;
    for (int i=0; i<k+d; i++) {
        inputVals.add(0);
    }
    double *inputs = new double[k+d];
    for (int i=0; i<k+d; i++) {
        inputs[i] = cycle.inputs[i];
        inputVals.set(i, cycle.inputs[i]);
    }
    // outputs are the target wav data to match with splines
    juce::Array<float> outputVals;
    for (int i=0; i<k+d; i++) {
        outputVals.add(0);
    }
    double *outputs = new double[k+d];
    outputs[0] = 0;
    outputs[k+d-1] = 0;
    for (int i=1; i<k+d-1; i++) {
        outputs[i] = interpFloat(inputs[i], samples);
        outputVals.set(i, interpFloat(inputVals[i], samples));
//        DBG("outputs[" << i << "] =  " << outputs[i]);
    }
    float incr = 1 / (float) k;
    juce::Array<float> knotVals;
    for (int i=0; i<k+2*d; i++) {
        knotVals.add((i-d) * incr);
    }
    double *knots = new double[k+2*d];
//    for (int i=0; i<k+d; i++) {
//        knots[i] = cycle.knots[i];
//    }
    fillinputs(d, k, inputs);
    fillknots(d, k, knots);
    double *trunccoeffs = new double[k+d];
    double *bcoeffs = new double[k+d];
    double *mat = new double[(k+d)*(k+d)];
    double *matinv = new double[(k+d)*(k+d)];
    double *bmat = new double[(k+d)*(k+d)];
    double *bmatinv = new double[(k+d)*(k+d)];
    double *temp = new double[(k+d)*(k+d)];
    fillmat(k, d, mat, inputs, knots);
    initbmat(k, d, bmat);
    fillbmat(k, d, bmat, knots);
    setidentity(k+d, matinv);
    setidentity(k+d, bmatinv);
    matcopy(d+k, mat, temp);
    gausselim(d+k, temp, matinv);
    matcopy(d+k, bmat, temp);
    gausselim(d+k, temp, bmatinv);
    gettrunccoeffs(k, d, trunccoeffs, matinv, outputs);
    getbcoeffs(k, d, bcoeffs, trunccoeffs, bmatinv);
    for (int i=0; i<k+d; i++) {
        float bcoeff = (float)bcoeffs[i];
        cycle.bcoeffs.set(i, bcoeff);
//        DBG("bcoeffs[" << i << "] =  " << bcoeffs[i]);
//        DBG("cycle.bcoeffs[" << i << "] =  " << cycle.bcoeffs[i]);
    }
    delete[] inputs;
    delete[] outputs;
    delete[] knots;
    delete[] bcoeffs;
    delete[] trunccoeffs;
    delete[] mat;
    delete[] matinv;
    delete[] bmat;
    delete[] bmatinv;
    delete[] temp;
}

void bsplinefit(int k, int d, unsigned n, int NumZeros, short *samples, short *samples2, short *samples3, double *inputs, double *outputs, double *knots, double *SPP, double *trunccoeffs, double *matinv, double *bcoeffsteps, double *bcoeffs, double *bmatinv, double *DeBoor, short *BSPsamples)
{
  // uses these function calls:
  // deboorfit(k, d, u0, SPP, bcoeffs, knots, samples3, DeBoor);
  // truncfit(k, d, u0, SPP[J], trunccoeffs, knots, samples2 );
  // getbcoeffs(k, d, bcoeffs, trunccoeffs, bmatinv);
  // getwaveoutputs(k, d, u0, inputs, outputs, SPP[J], samples);
  // gettrunccoeffs(k, d, trunccoeffs, matinv, outputs);
//  int J = 0;  // for loop on Periods
//  // period lengths = SPP[J]
//  // u0 = left endpoint, uk = right endpoint
//  // interval for interpolation is just SPP[J]*[0,1]
//  // so input values are just SPP[J]*inputs[i]
//  double u0 = 0.0f;
//
//  // loop on Periods found = NumZeros, starting with first Period J=1
//  // J^th sub-interval is [zeros[J-1],zeros[J]].
//  for (J=1; J<NumZeros; J++)
//  // for (J=1; J<41; J++)
//  {
//      // if (J > 12) { break; }
//
//// capture the wave file data for current cycle, using linear
//// interpolation bewteen samples, with left endpoint u0
//// outputs[] are normalized to be in [-1,1]
//
//   getwaveoutputs(k, d, u0, inputs, outputs, SPP[J], samples);
//
//// cout << endl << "outputs:" << endl << endl;
//// cout << i << ":  " << outputs[i] << endl;
//
//// coefficients for shifted power functions
//// given by matrix inverse times output values
//
//   gettrunccoeffs(k, d, trunccoeffs, matinv, outputs);
//
//   // coefficients for B-spline functions
//   // given by bmatinv matrix times trunccoeffs
//
//   getbcoeffs(k, d, bcoeffs, trunccoeffs, bmatinv);
//
//   // write the bcoeffs to BSPsamples array for output
//
//   writebcoeffs(J-1, k, d, bcoeffs, BSPsamples);
//
//// construct output wave data with truncated power spline coefficients
//// on interval from u0 to u0+SPP[J] and put in samples2[i]
//// for index values i in the above interval, from ceil(u0) to floor(u0+SPP[J])
//
//   truncfit(k, d, u0, SPP[J], trunccoeffs, knots, samples2 );
//
//// re-construct wave data with B-spline coefficients
//
//   deboorfit(k, d, u0, SPP[J], bcoeffs, knots, samples3, DeBoor);
//
//   // for checking output matches at zeros:
//   // if (J % 2 == 0)  // J even
//   // {
//   //    truncfit2(u0, SPP[J], samples2, samples3);
//   // }
//
//   // cout << J << "  u0 value is:  " << u0 << "  SPP[J]:  " << SPP[J] << endl;
//
//   u0 = u0 + SPP[J];

//  } // end first J loop on Periods

//  getbcoeffsteps(k, d, NumZeros, bcoeffsteps, bcoeffs);

//  for (J=41; J<NumZeros; J++)
//  {

//     getlinearbcoeffs(k, d, bcoeffsteps, bcoeffs);
//     writebcoeffs(J-1, k, d, bcoeffs, BSPsamples);
//     deboorfit(k, d, u0, SPP[J], bcoeffs, knots, samples3, DeBoor);

//     u0 = u0 + SPP[J];

//  } // end second J loop on Periods

  // the rest are zeros
//  int m = (int)floor(u0);
//  int r;
//  cout << "m value:  " << m << endl;
//
//  // for (r=0; r<n-1; r++)   //  write all as zero
//  for (r=m; r<(int)n; r++)      //  write samples past last interval as zero
//  {
//      // samples2[k] = samples2[k-98];
//      samples2[r] = 0;
//      samples3[r] = 0;
//  }
 
}    // end function bsplinefit

// write the bspline coefficients as shorts to array in col j
//void writebcoeffs(int j, int k, int d, float *bcoeffs, short *BSPsamples)
//{
//    int i=0;
//    for (i=0; i<k+d; i++)
//    {
//	BSPsamples[441*i + j] = (short) (32768.0f * bcoeffs[i]);
//    }
//}     // end function writebcoeffs
//
//void initBSPsamples(int k, int d, short *BSPsamples)
//{
//    int j=0;
//    for (j=0; j<(k+d)*441; j++)
//    {
//	BSPsamples[j] = 0;
//    }
//}


void fillinputs(int d, int k, double *inputs)
{

// This function puts k+d inputs values into the array
// inputs[].
// Set up the interval [0,1] broken into k subintervals 
// with the following scheme:
// Use the k subintervals broken at even steps of 1/k.
// Then we need to choose k+d points in these intervals
// for interpolation.  Call these values inputs[0],...,inputs[k+d-1].
// This array inputs[i] will include all of the interval end-points
// based on the regular steps of size 1/k, called "regular inputs",
// and also "extra insertions" at the half-way points in
// d-1 of these intervals, alternately at the left and
// right ends of the interval.  This yields for example:
// d=1: [0,1/k,2/k,...,(k-1)/k,1]  (no extra insertions)
// d=2: [0,1/2k,1/k,...,(k-1)/k,1] (1 extra insertion)
// d=3: [0,1/2k,1/k,...,(k-1)/k,1-1/2k,1]  (2 extra insertions)
// etc.
// This type of choice gives a unique interpolation function
// according to the Schoenberg-Whitney Theorem.
// (we could also do other interpolation schemes, including
//  derivatives, etc.)

  double delta = 1.0f / (2*(double)k);
  inputs[0]=0.0f; inputs[k+d-1]=1.0f;  // define endpoints

  // define values in pairs using half-steps delta
  // working alternately from the ends toward the middle

  int j = 1;

  while (j<d)  // there are d-1 extra inserts
  {
      if (j % 2 == 1)  // j odd, assign on the left
      {
	  inputs[j] = inputs[j-1]+delta;   // extra insert
	  inputs[j+1] = inputs[j]+delta;   // regular input
      }
      if (j % 2 == 0)  // j even, assign on the right
      {
	  inputs[k+d-1-(j-1)] = inputs[k+d-1-(j-2)]-delta;   // extra insert
	  inputs[k+d-1-j] = inputs[k+d-1-(j-1)]-delta;       // regular input
      }
      j++;
  }

  // now we have assigned inputs[0], inputs[d+k-1], and 2*(d-1) other
  // values, so there are k+d-2*d = k-d remaining.  These
  // are inputs[d] to inputs[k-1] for d odd, 
  // and inputs[d+1] to inputs[k] for d even.

  delta = delta*2;
  if (d % 2 == 1)  // d odd
  {
      for (j=d; j<k; j++)
      {
	  inputs[j] = inputs[j-1]+delta;
      }
  }
  if (d % 2 == 0)  // d even
  {
      for (j=d+1; j<k+1; j++)
      {
	  inputs[j] = inputs[j-1]+delta;
      }
  }

  // print array inputs[i]

  // for (j=0; j<d+k; j++)
  // {
  //    cout << "inputs[" << j << "] : " << inputs[j] << endl;
  // }

} // end function fillinputs


// compute uniform sequence of knot values
// t0=u0-d*delta, ti=t0+i*delta, i=1..N=k+2*d=k+d+d
// This will be used for B-spline basis
// Knot sequence for splines: t0,t1,...,tN where
// t0=u0-d*delta, ti=t0+i*delta,  td=u0, t(N-d)=uk, and N=k+2*d
// Truncated power function basis: (t-t0)^d_+,...,(t-t(N-d-1))^d_+,
// B-spline basis: B^d_0(t),...,B^d_(N-d-1)(t)
// Dimension: N-d = k+d.

void fillknots(int d, int k, double *knots)
{

  int i;
  double delta = 1.0f / (double)k;
  int N=k+2*d;  // N is last index of knot sequence
  
  knots[0]=-(double)d * delta;

  for (i=1; i< N+1; i++)
  {
      knots[i] = knots[i-1] + delta;
  //    cout << "knots[" << i << "] : " << knots[i] << endl;
  }

// print knots
// cout << endl << "knots:" << endl << endl;
// for (i=0; i<N+1; i++)
// {
//    cout << i << ":  " << knots[i] << endl;
// }

}  // end function fillknots()

//void getOutputs(int k, int d, float a, float *inputs, float *outputs, float numSamples, short *samples)
//{
//   int i;
//   float input;
//   for (i=0; i<d+k; i++)
//   {
//       input = a + inputs[i]*SPP;
//       outputs[i] = interp(input, samples);
//   }
//}

void getwaveoutputs(int k, int d, double u0, double *inputs, double *outputs, double SPP, short *samples)
{
//   int i;
//   double input;
//   for (i=0; i<d+k; i++)
//   {
//       input = u0 + inputs[i]*SPP;
//       outputs[i] = interp(input, samples);
//       //   cout << i << ":  " << outputs[i] << endl;
//   }
}

void gettrunccoeffs(int k, int d, double *trunccoeffs, double *matinv, double *outputs)
{
   int i,j;
   for (i=0; i<k+d; i++)
   {
       trunccoeffs[i] = 0;
       for (j=0; j<k+d; j++)
       {
           trunccoeffs[i] += matinv[(k+d)*i+j]*outputs[j];
       }
   }
}

// call this to use last calculated bcoeffs to produce step size
// for linear approx of further bcoeffs, to ramp down to zero
void getbcoeffsteps(int k, int d, int NumZeros, float *bcoeffsteps, float *bcoeffs)
{
    int i;
    int NumSteps = NumZeros-1-40;
    for (i=0; i<k+d; i++)
    {
	bcoeffsteps[i] = bcoeffs[i] / (float)NumSteps;
    }
}

void getlinearbcoeffs(int k, int d, float *bcoeffsteps, float *bcoeffs)
{
    int i;
    for (i=0; i<k+d; i++)
    {
	bcoeffs[i] = bcoeffs[i] - bcoeffsteps[i];
    }
}


void getbcoeffs(int k, int d, double *bcoeffs, double *trunccoeffs, double *bmatinv)
{
    int i,j;
    for (i=0; i<k+d; i++)
    {
        bcoeffs[i] = 0;
        for (j=0; j<k+d; j++)
        {
            bcoeffs[i] += bmatinv[(k+d)*i+j]*trunccoeffs[j];
        }
    }
}      

// write zeros to chosen periods instead of spline fit
// used for visualization and checking on ends of periods
void truncfit2(float u0, float SPP, short *samples2, short *samples3 )
{
  int h;
  int Imin = (int)ceil(u0);       // first sample index
  int Imax = (int)floor(u0+SPP);  // last sample index
  int M = (int)(Imax - Imin) + 1;   // number of samples
  if (M<1) return;
  for (h=Imin; h<Imax+1; h++)  // loop on samples
  {
      samples2[h] = 0;
      samples3[h] = 0;
  }
}   // end function truncfit2


void truncfit(int k, int d, float u0, float SPP, float *trunccoeffs, float *knots, short *samples2 )
{
  int i,h;
  float t, fval;
  int Imin = (int)ceil(u0);       // first sample index
  int Imax = (int)floor(u0+SPP);  // last sample index
  int M = (int)(Imax - Imin) + 1;   // number of samples
  if (M<1) return;
  for (h=Imin; h<Imax+1; h++)  // loop on samples
  {
      t = ((float)h-u0)/SPP;  // t value in interval [0,1]
      fval = 0.0f;
      for (i=0; i<d+1; i++)
      {
	  fval += trunccoeffs[i]*pow(t-knots[i],d);
      }
      for (i=d+1; i<d+k; i++)
      {
	  if (t > knots[i])
	  {
	     fval += trunccoeffs[i]*pow(t-knots[i],d);
	  }
      }
      samples2[h] = (short)(fval * 32768.0f);
  }
}   // end function truncfit

void deboorfit(int k, int d, float u0, float SPP, float *bcoeffs, float *knots, short *samples3, float *DeBoor)
{
  int i,K;  // use K for the usual J index in DeBoor Algorithm
  int N=k+2*d;
  int p;  // stage in DeBoor Algorithm
  // float t, fval;
  float denom, fac1, fac2;
  int h;
  float t, fval;
  int Imin = (int)ceil(u0);       // first sample index
  int Imax = (int)floor(u0+SPP);  // last sample index
  int M = (int)(Imax - Imin) + 1;   // number of samples
  if (M<1) return;
  for (h=Imin; h<Imax+1; h++)  // loop on samples
  {
      t = ((float)h-u0)/SPP;  // t value in interval [0,1]
      fval = 0.0f;
      // set K value
      for (i=1; i<N-d; i++)
      {
          if (t < knots[i])
          {
              K = i-1;
              break;
          }
      }
      for (i=0; i<d+k; i++)
      {
          DeBoor[(k+d)*i+0] = bcoeffs[i];
      }
      for (p=1; p<d+1; p++)
      {
          for (i=K-d+p; i<K+1; i++)
          {
              denom = (knots[i+d-(p-1)]-knots[i]);
              fac1 = (t-knots[i]) / denom;
              fac2 = (knots[i+d-(p-1)]-t) / denom;
              DeBoor[(k+d)*i+p] = fac1 * DeBoor[(k+d)*i+(p-1)]
              + fac2 * DeBoor[(k+d)*(i-1)+(p-1)];
          }
      }
      fval = DeBoor[(k+d)*K+d];
      samples3[h] = (short)(fval * 32768.0f);
  }
}   // end function deboorfit

void PrintHeader(fileheader* h, char* data, int i)
{
	 	  printf("\n");
	if(i) printf("riff_label:       %c%c%c%c.\n",h->riff_label[0],h->riff_label[1],h->riff_label[2],h->riff_label[3]);
		  printf("riff_size:        %u.\n",h->riff_size);
	if(i) printf("file_tag:         %c%c%c%c.\n",h->file_tag[0],h->file_tag[1],h->file_tag[2],h->file_tag[3]);
	if(i) printf("fmt_label:        %c%c%c%c.\n",h->fmt_label[0],h->fmt_label[1],h->fmt_label[2],h->fmt_label[3]);
	
		  printf("fmt_size:         %u.\n",h->fmt_size);
		  printf("audio_format:     %u.\n",h->audio_format);
		  printf("channel_count:    %u.\n",h->channel_count);
		  printf("sampling_rate:    %u.\n",h->sampling_rate);
		  
		  printf("bytes_per_second: %u.\n",h->bytes_per_second);
		  printf("bytes_per_sample: %u.\n",h->bytes_per_sample);
		  printf("bits_per_sample:  %u.\n",h->bits_per_sample);
	if(i) printf("data_label:       %c%c%c%c.\n",h->data_label[0],h->data_label[1],h->data_label[2],h->data_label[3]);
		  printf("data_size:        %u.\n",h->data_size);
	
	if(h->bits_per_sample == 16)
	{
		short* p = (short*)data;
		for(unsigned i=90;i<110;++i)
			printf("%d, ",p[i]);
		printf("\n");
	}
	else if(h->bits_per_sample == 8)
	{	
		char* p = (char*)data;
		for(unsigned i=0;i<4;++i)
			printf("%d, ",p[i]);
		printf("\n");
	}
}
