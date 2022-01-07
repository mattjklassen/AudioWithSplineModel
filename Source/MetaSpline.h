/*
  ==============================================================================

    MetaSpline.h
    Created: 27 Dec 2021 6:44:11pm
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once
#include "../JuceLibraryCode/JuceHeader.h"

class MetaSpline
{
public:
    // a meta-spline typically interpolates a particular bcoeff at all key cycles
    // refer to this bcoeff as bcoeffi, which is then an array of size numCycles
    // so it is used to generate missing values of that bcoeff at the non-key cycles
    int d;  // degree of (meta) B-splines, default d = 3
    // k=n-d  // number of subintervals, dim n=k+3 = # key cycles, so k=n-3
    int n;  // number of key cycles
    int m;  // indices of key cycles are j*m, j=0...n-1
    int numOutputs;  // should be equal to numCycles
    float a, b; // for interval [a,b], default to [0,1]
    Array<float> subintervals;  // default subinterval size = 1/k,
    // subinterval sequence:  [0,1/k,2/k,...,1] (length k+1)
    Array<float> inputs; // default: 0,1/(n-1),2/(n-1),...,1 (length n), determines target points
    Array<float> targets; // key cycle values of bcoeffs to match with meta-spline
    Array<float> outputs;  // output values computed from B-spline basis with DeBoor algorithm,
    // outputs to be computed for values between key cycle bcoeffs, ie. between indices 0, m, 2m, ...
    // total number of outputs should match total cycles, or (numKeyCycles-1)*mVal
    // one could also set the last values to equal to previous computed
    Array<float> knots; // default uniform knot sequence includes subintervals plus d more on each end
    Array<float> bcoeffs; // B-spline coefficients, n=d+k of them
    // The above data are enough to compute the bcoeffs

    MetaSpline (int _n, int _m);  // all other params are set to defaults

    MetaSpline (int _n, int _m, int _numOutputs);  // all other params are set to defaults
    
    MetaSpline();
    
    virtual ~MetaSpline ();
    
    void printData();
    
    void computeBcoeffs();
    
    float value(float t);
    
};
