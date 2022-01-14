/*
  ==============================================================================

    MetaSpline.cpp
    Created: 27 Dec 2021 6:44:11pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MetaSpline.h"
#include "../JuceLibraryCode/JuceHeader.h"

MetaSpline::MetaSpline()
{
    d = 3; n = 20; a = 0; b = 1; m=10;
}

MetaSpline::MetaSpline(int _n, int _m)
{
    // sets all parameters and arrays of MetaSpline except for targets, bcoeffs and outputs
    d = 3; n = _n; a = 0; b = 1; m = _m;
    int k = n - d;
    float incr = 1 / (float)k;
    //  n = dimension of splines, also number of inputs plus 2 (for natural cubic splines)
    //  inputs.size() = targets.size() = keys.size() = n-2
    int N = n + d;  // last index of knot sequence t_0,...,t_N
    
    subintervals.add(0);
    for (int i=1; i<k+1; i++)
    {
        subintervals.add(subintervals[i-1] + incr);
    }
    // now there are k subintervals for the interval [0,1]
    
    knots.add( - d * incr);
    for (int i=1; i<N+1; i++)
    {
        knots.add(knots[i-1] + incr);
    }
    // this is a simple uniform knot sequence
    
    // changing inputs to be suitable for natural cubic spline, ie. n-2=k+1 of them
    // inputs: 0,1/k,2/k,...,1
    inputs.add(0);
    for (int i=1; i<k+1; i++)
    {
        inputs.add(inputs[i-1] + incr);
    }
    // k+1=n-2 inputs are uniformly spaced in interval [0,1]
    
    for (int i=0; i<n; i++)
    {
        bcoeffs.add(0);
        targets.add(0);
    }
}

MetaSpline::MetaSpline(int _n, int _m, int _numOutputs)
{
    // sets all parameters and arrays of MetaSpline except for targets, bcoeffs and outputs
    d = 3; n = _n; a = 0; b = 1; m = _m; numOutputs = _numOutputs;
    int k = n - d;
    float incr = 1 / (float)k;
    //  n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
    
    subintervals.add(0);
    for (int i=1; i<k+1; i++)
    {
        subintervals.add(subintervals[i-1] + incr);
    }
    // now there are k subintervals inside interval [0,1]
    
    knots.add( - d * incr);
    for (int i=1; i<N+1; i++)
    {
        knots.add(knots[i-1] + incr);
    }
    // this is simple uniform knot sequence
    
    // changing inputs to be suitable for natural cubic spline, ie. n-2=k+1 of them
    // inputs: 0,1/k,2/k,...,1
    inputs.add(0);
    for (int i=1; i<k+1; i++)
    {
        inputs.add(inputs[i-1] + incr);
    }
    // k+1=n-2 inputs are uniformly spaced in interval [0,1]
    
    for (int i=0; i<n; i++)
    {
        bcoeffs.add(0);
        targets.add(0);
    }
    
    for (int i=0; i<numOutputs; i++)
    {
        outputs.add(0);
    }
}

MetaSpline:: ~MetaSpline ()
{
    
}

float MetaSpline::value(float t)  // assumes we have bcoeffs
{
    // assume t is in [0,1], 
    float output = 0;
    int J = 0;
    int k = n - d;  // n = dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
    juce::Array<float> controlCoeffs;
    for (int i=0; i<n*n; i++) {
        controlCoeffs.add(0);
    }
    for (int i=0; i<d+k; i++) {
        controlCoeffs.set((k+d)*i+0, bcoeffs[i]);
    }
    float incr = 1 / (float) k;
    juce::Array<float> knotVals;
    for (int i=0; i<N+1; i++) {
        knotVals.add((i-d) * incr);
    }
    if ((t > 0) && (t < 1)) {
        for (int i=1; i<N; i++)
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

        for (int p=1; p<d+1; p++)
        {
            for (int i=J-d+p; i<J+1; i++)
            {
              float denom = (knotVals[i+d-(p-1)]-knotVals[i]);
              float fac1 = (t-knotVals[i]) / denom;
              float fac2 = (knotVals[i+d-(p-1)]-t) / denom;
              controlCoeffs.set((k+d)*i+p, fac1 * controlCoeffs[(k+d)*i+(p-1)]
                  + fac2 * controlCoeffs[(k+d)*(i-1)+(p-1)]);
            }
        }
        output = controlCoeffs[(k+d)*J+d];
    }
    return output;
}


void MetaSpline::printData()
{
    DBG("degree d = " << d << ", k = " << n-d << ", a = " << a << ", b = " << b);
    DBG("numOutputs:  " << numOutputs);
//    DBG("subintervals: ");
//    for (int i=0; i<k+1; i++) {
//        DBG("u[" << i << "] = " << subintervals[i]);
//    }
//    DBG("inputs: ");
//    for (int i=0; i<k+d; i++) {
//        DBG("inputs[" << i << "] = " << inputs[i]);
//    }
//    DBG("diffs: ");
//    for (int i=1; i<k+d; i++) {
//        DBG("back diffs[" << i << "] = " << inputs[i] - inputs[i-1]);
//    }
    DBG("targets: ");
    for (int i=0; i<n; i++) {
        DBG("targets[" << i << "] = " << targets[i]);
    }
    DBG("knots: ");
    for (int i=0; i<n+d; i++) {
        DBG("knots[" << i << "] = " << knots[i]);
    }
    if (bcoeffs.size() > 0) {
        DBG("bcoeffs: ");
        for (int i=0; i<n; i++) {
            DBG("bcoeffs[" << i << "] = " << bcoeffs[i]);
        }
    }
    if (outputs.size() > 0) {
        DBG("outputs: ");
        for (int i=0; i<outputs.size(); i++) {
            DBG("outputs[" << i << "] = " << outputs[i]);
        }
    }
}
