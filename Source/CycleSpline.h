/*
  ==============================================================================

    CycleSpline.h
    Created: 17 Nov 2021 11:25:24am
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once
//#include "zeros.h"
//#include "bsp.h"
#include "../JuceLibraryCode/JuceHeader.h"

class CycleSpline
{
public:
    
    int d;  // degree of B-splines, default d = 3
    int k; // number of subintervals, can vary, use default 20
    float a, b; // for interval [a,b] on time axis, a = left end point, b = right end point
    Array<float> subintervals;  // default subinterval size = (b-a)/k, endpoints: u_0,...,u_k (k+1 values)
    Array<float> inputs;  // default at a, b, and subinterval breakpoints, and midpoints of outer two subintervals
    Array<float> knots;   // default uniform knot sequence includes subinterval endpoints, and d more on each end
    Array<float> bcoeffs; // B-spline coefficients, d+k of them
    // The above data are enough to compute the cycle values at samples between a and b.
    // If some data changes, the other data may need to be recomputed with reference to the original audio
    // data, using linear interpolation to compute new output values to be matched with splines.
    bool isKey;  // if true (default) then compute spline output values from above data
                 // if false, then compute values by interpolation using key cycles

    CycleSpline (int k, float a, float b);  // all other params are set to defaults

    virtual ~CycleSpline ();
    
};

