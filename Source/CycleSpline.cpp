
/*
  ==============================================================================

    CycleSpline.cpp
    Created: 17 Nov 2021 11:25:24am
    Author:  Matt Klassen

  ==============================================================================
*/

#include "CycleSpline.h"

CycleSpline::CycleSpline(int k, float a, float b)
{
    isKey = true;
    d = 3;
    float incr = (b - a) / k;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = k + 2*d;  // last index of knot sequence t_0,...,t_N
    subintervals = new float[k+1];
    subintervals[0] = a;
    for (int i=1; i<k+1; i++)
    {
        subintervals[i] = subintervals[i-1] + incr;
    }
    inputs = new float[n];
    inputs[0] = a; inputs[1] = a + incr / 2;
    inputs[n-1] = b; inputs[n-2] = b - incr / 2;
    for (int i=2; i<n-2; i++)
    {
        inputs[i] = inputs[i-1] + incr;
    }
    knots[0] = a - d * incr;
    for (int i=1; i<N+1; i++)
    {
        knots[i] = knots[i-1] + incr;
    }
    for (int i=0; i<k+d; i++)
    {
        bcoeffs[i] = 0;
    }
}

CycleSpline:: ~CycleSpline ()
{
    
}
