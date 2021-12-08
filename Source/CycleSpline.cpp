/*
  ==============================================================================

    CycleSpline.cpp
    Created: 17 Nov 2021 11:25:24am
    Author:  Matt Klassen

  ==============================================================================
*/

#include "CycleSpline.h"
#include "../JuceLibraryCode/JuceHeader.h"

CycleSpline::CycleSpline(int k, float a, float b)
{
    isKey = true;
    d = 3;
    float incr = (b - a) / k;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = k + 2*d;  // last index of knot sequence t_0,...,t_N
    
    subintervals.add(a);
    for (int i=1; i<k+1; i++)
    {
        subintervals.add(subintervals[i-1] + incr);
    }
    // now there are k subintervals inside interval [a,b] with a=subintervals[0], b=subintervals[k]
    
    inputs.add(a); inputs.add(a + incr / 2);
    for (int i=2; i<n-2; i++)
    {
        inputs.add(inputs[i-1] + incr);
    }
    inputs.add(b - incr / 2); inputs.add(b);
    
    knots.add(a - d * incr);
    for (int i=1; i<N+1; i++)
    {
        knots.add(knots[i-1] + incr);
    }
    
    for (int i=0; i<k+d; i++)
    {
        bcoeffs.add(0);
    }
}

CycleSpline:: ~CycleSpline ()
{
    
}


