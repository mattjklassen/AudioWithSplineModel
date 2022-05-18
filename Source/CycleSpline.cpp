/*
  ==============================================================================

    CycleSpline.cpp
    Created: 17 Nov 2021 11:25:24am
    Author:  Matt Klassen

  ==============================================================================
*/

#include "CycleSpline.h"
#include "../JuceLibraryCode/JuceHeader.h"

CycleSpline::CycleSpline()
{
    d = 3; k = 20; a = 0; b = 1; y0 = 0; y1 = 0;
}

// now this intialization is changed to use new knot sequence, inputs, and y0, y1
CycleSpline::CycleSpline(int _k, float _a, float _b)
{
    isKey = true;
    d = 3; k = _k; a = _a; b = _b;
    y0 = 0; y1 = 0;
    float incr = (b - a) / k;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
    
    subintervals.add(a);
    for (int i=1; i<k+1; i++)
    {
        subintervals.add(subintervals[i-1] + incr);
    }
    // now there are k subintervals inside interval [a,b] with a=subintervals[0], b=subintervals[k]
    
    incr = 1 / (float) k;
    inputs.add(0.5*incr);
    for (int i=1; i<k; i++) {
        inputs.add(i*incr);
    }
    inputs.add(1-0.5*incr);
    
    // New knot sequence: 0,0,0,0,1/k,2/k,...,(k-1)/k,1,1,1,1
    for (int i=0; i<d+1; i++) {
        knots.add(0);
    }
    for (int i=d+1; i<N-d; i++) {
        knots.add(knots[i-1]+incr);
    }
    for (int i=0; i<d+1; i++) {
        knots.add(1);
    }
    
    for (int i=1; i<k+d; i++)
    {
        bcoeffs.add(0);
    }
}

CycleSpline:: ~CycleSpline ()
{
    
}

float CycleSpline::value(float t)
{
    // assume t is in [0,1] and output is 0 at the ends
    float output = 0;
    int J = 0;
    int n = k + d;  // dimension of splines, also number of inputs
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

void CycleSpline::writeData()
{
    int n = k + d;
    auto file = File::getSpecialLocation(File::userHomeDirectory).getChildFile("cycle.txt");
    FileOutputStream output (file);
    if (output.openedOk())
    {
        output.setPosition (0);  // default would append
        output.truncate();
        for (int i=0; i<n; i++) {
            output.writeFloat(bcoeffs[i]);
        }
    }
    
//    DBG("degree d = " << d << ", k = " << k << ", a = " << a << ", b = " << b);
//    DBG("y0 = " << y0 << ", y1 = " << y1);
//    DBG("subintervals: ");
//    for (int i=0; i<k+1; i++) {
//        DBG("u[" << i << "] = " << subintervals[i]);
//    }
//    if (inputs.size() > 0) {
//        DBG("inputs: ");
//        for (int i=0; i<inputs.size(); i++) {
//            DBG("inputs[" << i << "] = " << inputs[i]);
//        }
//    }
//    DBG("diffs: ");
//    for (int i=1; i<k+d; i++) {
//        DBG("back diffs[" << i << "] = " << inputs[i] - inputs[i-1]);
//    }
//    if (targets.size() > 0) {
//        DBG("targets: ");
//        for (int i=0; i<targets.size(); i++) {
//            DBG("targets[" << i << "] = " << targets[i]);
//        }
//    }
//    if (knots.size() > 0) {
//        DBG("knots: ");
//        for (int i=0; i<k+2*d+1; i++) {
//            DBG("knots[" << i << "] = " << knots[i]);
//        }
//    }
//    if (bcoeffs.size() > 0) {
//        DBG("bcoeffs: ");
//        for (int i=0; i<k+d; i++) {
//            DBG("bcoeffs[" << i << "] = " << bcoeffs[i]);
//        }
//    }
//    if (outputs.size() > 0) {
//        DBG("outputs: ");
//        for (int i=0; i<outputs.size(); i++) {  // outputs.size()
//            DBG("outputs[" << i << "] = " << outputs[i]);
//        }
//    }
}

void CycleSpline::printData()
{
    DBG("degree d = " << d << ", k = " << k << ", a = " << a << ", b = " << b);
    DBG("y0 = " << y0 << ", y1 = " << y1);
//    DBG("subintervals: ");
//    for (int i=0; i<k+1; i++) {
//        DBG("u[" << i << "] = " << subintervals[i]);
//    }
//    if (inputs.size() > 0) {
//        DBG("inputs: ");
//        for (int i=0; i<inputs.size(); i++) {
//            DBG("inputs[" << i << "] = " << inputs[i]);
//        }
//    }
//    DBG("diffs: ");
//    for (int i=1; i<k+d; i++) {
//        DBG("back diffs[" << i << "] = " << inputs[i] - inputs[i-1]);
//    }
//    if (targets.size() > 0) {
//        DBG("targets: ");
//        for (int i=0; i<targets.size(); i++) {
//            DBG("targets[" << i << "] = " << targets[i]);
//        }
//    }
//    if (knots.size() > 0) {
//        DBG("knots: ");
//        for (int i=0; i<k+2*d+1; i++) {
//            DBG("knots[" << i << "] = " << knots[i]);
//        }
//    }
    if (bcoeffs.size() > 0) {
        DBG("bcoeffs: ");
        for (int i=0; i<k+d; i++) {
            DBG("bcoeffs[" << i << "] = " << bcoeffs[i]);
        }
    }
//    if (outputs.size() > 0) {
//        DBG("outputs: ");
//        for (int i=0; i<outputs.size(); i++) {  // outputs.size()
//            DBG("outputs[" << i << "] = " << outputs[i]);
//        }
//    }
}

// old version with simple knots sequence
//CycleSpline::CycleSpline(int _k, float _a, float _b)
//{
//    isKey = true;
//    d = 3; k = _k; a = _a; b = _b;
//    y0 = 0; y1 = 0;
//    float incr = (b - a) / k;
//    int n = k + d;  // dimension of splines, also number of inputs
//    int N = n + d;  // last index of knot sequence t_0,...,t_N
//
//    subintervals.add(a);
//    for (int i=1; i<k+1; i++)
//    {
//        subintervals.add(subintervals[i-1] + incr);
//    }
//    // now there are k subintervals inside interval [a,b] with a=subintervals[0], b=subintervals[k]
//
//    inputs.add(a);             // inputs[0]
//    inputs.add(a + incr / 2);  // inputs[1]
//    inputs.add(a + incr);      // inputs[2]
//    for (int i=3; i<n-2; i++)
//    {
//        inputs.add(inputs[i-1] + incr);  // inputs[3...(n-3)]
//    }
//    inputs.add(b - incr / 2);  // inputs[n-2]
//    inputs.add(b);             // inputs[n-1]
//
//    incr = 1 / (float) k;
//    knots.add(- d * incr);
//    for (int i=1; i<N+1; i++)
//    {
//        knots.add(knots[i-1] + incr);
//    }
//
//    for (int i=1; i<k+d; i++)
//    {
//        bcoeffs.add(0);
//    }
//}
