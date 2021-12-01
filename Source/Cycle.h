/*
  ==============================================================================

    Cycle.h
    Created: 28 Nov 2021 7:53:52pm
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once
#include "../JuceLibraryCode/JuceHeader.h"
#include "GraphComponent.h"

class Cycle
{
private:
    
    float a, b; // for interval [a,b] on time axis, a = left end point, b = right end point
//    int n;      // cycle number
    
public:
    
    Cycle (float a, float b); 

    virtual ~Cycle ();
    
    void shadeCycle(juce::Graphics& g, float leftEP, float rEP);
    
};
