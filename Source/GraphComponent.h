/*
  ==============================================================================

    GraphComponent.h
    Created: 18 Nov 2021 11:51:34am
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once

#include "zeros.h"
#include <string>
#include "../JuceLibraryCode/JuceHeader.h"
#include "CycleSpline.h"
#include "Cycle.h"
#include "bsp.h"

class GraphComponent    : public juce::Component
{
public:
    GraphComponent()
    {
        leftEP.setValue(0.0);
        rightEP.setValue(1200);
    }
    
    void setDataForGraph(short * _graphSamples, bool _audioLoaded, int _numSamples, float _magnify, float _leftEndPoint, float _rightEndPoint, unsigned _sampleCount);
    
    void setZerosForGraph( float * _zeros);
    
    short* graphSamples;
    bool audioLoaded = false;
    bool updateGraph = false;
    bool cyclesToGraph = false;
    bool callShadeCycles = false;
    bool hardLeft = true;
    float* zeros;
    int numSamples = 1000;
    float magnify = 1;
    float leftEndPoint;
    float rightEndPoint;
    float w, h;
    unsigned sampleCount;
    float addoffset = 0;
    float magfactor = 1;
    int startIndex = 0;
    int endIndex = 0;
    int highlightCycle = -1;  // positive number is cycle to highlight
    juce::Colour highlightColour;
    Array<juce::Colour> cycleColours;
//    Array<Cycle> cycles;
    Value leftEP;
    Value rightEP;
    juce::Point<int> doubleClick;
    
    juce::Point<float> signalToScreenCoords (juce::Point<float> P);
    juce::Point<float> screenToSignalCoords (juce::Point<float> Q);
    
    void setCycleColours();
    
    void paint (juce::Graphics& g) override;

    void resized() override;

private:
    
    void drawGraphBox(juce::Graphics& g, float w, float h);
    
    void mouseDown (const MouseEvent& event) override;
    
    void mouseDoubleClick (const MouseEvent& event) override;
    
    void graphSignal(juce::Graphics& g);

    void scaleInterval();

    int getMult(int numSamples);
    
    void drawDot(juce::Point<float> (P), juce::Graphics& g);
    
    void mouseWheelMove (const MouseEvent& event, const MouseWheelDetails & wheel) override;
    
    void mouseMagnify (const MouseEvent & event, float scaleFactor) override;
    
    void findCyclesToGraph();
    
    void shadeCycle(int n, juce::Graphics& g);
    
    void shadeCycles(juce::Graphics& g);
        
    int getCycleNumber(float t);
        
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (GraphComponent)
};



