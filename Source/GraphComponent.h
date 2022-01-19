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
#include "MetaSpline.h"
#include "bsp.h"


class GraphComponent    : public juce::Component
//                          public juce::Value::Listener
{
public:
    GraphComponent()
    {
        leftEP.setValue(0.0);
        rightEP.setValue(1200);
    }
    
    void setDataForGraph(AudioBuffer<float>& _floatBuffer, bool _audioLoaded, int _numSamples, float _magnify, float _leftEndPoint, float _rightEndPoint, unsigned _sampleCount, unsigned sampleRate, int kVal);
    
    void setZerosForGraph(Array<float>& _cycleZeros, Array<float>& _allZeros, Array<int> _samplesPerCycle, float _freqGuess, Array<float>& _maxSampleIndices, Array<float>& _maxSampleValues);
    
    void setMetaSplinesForGraph(Array<MetaSpline>& _metaSplineArray);
    void setModelForGraph(AudioBuffer<float>& _modelBuffer);
    void setkVal(int _kVal);
    void setmVal(int _mVal);
    void setAmplitudeFactor(float _amplitudeFactor);
    
    AudioBuffer<float> floatBuffer;
    AudioBuffer<float> modelBuffer;
    bool audioLoaded = false;
    bool updateGraph = false;
    bool modelLoaded = false;
    bool updateModelGraph = false;
    bool cyclesToGraph = false;
    bool callShadeCycles = false;
    bool graphSplineCycle = false;
    bool graphMetaSplinesOn = false;
    bool plotTargets = false;
    bool hardLeft = true;
    int numSamples = 1000;
    int kVal = 20;
    int mVal = 1;
    float freqGuess = 70;
    float magnify = 1;
    float leftEndPoint;
    float rightEndPoint;
    float w, h;
    unsigned sampleCount;
    unsigned sampleRate;
    float addoffset = 0;
    float magfactor = 1;
    float amplitudeFactor = 1;
    int startIndex = 0;
    int endIndex = 0;
    int metaSplineIndex = 0;
    int highlightCycle = -1;  // positive number is cycle to highlight
    juce::Colour highlightColour;
    Array<juce::Colour> cycleColours;
    Array<float> cycleZeros;    // zeros marking endpoints of cycles in audio sample
    Array<float> allZeros;      // all zeros in audio sample
    Array<int> samplesPerCycle;
    Array<float> maxSampleIndices;
    Array<float> maxSampleValues;
    Value leftEP;
    Value rightEP;
    Value dragged;
    juce::Point<int> doubleClick;
    juce::Path freeCurve;
    bool drawCurve = false;
    
    juce::Point<float> signalToScreenCoords (juce::Point<float> P);
    juce::Point<float> screenToSignalCoords (juce::Point<float> Q);
    
    void setCycleColours();
    
    void paint (juce::Graphics& g) override;

    void resized() override;
    
    void graphMetaSplines(juce::Graphics& g);
    
    void graphLinearMetaSplines(juce::Graphics& g);

    CycleSpline cycleToGraph = CycleSpline(20, 0, 1);
    
//    Array<CycleSpline> cyclesToPlay;
    
    Array<MetaSpline> metaSplineArray;
    
    CycleSpline cycles[10];
    
    void reset();
    
    float metaSplineFactor = 1;
    
private:
    
    void drawGraphBox(juce::Graphics& g, float w, float h);
    
    void mouseDown (const MouseEvent& event) override;
    
    void mouseDrag (const MouseEvent& event) override;
    
    void mouseDoubleClick (const MouseEvent& event) override;
    
    void graphSignal(juce::Graphics& g);
    
    void graphModel(juce::Graphics& g);

    void graphSpline(juce::Graphics& g, CycleSpline& cycle);
    
    void scaleInterval();

    int getMult(int numSamples);
    
    void drawDot(juce::Point<float> (P), juce::Graphics& g);
    
    void mouseWheelMove (const MouseEvent& event, const MouseWheelDetails & wheel) override;
    
    void mouseMagnify (const MouseEvent & event, float scaleFactor) override;
    
    void findCyclesToGraph();
    
    void shadeCycle(int n, juce::Graphics& g);
    
    void shadeCycles(juce::Graphics& g);
        
    int getCycleNumber(float t);
    
    void plotTargetPoints(juce::Graphics& g, CycleSpline& cycle);

    void plotMetaSplineTargets(juce::Graphics& g);

        
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (GraphComponent)
};



