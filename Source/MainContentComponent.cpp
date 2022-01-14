/*
  ==============================================================================

    AudioWithSplineModel.cpp
    Created: 15 Nov 2021 9:04:57am
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MainContentComponent.h"

MainContentComponent::MainContentComponent()
    : lEP(0.0), rEP(1200.0), signalScrollBar(false), state (Stopped)
{
    setButtonColours();
    addButtons();
    addSliders();
    addScrollbar();
          
    addAndMakeVisible(&graphView);
    
    lEP.referTo(graphView.leftEP);
    rEP.referTo(graphView.rightEP);
    
    setSize (1200, 800);

    formatManager.registerBasicFormats();
    transportSource.addChangeListener (this);
    
    lEP.addListener (this);
    rEP.addListener (this);

    setAudioChannels (0, 2);
    
    myADM.initialise(0, 2, nullptr, true);
    myADM.addAudioCallback(&player);
//    player.playTestSound();
}

void MainContentComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colours::lightgrey);
    g.setColour (juce::Colours::darkgrey);
    w = getWidth(); h = getHeight();
    drawScrollBox(g);
}

void MainContentComponent::changeListenerCallback (juce::ChangeBroadcaster* source)
{
    if (source == &transportSource)
    {
        if (transportSource.isPlaying())
            changeState (Playing);
        else if ((state == Stopping) || (state == Playing))
            changeState (Stopped);
        else if (Pausing == state)
            changeState (Paused);
    }
}

void MainContentComponent::resized()
{
    float w = getWidth();
    float h = getHeight();
    openButton.setBounds (10, 10, 70, 20);
    closeButton.setBounds (10, 70, 70, 20);
    playButton.setBounds (90, 10, 70, 20);
    stopButton.setBounds (170, 10, 70, 20);
    graphButton.setBounds (250, 10, 100, 20);
    shadeButton.setBounds (360, 10, 100, 20);
    targetButton.setBounds (470, 10, 100, 20);
    playCycleButton.setBounds (580, 10, 100, 20);
    computeModelButton.setBounds (690, 10, 100, 20);
    playModelButton.setBounds (800, 10, 100, 20);
    playCycleWithEnvButton.setBounds (910, 10, 100, 20);
    graphMetaSplinesButton.setBounds (1020, 10, 100, 20);
    signalScrollBar.setBounds (15, h-35, w-30, 20);
    graphView.setBounds (10, 100, w-20, h-145);
    auto freqGuessSliderLeft = 120;
    freqGuessSlider.setBounds (freqGuessSliderLeft, 40, 250, 20);
    freqGuessLabel.setFont(14.0f);
    auto frequencySliderLeft = 490;
    frequencySlider.setBounds (frequencySliderLeft, 40, 250, 20);
    frequencyLabel.setFont(14.0f);
    auto kValSliderLeft = 870;
    kValSlider.setBounds (kValSliderLeft, 40, 200, 20);
    kValLabel.setFont(14.0f);
    auto mValSliderLeft = 870;
    mValSlider.setBounds (mValSliderLeft, 70, 200, 20);
    mValLabel.setFont(14.0f);
    
}


void MainContentComponent::computeZeros()
{
    int freq = (int) freqGuess;     // = guess at frequency in Hz
    float lengthInSeconds = (float)sampleCount / (float)sampleRate;
    int periods = (int) (lengthInSeconds * (float)freqGuess) + 1;  // = # Cycles
    numCycles = periods;
    DBG("numCycles:  " << numCycles);
    int numKeyCycles = (int) ((float)numCycles / (float)mVal);
    DBG("numKeyCycles:  " << numKeyCycles);
    int numAllZeros = FindAllZerosFloat(sampleRate, periods, freq, floatBuffer, allZeros);
    float lastZero = allZeros[numAllZeros];
    int numZeros = FindZerosClosestToPeriods(sampleRate, periods, freq, cycleZeros, allZeros, samplesPerCycle, lastZero);
    lastSample = (int) lastZero;
    int actualLastSample = floatBuffer.getNumSamples() - 1;
    DBG("Sample before last zero: " << lastSample);
    DBG("Actual Last Sample: " << actualLastSample);
    DBG("Number of Periods: " << periods);
    DBG("Length in sec: " << lengthInSeconds);
    int totalSamples = 0;
    for (int i=0; i<samplesPerCycle.size(); i++) {
        totalSamples += samplesPerCycle[i];
    }
    DBG("Sum of Samples per Cycle: " << totalSamples);
    DBG("Last cycle zero: " << cycleZeros[cycleZeros.size()-1]);
    // this was testing setup of exp key pattern
//    keys.add(0); keys.add(1);
//    for (int j=0; j<7; j++) {
//        keys.add(keys[keys.size()-1] * 2);
//    }
//    keys.add(numCycles-1);
//    DBG("keys.size():  " << keys.size());
//    DBG("keys[size()-1]:  " << keys[keys.size()-1]);
}

void MainContentComponent::reComputeZeros()
{
    int freq = (int) freqGuess;     // = guess at frequency in Hz
    float lengthInSeconds = (float)sampleCount / (float)sampleRate;
    int periods = (int) (lengthInSeconds * freqGuess);  // = # Cycles
    int numAllZeros = allZeros.size() - 1;
    float lastZero = allZeros[numAllZeros];
    cycleZeros.clear();
    samplesPerCycle.clear();
    int numZeros = FindZerosClosestToPeriods(sampleRate, periods, freq, cycleZeros, allZeros, samplesPerCycle, lastZero);
    lastSample = (int) lastZero;
}

bool MainContentComponent::keyStateChanged(bool isKeyDown)
{
    if (KeyPress::isKeyCurrentlyDown(KeyPress::spaceKey))
    {
        DBG("SPACE DOWN");
        playButtonClicked();
    }
    KeyPress myKey = KeyPress('c', ModifierKeys::ctrlModifier, 0);
    KeyPress cmd1 = KeyPress('1', ModifierKeys::commandModifier, 0);
    KeyPress cmd3 = KeyPress('3', ModifierKeys::commandModifier, 0);
    if (cmd1.isCurrentlyDown())
    {
        DBG("cmd1 DOWN");
    }
    float incr = 0.001;
    if (KeyPress::isKeyCurrentlyDown(KeyPress::leftKey) || cmd1.isCurrentlyDown())
    {
        DBG("leftKey DOWN");
        graphView.magnify = 1.1 + incr;
        graphView.repaint();
    }
    if (KeyPress::isKeyCurrentlyDown(KeyPress::rightKey) || cmd3.isCurrentlyDown())
    {
        DBG("rightKey DOWN");
        graphView.magnify = 0.9 - incr;
        graphView.repaint();
    }
    if (KeyPress::isKeyCurrentlyDown(KeyPress::returnKey))
    {
        DBG("returnKey DOWN");
        int n = kVal + dVal;
        if (graphView.graphMetaSplinesOn) {
            graphView.metaSplineIndex += 1;
            DBG("graphing metaspline " << graphView.metaSplineIndex);
            if (graphView.metaSplineIndex == n) {
                graphView.metaSplineIndex = 0;
            }
            graphView.repaint();
        }
    }
    return true;
}

void MainContentComponent::valueChanged (Value &value)
{
    if ((value == lEP) || (value == rEP)) {
        // graph interval has changed, so need to change scrollBar size and/or position
        double newStart = ((double)lEP.getValue() / (double) sampleCount) * 1000;
        double newEnd =  ((double)rEP.getValue() / (double) sampleCount) * 1000;
        double newSize = newEnd - newStart;
        signalScrollBar.setCurrentRange(newStart, newSize, dontSendNotification);
        if ((int) newStart == 0) {
            graphView.hardLeft = true;
        }
        repaint();
    }
}



