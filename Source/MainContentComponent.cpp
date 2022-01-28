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
    graphView.setBounds (10, 100, w-40, h-145);
    signalScrollBar.setBounds (15, h-35, w-52, 20);
    auto ampSliderLeft = w - 25;
    amplitudeSlider.setBounds (ampSliderLeft, 100, 20, h-145);
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
    graphModelButton.setBounds (90, 70, 100, 20);
    interpSelector.setBounds (200, 70, 200, 20);
    normalizeCycleLengthButton.setBounds (410, 70, 200, 20);
    randomizeButton.setBounds (600, 70, 200, 20);
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
    // set normalized cycle lengths
    for (int i=0; i<cycleZeros.size(); i++) {
        normalizedCycleZeros.set(i, cycleZeros[1] * i);
    }
    samplesPerCycleGuess = (float) sampleRate / (float) freqGuess;
    float averageSamplesPerCycle = cycleZeros[cycleZeros.size()-1] / (float) cycleZeros.size();
    DBG("average samples per cycle:  " << averageSamplesPerCycle);
    DBG("number of zeros:  " << numAllZeros);
    DBG("samplesPerCycleGuess:  " << samplesPerCycleGuess);
//    samplesPerCycle.sort();
//    for (int i=0; i<samplesPerCycle.size()-1; i++) {
//        cout << samplesPerCycle[i] << " ";
//    }
}

void MainContentComponent::reComputeZeros()
{
    int freq = (int) freqGuess;     // = guess at frequency in Hz
    float lengthInSeconds = (float)sampleCount / (float)sampleRate;
    int periods = (int) (lengthInSeconds * freqGuess);  // = # Cycles
    numCycles = periods;
    DBG("numCycles:  " << numCycles);
    int numKeyCycles = (int) ((float)numCycles / (float)mVal);
    DBG("numKeyCycles:  " << numKeyCycles);
    int numAllZeros = allZeros.size() - 1;
    float lastZero = allZeros[numAllZeros];
    cycleZeros.clear();
    samplesPerCycle.clear();
    int numZeros = FindZerosClosestToPeriods(sampleRate, periods, freq, cycleZeros, allZeros, samplesPerCycle, lastZero);
    lastSample = (int) lastZero;
    samplesPerCycleGuess = (float) sampleRate / (float) freqGuess;
    DBG("samplesPerCycleGuess:  " << samplesPerCycleGuess);
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
    if (value == cycleRendered) {
//        DBG("value changing: cycleRendered = " << (int) cycleRendered.getValue());
//        randomizeCycleBcoeffs(cycleRendered.getValue());
    }
}

void MainContentComponent::setNew(CycleSpline& cycle)
{
    int k = cycle.k;
    int d = cycle.d;
    int n = k + d;
    int N = n + d;
    cycle.inputs.clear();
    cycle.targets.clear();
    cycle.knots.clear();
    
    float incr = 1 / (float) k;
    cycle.inputs.add(0.5*incr);
    for (int i=1; i<k; i++) {
        cycle.inputs.add(i*incr);
    }
    cycle.inputs.add(1-0.5*incr);
    
    // New knot sequence: 0,0,0,0,1/k,2/k,...,(k-1)/k,1,1,1,1
    for (int i=0; i<d+1; i++) {
        cycle.knots.add(0);
    }
    for (int i=1; i<k; i++) {
        cycle.knots.add(cycle.knots[i-1]+incr);
    }
    for (int i=0; i<d+1; i++) {
        cycle.knots.add(1);
    }
    DBG("New params set for CycleSpline knot sequence: ");
    DBG("N = " << N << "n = " << n << "k = " << k);
    for (int i=0; i<N+1; i++) {
        std::cout << cycle.knots[i] << " ";
    }
    DBG("");
}

void MainContentComponent::setParabolicTargets(float scale)
{
    int k = cycleParabola.k;
    int d = cycleParabola.d;
    int n = k + d;
    
    cycleParabola.inputs.clear();
    cycleParabola.targets.clear();
    cycleParabola.knots.clear();
    
    float incr = 1 / (float) k;
    cycleParabola.inputs.add(0);
    cycleParabola.inputs.add(0.5*incr);
    for (int i=1; i<k; i++) {
        cycleParabola.inputs.add(i*incr);
    }
    cycleParabola.inputs.add(1-0.5*incr);
    cycleParabola.inputs.add(1);
    
    for (int i=0; i<cycleParabola.inputs.size(); i++) {
        float t = cycleParabola.inputs[i];
//        float scale = 0.1;
        cycleParabola.targets.set(i, scale*4*t*(1-t));
    }
    cycleParabola.knots.add(-d*incr);
    for (int i=1; i<n; i++) {
        cycleParabola.knots.add(cycleParabola.knots[i-1]+incr);
    }
}

void MainContentComponent::computeParabolicBcoeffs()
{
    int k = cycleParabola.k, d = cycleParabola.d;
    int n = k + d;
    float val = 0;
    juce::Array<float> A;
    juce::Array<float> B;
    juce::Array<float> temp;
    for (int i=0; i<n*n; i++) {
        A.add(0);
        B.add(0);
        temp.add(0);
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            val = bSplineVal(k, j, cycleParabola.inputs[i]);  // A[i,j]
            A.set(i*n+j, val);
            temp.set(i*n+j, val);
            B.set(i*n+j, 0);
        }
    }
    for (int i=0; i<n; i++) {
        B.set(i*n+i, 1.0);
    }

    gaussElim(n, temp, B);
    Array<float> x;
    Array<float> y;
    for (int i=0; i<n; i++) {
        x.add(0);
        y.add(0);
    }
    y = multMatCol(n, B, cycleParabola.targets);
    for (int i=0; i<n; i++) {
        cycleParabola.bcoeffs.set(i, y[i]);
    }
}
