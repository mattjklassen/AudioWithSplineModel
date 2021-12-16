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
    
    addAndMakeVisible (&openButton);
    openButton.setButtonText ("Open");
    openButton.setColour (juce::TextButton::buttonColourId, buttonColours[0]);
    openButton.onClick = [this] { openButtonClicked(); };

    addAndMakeVisible (&playButton);
    playButton.setButtonText ("Play");
    playButton.onClick = [this] { playButtonClicked(); };
    playButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);
    playButton.setEnabled (false);

    addAndMakeVisible (&stopButton);
    stopButton.setButtonText ("Stop");
    stopButton.onClick = [this] { stopButtonClicked(); };
    stopButton.setColour (juce::TextButton::buttonColourId, buttonColours[2]);
    stopButton.setEnabled (false);
    
    addAndMakeVisible (&graphButton);
    graphButton.setButtonText ("Graph Signal");
    graphButton.onClick = [this] { graphButtonClicked(); };
    graphButton.setColour (juce::TextButton::buttonColourId, buttonColours[3]);
    
    addAndMakeVisible (&shadeButton);
    shadeButton.setButtonText ("Shade Cycles");
    shadeButton.onClick = [this] { shadeButtonClicked(); };
    shadeButton.setColour (juce::TextButton::buttonColourId, buttonColours[4]);
    
    addAndMakeVisible (&targetButton);
    targetButton.setButtonText ("Plot targets");
    targetButton.onClick = [this] { targetButtonClicked(); };
    targetButton.setColour (juce::TextButton::buttonColourId, buttonColours[0]);
    
    addAndMakeVisible (&playCycleButton);
    playCycleButton.setButtonText ("Play Cycle");
    playCycleButton.onClick = [this] { playCycleButtonClicked(); };
    playCycleButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);

    addAndMakeVisible (&signalScrollBar);
    signalScrollBar.setVisible(true);
    signalScrollBar.setRangeLimits(0, 1000, sendNotificationAsync);
    signalScrollBar.addListener(this);
    
    addAndMakeVisible (&freqGuessSlider);
    freqGuessSlider.setRange (20, 2000.0);
    freqGuessSlider.setTextValueSuffix (" Hz");
    freqGuessSlider.setNumDecimalPlacesToDisplay(2);
    freqGuessSlider.addListener (this);
    freqGuessSlider.setValue (220.0);
    freqGuessSlider.setSkewFactorFromMidPoint (100);
    
    addAndMakeVisible (&freqGuessLabel);
    freqGuessLabel.setText ("Frequency Guess", juce::dontSendNotification);
    freqGuessLabel.attachToComponent (&freqGuessSlider, true);
    
    addAndMakeVisible (&frequencySlider);
    frequencySlider.setRange (20, 2000.0);
    frequencySlider.setTextValueSuffix (" Hz");
    frequencySlider.setNumDecimalPlacesToDisplay(2);
    frequencySlider.addListener (this);
    frequencySlider.setValue (currentFrequency, juce::dontSendNotification);
    frequencySlider.setSkewFactorFromMidPoint (220);
    
    addAndMakeVisible (&frequencyLabel);
    frequencyLabel.setText ("Cycle Frequency", juce::dontSendNotification);
    frequencyLabel.attachToComponent (&frequencySlider, true);
    frequencySlider.onValueChange = [this]
    {
        if (currentSampleRate > 0.0)
            updateAngleDelta();
    };
    
    addAndMakeVisible (&kValSlider);
    kValSlider.setRange (5, 100);
    kValSlider.setNumDecimalPlacesToDisplay(0);
    kValSlider.addListener (this);
    kValSlider.setValue (20);
    kValSlider.setSkewFactorFromMidPoint (40);
    kValSlider.setTextBoxStyle (juce::Slider::TextBoxLeft, false, 40, kValSlider.getTextBoxHeight());
    
    addAndMakeVisible (&kValLabel);
    kValLabel.setText ("k = # subintervals", juce::dontSendNotification);
    kValLabel.attachToComponent (&kValSlider, true);
    
    addAndMakeVisible(&graphView);
    
    lEP.referTo(graphView.leftEP);
    rEP.referTo(graphView.rightEP);
    
    setSize (1200, 600);

    formatManager.registerBasicFormats();
    transportSource.addChangeListener (this);
    
    lEP.addListener (this);
    rEP.addListener (this);

    setAudioChannels (0, 2);
}

void MainContentComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colours::lightgrey);
    g.setColour (juce::Colours::blue);
    w = getWidth(); h = getHeight();
    
    drawScrollBox(g);
}

void MainContentComponent::setButtonColours()
{
    buttonColours.add (juce::Colour::fromHSV (0.5f, 0.5f, 0.7f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.7f, 0.7f, 0.7f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.6f, 0.5f, 0.5f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.8f, 0.5f, 0.7f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.6f, 0.3f, 0.4f, 0.6f));
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
    playButton.setBounds (90, 10, 70, 20);
    stopButton.setBounds (170, 10, 70, 20);
    graphButton.setBounds (250, 10, 100, 20);
    shadeButton.setBounds (360, 10, 100, 20);
    targetButton.setBounds (470, 10, 100, 20);
    playCycleButton.setBounds (580, 10, 100, 20);
    signalScrollBar.setBounds (15, h-35, w-30, 20);
    graphView.setBounds (10, 70, w-20, h-115);
    auto freqGuessSliderLeft = 120;
    freqGuessSlider.setBounds (freqGuessSliderLeft, 40, 250, 20);
    freqGuessLabel.setFont(14.0f);
    auto frequencySliderLeft = 490;
    frequencySlider.setBounds (frequencySliderLeft, 40, 250, 20);
    frequencyLabel.setFont(14.0f);
    auto kValSliderLeft = 870;
    kValSlider.setBounds (kValSliderLeft, 40, 200, 20);
    kValLabel.setFont(14.0f);

}
// getWidth() - sliderLeft - 10
void MainContentComponent::readAudioData2 (AudioFormatReader *reader) {
    
    sampleCount = (int) reader->lengthInSamples;
    sampleRate = (int) reader->sampleRate;
    floatBuffer.setSize ((int) reader->numChannels, (int) reader->lengthInSamples);
    reader->read(&floatBuffer,                     // juce AudioBuffer <float>
                 0,                               // start sample in buffer
                 (int) reader->lengthInSamples,   // number of samples in file data
                 0,                               // start sample to fill in buffer
                 true,                            // use Left channel (0)
                 false);                          // use Right channel (1)
}

void MainContentComponent::readAudioData (File file) {
//    old version
//    juce::FileInputStream inputStream (file);
//    inputStream.read(buffer, 44);
//    dataSize = *reinterpret_cast<unsigned*>(buffer+40);
//    sampleRate = *reinterpret_cast<unsigned*>(buffer+24);
//    sampleCount = dataSize/2;
//    data = new char[dataSize];
//    inputStream.read(data, dataSize);
//    samples = reinterpret_cast<short*>(data);
//    std::cout << "size of wav file data is:  " << dataSize << std::endl;
//    std::cout << "sample rate from wav file is:  " << sampleRate << std::endl;
//    std::cout << "number of data samples is:  " << sampleCount << std::endl;
}

void MainContentComponent::computeZeros()
{
    int freq = (int) freqGuess;     // = guess at frequency in Hz
    float lengthInSeconds = (float)sampleCount / (float)sampleRate;
    int periods = (int) (lengthInSeconds * freqGuess);  // = # Cycles
    int numAllZeros = FindAllZerosFloat(sampleRate, periods, freq, floatBuffer, allZeros);
    float lastZero = allZeros[numAllZeros];
    int numZeros = FindZerosClosestToPeriods(sampleRate, periods, freq, cycleZeros, allZeros, samplesPerCycle, lastZero);
}

void MainContentComponent::reComputeZeros()
{
    int freq = (int) freqGuess;     // = guess at frequency in Hz
    float lengthInSeconds = (float)sampleCount / (float)sampleRate;
    int periods = (int) (lengthInSeconds * freqGuess);  // = # Cycles
    int numAllZeros = allZeros.size() - 1;
    float lastZero = allZeros[numAllZeros];
    cycleZeros = Array<float>();
    samplesPerCycle = Array<int>();
    int numZeros = FindZerosClosestToPeriods(sampleRate, periods, freq, cycleZeros, allZeros, samplesPerCycle, lastZero);
}

void MainContentComponent::shadeButtonClicked()
{
    if (graphView.audioLoaded)
    {
        graphView.updateGraph = true;
        graphView.setCycleColours();
        
        if (graphView.callShadeCycles) {
            graphView.callShadeCycles = false;
        } else {
            graphView.callShadeCycles = true;
        }
    }
    graphView.repaint();
}

void MainContentComponent::targetButtonClicked()
{
    int n = graphView.highlightCycle;
    if (n > -1) {
        if (graphView.plotTargets) {
            graphView.plotTargets = false;
        } else {
            graphView.plotTargets = true;
        }
    }
    graphView.repaint();
}

void MainContentComponent::graphButtonClicked()
{
    if (graphView.audioLoaded) {
        if (graphView.updateGraph) {
            graphView.updateGraph = false;
        } else {
            graphView.updateGraph = true;
        }
    }
    graphView.repaint();
}

void MainContentComponent::stopButtonClicked()
{
    if (state == Paused)
        changeState (Stopped);
    else
        changeState (Stopping);
}

void MainContentComponent::playButtonClicked()
{
//        std::cout << "play button clicked" << std::endl;
    if ((state == Stopped) || (state == Paused))
        changeState (Starting);
    else if (state == Playing)
        changeState (Pausing);
}

void MainContentComponent::openButtonClicked()
{
//        std::cout << "open button clicked" << std::endl;
    chooser = std::make_unique<juce::FileChooser> ("Select a Wave file to play...",
                                                   juce::File{},
                                                   "*.wav");
    auto chooserFlags = juce::FileBrowserComponent::openMode
                      | juce::FileBrowserComponent::canSelectFiles;

    chooser->launchAsync (chooserFlags, [this] (const FileChooser& fc)
    {
        auto file = fc.getResult();

        if (file != File{})                                                      
        {
            reader = formatManager.createReaderFor (file);
            readAudioData2(reader);
            audioDataLoaded = true;
            initSplineArrays();
            computeZeros();
            graphView.setDataForGraph(floatBuffer, audioDataLoaded, numSamples, magnify,
                                      leftEndPoint, rightEndPoint, sampleCount, sampleRate, kVal);
            graphView.setZerosForGraph(cycleZeros, allZeros, samplesPerCycle, freqGuess);
        
            // setting up to play back audio file (as in JUCE tutorial PlayingSoundFilesTutorial)
            if (reader != nullptr)
            {
                auto newSource = std::make_unique<juce::AudioFormatReaderSource> (reader, true);   // [11]
                transportSource.setSource (newSource.get(), 0, nullptr, reader->sampleRate);       // [12]
                playButton.setEnabled (true);                                                      // [13]
                readerSource.reset (newSource.release());                                          // [14]
            }
        }
    });

}

int MainContentComponent::getCycleNumber(float t)
{
    int i = 0;
    while (cycleZeros[i] < t) {
        i++;
    }
    return i-1;
}

void MainContentComponent::drawScrollBox(juce::Graphics& g)
{
    juce::Path path;
    path.startNewSubPath (juce::Point<float> (10, h-40));
    path.lineTo (juce::Point<float> (w-10, h-40));
    path.lineTo (juce::Point<float> (w-10, h-10));
    path.lineTo (juce::Point<float> (10, h-10));
    path.lineTo (juce::Point<float> (10, h-40));
    float myThickness = 1.0;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.strokePath (path, myType);
}

void MainContentComponent::scrollBarMoved (ScrollBar* scrollBarThatHasMoved ,
                                 double newRangeStart)
{
    double newRangeL = scrollBarThatHasMoved->getCurrentRange().getStart();
    graphView.leftEndPoint = (float) newRangeL / 1000 * (float) graphView.sampleCount;
    graphView.rightEndPoint = graphView.leftEndPoint + (float) graphView.numSamples;
    if (graphView.rightEndPoint > graphView.sampleCount)
    {
        graphView.rightEndPoint = (float) graphView.sampleCount;
        graphView.leftEndPoint = graphView.rightEndPoint - (float) graphView.numSamples;
    }
    if (newRangeL == 0) {
        graphView.hardLeft = true;
        DBG("set hardLeft true");
    }
    graphView.repaint();
}

void MainContentComponent::sliderValueChanged (juce::Slider* slider)
{
    if (! audioDataLoaded) {
        return;
    }
    if (slider == &freqGuessSlider)
    {
        freqGuess = freqGuessSlider.getValue();
        reComputeZeros();
        graphView.setZerosForGraph(cycleZeros, allZeros, samplesPerCycle, freqGuess);
        graphView.callShadeCycles = false;
    }
    if (slider == &kValSlider)
    {
        kVal = kValSlider.getValue();
        graphView.setkVal(kVal);
        graphView.repaint();
    }
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

void MainContentComponent::changeState (TransportState newState)
{
    if (state != newState)
    {
        state = newState;

        switch (state)
        {
            case Stopped:
                playButton.setButtonText ("Play");
                stopButton.setButtonText ("Stop");
                stopButton.setEnabled (false);
                transportSource.setPosition (0.0);
                break;

            case Starting:
                transportSource.start();
                break;

            case Playing:
                playButton.setButtonText ("Pause");
                stopButton.setButtonText ("Stop");
                stopButton.setEnabled (true);
                break;

            case Pausing:
                transportSource.stop();
                break;

            case Paused:
                playButton.setButtonText ("Resume");
                stopButton.setButtonText ("Return to Zero");
                break;

            case Stopping:
                transportSource.stop();
                break;
        }
    }
}

void MainContentComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    if (playCycleOn) {
        currentSampleRate = sampleRate;
        updateAngleDelta();
    } else {
        transportSource.prepareToPlay (samplesPerBlockExpected, sampleRate);
    }
}

void MainContentComponent::updateAngleDelta()
{
    auto cyclesPerSample = frequencySlider.getValue() / currentSampleRate;
    // cycles/sample = fractions of sampling frequency
    angleDelta = cyclesPerSample * 2.0 * juce::MathConstants<double>::pi; 
    // angelDelta = frequency in radians/sample
}

void MainContentComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    if (playCycleOn) {
        auto* leftBuffer  = bufferToFill.buffer->getWritePointer (0, bufferToFill.startSample);
        auto* rightBuffer = bufferToFill.buffer->getWritePointer (1, bufferToFill.startSample);
        auto Pi = juce::MathConstants<double>::pi;
        auto localTargetFrequency = targetFrequency;
        auto frequencyIncrement = (localTargetFrequency - currentFrequency) / bufferToFill.numSamples;
        for (auto sample = 0; sample < bufferToFill.numSamples; ++sample)
        {
            // compute currentSample value with spline of selected cycle:
            if (currentAngle > 2.0 * Pi) {
                currentAngle -= 2.0 * Pi;
            }
            float currentSample = computeSpline(currentAngle / (2 * Pi));
            currentFrequency += frequencyIncrement;
            currentAngle += angleDelta;
            leftBuffer[sample]  = currentSample;
            rightBuffer[sample] = currentSample;
        }
    } else {
        if (readerSource.get() == nullptr)
        {
            bufferToFill.clearActiveBufferRegion();
            return;
        }
        transportSource.getNextAudioBlock (bufferToFill);
    }
}

void MainContentComponent::playCycleButtonClicked()
{
    setSplineArrays();
    // toggle playing of selected cycle on/off
    if (playCycleOn) {
        playCycleOn = false;
    } else {
        playCycleOn = true;
    }
    prepareToPlay(512, sampleRate);
}

void MainContentComponent::setSplineArrays()
{
    int d = 3;
    int k = graphView.cycleToGraph.k;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
//    controlCoeffs.clear();
//    knotVals.clear();
    for (int i=0; i<n*n; i++) {
        controlCoeffs.set(i, 0);
    }
    float incr = 1 / (float) k;
    for (int i=0; i<N+1; i++) {
        knotVals.set(i, (i-d) * incr);
//        DBG("knot vals [" << i << "] = " << knotVals[i]);
    }
    for (int i=0; i<n; i++) {
        controlCoeffs.set(n*i, graphView.cycleToGraph.bcoeffs[i]);
//        DBG("control coeffs [" << n*i << "] = " << controlCoeffs[n*i]);
    }
}

float MainContentComponent::computeSpline(float t)
{
    // assume t is in [0,1] and output is 0 at the ends
    int d = 3;
//    int k = graphView.cycleToGraph.k;
    int k = kVal;
    float output = 0;
    int J = 0;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
//    juce::Array<float> controlCoeffs;
//    for (int i=0; i<n*n; i++) {
//        controlCoeffs.add(0);
//    }
//    for (int i=0; i<d+k; i++) {
//        controlCoeffs.set((k+d)*i+0, graphView.cycleToGraph.bcoeffs[i]);
//    }
//    float incr = 1 / (float) k;
//    juce::Array<float> knotVals;
//    for (int i=0; i<N+1; i++) {
//        knotVals.add((i-d) * incr);
//    }
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

void MainContentComponent::initSplineArrays() {
    
    
}
