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

    addAndMakeVisible (&closeButton);
    closeButton.setButtonText ("Close");
    closeButton.setColour (juce::TextButton::buttonColourId, buttonColours[3]);
    closeButton.onClick = [this] { closeButtonClicked(); };
    
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

    addAndMakeVisible (&playModelButton);
    playModelButton.setButtonText ("Play Model");
    playModelButton.onClick = [this] { playModelButtonClicked(); };
    playModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[3]);
    
    addAndMakeVisible (&computeModelButton);
    computeModelButton.setButtonText ("Compute Model");
    computeModelButton.onClick = [this] { computeModelButtonClicked(); };
    computeModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[2]);
    
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
    kValSlider.setValue (10);
    kValSlider.setSkewFactorFromMidPoint (40);
    kValSlider.setTextBoxStyle (juce::Slider::TextBoxLeft, false, 40, kValSlider.getTextBoxHeight());
    
    addAndMakeVisible (&kValLabel);
    kValLabel.setText ("k = # subintervals", juce::dontSendNotification);
    kValLabel.attachToComponent (&kValSlider, true);
    
    addAndMakeVisible (&mValSlider);
    mValSlider.setRange (1, 100);
    mValSlider.setNumDecimalPlacesToDisplay(0);
    mValSlider.addListener (this);
    mValSlider.setValue (20);
    mValSlider.setSkewFactorFromMidPoint (40);
    mValSlider.setTextBoxStyle (juce::Slider::TextBoxLeft, false, 40, mValSlider.getTextBoxHeight());
    
    addAndMakeVisible (&mValLabel);
    mValLabel.setText ("m = key cycle mult", juce::dontSendNotification);
    mValLabel.attachToComponent (&mValSlider, true);
    
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
    closeButton.setBounds (10, 70, 70, 20);
    playButton.setBounds (90, 10, 70, 20);
    stopButton.setBounds (170, 10, 70, 20);
    graphButton.setBounds (250, 10, 100, 20);
    shadeButton.setBounds (360, 10, 100, 20);
    targetButton.setBounds (470, 10, 100, 20);
    playCycleButton.setBounds (580, 10, 100, 20);
    computeModelButton.setBounds (690, 10, 100, 20);
    playModelButton.setBounds (800, 10, 100, 20);
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
    int periods = (int) (lengthInSeconds * (float)freqGuess) + 1;  // = # Cycles
    numCycles = periods;
    DBG("numCycles:  " << numCycles);
    int numKeyCycles = (int) ((float)numCycles / (float)mVal);
    DBG("numKeyCycles:  " << numKeyCycles);
    int numAllZeros = FindAllZerosFloat(sampleRate, periods, freq, floatBuffer, allZeros);
    float lastZero = allZeros[numAllZeros];
    int numZeros = FindZerosClosestToPeriods(sampleRate, periods, freq, cycleZeros, allZeros, samplesPerCycle, lastZero);
    lastSample = (int) lastZero;
    float actualLastSample = floatBuffer.getNumSamples();
    DBG("Sample before last zero: " << lastSample);
    DBG("Actual Last Sample: " << actualLastSample);
    DBG("Number of Periods: " << periods);
    DBG("Length in sec: " << lengthInSeconds);
    int totalSamples = 0;
    for (int i=0; i<samplesPerCycle.size(); i++) {
        totalSamples += samplesPerCycle[i];
    }
    DBG("Sum of Samples per Cycle: " << totalSamples);
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
    playCycleOn = false;
    playModelOn = false;
    playWavFileOn = true;
    if ((state == Stopped) || (state == Paused))
        changeState (Starting);
    else if (state == Playing)
        changeState (Pausing);
}

void MainContentComponent::closeButtonClicked()
{
    cycleZeros = Array<float>();
    allZeros = Array<float>();
    graphView.updateGraph = false;
    graphView.audioLoaded = false;
    graphView.callShadeCycles = false;
    graphView.graphSplineCycle = false;
    graphView.plotTargets = false;
    graphView.highlightCycle = -1;
    playCycleOn = false;
    playModelOn = false;
    playWavFileOn = false;
    graphView.repaint();
}

void MainContentComponent::openButtonClicked()
{
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
            reader = formatManager.createReaderFor (file);  // file is InputStream
            readAudioData2(reader);
            audioDataLoaded = true;
            initSplineArrays();
            cycleZeros.clear();
            allZeros.clear();
            samplesPerCycle.clear();
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
                playWavFileOn = true;
                prepareToPlay(512, sampleRate);
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
        DBG("new kVal:  " << kVal);
        graphView.setkVal(kVal);
        graphView.repaint();
    }
    if (slider == &mValSlider)
    {
        mVal = mValSlider.getValue();
        DBG("new mVal:  " << mVal);
        graphView.setmVal(mVal);
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
    }
    if (playModelOn) {
        currentSampleRate = sampleRate;
    }
    if (playWavFileOn) {
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

// Need to do other presets for cycle interpolation, like exponential:
// 0,1,2,4,8,16,32,64,128,256,...
// Need to break up next function into smaller pieces:
// 1. prep key cycles
// 2. write key cycles to buffer
// 3. write key cycle bcoeffs to array
// 4. prep meta-splines
// 5. prep non-key cycles
// 6. write non-key cycles to buffer

// Write model to buffer with uniform (or other) cycle interpolation
void MainContentComponent::writeModelToBuffer()
{
    if (mVal == 1) {
        if (expCycleInterp) {
            
        } else {
            writeBuffer.clear();
            writeBuffer.setSize(1, lastSample+1);  // channels, samples
            float a = 0, b = 1;
            CycleSpline cycle = CycleSpline(kVal, a, b);
            // compute each cycle and write out the computed samples to writeBuffer
            auto* leftBuffer = writeBuffer.getWritePointer (0, 0);
            int currentSampleIndex = 0;
            for (int i=0; i<numCycles; i++)
            {
                a = cycleZeros[i];
                b = cycleZeros[i+1];
                cycle = CycleSpline(kVal, a, b);
                computeCycleBcoeffs(cycle, floatBuffer);
                computeCycleSplineOutputs(cycle);
                int M = cycle.outputs.size();
                for (int j=0; j<M; ++j) {
                    leftBuffer[currentSampleIndex] = cycle.outputs[j];
                    currentSampleIndex += 1;
                }
            }
        }
    }
    if (mVal > 1) {
        // Do cycle interpolation with key cycle indices 0,m,2m,...
        writeBuffer.clear();
        writeBuffer.setSize(1, lastSample+1);  // channels, samples
        Array<float> keyBcoeffs;
        int numKeyCycles = numCycles / mVal;  // dim of meta splines
        float a = 0, b = 1;
        CycleSpline cycle = CycleSpline(kVal, a, b);
        auto* leftBuffer  = writeBuffer.getWritePointer (0, 0);
        for (int t=0; t<lastSample+1; t++) {
//            writeBuffer.setSample(0, t, 0);
            leftBuffer[t] = 0;
        }
        int i = 0;
        int k = kVal;
        int d = dVal;
        int n = k + d;  // dim of cycle splines
        // compute key cycles and write outputs to buffer
        while (i < numKeyCycles)
        {
            a = cycleZeros[mVal*i];
            b = cycleZeros[mVal*i+1];
            cycle = CycleSpline(kVal, a, b);
            // compute each key cycle bcoeffs and store in array of floats
            computeCycleBcoeffs(cycle, floatBuffer);
            for (int p=0; p<n; p++) {
                keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
            }
            computeCycleSplineOutputs(cycle);
            int M = cycle.outputs.size();
            int L = (int) a + 1;  // first sample index to write
            for (int j=0; j<M; j++) {
                writeBuffer.setSample(0, j+L, cycle.outputs[j]);
            }
            i += 1;
        }
        // Now key cycles are written to writeBuffer, and key cycle bcoeffs are in array keyBcoeffs
        // next compute bcoeffs for intermediate cycles (from meta-splines) and write outputs
        // meta-splines use the keyBcoeffs as targets
        MetaSpline spline = MetaSpline(numKeyCycles, mVal);
        Array<MetaSpline> metaSplineArray;
        for (int i=0; i<n; i++) {  // loop on each bcoeff of cycles with same fixed kVal
            spline = MetaSpline(numKeyCycles, mVal,numCycles);   // n = numKeyCycles (targets)
            for (int j=0; j<numKeyCycles; j++) {
                spline.targets.set(j, keyBcoeffs[j*n+i]);
            }
            computeMetaSplineBcoeffs(spline);
//            spline.printData();
            computeMetaSplineOutputs(spline);
//            spline.printData();
            metaSplineArray.add(spline);
//            DBG("size of metaSplineArray:  " << metaSplineArray.size());
//            spline.printData();
//            DBG("mVal :  " << mVal << "  kVal :  " << kVal);
//            DBG("numKeyCycles: " << numKeyCycles);
        }  // now we have n meta-splines, one for each of n bcoeffs
        
        // compute bcoeffs and write outputs for intermediate cycles
        i = 0;  // i is cycle number, now compute all cycles up to the last key cycle
        n = kVal + dVal;
        while (i < mVal * numKeyCycles) {
            if (i % mVal == 0) { // key cycles already written
            } else {
                // compute bcoeffs for cycle i from meta-spline and write outputs
                // cycle i starting sample is (int) cycle.a + 1
                a = cycleZeros[i];
                b = cycleZeros[i+1];
                cycle = CycleSpline(kVal, a, b);  // this is cycle i
                for (int j=0; j<n; j++) {
                    cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
                }
                computeCycleSplineOutputs(cycle);
//                cycle.printData();
//                metaSplineArray[5].printData();
                int M = cycle.outputs.size();
                int L = (int) a + 1;  // first sample index to write
                for (int j=0; j<M; j++) {
                    writeBuffer.setSample(0, j+L, cycle.outputs[j]);
                }
            }
            i += 1;
        }
    }
}

void MainContentComponent::writeWavFile()
{
    if (sampleRate > 0)
    {
        // Create an OutputStream to write to our destination file...
        outputFile.deleteFile();
        if (auto fileStream = std::unique_ptr<FileOutputStream> (outputFile.createOutputStream()))
        {
            // Now create a WAV writer object that writes to our output stream...
            WavAudioFormat wavFormat;
            if (auto writer = wavFormat.createWriterFor (fileStream.get(), sampleRate, 1, 16, {}, 0))
            {
                fileStream.release(); // (passes responsibility for deleting the stream to the writer object that is now using it)
                if (writer != nullptr) {
                        writer->writeFromAudioSampleBuffer (writeBuffer, 0, writeBuffer.getNumSamples());
                }
                delete writer;
            }
        }
    }
}

void MainContentComponent::computeModelButtonClicked()
{
    writeModelToBuffer();
    DBG("spline model computed");
    writeWavFile();
    DBG("output.wav written");
}

void MainContentComponent::playModelButtonClicked()
{
    player.play(&writeBuffer, false, true);
//    player.play(&floatBuffer, false, true);
//    playModelOn = true;
//    playCycleOn = false;
//    playWavFileOn = false;
//    prepareToPlay(512, sampleRate);
}

void MainContentComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
//    if (playModelOn) {
//    }
    if (playCycleOn) {
        auto* leftBuffer  = bufferToFill.buffer->getWritePointer (0, bufferToFill.startSample);
        auto* rightBuffer = bufferToFill.buffer->getWritePointer (1, bufferToFill.startSample);
        int control = 0;
        auto Pi = juce::MathConstants<double>::pi;
        auto localTargetFrequency = targetFrequency;
        auto frequencyIncrement = (localTargetFrequency - currentFrequency) / bufferToFill.numSamples;
        for (auto sample = 0; sample < bufferToFill.numSamples; ++sample)
        {
            // compute currentSample value with spline of selected cycle:
            if (currentAngle > 2.0 * Pi) {
                currentAngle -= 2.0 * Pi;
                control += 1;
                if (control == 3) {
                    control = 0;
                }
            }
            float currentSample = computeSpline(control, currentAngle / (2 * Pi));
            currentFrequency += frequencyIncrement;
            currentAngle += angleDelta;
            leftBuffer[sample]  = currentSample;
            rightBuffer[sample] = currentSample;
        }
    }
    if (playWavFileOn)
    {
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
    playWavFileOn = false;
    playModelOn = false;
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
    // changing the array of control coeffs to be 3*n rows and 4 cols
    // so c_i,j is now c[4*(i+m*n)+j], for cols j=0,1,2,3 and 3 sets of rows for m=0,1,2
    for (int i=0; i<12*n; i++) {  // (3*n)*4 = rows*cols 
        controlCoeffs.set(i, 0);
    }
    float incr = 1 / (float) k;
    for (int i=0; i<N+1; i++) {
        knotVals.set(i, (i-d) * incr);
//        DBG("knot vals [" << i << "] = " << knotVals[i]);
    }
    for (int i=0; i<n; i++) {
        controlCoeffs.set(4*i, graphView.cycleToGraph.bcoeffs[i]);
//        DBG("control coeffs [" << n*i << "] = " << controlCoeffs[n*i]);
    }
    for (int m=0; m<3; m++) {
        for (int i=0; i<n; i++) {
            controlCoeffs.set(4*(i+m*n), graphView.cycles[m].bcoeffs[i]);
        }
    }
}

void MainContentComponent::fillBcoeffs()
{
//    if (cycleRendered == 0) {
//        // fill controlCoeffs with bcoeffs in first column
//    }
}

float MainContentComponent::computeSpline(int control, float t)
{
    // assume t is in [0,1] and output is 0 at the ends
    int m = control;
    m = 0;  // set this to ignore the other two cycles in controlCoeffs
    int d = 3;
    int k = kVal;
    float output = 0;
    int J = 0;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
//    juce::Array<float> controlCoeffs and knotVals are created and updated outside this function
//    in order to keep those out of the AudioCallback, increasing efficiency dramatically.
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
        // controlCoeffs[i,j] is controlCoeffs[4*i+j], row=i,col=j
        // m=0,1,2 shifts controlCoeffs by m*n
        for (int p=1; p<d+1; p++)
        {
            for (int i=J-d+p; i<J+1; i++)
            {
              float denom = (knotVals[i+d-(p-1)]-knotVals[i]);
              float fac1 = (t-knotVals[i]) / denom;
              float fac2 = (knotVals[i+d-(p-1)]-t) / denom;
              controlCoeffs.set(4*(i+m*n)+p, fac1 * controlCoeffs[4*(i+m*n)+(p-1)]
                  + fac2 * controlCoeffs[4*(i+m*n-1)+(p-1)]);
            }
        }
        output = controlCoeffs[4*(J+m*n)+d];
    }
    return output;
}

void MainContentComponent::initSplineArrays() {
    // create the spline arrays juce::Array<float> controlCoeffs and knotVals here
    // controlCoeffs can be about 200 rows, but only needs four columns if d = 3.
}

// Cycle Interpolation ToDo List:

// The cycle interpolation model should have chosen key cycles and meta data which shows how
// to construct the interpolating cycles, including meta-splines used to generate the B-spline
// coefficients of non-key cycles, and lengths of intervals for each cycle (if changing).

// If we create the simplest model first, it should be to use the cycle lengths and samples per cycle
// from a wav file, so the model just follows the original data but creates a spline for each cycle
// and then generates sample data from those splines.  This method is used if m=1 to compute model.

// If m>1 then we are doing cycle interpolation with key cycles j*m, j=0,1,2,...

// Suppose for now that we are going to use constant cycle length, and only the B-spline
// coefficients are changing from cycle to cycle, not the cycle lengths or subinterval sizes.

// To deal with the realtime rendering of spline model cycle by cycle:
// I think I will have three sets of arrays which are initialized with the B-spline coeffs
// for the first three cycles.  When the first cycle is finished rendering, we can start to
// write the coefficients of the 4th cycle over those for the first cycle.  The AudioCallBack
// can then pass an int variable outside to say which cycle has finished rendering.  This
// could be a Value type which has a listener that responds by writing the next B-spline
// coefficients.  So cycles are organized mod 3. When cycle i finishes rendering, overwrite
// its B-spline coefficients with those for cycle i + 3.
// write cycleBcoeffs for i=0,1,2
// render 0, overwrite with 3
// render 1, overwrite with 4
// render 2, overwrite with 5 etc.
// The time it takes to write one cycle's bcoeffs should not be longer than the render/playback
// of one cycle.  This way the new coeffs will be in place before needed.
// If new cycle bcoeffs are being computed from a meta-spline, these should still be computable
// in realtime in the above scheme.

// The selection of key cycles should be done with some preset schemes.
// 1. (Simple regular) Select every m^th cycle as key, and interpolate all others
// - this method will kill oscillating patterns between key cycles, and possibly also subharmonics
// * meta-splines should be used to generate B-spline coefficients of non-key cycles.
// * meta-splines should be applied in groups of cycles of reasonable length.
// ex: 400 cycles in 1 sec, use every 10th cycle as key, so 40 key cycles.
// - suppose sample rate 48 kHz, f_0 guess 400, so 120 samples per cycle
// - suppose each key cycle has n=30 interpolation points (targets), so there are 30 meta-splines
// which describe the coefficients as they evolve through the non-key cycles.
// The number of data points in this example is 1) the B-spline coeffs of key cycles: 40*30=1200
// 2) the meta-spline coeffs for each of 30 meta-splines, with 40 targets (from key cycles): 1200
// so 2400 float values, compared with 48000 original samples, or 1/20 = 5% of the original data.

// 2. (Exponential scheme) Select cycles with indices 0,1,2,4,8,16,32...
// - meta-splines could still be used even though the separation of targets is nonuniform.
// - in fact, one could even still use uniform separation for the splines but then apply them
//   to the nonuniform intervals.
// ex. as above, now with only 8 key cycles: 0,1,2,4,8,16,32,40 (include last one).
// So data is now 1) 8*30=240 2) also 240, so 480 floats, or data reduction of 480/48000=.01
// for a reduction to 1% of original data.

// An important point to consider in the above examples is that one can make modifications to
// the meta-splines without changing the reduction of data, for instance to give some oscillations
// to the pattern of non-key cycles, or to experiment with the degrees of the splines, etc.

// A meta-spline can simply be defined on the interval [0,1] with equally spaced intermediate
// points defining subintervals as before.  If we are working on the spline for bcoeff b_0, and
// we have k+1 key cycles, then we can define the subintervals, knot sequence, and target points
// similarly to how we did for the cycle splines.  The target values are now B-spline coeffs at each of
// of the key cycles.

// Another note regarding computation: we can choose a smaller cycle length in order to lower the
// size of linear systems n=k+d but then create a sequence of these smaller cycles in chunks to
// perform the role of key cycles.  For cycle interpolation, however, this would only work if the
// chunk of cycles is still part of a sort of periodic pattern.






