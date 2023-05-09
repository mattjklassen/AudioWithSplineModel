/*
  ==============================================================================

    Buttons.cpp
    Created: 8 Jan 2022 3:02:15pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MainContentComponent.h"

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
            cycleZeros.clear();
            allZeros.clear();
            samplesPerCycle.clear();
            computeZeros();
            findMaxValuesPerCycle(maxSampleIndices, maxSampleValues, cycleZeros, floatBuffer);
            graphView.setDataForGraph(floatBuffer, audioDataLoaded, numSamples, magnify,
                                      leftEndPoint, rightEndPoint, sampleCount, sampleRate, kVal);
            graphView.setZerosForGraph(cycleZeros, allZeros, samplesPerCycle, freqGuess, maxSampleIndices, maxSampleValues);
        
            DBG("sample rate is:  " << (int)sampleRate);
            
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


void MainContentComponent::openOutput()
{
        auto file = outputFile;

        if (file != File{})                                                      
        {
            reader = formatManager.createReaderFor (file);  // file is InputStream
            readAudioData2(reader);
            audioDataLoaded = true;
            cycleZeros.clear();
            allZeros.clear();
            samplesPerCycle.clear();
            computeZeros();
            findMaxValuesPerCycle(maxSampleIndices, maxSampleValues, cycleZeros, floatBuffer);
            graphView.setDataForGraph(floatBuffer, audioDataLoaded, numSamples, magnify,
                                      leftEndPoint, rightEndPoint, sampleCount, sampleRate, kVal);
            graphView.setZerosForGraph(cycleZeros, allZeros, samplesPerCycle, freqGuess, maxSampleIndices, maxSampleValues);
        
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
}

void MainContentComponent::addButtons()
{
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
    
    addAndMakeVisible (&graphMetaSplinesButton);
    graphMetaSplinesButton.setButtonText ("Graph MetaSplines");
    graphMetaSplinesButton.onClick = [this] { graphMetaSplinesButtonClicked(); };
    graphMetaSplinesButton.setColour (juce::TextButton::buttonColourId, buttonColours[0]);
    
    addAndMakeVisible (&nextRandomButton);
    nextRandomButton.setButtonText ("Next Random");
    nextRandomButton.onClick = [this] { nextRandomButtonClicked(); };
    nextRandomButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);
    
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

    addAndMakeVisible (&playCycleWithEnvButton);
    playCycleWithEnvButton.setButtonText ("Play Cycle Env");
    playCycleWithEnvButton.onClick = [this] { playCycleWithEnvButtonClicked(); };
    playCycleWithEnvButton.setColour (juce::TextButton::buttonColourId, buttonColours[4]);

    addAndMakeVisible (&playModelButton);
    playModelButton.setButtonText ("Play Model");
    playModelButton.onClick = [this] { playModelButtonClicked(); };
    playModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[3]);
    
    addAndMakeVisible (&computeModelButton);
    computeModelButton.setButtonText ("Compute Model");
    computeModelButton.onClick = [this] { computeModelButtonClicked(); };
    computeModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[2]);

    addAndMakeVisible (&graphModelButton);
    graphModelButton.setButtonText ("Graph Model");
    graphModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[0]);
    graphModelButton.onClick = [this] { graphModelButtonClicked(); };

    addAndMakeVisible (&normalizeCycleLengthButton);
    normalizeCycleLengthButton.setButtonText ("Normalize Cycle Length");
    normalizeCycleLengthButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);
    normalizeCycleLengthButton.onClick = [this] { normalizeCycleLengthButtonClicked(); };
    normalizeCycleLengthButton.changeWidthToFitText();
    
    addAndMakeVisible (&randomizeButton);
    randomizeButton.setButtonText ("Randomize");
    randomizeButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);
    randomizeButton.onClick = [this] { randomizeButtonClicked(); };
    randomizeButton.changeWidthToFitText();
    
    addAndMakeVisible (&computeFFTButton);
    computeFFTButton.setButtonText ("Compute FFT");
    computeFFTButton.setColour (juce::TextButton::buttonColourId, buttonColours[3]);
    computeFFTButton.onClick = [this] { fftButtonClicked(); };
    computeFFTButton.changeWidthToFitText();
    
    addAndMakeVisible (&CAmodelButton);
    CAmodelButton.setButtonText ("CA model");
    CAmodelButton.setColour (juce::TextButton::buttonColourId, buttonColours[4]);
    CAmodelButton.onClick = [this] { CAmodelButtonClicked(); };
    CAmodelButton.changeWidthToFitText();
    
    addAndMakeVisible (&loadBcoeffsButton);
    loadBcoeffsButton.setButtonText ("load Bcoeffs");
    loadBcoeffsButton.setColour (juce::TextButton::buttonColourId, buttonColours[3]);
    loadBcoeffsButton.onClick = [this] { loadBcoeffsButtonClicked(); };
    loadBcoeffsButton.changeWidthToFitText();
    
    addAndMakeVisible (&averageBcoeffsButton);
    averageBcoeffsButton.setButtonText ("average Bcoeffs");
    averageBcoeffsButton.setColour (juce::TextButton::buttonColourId, buttonColours[4]);
    averageBcoeffsButton.onClick = [this] { averageBcoeffsButtonClicked(); };
    averageBcoeffsButton.changeWidthToFitText();
    
    addAndMakeVisible (&nextBsplineButton);
    nextBsplineButton.setButtonText ("next Bspline");
    nextBsplineButton.setColour (juce::TextButton::buttonColourId, buttonColours[2]);
    nextBsplineButton.onClick = [this] { nextBsplineButtonClicked(); };
    nextBsplineButton.changeWidthToFitText();
    
    addAndMakeVisible (&useDeltaModelButton);
    useDeltaModelButton.setButtonText ("use Delta Model");
    useDeltaModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);
    useDeltaModelButton.onClick = [this] { useDeltaModelButtonClicked(); };
    useDeltaModelButton.changeWidthToFitText();
    
    addAndMakeVisible (&useModelButton);
    useModelButton.setButtonText ("use Model");
    useModelButton.setColour (juce::TextButton::buttonColourId, buttonColours[1]);
    useModelButton.onClick = [this] { useModelButtonClicked(); };
    useModelButton.changeWidthToFitText();
    
    addAndMakeVisible (interpSelector);
    interpSelector.addItem ("Regular Cycle Interp",  1);
    interpSelector.addItem ("Exponential Cycle Interp",   2);
    interpSelector.addItem ("No Cycle Interp", 3);
    interpSelector.addItem ("Fibonacci Cycle Interp", 4);
    interpSelector.addItem ("Ends Only Cycle Interp", 5);
    interpSelector.addItem ("Other Cycle Interp", 6);
    interpSelector.setColour (juce::ComboBox::backgroundColourId, buttonColours[3]);
    interpSelector.onChange = [this] { interpSelectionChanged(); };
    interpSelector.setSelectedId (1);
    
    setSize (400, 200);
}

void MainContentComponent::setInterpSelectionsFalse()
{
    regularCycleInterp = false;
    expCycleInterp = false;
    noCycleInterp = false;
    fibonacciCycleInterp = false;
    endsOnlyCycleInterp = false;
    otherCycleInterp = false;
}

void MainContentComponent::interpSelectionChanged()
{
    setInterpSelectionsFalse();
    switch (interpSelector.getSelectedId())
    {
        case 1: regularCycleInterp = true;
//            DBG("regularCycleInterp = true");
                break;
        case 2: expCycleInterp = true;
                break;
        case 3: noCycleInterp = true;
                break;
        case 4: fibonacciCycleInterp = true;
                break;
        case 5: endsOnlyCycleInterp = true;
                break;
        case 6: otherCycleInterp = true;
                break;
        default: otherCycleInterp = true;
                break;
    }
}

void MainContentComponent::setButtonColours()
{
    buttonColours.add (juce::Colour::fromHSV (0.5f, 0.5f, 0.7f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.7f, 0.7f, 0.7f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.6f, 0.5f, 0.5f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.8f, 0.5f, 0.7f, 0.6f));
    buttonColours.add (juce::Colour::fromHSV (0.6f, 0.3f, 0.4f, 0.6f));
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
    DBG("highlightCycle = " << n);
    if ((n > -1) || graphView.graphMetaSplinesOn) {
        if (graphView.plotTargets) {
            graphView.plotTargets = false;
            DBG("setting plotTargets = false");
        } else {
            graphView.plotTargets = true;
            DBG("setting plotTargets = true");
        }
    }
    graphView.repaint();
}

void MainContentComponent::graphButtonClicked()
{
    if (graphView.graphMetaSplinesOn) {
        graphView.graphMetaSplinesOn = false;
        graphView.updateGraph = false;
        graphView.reset();
    }
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
    oneCycleWithEnvelope = false;
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
    oneCycleWithEnvelope = false;
    graphView.graphMetaSplinesOn = false;
    graphView.repaint();
}

void MainContentComponent::computeModelButtonClicked()
{
    writeModelToBuffer();
    DBG("spline model computed");
    writeWavFile();
    DBG("output.wav written");
    graphView.keysToAdd.clear();
    graphView.keysToRemove.clear();
    graphView.modelLoaded = true;
    graphView.setModelForGraph(writeBuffer);
    graphView.setAllCycleArray(allCycleArray);
    graphView.setDeltaCycleArray(deltaCycleArray);
}


void MainContentComponent::CAmodelButtonClicked()
{
    closeButtonClicked();
    modelWithCA = true;
    initECA();
    lengthInSecondsCA = 1;
//    freqGuess = 261.6;  // middle C
    freqGuess = 264;
    numCycles = (int)(freqGuess * lengthInSecondsCA);
    DBG("freqGuess: " << freqGuess << " numCycles: " << numCycles);
    keysCA = 5;
    mVal = 5;
    useADSR = true;
    writeScaleToBuffer();
    writeWavFile();
    DBG("output.wav written");
    graphView.keysToAdd.clear();
    graphView.keysToRemove.clear();
    graphView.modelLoaded = true;
    graphView.setModelForGraph(writeBuffer);
    openOutput();
    graphButtonClicked();
}

void MainContentComponent::loadBcoeffs(String filename, int d, int k, int numKeys, Array<float>& Bcoeffs)
{
    DBG("loading Bcoeffs array for: " << filename);
    int n = k + d;  // dimension or number of bcoeffs for instrument models
    auto file = File::getSpecialLocation(File::userHomeDirectory).getChildFile(filename);
    FileInputStream input (file);
    float val = 0;
    if (input.openedOk())
    {
        int size = (int) input.getTotalLength();
        // new comment
        char data[size];
        input.read(data, size);
        float *bcoeffs = reinterpret_cast<float*>(data);
        float max = 0;
        for (int i=0; i<numKeys; i++) {
            for (int p=0; p<n; p++) {
                val = bcoeffs[i*n+p];
                if (abs(val) > max) {
                    max = abs(val);
                }
                if (abs(val) > 0.5) {
//                    DBG("bcoeff[" << i*n+p << "]: " << val);
                }
                Bcoeffs.set(i*n+p, val);
            }
        }
//        DBG("max bcoeff: " << max);
    }
}

void MainContentComponent::averageBcoeffsButtonClicked()
{
    int n = 33, numKeys = 18;
    float val1 = 0, val2 = 0, val3 = 0;
    for (int i=0; i<numKeys; i++) {
        for (int p=0; p<n; p++) {
            val1 = guitarKeyBcoeffs[i*n+p];
            val2 = fluteKeyBcoeffs[i*n+p];
            val3 = celloKeyBcoeffs[i*n+p];
            keyBcoeffs.set(i*n+p, (1*val1 + 1*val2 + 1*val3)/3);
        }
    }
    modelWithSmall = true;
    modelWithDelta = true;
    loadSmallModel();
}

// Other instruments and timbres:
// violin, G scale from G string, includes middle C, also pizz and vibrato
// cello, redo, include middle C and vibrato/pizz
// splinusoid plain cubic
// alternating up/down on splinusoid coeffs
// B-spline coeff single central one, say number 16 out of 0-32
// several central B-coeffs say -.5,+.7,-.5 on 15,16,17
// maybe accent the octave by putting in half period repeat inside one cycle with bcoeffs
// this could be 7,8 then 24,25 with .7,-.7
// guitar redo, middle C and up, also harmonics, pizz, and vibrato

void MainContentComponent::loadBcoeffsButtonClicked()
{
    numInstr = 0;
    loadBcoeffs("KeyBcoeffs/guitar.keyBcoeffs", 3, 30, 18, guitarKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/guitarpizz.keyBcoeffs", 3, 30, 18, guitarpizzKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/guitarpont.keyBcoeffs", 3, 30, 18, guitarpontKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/guitarharm.keyBcoeffs", 3, 30, 18, guitarharmKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/theorbo.keyBcoeffs", 3, 30, 18, theorboKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/theorbopont.keyBcoeffs", 3, 30, 18, theorbopontKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/flute.keyBcoeffs", 3, 30, 18, fluteKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/cello.keyBcoeffs", 3, 30, 18, celloKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/cellopizz.keyBcoeffs", 3, 30, 18, cellopizzKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/cellopont.keyBcoeffs", 3, 30, 18, cellopontKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/marimba.keyBcoeffs", 3, 30, 18, marimbaKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/clarinet.keyBcoeffs", 3, 30, 18, clarinetKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/strings.keyBcoeffs", 3, 30, 18, stringsKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/trumpet.keyBcoeffs", 3, 30, 18, trumpetKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/handbell.keyBcoeffs", 3, 30, 18, handbellKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/piano.keyBcoeffs", 3, 30, 18, pianoKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/violin.keyBcoeffs", 3, 30, 18, violinKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/violinpizz.keyBcoeffs", 3, 30, 18, violinpizzKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/violinpont.keyBcoeffs", 3, 30, 18, violinpontKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/bassoon.keyBcoeffs", 3, 30, 18, bassoonKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/celeste.keyBcoeffs", 3, 30, 18, celesteKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/cimbalom.keyBcoeffs", 3, 30, 18, cimbalomKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/enghorn.keyBcoeffs", 3, 30, 18, enghornKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/frhorn.keyBcoeffs", 3, 30, 18, frhornKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/harp.keyBcoeffs", 3, 30, 18, harpKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/oboe.keyBcoeffs", 3, 30, 18, oboeKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/organ.keyBcoeffs", 3, 30, 18, organKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/vibraphone.keyBcoeffs", 3, 30, 18, vibraphoneKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/wurli.keyBcoeffs", 3, 30, 18, wurliKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/splinusoid.keyBcoeffs", 3, 30, 18, splinusoidKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/splinufuzz.keyBcoeffs", 3, 30, 18, splinufuzzKeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/splinumid1.keyBcoeffs", 3, 30, 18, splinumid1KeyBcoeffs); numInstr += 1;
    loadBcoeffs("KeyBcoeffs/splinumid2.keyBcoeffs", 3, 30, 18, splinumid2KeyBcoeffs); numInstr += 1;
    
    DBG("loaded " << numInstr << " instruments");
//    keyBcoeffs = guitarKeyBcoeffs;

//    int n = 33, numKeys = 18;
//    for (int i=0; i<numKeys; i++) {
//        for (int p=0; p<n; p++) {
//            DBG("" << keyBcoeffs[i*n+p]);
//        }
//    }
    
    modelWithSmall = true;
    modelWithDelta = true;
    loadSmallModel();
}

void MainContentComponent::fftButtonClicked()
{
    // first block is computing ECA examples, hijacking this button
    // placeholder hijack to get an int r to use for ECA
    int r = frequencySlider.getValue();
    Array<int> input;
    Array<int> output;
    for (int i=0; i<41; i++) {
        input.add(0);
    }
    input.set(20, 1);
    int val = 0;
    for (int k=0; k<21; k++) {
        for (int i=0; i<41; i++) {
            val = input[i];
            if (val == 0) {
//                cout << "  ";
                cout << "X";
            }
            if (val == 1) {
                cout << "1";
            }
        }
        cout << endl;
        input = computeCA(r, input);
    }
    cout << endl;
    return;
    
    // next block is unrelated to previous, doing fft
    nextFFTBlockReady = false;
    if (floatBuffer.getNumSamples() > fftSize) {
        for (int i=0; i < fftSize; i++) {
            pushNextSampleIntoFifo(floatBuffer.getSample(0, i));
        }
    } else {
        DBG("floatBuffer not full");
    }
    if (nextFFTBlockReady) {
        forwardFFT.performFrequencyOnlyForwardTransform (fftData.data());
        for (int j=0; j < 60; j++) {
            DBG("fftData[" << j << "]: " << fftData[j]);
        }
    }
    nextFFTBlockReady = false;
    if (floatBuffer.getNumSamples() > fftSize) {
        for (int i=0; i < fftSize; i++) {
            pushNextSampleIntoFifo(floatBuffer.getSample(0, i));
        }
    } else {
        DBG("floatBuffer not full");
    }
    if (nextFFTBlockReady) {
        forwardFFT.performRealOnlyForwardTransform (fftData.data());
        int k = 0;
        float real = 0, complex = 0, magnitude = 0;
        for (int j=0; j < 60; j++) {
            k = 2*j;
            real = fftData[k];
            complex = fftData[k+1];
            magnitude = real*real + complex*complex;
            magnitude = sqrt(magnitude);
            DBG("magnitude[" << j << "]: " << magnitude);
        }
    }
}

void MainContentComponent::playModelButtonClicked()
{
    oneCycleWithEnvelope = false;
    player.play(&writeBuffer, false, true);
}

void MainContentComponent::playCycleButtonClicked()
{
    if (graphView.highlightCycle == -1) {
        return;
    }
    playWavFileOn = false;
    playModelOn = false;
    oneCycleWithEnvelope = false;
    setSplineArrays();
    if (graphView.graphNewSplineCycle) {
        samplesPerSelectedCycle = graphView.cycleNew.outputs.size();
    } else {
        samplesPerSelectedCycle = graphView.cycleToGraph.outputs.size();
    }
    // toggle playing of selected cycle on/off
    if (playCycleOn) {
        playCycleOn = false;
        DBG("set playCycleOn = false");
    } else {
        playCycleOn = true;
        DBG("set playCycleOn = true");
    }
    prepareToPlay(512, sampleRate);
}

void MainContentComponent::playCycleWithEnvButtonClicked()
{
    oneCycleWithEnvelope = true;
    // use cycleToGraph to generate audio buffer with envelope
    computeCycleWithEnv();
    DBG("spline model computed");
    writeWavFile();
    DBG("output.wav written");
    player.play(&writeBuffer, false, true);
}

void MainContentComponent::graphMetaSplinesButtonClicked()
{
    graphView.setMetaSplinesForGraph(metaSplineArray);
    graphView.graphMetaSplinesOn = true;
    graphView.plotTargets = true;
    DBG("setting plotTargets to true");
    graphView.repaint();
}

void MainContentComponent::graphModelButtonClicked()
{
    graphView.graphMetaSplinesOn = false;
    if (graphView.modelLoaded) {
        if (graphView.updateModelGraph) {
            graphView.updateModelGraph = false;
        } else {
            graphView.updateModelGraph = true;
        }
    }
//    if 
    graphView.setKeys(keys);
    graphView.repaint();
}

void MainContentComponent::normalizeCycleLengthButtonClicked()
{
    if (normalizeCycleLengthButton.getToggleState()) {
        normalizeCycleLength = true;
        DBG("normalizeCycleLength is now: true");
    } else {
        normalizeCycleLength = false;
        DBG("normalizeCycleLength is now: false");
    }
    graphView.repaint();
}

void MainContentComponent::randomizeButtonClicked()
{
    if (randomizeButton.getToggleState()) {
        randomizeBcoeffs = true;
        DBG("randomizeBcoeffs is now: true");
    } else {
        randomizeBcoeffs = false;
        DBG("randomizeBcoeffs is now: false");
    }
}

void MainContentComponent::useDeltaModelButtonClicked()
{
    if (useDeltaModelButton.getToggleState()) {
        modelWithDelta = true;
        graphView.graphDeltaModel = true;
        DBG("useDeltaModel is now: true");
    } else {
        modelWithDelta = false;
        graphView.graphDeltaModel = false;
        DBG("useDeltaModel is now: false");
    }
}

void MainContentComponent::useModelButtonClicked()
{
    if (useModelButton.getToggleState()) {
        modelWithoutDelta = true;
        DBG("useModel is now: true");
    } else {
        modelWithoutDelta = false;
        DBG("useModel is now: false");
    }
}

void MainContentComponent::nextRandomButtonClicked()
{
    if (graphView.graphNewSplineCycle) {
        graphView.randomizeNewCycle();
        int n = kVal + 3;
        for (int i=0; i<n; i++) {
            if (graphView.graphNewSplineCycle) {
                controlCoeffs.set(4*i, graphView.cycleNew.bcoeffs[i]);
            }
        }
    }

}

void MainContentComponent::nextBsplineButtonClicked()
{
    graphView.nextBspline += 1;
    if (graphView.nextBspline == kVal + 3) {
        graphView.nextBspline = 0;
    }
    graphView.graphNextBspline();
}
