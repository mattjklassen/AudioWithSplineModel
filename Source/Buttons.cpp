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
    playCycleWithEnvButton.setButtonText ("Play Cycle with Env");
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
    if ((n > -1) || graphView.graphMetaSplinesOn) {
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
    if (mVal > 1) {
        noCycleInterp = false;
        regularCycleInterp = true;
    }
    writeModelToBuffer();
//    writeCycleInterpModelToBuffer();
    DBG("spline model computed");
    writeWavFile();
    DBG("output.wav written");
}

void MainContentComponent::playModelButtonClicked()
{
    oneCycleWithEnvelope = false;
    player.play(&writeBuffer, false, true);
}

void MainContentComponent::playCycleButtonClicked()
{
    playWavFileOn = false;
    playModelOn = false;
    oneCycleWithEnvelope = false;
    setSplineArrays();
    // toggle playing of selected cycle on/off
    if (playCycleOn) {
        playCycleOn = false;
//        DBG("set playCycleOn = false");
    } else {
        playCycleOn = true;
//        DBG("set playCycleOn = true");
//        DBG("controlCoeffs: ");
//        for (int i=0; i<controlCoeffs.size(); i++) {
//            DBG(controlCoeffs[i]);
//        }
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
    graphView.repaint();
}
