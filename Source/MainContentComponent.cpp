/*
  ==============================================================================

    MainContentComponent.cpp
    Created: 29 Nov 2021 11:39:48am
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MainContentComponent.h"

MainContentComponent::MainContentComponent()
    :  lEP(0.0), rEP(1200.0), signalScrollBar(false), state (Stopped)
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
    
    addAndMakeVisible (&signalScrollBar);
    signalScrollBar.setVisible(true);
    signalScrollBar.setRangeLimits(0, 1000, sendNotificationAsync);
    signalScrollBar.addListener(this);
    
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
    signalScrollBar.setBounds (15, h-35, w-30, 20);
    graphView.setBounds (10, 70, w-20, h-115);
}

void MainContentComponent::readAudioData (File file) {
    juce::FileInputStream inputStream (file);
    inputStream.read(buffer, 44);
    dataSize = *reinterpret_cast<unsigned*>(buffer+40);
    sampleRate = *reinterpret_cast<unsigned*>(buffer+24);
    sampleCount = dataSize/2;
    data = new char[dataSize];
    inputStream.read(data, dataSize);
    samples = reinterpret_cast<short*>(data);
    std::cout << "size of wav file data is:  " << dataSize << std::endl;
    std::cout << "sample rate from wav file is:  " << sampleRate << std::endl;
    std::cout << "number of data samples is:  " << sampleCount << std::endl;

}

void MainContentComponent::computeZeros()
{
    int freq = 225;
    int periods = 300;  // = Cycles
    float *allzeros = new float[10 * periods];
    float *zeros = new float[2*periods];
    float *SPP = new float[2*periods];
    int numAllZeros = FindAllZeros(periods, freq, data, allzeros);
//        std::cout << "number of (all) zeros:  " << numAllZeros << std::endl;
    float lastZero = allzeros[numAllZeros];
//        std::cout << "last zero:  " << lastZero << std::endl;
    int numZeros = FindZerosClosestToPeriods(periods, freq, zeros, allzeros, SPP, lastZero);
    for (int i=0; i<numZeros; i++) {
        cycleZeros.add (zeros[i]);
    }
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
                                                   "*.wav");                     // [7]
    auto chooserFlags = juce::FileBrowserComponent::openMode
                      | juce::FileBrowserComponent::canSelectFiles;

    chooser->launchAsync (chooserFlags, [this] (const FileChooser& fc)           // [8]
    {
        auto file = fc.getResult();

        if (file != File{})                                                      // [9]
        {
            auto* reader = formatManager.createReaderFor (file);                 // [10]

            readAudioData (file);
            audioDataLoaded = true;
            computeZeros();
            graphView.setDataForGraph(samples, audioDataLoaded, numSamples, magnify, leftEndPoint, rightEndPoint, sampleCount);
            graphView.setZerosForGraph(cycleZeros.getRawDataPointer());
        
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

void MainContentComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    if (readerSource.get() == nullptr)
    {
        bufferToFill.clearActiveBufferRegion();
        return;
    }
    transportSource.getNextAudioBlock (bufferToFill);
}

