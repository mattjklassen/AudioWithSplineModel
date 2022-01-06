/*
  ==============================================================================

   This file is part of the JUCE tutorials.
   Copyright (c) 2020 - Raw Material Software Limited

   The code included in this file is provided under the terms of the ISC license
   http://www.isc.org/downloads/software-support-policy/isc-license. Permission
   To use, copy, modify, and/or distribute this software for any purpose with or
   without fee is hereby granted provided that the above copyright notice and
   this permission notice appear in all copies.

   THE SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY, AND ALL WARRANTIES,
   WHETHER EXPRESSED OR IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR
   PURPOSE, ARE DISCLAIMED.

  ==============================================================================
*/

/*******************************************************************************
 The block below describes the properties of this PIP. A PIP is a short snippet
 of code that can be read by the Projucer and used to generate a JUCE project.

 BEGIN_JUCE_PIP_METADATA

 name:             PlayingSoundFilesTutorial
 version:          1.0.0
 vendor:           JUCE
 website:          http://juce.com
 description:      Plays audio files.

 dependencies:     juce_audio_basics, juce_audio_devices, juce_audio_formats,
                   juce_audio_processors, juce_audio_utils, juce_core,
                   juce_data_structures, juce_events, juce_graphics,
                   juce_gui_basics, juce_gui_extra
 exporters:        xcode_mac, vs2019, linux_make

 type:             Component
 mainClass:        MainContentComponent

 useLocalCopy:     1

 END_JUCE_PIP_METADATA

*******************************************************************************/


#pragma once

#include "zeros.h"
#include <string>
#include "../JuceLibraryCode/JuceHeader.h"
#include "CycleSpline.h"
#include "bsp.h"
#include "GraphComponent.h"



//==============================================================================
class MainContentComponent   : public juce::AudioAppComponent,
                               public juce::ChangeListener,
                               public juce::ScrollBar::Listener,
                               public juce::Value::Listener,
                               public juce::Slider::Listener
{
public:
    
    MainContentComponent();

    ~MainContentComponent() override
    {
        shutdownAudio();
    }
    
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;

    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override;

    void releaseResources() override
    {
        transportSource.releaseResources();
    }

    void resized() override;

    void changeListenerCallback (juce::ChangeBroadcaster* source) override;
    
    float computeSpline(int control, float t);
    
    bool audioDataLoaded = false;

    Value lEP;  // left End Point for graphing interval in graphView
    Value rEP;  // right End Point for graphing interval in graphView
//    Value cycleRendered;  // set to 1 or true when cycle is rendered

    
private:
    
    enum TransportState
    {
        Stopped, Starting, Playing, Pausing, Paused, Stopping
    };

    void scrollBarMoved (ScrollBar* scrollBarThatHasMoved, double newRangeStart) override;
    
    void sliderValueChanged (juce::Slider* slider) override;
    
    void valueChanged (Value &value) override;

    void changeState (TransportState newState);

    void paint (juce::Graphics& g) override;

    int getMult(int numSamples);
    
    void setCycleColours();
    
    void setButtonColours();
    
    void computeZeros();
    
    void reComputeZeros();
    
    void drawScrollBox(juce::Graphics& g);
    
    void mouseDown (const MouseEvent& event) override
    {
        DBG ("Clicked at: " << event.getPosition().toString());
        std::cout << "lEP: " << lEP.toString() << "  rEP: " << rEP.toString() << std::endl;
    }
    
    void mouseDrag (const MouseEvent& event) override
    {
//        DBG ("Dragged at: " << event.getPosition().toString());
    }
  
    int getCycleNumber(float t);
    
    void openButtonClicked();

    void closeButtonClicked();
    
    void playButtonClicked();
    
    void stopButtonClicked();

    void graphButtonClicked();

    void shadeButtonClicked();
   
    void targetButtonClicked();
    
    void playCycleButtonClicked();
    
    void computeModelButtonClicked();
    
    void playModelButtonClicked();
    
    void readAudioData (File file);

    void readAudioData2 (AudioFormatReader *reader);
    
    void updateAngleDelta();
    
    void setSplineArrays();
    
    void initSplineArrays();
    
    void fillBcoeffs();
    
    void writeWavFile();
    
    void writeModelToBuffer();
    
    struct fileheader
    {
        char riff_label[4]; // (00) = {'R','I','F','F'}
        unsigned riff_size; // (04) = 36 + data_size
        char file_tag[4]; // (08) = {'W','A','V','E'}
        char fmt_label[4]; // (12) = {'f','m','t',' '}
        unsigned fmt_size; // (16) = 16
        unsigned short audio_format; // (20) = 1
        unsigned short channel_count; // (22) = 1 or 2
        unsigned sampling_rate; // (24) = (anything)
        unsigned bytes_per_second; // (28) = (see above)
        unsigned short bytes_per_sample; // (32) = (see above)
        unsigned short bits_per_sample; // (34) = 8 or 16
        char data_label[4]; // (36) = {'d','a','t','a'}
        unsigned data_size; // (40) = # bytes of data
    };
    
    //==========================================================================

    juce::TextButton openButton;
    juce::TextButton closeButton;
    juce::TextButton playButton;
    juce::TextButton stopButton;
    juce::TextButton graphButton;
    juce::TextButton shadeButton;
    juce::TextButton targetButton;
    juce::TextButton selectCycleButton;
    juce::TextButton playCycleButton;
    juce::TextButton computeModelButton;
    juce::TextButton playModelButton;
    juce::ScrollBar signalScrollBar;
    juce::Slider freqGuessSlider;
    juce::Label  freqGuessLabel;
    juce::Slider frequencySlider;
    juce::Label  frequencyLabel;
    juce::Slider kValSlider;
    juce::Label  kValLabel;
    juce::Slider mValSlider;
    juce::Label  mValLabel;
    
    std::unique_ptr<juce::FileChooser> chooser;

    juce::AudioFormatManager formatManager;
    std::unique_ptr<juce::AudioFormatReaderSource> readerSource;
    juce::AudioTransportSource transportSource;
    TransportState state;
    juce::AudioBuffer<float> floatBuffer;
    juce::AudioBuffer<float> writeBuffer;
    juce::AudioFormatReader *reader;
    juce::SoundPlayer player;
    juce::AudioDeviceManager myADM;
    
    unsigned dataSize,          // size of data in bytes
             sampleRate,        // sample rate
             sampleCount;       // number of data samples
    float freqGuess = 220;
    char buffer[44];            // for wav file header
    float w = 0, h = 0;         // width and height of display for MainContentComponent
    float leftEndPoint = 0;     // left endpoint in (float)samples of Interval to graph
    float rightEndPoint = 2000; // right endpoint in (float)samples of Interval to graph
    float addoffset = 0;        // accumulate for left-right Interval shifting
    float magnify = 1;          // to scale the size of Interval
    float magfactor = 1;        // accumulate to give magnification
    int numSamples = 2000;      // for width of (displayed) graph time interval in samples
    int kVal = 20;              // k = number of subintervals for splines
    int mVal = 20;               // m = multiple for key cycles (simple regular model)
    int dVal = 3;               // d = degree for splines, default 3
    // zeros in audio data are computed based on linear interpolation between audio samples
    Array<float> cycleZeros;    // zeros marking endpoints of cycles in audio sample
    Array<float> allZeros;      // all zeros in audio sample
    Array<int> samplesPerCycle; // number of samples in each cycle of audio sample until last zero
    int lastSample;             // index of last sample = (int) lastZero
    int numCycles;              // number of cycles computed, (int) (lengthInSeconds * freqGuess)
    int startIndex = 0;         // index of first cycle to graph
    int endIndex = 0;           // index of last cycle to graph
    Array<juce::Colour> buttonColours;
    juce::Point<int> doubleClick;
    GraphComponent graphView;
    bool playCycleOn = false;
    bool playModelOn = false;
    bool playWavFileOn = false;
    bool expCycleInterp = false;
    float currentSampleRate = 0.0, currentAngle = 0.0, angleDelta = 0.0;
    double currentFrequency = 220.0, targetFrequency = 220.0;
    // set some other arrays and variables to be used in computing cycleToGraph in getNextAudioBlock
    juce::Array<float> controlCoeffs;   // need to set this to size 4*n where n is max number of bcoeffs
    juce::Array<float> knotVals;        // this only depends on k and d
    File outputFile = File("~/output.wav");
    
//    juce::Array<float> controlCoeffs0;  // These are for triple-buffering so that AudioCallBack
//    juce::Array<float> controlCoeffs1;  // can cycle through different sets of bcoeffs
//    juce::Array<float> controlCoeffs2;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainContentComponent)
};

