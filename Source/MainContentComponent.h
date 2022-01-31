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
                    
    bool keyStateChanged(bool isKeyDown) override;
    
    void resized() override;

    void changeListenerCallback (juce::ChangeBroadcaster* source) override;
    
    float computeSpline(int control, float t);
    
    void randBcoeff(int n);
    
    bool audioDataLoaded = false;

    Value lEP;  // left End Point for graphing interval in graphView
    Value rEP;  // right End Point for graphing interval in graphView
    Value cycleRendered;  // set to i when cycle[i] is rendered by audio callback

    
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
//        DBG ("Clicked at: " << event.getPosition().toString());
//        std::cout << "lEP: " << lEP.toString() << "  rEP: " << rEP.toString() << std::endl;

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
    
    void nextRandomButtonClicked();
    
    void playCycleWithEnvButtonClicked();
    
    void graphMetaSplinesButtonClicked();
    
    void normalizeCycleLengthButtonClicked();
    
    void randomizeCycleBcoeffs(int cycleRendered);
    
    void interpSelectionChanged();
    
    void graphModelButtonClicked();
    
    void readAudioData (File file);

    void readAudioData2 (AudioFormatReader *reader);
    
    void updateAngleDelta();
    
    void setSplineArrays();

    void fillBcoeffs();
    
    void writeWavFile();
    
    void writeModelToBuffer();
    
    void writeCycleInterpModelToBuffer();
    
    void modelWithoutCycleInterp();
    
    void writeKeyCyclesToBuffer();
    
    void writeCycleWithEnvToBuffer();
    
    void computeCycleWithEnv();
    
    void computeMetaSplines();
    
    void setKeyCycleIndices();
    
    void computeKeyCycles();
    
    void writeNonKeyCyclesToBuffer();
    
    void writeNormalizedCyclesToBuffer();
    
    void computeAllCycles();
    
    bool isKey(int i);
    
    int keyIndex(int i);
        
    void findMaxValuesPerCycle(Array<float>&  maxSampleIndices, Array<float>&  maxSampleValues, Array<float>& cycleZeros, AudioBuffer<float>& samples);
   
    void addButtons();
    
    void addSliders();
    
    void addScrollbar();
    
    void randomizeButtonClicked();
    
    void setInterpSelectionsFalse();
    
    void setParabolicTargets(float scale);
    
    void computeParabolicBcoeffs();
    
    CycleSpline cycleParabola = CycleSpline(20, 0, 1);
    
    void setNew(CycleSpline& cycle);
    
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
    juce::TextButton graphModelButton;
    juce::TextButton expCycleInterpButton;
    juce::TextButton regCycleInterpButton;
    juce::TextButton shadeButton;
    juce::TextButton targetButton;
    juce::TextButton selectCycleButton;
    juce::TextButton playCycleButton;
    juce::TextButton computeModelButton;
    juce::TextButton playModelButton;
    juce::TextButton playCycleWithEnvButton;
    juce::TextButton graphMetaSplinesButton;
    juce::ToggleButton normalizeCycleLengthButton;
    juce::ToggleButton randomizeButton;
    juce::ScrollBar signalScrollBar;
    juce::Slider freqGuessSlider;
    juce::Label  freqGuessLabel;
    juce::Slider frequencySlider;
    juce::Label  frequencyLabel;
    juce::Slider kValSlider;
    juce::Label  kValLabel;
    juce::Slider mValSlider;
    juce::Slider amplitudeSlider;
    juce::Label  mValLabel;
    juce::ComboBox interpSelector;
    juce::Random random;
    juce::TextButton nextRandomButton;
    
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
    float freqGuess = 220;      // guess to determine cycle lengths
    float samplesPerCycleGuess = 200.5;  // will calculate based on sampleRate/freqGuess
    char buffer[44];            // for wav file header
    float w = 0, h = 0;         // width and height of display for MainContentComponent
    float leftEndPoint = 0;     // left endpoint in (float)samples of Interval to graph
    float rightEndPoint = 2000; // right endpoint in (float)samples of Interval to graph
    float addoffset = 0;        // accumulate for left-right Interval shifting
    float magnify = 1;          // to scale the size of Interval
    float magfactor = 1;        // accumulate to give magnification
    int numSamples = 2000;      // for width of (displayed) graph time interval in samples
    int kVal = 50;              // k = number of subintervals for splines
    int mVal = 20;               // m = multiple for key cycles (simple regular model)
    int dVal = 3;                // d = degree for splines, default 3
    int sampleRendered = 0;
    int samplesPerSelectedCycle = 0;
    float amplitudeFactor = 1;
    // zeros in audio data are computed based on linear interpolation between audio samples
    Array<float> cycleZeros;    // zeros marking endpoints of cycles in audio sample
    Array<float> normalizedCycleZeros;  // make all cycles same length
    Array<float> allZeros;      // all zeros in audio sample
    Array<float> keyBcoeffs;
    Array<int> samplesPerCycle; // number of samples in each cycle of audio sample until last zero
    Array<int> keys;
    Array<MetaSpline> metaSplineArray;
    Array<CycleSpline> keyCycleArray;
    Array<float> maxSampleIndices;  // per cycle, sample index for max abs value
    Array<float> maxSampleValues;   // per cycle, abs max value at sample
    int lastSample;             // index of last sample = (int) lastZero
    int numCycles;              // number of cycles computed, (int) (lengthInSeconds * freqGuess)
    int startIndex = 0;         // index of first cycle to graph
    int endIndex = 0;           // index of last cycle to graph
    int control = 0;            // for not randomizing cycles
    int cycleInUse = 0;         // in use by audio callback
    int bcoeffRandomized = 0;   // changed by audio callback
    Array<juce::Colour> buttonColours;
    juce::Point<int> doubleClick;
    GraphComponent graphView;
    bool playCycleOn = false;
    bool playModelOn = false;
    bool playWavFileOn = false;
    bool expCycleInterp = false;
    bool fibonacciCycleInterp = false;
    bool noCycleInterp = false;
    bool otherCycleInterp = false;
    bool regularCycleInterp = true;
    bool endsOnlyCycleInterp = false;
    bool normalizeCycleLength = false;
    bool addSubharmonic = false;
    bool useEnvelope = true;
    bool randomizeBcoeffs = false;
    bool keysWithoutCycleInterp = false;
    bool oneCycleWithEnvelope = false;
    float currentSampleRate = 0.0, currentAngle = 0.0, angleDelta = 0.0;
    double currentFrequency = 220.0, targetFrequency = 220.0;
    // set some other arrays and variables to be used in computing cycleToGraph in getNextAudioBlock
    juce::Array<float> controlCoeffs;   // need to set this to size 4*n where n is max number of bcoeffs
    juce::Array<float> knotVals;        // this only depends on k and d
    File outputFile = File::getSpecialLocation(File::userHomeDirectory).getChildFile("output.wav");
    
//    juce::Array<float> controlCoeffs0;  // These are for triple-buffering so that AudioCallBack
//    juce::Array<float> controlCoeffs1;  // can cycle through different sets of bcoeffs
//    juce::Array<float> controlCoeffs2;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainContentComponent)
};

