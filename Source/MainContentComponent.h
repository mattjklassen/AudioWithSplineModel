/*
  ==============================================================================

    MainContentComponent.h
    Created: 29 Nov 2021 11:39:48am
    Author:  Matt Klassen

  ==============================================================================
*/



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
                               public juce::Value::Listener
{
public:
    
    MainContentComponent();

    ~MainContentComponent() override
    {
        shutdownAudio();
    }
    
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override
    {
        transportSource.prepareToPlay (samplesPerBlockExpected, sampleRate);
    }

    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override;

    void releaseResources() override
    {
        transportSource.releaseResources();
    }

    void resized() override;

    void changeListenerCallback (juce::ChangeBroadcaster* source) override;
    
    bool audioDataLoaded = false;

    Value lEP;  // left End Point for graphing interval in graphView
    Value rEP;  // right End Point for graphing interval in graphView
 
private:
    
    enum TransportState
    {
        Stopped, Starting, Playing, Pausing, Paused, Stopping
    };

    void scrollBarMoved (ScrollBar* scrollBarThatHasMoved, double newRangeStart) override;
    
    void valueChanged (Value &value) override;

    void changeState (TransportState newState);

    void paint (juce::Graphics& g) override;

    int getMult(int numSamples);
    
    void setCycleColours();
    
    void setButtonColours();
    
    void computeZeros();
    
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

    void playButtonClicked();
    
    void stopButtonClicked();

    void graphButtonClicked();

    void shadeButtonClicked();
   
    void readAudioData (File file);
    
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
    juce::TextButton playButton;
    juce::TextButton stopButton;
    juce::TextButton graphButton;
    juce::TextButton shadeButton;
    juce::TextButton selectCycleButton;
    juce::ScrollBar signalScrollBar;
    
    std::unique_ptr<juce::FileChooser> chooser;

//! [members]
    juce::AudioFormatManager formatManager;
    std::unique_ptr<juce::AudioFormatReaderSource> readerSource;
    juce::AudioTransportSource transportSource;
    TransportState state;
    
    unsigned dataSize,              // size of data in bytes
             sampleRate,              // sample rate
             sampleCount;             // number of data samples
    char buffer[44];            // for wav file header
    char *data = 0;             // pointer to array of wav file data in bytes
    short *samples = 0;         // pointer to array of wav file data in samples
//    bool updateGraph = false;
    float w = 0, h = 0;
    float leftEndPoint = 0;     // left endpoint (in (float)samples) of Interval to graph
    float rightEndPoint = 2000; // right endpoint (in (float)samples) of Interval to graph
    float addoffset = 0;        // accumulate for left-right Interval shifting
    float magnify = 1;          // to scale the size of Interval
    float magfactor = 1;        // accumulate to give magnification
    int numSamples = 2000;      // for width of (displayed) graph time interval in samples
//    bool callShadeCycles = false;
//    bool cyclesToGraph = false;
    Array<float> cycleZeros;    // dynamically resizable array
    int startIndex = 0;
    int endIndex = 0;
//    Array<juce::Colour> cycleColours;
    Array<juce::Colour> buttonColours;
    juce::Point<int> doubleClick;
    GraphComponent graphView;
//    std::unique_ptr<juce::ScrollBar> scrollBar;
//! [members]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainContentComponent)
};

