/*
  ==============================================================================

    Sliders.cpp
    Created: 8 Jan 2022 3:14:33pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MainContentComponent.h"

void MainContentComponent::addSliders()
{
    
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
    kValSlider.setValue (kVal);
    kValSlider.setSkewFactorFromMidPoint (40);
    kValSlider.setTextBoxStyle (juce::Slider::TextBoxLeft, false, 40, kValSlider.getTextBoxHeight());
    
    addAndMakeVisible (&kValLabel);
    kValLabel.setText ("k = # subintervals", juce::dontSendNotification);
    kValLabel.attachToComponent (&kValSlider, true);
    
    addAndMakeVisible (&mValSlider);
    mValSlider.setRange (1, 100);
    mValSlider.setNumDecimalPlacesToDisplay(0);
    mValSlider.addListener (this);
    mValSlider.setValue (mVal);
    mValSlider.setSkewFactorFromMidPoint (40);
    mValSlider.setTextBoxStyle (juce::Slider::TextBoxLeft, false, 40, mValSlider.getTextBoxHeight());
    
    addAndMakeVisible (&mValLabel);
    mValLabel.setText ("m = key cycle mult", juce::dontSendNotification);
    mValLabel.attachToComponent (&mValSlider, true);
  
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
        findMaxValuesPerCycle(maxSampleIndices, maxSampleValues, cycleZeros, floatBuffer);
        graphView.setZerosForGraph(cycleZeros, allZeros, samplesPerCycle, freqGuess, maxSampleIndices, maxSampleValues);
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
        if (mVal == 1) {
            noCycleInterp = true;
        }
        DBG("new mVal:  " << mVal);
        graphView.setmVal(mVal);
        graphView.repaint();
    }
}
