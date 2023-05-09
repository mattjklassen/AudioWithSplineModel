/*
  ==============================================================================

    ScrollBars.cpp
    Created: 8 Jan 2022 3:40:05pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "MainContentComponent.h"

void MainContentComponent::addScrollbar()
{
    addAndMakeVisible (&signalScrollBar);
    signalScrollBar.setVisible(true);
    signalScrollBar.setRangeLimits(0, 1000, sendNotificationAsync);
    signalScrollBar.addListener(this);
}

void MainContentComponent::drawScrollBox(juce::Graphics& g)
{
    juce::Path path;
    path.startNewSubPath (juce::Point<float> (10, h-40));
    path.lineTo (juce::Point<float> (w-31, h-40));
    path.lineTo (juce::Point<float> (w-31, h-10));
    path.lineTo (juce::Point<float> (10, h-10));
    path.lineTo (juce::Point<float> (10, h-40));
    float myThickness = 1.0;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.strokePath (path, myType);
}

// don't allow graph window to have rightEndPoint past end of last full cycle
void MainContentComponent::scrollBarMoved (ScrollBar* scrollBarThatHasMoved ,
                                 double newRangeStart)
{
    int n = graphView.cycleZeros.size();
    double newRangeL = scrollBarThatHasMoved->getCurrentRange().getStart();
    graphView.leftEndPoint = (float) newRangeL / 1000 * (float) graphView.sampleCount;
    graphView.rightEndPoint = graphView.leftEndPoint + (float) graphView.numSamples;
    if (graphView.rightEndPoint > graphView.cycleZeros[n-1]) {
        graphView.rightEndPoint = (float) graphView.cycleZeros[n-1];
        graphView.leftEndPoint = graphView.rightEndPoint - (float) graphView.numSamples;
        graphView.hardRight = true;
        DBG("set hardRight true");
    }
    if (newRangeL == 0) {
        graphView.hardLeft = true;
//        DBG("set hardLeft true");
    }
    graphView.repaint();
}

