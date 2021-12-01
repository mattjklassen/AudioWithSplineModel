/*
  ==============================================================================

    GraphComponent.cpp
    Created: 18 Nov 2021 11:51:34am
    Author:  Matt Klassen

  ==============================================================================
*/


#include "GraphComponent.h"
#include "Cycle.h"

void GraphComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colours::lightgrey);
    g.setColour (juce::Colours::blue);
    w = getWidth();
    h = getHeight();
    
    drawGraphBox(g, w, h);
    

    if (updateGraph) {
        graphSignal(g);
    }
    if (audioLoaded) {
        findCyclesToGraph();
    }
    if (callShadeCycles) {
        shadeCycles(g);
    }
}

void GraphComponent::drawGraphBox(juce::Graphics& g, float w, float h)
{
    juce::Path path;
    path.startNewSubPath (juce::Point<float> (1, 1));
    path.lineTo (juce::Point<float> (w-1, 1));
    path.lineTo (juce::Point<float> (w-1, h-1));
    path.lineTo (juce::Point<float> (1, h-1));
    path.lineTo (juce::Point<float> (1, 1));
    // line through middle for t or x axis
    path.startNewSubPath (juce::Point<float> (0, h/2));
    path.lineTo (juce::Point<float> (w, h/2));
    float myThickness = 1.0;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.strokePath (path, myType);
}

juce::Point<float> GraphComponent::signalToScreenCoords (juce::Point<float> P)
{
    // signal Coords are:
    // t ranges from leftEndPoint to rightEndPoint in (float) samples
    // s ranges from -1 to 1 (if s = (float) samples[j] / 32768.0;)
    // screen Coords are:
    // x ranges from 0 to w, y from 0 down to h
    // y = (1 - s) * h/2 or s = 1 - (2/h) * y
    juce::Point<float> Q;
    float xincr = float(w)/float(numSamples);  // one sample in screen coords
    float tval = P.getX() - leftEndPoint;
    float sval = P.getY();
    Q.setX(tval * xincr);
    Q.setY((1 - sval) * h/2);
    return Q;
}

juce::Point<float> GraphComponent::screenToSignalCoords (juce::Point<float> Q)
{
    juce::Point<float> P;
    float y = Q.getY();
    float x = Q.getX();
    float xincr = float(w)/float(numSamples);  // one sample in screen coords
    P.setX(leftEndPoint + x / xincr);          // t
    P.setY(1 - 2/h * y);                       // s
    return P;
}

void GraphComponent::resized() {}


void GraphComponent::setDataForGraph(short * _graphSamples, bool _audioLoaded, int _numSamples, float _magnify, float _leftEndPoint, float _rightEndPoint, unsigned _sampleCount)
{
    graphSamples = _graphSamples;
    audioLoaded = _audioLoaded;
    numSamples = _numSamples;
    magnify = _magnify;
    leftEndPoint = _leftEndPoint;
    rightEndPoint = _rightEndPoint;
    sampleCount = _sampleCount;
}

void GraphComponent::setZerosForGraph(float * _zeros)
{
    zeros = _zeros;
}


void GraphComponent::setCycleColours()
{
    cycleColours.add (juce::Colour::fromHSV (0.5f, 0.5f, 0.7f, 0.2f));
    cycleColours.add (juce::Colour::fromHSV (0.7f, 0.7f, 0.7f, 0.2f));
    cycleColours.add (juce::Colour::fromHSV (0.6f, 0.5f, 0.5f, 0.3f));
    cycleColours.add (juce::Colour::fromHSV (0.8f, 0.5f, 0.7f, 0.2f));
    cycleColours.add (juce::Colour::fromHSV (0.6f, 0.3f, 0.4f, 0.2f));
    highlightColour = juce::Colour::fromHSV (0.3f, 0.5f, 0.8f, 0.2f);
}

void GraphComponent::findCyclesToGraph ()
{
    int i = 0;
    startIndex = 0;
    while (leftEndPoint > zeros[i]) {
        i++;
    }
    if (i > 0) {
        startIndex = i-1;
    }
    while (rightEndPoint > zeros[i]) {
        i++;
        if (i > 300) {
            break;
        }
    }
    endIndex = i-1;
    cyclesToGraph = true;
    std::cout << "cycles to graph:  " << startIndex << " to " << endIndex << std::endl;
}

// to graph signal: leftEndPoint and rightEndPoint determine the graph window
// audio interval is [leftEndPoint, rightEndPoint] = [a,b]
// after scaling by m, [a,b] becomes interval of length (b-a)/m centered at (a+b)/2,
// which is: [(a+b)/2-(b-a)/2m,(a+b)/2+(b-a)/2m]
// graph y=f(x) with x=0 to some number (samples), y=-1 to 1
// in pixels: x goes from 0 to w, and y from 0 (down) to h

void GraphComponent::graphSignal(juce::Graphics& g)
{
    juce::Path graph;
    g.setColour (juce::Colours::black);
    scaleInterval();
    // set Value objects which will be used to set scrollBar
    leftEP.setValue(leftEndPoint);
    rightEP.setValue(rightEndPoint);
    int startSample = (int) leftEndPoint;
    float x = 0;
    float s = (float) graphSamples[startSample] / 32768.0;
    float y = (1 - s) * h/2;
    graph.startNewSubPath (juce::Point<float> (x, y));

    int mult = getMult(numSamples);  // to graph fewer samples as numSamples grows
    float xincr = w / float(numSamples); // convert 1 sample to screen coords
    int j = 1;
    for (int i=1; i < numSamples; i++) {
        if (i * mult > numSamples) {
            break;
        }
        x = x + mult * xincr;
        j = startSample + i * mult;
        s = (float) graphSamples[j] / 32768.0;
        y = (1 - s) * h/2;
        graph.lineTo (juce::Point<float> (x,y));
        g.setColour (juce::Colours::green);
        if (numSamples < 400) {
            drawDot(juce::Point<float> (x,y), g);
        }
    }
    float myThickness = 1;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.setColour (juce::Colours::darkgrey);
    g.strokePath (graph, myType);
}

int GraphComponent::getMult(int numSamples) {
    int mult = 1;
    if (numSamples > 4000) {
        mult = 2;
        if (numSamples > 8000) {
            mult = 4;
            if (numSamples > 16000) {
                mult = 8;
                if (numSamples > 32000) {
                    mult = 16;
                    if (numSamples > 64000) {
                        mult = 32;
                        if (numSamples > 128000) {
                            mult = 64;
                        }
                    }
                }
            }
        }
    }
    return mult;
}

void GraphComponent::drawDot(juce::Point<float> (P), juce::Graphics& g)
{
    g.fillEllipse (P.getX()-1.5, P.getY()-1.5, 3, 3);
}

// Need to process cycles in view one by one, to shade, then also to select
// We should have a cycle object which has the following data and methods:
// float a, b - interval endpoints of cycle
// int cycle number
// color can be based on cycle number
// shade method, takes graphics context leftEP and rightEP

void GraphComponent::shadeCycle(int n, juce::Graphics& g)
{
    float xincr = float(w) / float(numSamples);  // one sample in screen coords
    float tval1 = zeros[n] - leftEndPoint;
    juce::Point<float> P (1, 1);  // upper left pt of rectangle to shade
    if (tval1 > 0) {
        P.setX(tval1 * xincr);
    }
    float tval2 = rightEndPoint - zeros[n+1];
    juce::Point<float> Q (w - 1, h - 1);  // lower right pt of rectangle to shade
    if (tval2 > 0) {
        Q.setX(w + 1 - tval2 * xincr);
    }
    juce::Rectangle<float> toShade (P, Q);
    g.setColour (cycleColours[n % 5]);
    if (n == highlightCycle) {
        g.setColour(highlightColour);
    }
    g.fillRect(toShade);
    g.setFont (10.0f);
    float midpt = (zeros[n] + zeros[n+1])/2 - leftEndPoint;
    float x = midpt * xincr;
    auto cycleNumber = std::to_string(n);
    juce::Point<float> ptL (x-20, 100);
    juce::Point<float> ptR (x+20, 120);
    juce::Rectangle<float> Rect (ptL, ptR);
    g.setColour(juce::Colours::darkgrey);
    g.drawText (cycleNumber, Rect, juce::Justification::centred, true);
}

void GraphComponent::shadeCycles(juce::Graphics& g)
{
    for (int i = startIndex; i < endIndex + 1; i++)
    {
        shadeCycle(i, g);
    }
}

// trackpad 2-finger gesture is mouseWheelMove in JUCE
void GraphComponent::mouseWheelMove (const MouseEvent& event, const MouseWheelDetails & wheel)
{
    addoffset -= wheel.deltaX;
    if ((addoffset > 0.1) || (addoffset < -0.1))
    {
        if (leftEndPoint + 100 * addoffset < 0) {
            addoffset = 0;
            hardLeft = true;
            return;
        }
        hardLeft = false;
        leftEndPoint += 100 * addoffset;
        rightEndPoint += 100 * addoffset;
        addoffset = 0;
        if (leftEndPoint < 0) {
            leftEndPoint = 0;
        }
        if (rightEndPoint > sampleCount) {
            rightEndPoint = sampleCount;
            leftEndPoint = rightEndPoint - (float) numSamples;
        }
        repaint();
    }
}

void GraphComponent::mouseMagnify (const MouseEvent & event, float scaleFactor)
{
    magfactor *= scaleFactor;
    if ((magfactor > 1.1) || (magfactor < 0.9)) {
//            DBG ("magfactor is: " << magfactor);
        magnify = magfactor;
        repaint();
        magfactor = 1;
    }
}

void GraphComponent::scaleInterval()
{
    if ((magnify > 1.1) || (magnify < 0.9)) {
        // change interval [a,b] to [(a+b)/2-(b-a)/2m,(a+b)/2+(b-a)/2m]
        float a = leftEndPoint, b = rightEndPoint, m = magnify;
//            std::cout "a: " << a << " b: " << b << " m: " << m << std::endl;
        if (a>0) {
            leftEndPoint = (a+b)/2 - (b-a)/(2*m);
            rightEndPoint = (a+b)/2 + (b-a)/(2*m);
        } else {
            rightEndPoint = b / m;
        }
        if (leftEndPoint < 0) {
            leftEndPoint = 0;
        }
        if (rightEndPoint > sampleCount) {
            rightEndPoint = sampleCount;
            leftEndPoint = rightEndPoint - (float) numSamples;
        }
        numSamples = (int) rightEndPoint - leftEndPoint;
        magnify = 1;
    }
}

int GraphComponent::getCycleNumber(float t)
{
    int i = 0;
    while (zeros[i] < t) {
        i++;
    }
    return i-1;
}
            
void GraphComponent::mouseDoubleClick (const MouseEvent& event)
{
    juce::Point<int> P = event.getPosition();
//    DBG ("Double-Clicked at: " << P.toString());
    juce::Point<float> Q (P.getX(),P.getY());
    juce::Point<float> S = screenToSignalCoords(Q);
    int n = getCycleNumber(S.getX());
    DBG("selected cycle: " << n);
    if (highlightCycle == n) {
        highlightCycle = -1;
    } else {
        highlightCycle = n;
    }
    repaint();
}

void GraphComponent::mouseDown (const MouseEvent& event)
{
//    DBG ("Clicked at: " << event.getPosition().toString());
}