/*
  ==============================================================================

    GraphComponent.cpp
    Created: 18 Nov 2021 11:51:34am
    Author:  Matt Klassen

  ==============================================================================
*/


#include "GraphComponent.h"
#include "CycleSpline.h"
#include "zeros.h"
#include "bsp.h"

void GraphComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colours::lightgrey);
    w = getWidth();
    h = getHeight();
        
    g.setColour (juce::Colours::darkgrey);
    drawGraphBox(g, w, h);
    g.setColour (juce::Colours::blue);
    if (graphMetaSplinesOn) {
//        graphMetaSplines(g);
        graphLinearMetaSplines(g);
        if (plotTargets) {
            plotMetaSplineTargets(g);
        }
        return;
    }
    if (updateGraph) {
        graphSignal(g);
    }
    if (audioLoaded) {
        findCyclesToGraph();
    }
    if (callShadeCycles) {
        shadeCycles(g);
    }
    if (graphSplineCycle) {
        graphSpline(g, cycleToGraph);
    }
    if (graphNewSplineCycle) {
        graphSpline(g, cycleNew);
    }
    if (graphParabolicSpline) {
        graphSpline(g, cycleParabola);
    }
    if (plotTargets) {
        plotTargetPoints(g, cycleToGraph);
    }
    if (mouseOverCycleZero) {
        g.setColour (juce::Colours::yellow);
        drawLargeDot(zeroPoint, g);
    }
}

void GraphComponent::reset()
{
    leftEndPoint = 0;
    rightEndPoint = 1000;
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
    float xincr = w / (rightEndPoint - leftEndPoint);  // one sample in screen coords
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
    float xincr = w / (rightEndPoint - leftEndPoint);  // one sample in screen coords
    P.setX(leftEndPoint + x / xincr);          // t
    P.setY(1 - 2/h * y);                       // s
    return P;
}

void GraphComponent::resized() {}


void GraphComponent::setDataForGraph(AudioBuffer<float>& _floatBuffer, bool _audioLoaded, int _numSamples, float _magnify, float _leftEndPoint, float _rightEndPoint, unsigned _sampleCount, unsigned _sampleRate, int _kVal)
{
    floatBuffer = _floatBuffer;
    audioLoaded = _audioLoaded;
    numSamples = _numSamples;
    magnify = _magnify;
    leftEndPoint = _leftEndPoint;
    rightEndPoint = _rightEndPoint;
    sampleCount = _sampleCount;
    sampleRate = _sampleRate;
    kVal = _kVal;
}

void GraphComponent::setAmplitudeFactor(float _amplitudeFactor)
{
    amplitudeFactor = _amplitudeFactor;
}

void GraphComponent::setModelForGraph(AudioBuffer<float>& _modelBuffer)
{
    modelBuffer = _modelBuffer;
}

void GraphComponent::setMetaSplinesForGraph(Array<MetaSpline>& _metaSplineArray)
{
    metaSplineArray = _metaSplineArray;
}

void GraphComponent::setkVal(int _kVal)
{
    kVal = _kVal;
}

void GraphComponent::setmVal(int _mVal)
{
    mVal = _mVal;
}

void GraphComponent::setKeys(Array<int>& _keys)
{
    keys = _keys;
}

// returns true if i is the value of some key cycle
bool GraphComponent::isKey(int i)
{
    for (int j=0; j<keys.size(); j++) {
        if (keys[j] == i) {
            return true;
        }
    }
    return false;
}

void GraphComponent::setZerosForGraph(Array<float>& _cycleZeros, Array<float>& _allZeros, Array<int> _samplesPerCycle, float _freqGuess, Array<float>& _maxSampleIndices, Array<float>& _maxSampleValues)
{
    cycleZeros = _cycleZeros;
    allZeros = _allZeros;
    samplesPerCycle = _samplesPerCycle;
    freqGuess = _freqGuess;
    maxSampleIndices = _maxSampleIndices;
    maxSampleValues = _maxSampleValues;
    
//    DBG("zeros are set");
//    for (int i = 0; i<10; i++) {
//        DBG("allZeros[" << i << "]:  " << allZeros[i]);
//    }
//    for (int i = 0; i<10; i++) {
//        DBG("cycleZeros[" << i << "]:  " << cycleZeros[i]);
//    }
//    DBG("number of cycles:  " << cycleZeros.size());
//    DBG("number of zeros:  " << allZeros.size());
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
    while (leftEndPoint > cycleZeros[i]) {
        i++;
    }
    if (i > 0) {
        startIndex = i-1;
    }
    while (rightEndPoint > cycleZeros[i]) {
        i++;
        if (i > cycleZeros.size()-1) {
            break;
        }
    }
    endIndex = i-1;
    cyclesToGraph = true;
//    std::cout << "cycles to graph:  " << startIndex << " to " << endIndex << std::endl;
}

// to graph signal: leftEndPoint and rightEndPoint determine the graph window
// audio interval is [leftEndPoint, rightEndPoint] = [a,b]
// after scaling by m, [a,b] becomes interval of length (b-a)/m centered at (a+b)/2,
// which is: [(a+b)/2-(b-a)/2m,(a+b)/2+(b-a)/2m]
// graph y=f(x) with x=0 to some number (samples), y=-1 to 1
// in pixels: x goes from 0 to w, and y from 0 (down) to h

void GraphComponent::graphSignal(juce::Graphics& g)
{
    juce::Path signalGraph;
    juce::Path modelGraph;
    scaleInterval();
    // set Value objects which will be used to set scrollBar
    leftEP.setValue(leftEndPoint);
    rightEP.setValue(rightEndPoint);
    int startSample = (int) leftEndPoint;
    float x = 0;
    float s = floatBuffer.getSample(0, startSample);
    float y = (1 - s) * h/2;
    signalGraph.startNewSubPath (juce::Point<float> (x, y));
    if (updateModelGraph) {
        modelGraph.startNewSubPath (juce::Point<float> (x, y));
    }
    int mult = getMult(numSamples);  // to graph fewer samples as numSamples grows
    float xincr = w / (rightEndPoint - leftEndPoint); // convert 1 sample to screen coords
    int j = 1;
    for (int i=1; i < numSamples; i++) {
        if (i * mult > numSamples) {
            break;
        }
        x = x + mult * xincr;
        j = startSample + i * mult;
        if (j > sampleCount-1) {
            j = sampleCount-1;
            i = numSamples;
        }
        s = floatBuffer.getSample(0, j);
        s *= amplitudeFactor;
        y = (1 - s) * h/2;
        signalGraph.lineTo (juce::Point<float> (x,y));
        g.setColour (juce::Colours::green);
        if (numSamples < 400) {
            drawDot(juce::Point<float> (x,y), g);
        }
        if (updateModelGraph) {
            s = modelBuffer.getSample(0, j);
            s *= amplitudeFactor;
            y = (1 - s) * h/2;
            modelGraph.lineTo (juce::Point<float> (x,y));
        }
    }
    float myThickness = 1;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.setColour (juce::Colours::darkgrey);
    g.strokePath (signalGraph, myType);
    if (updateModelGraph) {
        g.setColour (juce::Colours::blue);
        g.strokePath (modelGraph, myType);
    }
}

void GraphComponent::graphLinearMetaSplines(juce::Graphics& g)
{
    MetaSpline spline = metaSplineArray[metaSplineIndex];
    float incr = 1 / (float)spline.outputs.size();
    leftEndPoint = 0;
    rightEndPoint = 1;

    // graph meta-spline from 0 to 1 using piecewise linear path between outputs
    juce::Path graph;
    juce::Point<float> P_0(0,spline.outputs[0]);
    P_0 = signalToScreenCoords(P_0);
    graph.startNewSubPath (P_0);
    float x, y;
    for (int i=1; i<spline.outputs.size(); i++) {
        x = (float)i * incr;
        y = spline.outputs[i];
        y *= amplitudeFactor;
        juce::Point<float> P(x,y);
        P = signalToScreenCoords(P);
        graph.lineTo(P);
    }
//    juce::Point<float> P_1(1,spline.value(1));
//    P_1 = signalToScreenCoords(P_1);
//    graph.lineTo(P_1);
    float myThickness = 1;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.setColour (juce::Colours::blue);
    g.strokePath (graph, myType);
}

void GraphComponent::graphMetaSplines(juce::Graphics& g)
{
    MetaSpline spline = metaSplineArray[metaSplineIndex];
//    int n = spline.n;
    float incr = 0.001;
    leftEndPoint = 0;
    rightEndPoint = 1;
    // find max value of spline to normalize
    float max = 0;
    for (int j=0; j<1001; j++) {
        float val = abs(spline.value(j * incr));
        if (val > max) {
            max = val;
        }
    }
    float scale = 1;
    if (max > 1) {
        scale = 1 / max;
    }
    metaSplineFactor = scale;
    // graph meta-spline from 0 to 1 using spline.value
    juce::Path graph;
    juce::Point<float> P_0(0,spline.value(0));
    P_0 = signalToScreenCoords(P_0);
    graph.startNewSubPath (P_0);
    float x, y;
    for (int i=0; i<1001; i++) {
        x = (float)i * incr;
        y = scale * spline.value(x);
        juce::Point<float> P(x,y);
        P = signalToScreenCoords(P);
        graph.lineTo(P);
    }
    juce::Point<float> P_1(1,spline.value(1));
    P_1 = signalToScreenCoords(P_1);
    graph.lineTo(P_1);
    float myThickness = 1;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.setColour (juce::Colours::blue);
    g.strokePath (graph, myType);
}

void GraphComponent::graphSpline (juce::Graphics& g, CycleSpline& cycle)
{
    // graph cycle spline from a to b using outputs at each sample
    juce::Path graph;
    juce::Point<float> P_a(cycle.a,0);
    P_a = signalToScreenCoords(P_a);
    graph.startNewSubPath (P_a);
    float x, y;
    int I = (int) cycle.a + 1;
    int J = (int) cycle.b + 1;
    for (int i=0; i<J-I; i++) {
        x = (float) i+I;
        y = cycle.outputs[i];
        y *= amplitudeFactor;
        juce::Point<float> P(x,y);
        P = signalToScreenCoords(P);
        graph.lineTo(P);
    }
    juce::Point<float> P_b(cycle.b,0);
    P_b = signalToScreenCoords(P_b);
    graph.lineTo(P_b);
    float myThickness = 1;
    juce::PathStrokeType myType = PathStrokeType(myThickness);
    g.setColour (juce::Colours::blue);
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

void GraphComponent::drawLargeDot(juce::Point<float> (P), juce::Graphics& g)
{
    g.fillEllipse (P.getX()-2.5, P.getY()-2.5, 5, 5);
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
    float tval1 = cycleZeros[n] - leftEndPoint;
    juce::Point<float> P (1, 1);  // upper left pt of rectangle to shade
    if (tval1 > 0) {
        P.setX(tval1 * xincr);
    }
    float tval2 = rightEndPoint - cycleZeros[n+1];
    juce::Point<float> Q (w - 1, h - 1);  // lower right pt of rectangle to shade
    if (tval2 > 0) {
        Q.setX(w + 1 - tval2 * xincr);
    }
    juce::Rectangle<float> toShade (P, Q);
    g.setColour (cycleColours[n % 5]);
    if (isKey(n)) {
        g.setColour(highlightColour);
    }
    if (n == highlightCycle) {
        g.setColour(highlightColour);
    }
    g.fillRect(toShade);
    g.setFont (10.0f);
    float midpt = (cycleZeros[n] + cycleZeros[n+1])/2 - leftEndPoint;
    float x = midpt * xincr;
    auto cycleNumber = std::to_string(n);
    juce::Point<float> ptL (x-20, 40);
    juce::Point<float> ptR (x+20, 60);
    juce::Rectangle<float> Rect (ptL, ptR);
    g.setColour(juce::Colours::darkgrey);
    g.drawText (cycleNumber, Rect, juce::Justification::centred, true);
    auto samplesPC = std::to_string(samplesPerCycle[n]);
    juce::Point<float> ptL2 (x-20, h-50);
    juce::Point<float> ptR2 (x+20, h-30);
    juce::Rectangle<float> Rect2 (ptL2, ptR2);
    g.setColour(juce::Colours::darkgrey);
    g.drawText (samplesPC, Rect2, juce::Justification::centred, true);
    juce::Point<float> cycleZeroLeft (cycleZeros[n], 0);
    cycleZeroLeft = signalToScreenCoords(cycleZeroLeft);
    g.setColour (juce::Colours::darkgrey);
    drawLargeDot(cycleZeroLeft, g);
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
        if (rightEndPoint + 100 * addoffset > sampleCount-1) {
            addoffset = 0;
            hardRight = true;
            return;
        }
        hardLeft = false;
        hardRight = false;
        leftEndPoint += 100 * addoffset;
        rightEndPoint += 100 * addoffset;
        addoffset = 0;
        if (leftEndPoint < 0) {
            leftEndPoint = 0;
        }
        if (rightEndPoint > sampleCount-1) {
            rightEndPoint = sampleCount-1;
//            leftEndPoint = rightEndPoint - (float) numSamples;
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
    while (cycleZeros[i] < t) {
        i++;
    }
    return i-1;
}
            
void GraphComponent::mouseDoubleClick (const MouseEvent& event)
{
    if (! audioLoaded) {
        return;
    }
    juce::Point<int> P = event.getPosition();
    juce::Point<float> Q (P.getX(),P.getY());
    juce::Point<float> S = screenToSignalCoords(Q);
    int n = getCycleNumber(S.getX());
    if (highlightCycle == n) {
        highlightCycle = -1;
        graphSplineCycle = false;
        plotTargets = false;
    } else {
        highlightCycle = n;
    }
    if (highlightCycle == n) {
        float a = cycleZeros[n];
        float b = cycleZeros[n+1];
        cycleToGraph = CycleSpline(kVal, a, b);
//        DBG("selected cycle: " << n << "  a = " << a << "  b = " << b << " k = " << kVal);
        computeCycleBcoeffs(cycleToGraph, floatBuffer);
        computeCycleSplineOutputs(cycleToGraph);
//        cycleToGraph.printData();
        graphSplineCycle = true;

        // How many cycles do we need in this array in order to render one chunk
        // while modifying others?  One call to audio callback could require 5 or
        // more of these cycles if they are around 100 samples each.  So maybe 10?
        
        cycles[0] = cycleToGraph;
    }
    repaint();
}


void GraphComponent::plotMetaSplineTargets(juce::Graphics& g)
{
    MetaSpline spline = metaSplineArray[metaSplineIndex];
    leftEndPoint = 0;
    rightEndPoint = 1;
    juce::Point<float> P;
    juce::Point<float> Q;
    g.setColour (juce::Colours::red);
    int n = spline.n;
    for (int i=0; i<n; i++) {
        P.setX(spline.inputs[i]);
        P.setY(spline.targets[i] * metaSplineFactor);
        Q = signalToScreenCoords(P);
        drawDot(Q, g);
    }
    repaint();
}

void GraphComponent::plotTargetPoints(juce::Graphics& g, CycleSpline& cycle)
{
    juce::Point<float> P;
    juce::Point<float> Q;
    g.setColour (juce::Colours::red);
    int n = cycle.inputs.size();
    for (int i=0; i<n; i++) {
        P.setX(cycle.inputs[i]);
        float y = cycle.targets[i];
        y *= amplitudeFactor;
        P.setY(y);
        Q = signalToScreenCoords(P);
        drawDot(Q, g);
    }
    repaint();
}

void GraphComponent::setParabolicTargets(float scale)
{
    int k = cycleParabola.k;
    int d = cycleParabola.d;
    int n = k + d;
    int N = n + d;
    cycleParabola.inputs.clear();
    cycleParabola.targets.clear();
    cycleParabola.knots.clear();
    float incr = 1 / (float) k;
    cycleParabola.inputs.add(0);
    cycleParabola.inputs.add(0.5*incr);
    for (int i=1; i<k; i++) {
        cycleParabola.inputs.add(i*incr);
    }
    cycleParabola.inputs.add(1-0.5*incr);
    cycleParabola.inputs.add(1);
    for (int i=0; i<cycleParabola.inputs.size(); i++) {
        float t = cycleParabola.inputs[i];
        cycleParabola.targets.set(i, scale*4*t*(1-t));
    }
    cycleParabola.knots.add(-d*incr);  // knot t_0
    for (int i=1; i<N+1; i++) {
        cycleParabola.knots.add(cycleParabola.knots[i-1]+incr);
    }
}

void GraphComponent::mouseDown (const MouseEvent& event)
{
    // need to capture mouse hover when near a cycle zero on the timeline and highlight zero
    // then click and drag should move cycle zero to next zero (AllZeros) left or right

    if (!mouseOverCycleZero)
    {
        PopupMenu m;
        m.addItem (1, "add key cycle");
        m.addItem (2, "remove key cycle");
        m.addItem (3, "graph cycle spline");
        m.addItem (4, "graph new cycle spline");
        
        const int result = m.show();
        juce::Point<int> P = event.getPosition();
        juce::Point<float> Q (P.getX(),P.getY());
        juce::Point<float> S = screenToSignalCoords(Q);
        int n = getCycleNumber(S.getX());
        DBG("selected cycle: " << n);
        
        if (result == 0)
        {
            // user dismissed the menu without picking anything
        }
        else if (result == 1)
        {
            keysToAdd.add(n);
            DBG("adding key cycle " << n);
        }
        else if (result == 2)
        {
            keysToRemove.add(n);
            DBG("removing key cycle " << n);
        }
        else if (result == 3)
        {
            if (! audioLoaded) {
                return;
            }
            DBG("graphing cycle spline " << n);
            juce::Point<int> P = event.getPosition();
            juce::Point<float> Q (P.getX(),P.getY());
            juce::Point<float> S = screenToSignalCoords(Q);
            int n = getCycleNumber(S.getX());
            if (highlightCycle == n) {
                highlightCycle = -1;
                graphSplineCycle = false;
                plotTargets = false;
            } else {
                highlightCycle = n;
            }
            if (highlightCycle == n) {
                float a = cycleZeros[n];
                float b = cycleZeros[n+1];
                cycleToGraph = CycleSpline(kVal, a, b);
                computeCycleBcoeffs(cycleToGraph, floatBuffer);
                computeCycleSplineOutputs(cycleToGraph);
                graphSplineCycle = true;
                graphNewSplineCycle = false;
                cycleToGraph.printData();
                cycles[0] = cycleToGraph;
            }
            repaint();
        }
        else if (result == 4) {
            float a = cycleZeros[n];
            float b = cycleZeros[n+1];
//            cycleParabola = CycleSpline(kVal, a, b);
//            setParabolicTargets(0.01);
//            computeParabolicBcoeffs();
//            computeCycleSplineOutputs(cycleParabola);
//            graphParabolicSpline = true;
//            repaint();
            // new cycleSpline has knot sequence: 0,0,0,0,1/k,...,1,1,1,1
            DBG("graphing new cycle spline " << n);
            if (highlightCycle == n) {
                highlightCycle = -1;
                graphNewSplineCycle = false;
                plotTargets = false;
            } else {
                highlightCycle = n;
            }
            if (highlightCycle == n) {
                cycleNew = CycleSpline(kVal, a, b);
                setNewTargets(floatBuffer);   // initialize data for cycleNew
                computeNewBcoeffs(floatBuffer);
                resetNewBcoeffs();
                computeCycleSplineOutputs(cycleNew);
                graphNewSplineCycle = true;
                graphSplineCycle = false;
                cycles[0] = cycleNew;
                cycleNew.printData();
                repaint();
            }
        }
    }
}

void GraphComponent::resetNewBcoeffs()
{
    int n = kVal + 3;
    for (int i=1; i<n-1; i++) {
        cycleNew.bcoeffs.set(i, 0);
    }
    int middle = n / 2;
    cycleNew.bcoeffs.set(middle, 0.5);
}

void GraphComponent::iterateCA() {
    
    int n = kVal + 3;
//    for (int i=2; i<n-2; i++) {
//        cycleNew.bcoeffs.set(i, 0);
//    }
    Array<float> temp;
    for (int i=0; i<n; i++) {
        temp.set(i, cycleNew.bcoeffs[i]);
    }
    for (int i=2; i<n-2; i++) {
//        float r = juce::Random::getSystemRandom().nextFloat();
//        float s = juce::Random::getSystemRandom().nextFloat();
//        s = 2 * s - 1; // random float in [-1,1]
//        s *= 0.1;
        float a1 = temp[i-1];
        float a2 = temp[i];
        float a3 = temp[i+1];
        float h = 1.5; // / (float) kVal;
        float val = (a1 + a3) / h;
//        float val = (a1 + a2 + a3) / h;
//        if (r < 0.1) {
//            val += s;
//        }
        cycleNew.bcoeffs.set(i, val);
    }
}

void GraphComponent::randomizeNewCycle()
{
    int n = kVal + 3;
    iterateCA();
    
    // Paul, below is the version that I first showed you, just adding random amounts

//    for (int i=2; i<n-2; i++) {
//        float r = juce::Random::getSystemRandom().nextFloat();
//        r = 2 * r - 1; // random float in [-1,1]
//        r *= 0.1;
//        float val = cycleNew.bcoeffs[i] + r;
//        if (abs(val) < 0.8) {
//            cycleNew.bcoeffs.set(i, val);
//        }
//    }
    
    // this is one is basically a low pass filter
//    for (int i=2; i<n-2; i++) {
//        float r = juce::Random::getSystemRandom().nextFloat();
//        r = 2 * r - 1; // random float in [-1,1]
//        r *= 0.1;
//        float val = (cycleNew.bcoeffs[i] + cycleNew.bcoeffs[i-1])/2;
//        if (abs(val) < 0.9) {
//            cycleNew.bcoeffs.set(i, val);
//            DBG("bcoeffs[" << i << "]:  " << cycleNew.bcoeffs[i]);
//        }
//    }
//    cycleNew.printData();
    
    // Paul, next is the final example we did together this morning

//    Array<float> temp;
//    for (int i=0; i<n; i++) {
//        temp.set(i, cycleNew.bcoeffs[i]);
//    }
//    for (int i=2; i<n-2; i++) {
//
//        float r = juce::Random::getSystemRandom().nextFloat();
//        float s = juce::Random::getSystemRandom().nextFloat();
//        s = 2 * s - 1; // random float in [-1,1]
//        s *= 0.1;
//        float a1 = temp[i-1];
//        float a2 = temp[i];
//        float a3 = temp[i+1];
//        float h = 4; // / (float) kVal;
//        float val = (a3 - 2*a2 + a1) / h;
//        if (r < 0.1) {
//            val += s;
//        }
//        cycleNew.bcoeffs.set(i, val);
//    }
    computeCycleSplineOutputs(cycleNew);
    graphNewSplineCycle = true;
    graphSplineCycle = false;
    repaint();
}

void GraphComponent::computeNewBcoeffs(AudioBuffer<float>& floatBuffer)
{
    // assume bcoeffs[0] = bcoeffs[n-1] = 0, and compute the other n-2
    // this means B-spline graph will hit target 0 at the ends
    // also bcoeffs[1] and bcoeffs[n-2] control the derivatives at the ends
    int k = cycleNew.k, d = cycleNew.d;
    int n = k + d;  // n is full dimension
    float val = 0;
    juce::Array<float> A;
    juce::Array<float> B;
    juce::Array<float> temp;
    for (int i=0; i<(n-2)*(n-2); i++) {  // A is n-2 x n-2
        A.add(0);
        B.add(0);
        temp.add(0);
    }
    // linear system rows i, columns j to solve for c_1 ... c_{n-2}
    // in system Ax=b these are indexed 0 ... n-3
    // the entry A[i,j] should be B^3_j(s_i) for input s_i, but B-splines
    // are shifted forward by one index, so B^3_{j+1}(s_i)
    for (int i=0; i<n-2; i++) {
        for (int j=0; j<n-2; j++) {
            // shift j up by one since newBsplineVal works with nxn system
            // and outer inputs are 0 and 1 with targets 0 already achieved
            val = newBsplineVal(k, j+1, cycleNew.inputs[i]);  // A[i,j]
            A.set(i*(n-2)+j, val);
            temp.set(i*(n-2)+j, val);
            B.set(i*(n-2)+j, 0);
        }
    }
    for (int i=0; i<n-2; i++) {
        B.set(i*(n-2)+i, 1.0);
    }

    gaussElim(n-2, temp, B);
    Array<float> x;
    Array<float> y;
    for (int i=0; i<n-2; i++) {
        x.add(0);
        y.add(0);
    }
    y = multMatCol(n-2, B, cycleNew.targets);
    for (int i=1; i<n-1; i++) {
        cycleNew.bcoeffs.set(i, y[i-1]);
    }
}

void GraphComponent::computeParabolicBcoeffs()
{
    int k = cycleParabola.k, d = cycleParabola.d;
    int n = k + d;
    float val = 0;
    juce::Array<float> A;
    juce::Array<float> B;
    juce::Array<float> temp;
    for (int i=0; i<n*n; i++) {
        A.add(0);
        B.add(0);
        temp.add(0);
    }
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            val = bSplineVal(k, j, cycleParabola.inputs[i]);  // A[i,j]
            A.set(i*n+j, val);
            temp.set(i*n+j, val);
            B.set(i*n+j, 0);
        }
    }
    for (int i=0; i<n; i++) {
        B.set(i*n+i, 1.0);
    }

    gaussElim(n, temp, B);
    Array<float> x;
    Array<float> y;
    for (int i=0; i<n; i++) {
        x.add(0);
        y.add(0);
    }
    y = multMatCol(n, B, cycleParabola.targets);
    for (int i=0; i<n; i++) {
        cycleParabola.bcoeffs.set(i, y[i]);
    }
}

void GraphComponent::mouseDrag (const MouseEvent& event)
{
    DBG ("Dragging at: " << event.getPosition().toString());
}

void GraphComponent::mouseMove (const MouseEvent& event)
{
    mouseOverCycleZero = false;
    if (updateGraph)
    {
        juce::Point<int> P = event.getPosition();
        juce::Point<float> Q (P.getX(),P.getY());
        juce::Point<float> S = screenToSignalCoords(Q);
        if (abs(S.getY()) < 0.01) {
            int i = startIndex;
            while (i < endIndex + 1) {
                float t = cycleZeros[i];
                if (abs(t - S.getX()) < 10) {
                    zeroPoint.setX(t);
                    zeroPoint.setY(0);
                    zeroPoint = signalToScreenCoords(zeroPoint);
                    mouseOverCycleZero = true;
                    repaint();
                }
                i += 1;
            }
        } else {
            mouseOverCycleZero = false;
        }
        repaint();
    }
}

void GraphComponent::setNewTargets(AudioBuffer<float>& floatBuffer)
{
    int k = cycleNew.k;
    int d = cycleNew.d;
    int n = k + d;
    int N = n + d;
    cycleNew.inputs.clear();
    cycleNew.targets.clear();
    cycleNew.knots.clear();
    cycleNew.bcoeffs.clear();
    DBG("printData after clear:");
    cycleNew.printData();
    
    float incr = 1 / (float) k;
    cycleNew.inputs.add(0.5*incr);
    for (int i=1; i<k; i++) {
        cycleNew.inputs.add(i*incr);
    }
    cycleNew.inputs.add(1-0.5*incr);
    
    // New knot sequence: 0,0,0,0,1/k,2/k,...,(k-1)/k,1,1,1,1
    for (int i=0; i<d+1; i++) {
        cycleNew.knots.add(0);
    }
    for (int i=4; i<N-d; i++) {
        cycleNew.knots.add(cycleNew.knots[i-1]+incr);
    }
    for (int i=0; i<d+1; i++) {
        cycleNew.knots.add(1);
    }
    
    // set bcoeffs to 0 then will compute n-2 of them c_1,...,c_{n-1}
    // in function computeNewBcoeffs
    for (int i=0; i<n; i++) {
        cycleNew.bcoeffs.add(0);
    }
    
    float a = cycleNew.a;
    float b = cycleNew.b;
    for (int i=0; i<n-2; i++) {
        float t = cycleNew.inputs[i];
        float output = interpFloat(a + t * (b-a), floatBuffer);
        cycleNew.targets.set(i, output);
    }
    
    DBG("New params set for CycleSpline knot sequence: ");
    DBG("N = " << N << " n = " << n << " k = " << k);
    for (int i=0; i<N+1; i++) {
        std::cout << cycleNew.knots[i] << " ";
    }
    DBG("");
}
