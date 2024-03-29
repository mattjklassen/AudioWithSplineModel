/*
  ==============================================================================

    CycleInterp.cpp
    Created: 8 Jan 2022 3:50:44pm
    Author:  Matt Klassen

  ==============================================================================
*/

// The standard method has become as of June 7 (SMC week) 2022:
// First compute delta model which gives smoothest breakdown into cycles allowing for y0 and y1 nonzero.
// Then can use other interpolation methods based on delta model, ie. use cycles from delta as key cycles
// and interpolate the rest with any interpolation method.  The interpolated model can either use the
// delta values y0 and y1 or not.  Not using them gives a smaller model, based on less data.  When using
// the delta values, one simply interpolates B-spline coefficients and constructs intermediate cycles with
// those plus delta values.  The delta model cycles are now stored in deltaCycleArray, and the interpolation
// model cycles are stored in allCycleArray.

// To compute delta model, select "No Cycle Interp" and "Use Delta Model" and choose k value,
// then "Compute Model".
// To interpolate with key cycles from delta model, leave "use Delta Model" selected, choose m,
// and choose interpolation method, then "Compute Model".  In this case a, b, y0 and y1 are used
// from the delta model for all cycles.
// In both cases the the model is written to writeBuffer and this is written to file output.wav

#include "MainContentComponent.h"

void MainContentComponent::initTimbreCoeffs()
{
    for (int i=0; i<33; i++) {
        timbreCoeffs.add(0);
    }
}

void MainContentComponent::getSplinusoidKeyBcoeffs()
{
    int n = 33;
    computeSplinusoidBcoeffs();
    auto file = File::getSpecialLocation(File::userHomeDirectory).getChildFile("keyBcoeffs");
    FileOutputStream output (file);
    if (output.openedOk())
    {
        output.setPosition (1);  // default would append
        output.truncate();
    }
    // write key cycle bcoeffs from splinusoid to keyBcoeffs array
    for (int i=0; i<keys.size(); i++) {
        for (int p=0; p<n; p++) {
            keyBcoeffs.set(i*n+p, splinusoid.bcoeffs[p]);
            output.writeFloat(keyBcoeffs[i*n+p]);
        }
    }
}

void MainContentComponent::writeScaleToBuffer()
{
    int n = dVal + kVal;
    // use a set of key cycles with bcoeffs in small model of instrument
    DBG("writing scale to buffer");
//    if (sampleRate < 44100) {  // ie. no audio loaded, so will use spline generator
        sampleRate = 44100;
        DBG("setting sampleRate to 44100");
//    }
//    numCycles = (int)(freqGuess * lengthInSecondsCA);
    DBG("numCycles: " << numCycles << " sampleRate: " << (int)sampleRate << " freqGuess: " << freqGuess);
    samplesPerCycleCA = (int)(((float)sampleRate) / freqGuess);
    DBG("samplesPerCycleCA is computed in writeScaleToBuffer: " << samplesPerCycleCA);
    samplesPerCycleGuess = (float)sampleRate / freqGuess;
    DBG("samplesPerCycleGuess is computed in writeScaleToBuffer: " << samplesPerCycleGuess);
    DBG("period length should be samplePerCycleGuess");
    writeBuffer.clear();
    writeBuffer.setSize(1, sampleRate * lengthInSecondsCA * 13);  // channels, samples
//    writeBuffer.setSize(1, sampleRate+1);
    DBG("writeBuffer.size(): " << writeBuffer.getNumSamples());
    uint32 time1 = Time::getMillisecondCounter();
    initTimbreCoeffs();
//    timbreCoeffs.set(31, 1);
//    initECA();
//    vectorECA.set(31, 1);
    // Now compute 1 sec of middle C, then multiply freq and compute 12 steps
    for (int i=0; i<13; i++) {
//        if (i > 0) {
//            vectorECA = computeCA(30, vectorECA);
//            for (int j=0; j<n; j++) {
//                timbreCoeffs.set(j, (float)vectorECA[j]);
//            }
//        }
//        printECA();
        float b = (float)i / (float)12;
        float a = 1 - b;
        // Need to set timbreCoeffs here
        timbreCoeffs.set(0, a);
        timbreCoeffs.set(20, b);
        mixKeyBcoeffs();
//        keyBcoeffs = splinumid1KeyBcoeffs;
        setKeyCycleIndices();
//        getSplinusoidKeyBcoeffs();
        computeMetaSplines();
        writeNextCycleToBuffer(i);
        freqGuess *= 1.05946309436;  // = 2^(1/12)
        // this advances fund freq by one semitone
        numCycles = (int)(freqGuess * lengthInSecondsCA);
        samplesPerCycleGuess = sampleRate / freqGuess;
        samplesPerCycleCA = (int) (sampleRate / freqGuess);
//        samplesPerCycleCA = (int) samplesPerCycleGuess;
    }
    graphView.graphCAmodel = true;
    graphView.repaint();
    uint32 time2 = Time::getMillisecondCounter();
    float timeVal = (float(time2)-float(time1)) * 0.001;
    DBG("time to compute model (sec):  " << timeVal);
}


// this version is basically same as writeNormalizedCyclesToBuffer() but takes input J
void MainContentComponent::writeNextCycleToBuffer(int J)
{
    // J = index of next sample to write in buffer, writing in 1 sec chunks use J*sampleRate + ...
    // write all cycles using constant cycle length averageSamplesPerCycle
    DBG("computing and writing normalized cycles");
    DBG("freqGuess: " << freqGuess);
    DBG("samplesPerCycleCA: " << samplesPerCycleCA);
    DBG("samplesPerCycleGuess: " << samplesPerCycleGuess);
    // choose number of cycles for Attack part A to be keysCA * mVal
    // for example, 5*5=25 works well for A450
    // Note: we are not using keysCA other than to compute A if modelWithDelta is true
    float A = (keysCA * mVal) * samplesPerCycleCA;
    const float D = A / 2;
    const float S = (freqGuess * samplesPerCycleCA - A - D) / 4 ;
    const float R = 3 * S;
//    DBG("A: " << A << "  D: " << D << "  S: " << S << "  R: " << R);
    int i = 0;
    float a = 0, b=1, ratio = 1;
    CycleSpline cycle;
    int n = kVal + dVal;
//    float total = (float) numCycles;
    while (i < numCycles) {  // i is cycle number whether key or non-key
        // the next block forces each cycle to have length samplesPerCycleGuess
        // as computed by computeCycleSplineOutputs(cycle) using only a and b
        // in order to achieve correctly the fundamental frequency freqGuess
        // we also need to
        if (i == 0) {
            a = 0;
//            b = (float)((int)samplesPerCycleCA) - 0.5;
            b = samplesPerCycleGuess;
        } else {
            a = b;
//            b = a + (int)samplesPerCycleCA;
            b = a + samplesPerCycleGuess;
        }
        if (isKey(i)) {   // use key cycle i
            int j = keyIndex(i);
            if (modelWithSmall) {
                cycle = CycleSpline(kVal, a, b);  // this is key cycle j
                for (int p=0; p<n; p++) {
                    cycle.bcoeffs.set(p, keyBcoeffs[j*n+p]);
                }
            } else {
                if (modelWithDelta) {
                    cycle = deltaCycleArray[i];
                } else {
                    cycle = keyCycleArray[j];
                }
            }
            cycle.a = a; cycle.b = b;   // reset to normalized a and b
            cycle.y0 = 0; cycle.y1 = 0; // reset to avoid Delta floating spline
//            DBG("writing key cycle: " << i);
        } else {
            // set bcoeffs for (non-key) cycle i using meta-spline when i < index of last key cycle
            if (i < keys[keys.size()-1]) {
                cycle = CycleSpline(kVal, a, b);  // this is cycle i
                for (int j=0; j<n; j++) {
                    cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
                }
            } else {  // will not use this block if modelWithDelta true
                cycle = CycleSpline(kVal, a, b);
                for (int j=0; j<n; j++) {
                    // use splinusoid
                    cycle.bcoeffs.set(j, splinusoid.bcoeffs[j]);
                }
            }
        }
        computeCycleSplineOutputs(cycle);
        ratio = 1;  // this overrides the previous block
        int M = cycle.outputs.size();
//        DBG("cycle " << i << " numSamples M: " << M << " a: " << a << " b: " << b << " b-a: " << b-a);
        int L = ((int)a) + 1;  // first sample index to write
        if (i == 0) {
            L = 0;
        }
        for (int j=0; j<M; j++) {  // writing samples loop
            if (useADSR) {
                ratio = 0.5;  // Sustain value 1/2
                if (j+L < A) {
                    ratio = (1/A) * (j+L);  // Attack value goes 0 to 1
                } else {
                    if (j+L < A+D) {
                        ratio = 1 - (1/(2*D))*(j+L-A);  // Decay value, 1 to 1/2
                    } else {
                        if (j+L > A+D+S) {
                            ratio = -(1/(2*R))*(j+L-(A+D+S+R));  // Release value, 1/2 to 0
                        }
                    }
                }
            }
            // here we restrict envelope only to the release part
            ratio = 1;
            if (j+L > A+D+S) {
                ratio = -(1/R)*(j+L-(A+D+S+R));  // Release ratio, 1 to 0
            }
            writeBuffer.setSample(0, J*sampleRate + j+L, ratio * cycle.outputs[j]);
        }
        i += 1;
    }  // end while (i < numCycles)
}


void MainContentComponent::mixKeyBcoeffs()
{
    int n = 33, numKeys = 18;
    Array<float> instrVals;  // there are n=33 instr
    float denom = 0;
//    int numInstr = instrVals.size();
    for (int j=0; j<n; j++) {
        denom += timbreCoeffs[j];
    }
    if (denom > 0) {
        // OK
    } else {
        DBG("timbre Coeffs are all zero");
        return;
    }
    for (int i=0; i<numKeys; i++) {
        for (int p=0; p<n; p++) {
            instrVals.clear();
            // plucked
            instrVals.add(guitarKeyBcoeffs[i*n+p]);       // 0
            instrVals.add(guitarpizzKeyBcoeffs[i*n+p]);   // 1
            instrVals.add(guitarpontKeyBcoeffs[i*n+p]);   // 2
            instrVals.add(guitarharmKeyBcoeffs[i*n+p]);   // 3
            instrVals.add(theorboKeyBcoeffs[i*n+p]);      // 4
            instrVals.add(theorbopontKeyBcoeffs[i*n+p]);  // 5
            instrVals.add(harpKeyBcoeffs[i*n+p]);         // 6
            instrVals.add(violinpizzKeyBcoeffs[i*n+p]);   // 7
            instrVals.add(cellopizzKeyBcoeffs[i*n+p]);    // 8
            
            // bowed
            instrVals.add(violinKeyBcoeffs[i*n+p]);       // 9
            instrVals.add(violinpontKeyBcoeffs[i*n+p]);   // 10
            instrVals.add(celloKeyBcoeffs[i*n+p]);        // 11
            instrVals.add(cellopontKeyBcoeffs[i*n+p]);    // 12
            instrVals.add(stringsKeyBcoeffs[i*n+p]);      // 13
            
            // hammered or bells
            instrVals.add(pianoKeyBcoeffs[i*n+p]);        // 14
            instrVals.add(marimbaKeyBcoeffs[i*n+p]);      // 15
            instrVals.add(vibraphoneKeyBcoeffs[i*n+p]);   // 16
            instrVals.add(celesteKeyBcoeffs[i*n+p]);      // 17
            instrVals.add(cimbalomKeyBcoeffs[i*n+p]);     // 18
            instrVals.add(handbellKeyBcoeffs[i*n+p]);     // 19

            // winds
            instrVals.add(fluteKeyBcoeffs[i*n+p]);        // 20
            instrVals.add(bassoonKeyBcoeffs[i*n+p]);      // 21
            instrVals.add(clarinetKeyBcoeffs[i*n+p]);     // 22
            instrVals.add(enghornKeyBcoeffs[i*n+p]);      // 23
            instrVals.add(oboeKeyBcoeffs[i*n+p]);         // 24
            
            // brass
            instrVals.add(trumpetKeyBcoeffs[i*n+p]);      // 25
            instrVals.add(frhornKeyBcoeffs[i*n+p]);       // 26
            
            // organ
            instrVals.add(organKeyBcoeffs[i*n+p]);        // 27
            instrVals.add(wurliKeyBcoeffs[i*n+p]);        // 28
            
            // synth
            instrVals.add(splinusoidKeyBcoeffs[i*n+p]);   // 29
            instrVals.add(splinufuzzKeyBcoeffs[i*n+p]);   // 30
            instrVals.add(splinumid1KeyBcoeffs[i*n+p]);   // 31
            instrVals.add(splinumid2KeyBcoeffs[i*n+p]);   // 32

            float newVal = 0;
            for (int j=0; j<n; j++) {
                newVal += timbreCoeffs[j] * instrVals[j];
            }
            keyBcoeffs.set(i*n+p, newVal / denom);
        }
    }
}

void MainContentComponent::writeModelToBuffer()
{
    if (otherCycleInterp) {
        if (writeShortFade) {
            // use this to write 1 sec of original audio with fade out
            numSamples = 44100;
            writeBuffer.clear();
            writeBuffer.setSize(1, numSamples);  // channels, samples
            float val = 0;
            float denom = 29400;
            float ramp = 0;
            for (int i=0; i<numSamples; i++) {
                val = floatBuffer.getSample(0, i);
                if (i > 14700) {
                    ramp = (float)(44100-i) / denom;
                    val *= ramp;
                }
                writeBuffer.setSample(0, i, val);
            }
            writeWavFile();
            return;
        }
        if (audioDataLoaded) {
            // move on for case of added or removed cycles
        } else {  // use this area to do other models with no audio loaded
            // here we build a cubic spline model with no initial audio data
            // next: need to write splinusoid bcoeffs to file
            // do 1 second of data at 262 Hz
            if (sampleRate == 0) {
                sampleRate = 44100;
            }
            numCycles = 262;
            int samplesPerCycle = (int) ((float)sampleRate / (float)262);
            int numSamples = numCycles * samplesPerCycle;
            numSamples = 44100;
            writeBuffer.clear();
            writeBuffer.setSize(1, numSamples);  // channels, samples
            computeSplinusoidBcoeffs();
            initECA();
            int k = kVal;
            int n = k + dVal;
            CycleSpline cycle = CycleSpline(k, 0, 1);
            rVal = 30;           // to use ECA(30)
            float a = 0, b = 1;
            int i = 0;
            while (i < numCycles) {
                // the next block forces each cycle to have samplesPerCycle samples
                // as computed by computeCycleSplineOutputs(cycle) using only a and b
                if (i == 0) {
                    a = -0.5;
                    b = (float)(samplesPerCycle) - 0.5;
                } else {
                    a = b;
                    b = a + samplesPerCycle;
                }
                splinusoid.a = a; splinusoid.b = b;
                // recompute cycle using splinusoid bcoeffs modified by ECA(30):
//                if ((i < 4*16) && (i % 4 == 0)) {
//                    cycle.a = a; cycle.b = b;
//                    for (int p=0; p<n; p++) {
//                        float mult = 0.8 + 0.3 * vectorECA[p];
//                        cycle.bcoeffs.set(p, splinusoid.bcoeffs[p] * mult);
//                        keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
//                    }
//                    vectorECA = computeCA(rVal, vectorECA);  // iterate next row for ECA(30)
//                    computeCycleSplineOutputs(cycle);
//                    int M = cycle.outputs.size();
//                    int L = ((int)a) + 1;  // first sample index to write
//                    if (i == 0) {
//                        L = 0;
//                    }
//                    for (int j=0; j<M; j++) {
//                        writeBuffer.setSample(0, j+L, cycle.outputs[j]);
//                    }
//                    i += 1;
//                } else {  // otherwise just use splinusoid
//                if (i > 0) {
//                    float r = juce::Random::getSystemRandom().nextFloat();
//                    r = 2 * r - 1;
//                    r *= 0.3;  // random float in [-.1,.1]
//                    splinusoid.y0 = splinusoid.y1;
//                    splinusoid.y1 = r;
//                }
                computeCycleSplineOutputs(splinusoid);
                int M = splinusoid.outputs.size();
                int L = ((int)a) + 1;  // first sample index to write
                if (i == 0) {
                    L = 0;
                }
                for (int j=0; j<M; j++) {
                    writeBuffer.setSample(0, j+L, splinusoid.outputs[j]);
                }
                i += 1;
            }  // end while i < numCycles
            return;
        }  // end if no audio data
    }  // end if otherCycleInterp which can be used for different purposes
    if (modelWithCA) {
        // this uses one (audio data) cycle to extract bcoeffs from,
        // use selected cycleNew from graphView
        // or extract one cycle from beginning of audio data
        // or use splinusoid
        // or use a set of key cycles from Delta Model in allCycleArray
        DBG("modeling with CA");
        if (sampleRate == 0) {  // ie. no audio loaded, so will use spline generator
            sampleRate = 44100;
            DBG("setting sampleRate to 44100");
        }
        numCycles = freqGuess * lengthInSecondsCA;
        DBG("numCycles set to: " << numCycles);
        samplesPerCycleCA = sampleRate / freqGuess;
        writeBuffer.clear();
        writeBuffer.setSize(1, sampleRate * lengthInSecondsCA);  // channels, samples
        DBG("writeBuffer.size(): " << writeBuffer.getNumSamples());
        uint32 time1 = Time::getMillisecondCounter();
        setKeyCycleIndices();
        modelWithCAplusCycleInterp();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
//        graphView.graphCAmodel = true;
        return;
    }
    if (randomizeBcoeffs) {
        uint32 time1 = Time::getMillisecondCounter();
        setKeyCycleIndices();
        modelWithCAplusCycleInterp();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
        return;
    }
    graphView.graphCAmodel = false;
    if (noCycleInterp) {
        uint32 time1 = Time::getMillisecondCounter();
        if (modelWithDelta) {
            computeModelWithDelta();
        } else {
            modelWithoutCycleInterp();
        }
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
        return;
    }
    if (normalizeCycleLength) {
        uint32 time1 = Time::getMillisecondCounter();
        setKeyCycleIndices();
        computeKeyCycles();
    //  writeKeyCyclesToBuffer();
        computeMetaSplines();
    //  computeAllCycles();
    //  writeNonKeyCyclesToBuffer();
        writeNormalizedCyclesToBuffer();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
    } else {
        // this covers the main cases of cycle interp like regular etc. --MAIN CASE--
        // here we do not want to graph Delta Model, only possibly use its key cycles
        // when modelWithDelta = true.  So we set graphView.graphDeltaModel = false
        DBG("Main Case of Cycle Interpolation");
        graphView.graphDeltaModel = false;
        uint32 time1 = Time::getMillisecondCounter();
        setKeyCycleIndices();
        // if modelWithDelta = true, or modelWithoutDelta = true
        // then the next function changes to use already computed key cycles
        // if they are both false, then keyCycles are computed
        writeKeyCyclesToBuffer();
        // metasplines are computed based on B-spline coeffs only, not y0 and y1
        computeMetaSplines();
        // non-key cycles are updated to use the y0 and y1 values from the delta model
        writeNonKeyCyclesToBuffer();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
    }
}

// compute mean square error on cycle, comparing current cycle spline to audio data
float MainContentComponent::computeCycleSplineError(CycleSpline& cycle, AudioBuffer<float>& samples)
{
    int N = cycle.outputs.size();
    float output = 0, sample = 0, val = 0;
    float a = cycle.a;
    int j = (int) a + 1;
    for (int i=0; i<N; i++) {
        sample = samples.getSample(0, i+j);
        val = sample - cycle.outputs[i];
        output += val * val;
    }
    return output / (float) N;
}

// use bcoeffs from cycle1 to compute error on cycle2
float MainContentComponent::compareCycleSplineError(CycleSpline& cycle1,  CycleSpline& cycle2, AudioBuffer<float>& samples)
{
    int N = cycle2.outputs.size();
    float output = 0, sample = 0, val = 0;
    float a = cycle2.a;
    int j = (int) a + 1;
    // now set bcoeffs of cycle2 to those of cycle1 and compute outputs
    Array<float> temp;
    int k = cycle2.k, d = cycle2.d, n = k+d;
    for (int p=0; p<n; p++) {
        temp.add(cycle2.bcoeffs[p]);
        cycle2.bcoeffs.set(p, cycle1.bcoeffs[p]);
    }
    computeCycleSplineOutputs(cycle2);
    // compute error for cycle2 using cycle1 bcoeffs
    for (int i=0; i<N; i++) {
        sample = samples.getSample(0, i+j);
        val = sample - cycle2.outputs[i];
        output += val * val;
    }
    // reset cycle2 bcoeffs and outputs back
    for (int p=0; p<n; p++) {
        cycle2.bcoeffs.set(p, temp[p]);
    }
    computeCycleSplineOutputs(cycle2);
    return output / (float) N;
}

// cubic polynomial p(t) = t*(t-1/2)*(t-1) has B-spline coefficients
// computed from its polar form F[x,y,z] = x*y*z -(3/2)(x*y + x*z + y*z) + (1/2)*(x+y+z)
void MainContentComponent::computeCubicSinusoidBcoeffs()
{
    int k = kVal;
    int d = dVal;
    int n = k + d;
    int N = n + d;
    cubicSinusoidSpline.k = k;
    cubicSinusoidSpline.d = d;
    cubicSinusoidSpline.inputs.clear();
    cubicSinusoidSpline.targets.clear();
    cubicSinusoidSpline.knots.clear();
    cubicSinusoidSpline.bcoeffs.clear();
    float incr = 1 / (float) k;
    
    cubicSinusoidSpline.knots.add(0);  // knot t_0
    cubicSinusoidSpline.knots.add(0);  // knot t_1
    cubicSinusoidSpline.knots.add(0);  // knot t_2
    cubicSinusoidSpline.knots.add(0);  // knot t_3
    for (int i=4; i<N-3; i++) {
        cubicSinusoidSpline.knots.add(cubicSinusoidSpline.knots[i-1]+incr);
    }
    cubicSinusoidSpline.knots.add(1);  // knot t_{N-3}
    cubicSinusoidSpline.knots.add(1);  // knot t_{N-2}
    cubicSinusoidSpline.knots.add(1);  // knot t_{N-1}
    cubicSinusoidSpline.knots.add(1);  // knot t_N
    
    float val = 0, a = 0, b = 0, c = 0;
    for (int i=0; i<n; i++) {
        a = cubicSinusoidSpline.knots[i+1];
        b = cubicSinusoidSpline.knots[i+2];
        c = cubicSinusoidSpline.knots[i+3];
        val = a*b*c - (0.5)*(a*b + a*c + b*c) + (a+b+c) / (float)6;
        val *= 3;
        if (i > 0 && i < n-1) {
            if (mVal > 1) {
                if (i % 2 == 0)
                    val += (float(mVal) / (float)1000);
                else {
                    val -= (float(mVal) / (float)1000);
                }
            }
        }
        cubicSinusoidSpline.bcoeffs.set(i, val);
    }
}

// Define here the "splinusoid" which looks like one period of sin(2Pi*t) on [0,1], but
// with slope at t=0 equal to 6 (approx 2*Pi).
// splinusoid = cubic spline based on p(t) = 6*t - 32*t^3, defined as f(t) =
// p(t) on [0,1/4], -p(t-1/2) on [1/4,3/4], and p(t-1) on [3/4,1].
// f(t) has B-spline coefficients for knot sequence 0,0,0,0,1/4,1/2,3/4,1,1,1,1 given as:
// c_0=0, c_1=1/2, c_2=3/2, c_3=0, c_4=-3/2, c_5=-1/2, c_6=0, but this is only k=4.
// Below we compute its B-spline coefficients for standard knot sequence 0,0,0,0,1/k,...,(k-1)/k,1,1,1,1
// computed from the polar form for p(t): F[x,y,z] = 2(x+y+z)-32xyz, -p(t-1/2): -F[x-1/2,y-1/2,z-1/2],
// and p(t-1): F[x-1,y-1,z-1].  For simplicity, suppose k=4L, so L/k=1/4, 2L/k=1/2, 3L/k=3/4.
// Also note t_4=1/k,...,t_{N-4}=t_{k+2}=(k-1)/k, so t_{i+3}=i/k, for 3<i<N-3.
// The polar form switches when t_J reaches 1/4=L/k (so J=L+3)
// and then again when t_J reaches 3/4=3L/k (so J=3L+3).
// So for J=0,...,n-1, we compute Bcoeffs c_J as:
// c_J = F[t_{J+1},t_{J+2},t_{J+3}] for J=0,...,L+2,
// c_J = -F[t_{J+1}-1/2,t_{J+2}-1/2,t_{J+3}-1/2] for J=L+3,...,3L+2,
// c_J = F[t_{J+1}-1,t_{J+2}-1,t_{J+3}-1] for J=3L+3,...,n-1=k+2=N-d-1.

void MainContentComponent::computeSplinusoidBcoeffs()
{
    int k = kVal;
    int d = dVal;
    int n = k + d;
    int N = n + d;
    int L = k / 4;
    float gain = 0.7;
    splinusoid.k = k;
    splinusoid.d = d;
    splinusoid.inputs.clear();
    splinusoid.targets.clear();
    splinusoid.knots.clear();
    splinusoid.bcoeffs.clear();
    float incr = 1 / (float) k;
    
    splinusoid.knots.add(0);  // knot t_0
    splinusoid.knots.add(0);  // knot t_1
    splinusoid.knots.add(0);  // knot t_2
    splinusoid.knots.add(0);  // knot t_3
    for (int i=4; i<N-3; i++) {
        splinusoid.knots.add(splinusoid.knots[i-1]+incr);
    }
    splinusoid.knots.add(1);  // knot t_{N-3}
    splinusoid.knots.add(1);  // knot t_{N-2}
    splinusoid.knots.add(1);  // knot t_{N-1}
    splinusoid.knots.add(1);  // knot t_N
    
    float val = 0, a = 0, b = 0, c = 0;
    for (int i=0; i<n; i++) {
        a = splinusoid.knots[i+1];
        b = splinusoid.knots[i+2];
        c = splinusoid.knots[i+3];
        if (i < L+3) {
            val = 2*(a+b+c) - 32*a*b*c;
        } else {
            if (i < 3*L+3) {
                a -= 0.5; b -= 0.5; c -= 0.5;
                val = -2*(a+b+c) + 32*a*b*c;
            } else {
                a -= 1; b -= 1; c -= 1;
                val = 2*(a+b+c) - 32*a*b*c;
            }
        }
        // oscillation pattern for bcoeffs, use in instrument splinufuzz:
//        if ((i>1) && (i<n-2)) {
//            val += (i % 2) * 0.2;
//            val -= ((i+1) % 2) * 0.2;
//        }
        // for instrument splinumid1 we use only: if (i == 16) { val = 0.8; }
        // for instrument splinumid2 we use all i=12,...,20 below:
        // these are no longer related to the splinusoid, but are simply setting some fixed mid-values
//        val = 0;
//        if (i == 12) { val = 0.4; }
//        if (i == 13) { val = -0.5; }
//        if (i == 14) { val = 0.6; }
//        if (i == 15) { val = -0.7; }
//        if (i == 16) { val = 0.8; }
//        if (i == 17) { val = -0.7; }
//        if (i == 18) { val = 0.6; }
//        if (i == 19) { val = -0.5; }
//        if (i == 20) { val = 0.4; }
//        gain = 1;
        // for square wave spline approximation set coeffs to 1/2 if i < 16, 0 if i = 16, -1/2 if i > 16
        if (i < 16) { val = 0.5; }
        if (i == 16) { val = 0; }
        if (i > 16) { val = -0.5; }
        gain = 1;
        
        // set gain
        val *= gain;
        // write bcoeff
        splinusoid.bcoeffs.set(i, val);
    }
}

// Define "splinusoid2" based on splinusoid but extended to interval [0,3/2] by
// adding on another splinusoid of half amplitude and half period length in [1,3/2].
// Then we scale back to the interval [0,1] again, and call it h(t).
// Suppose f(t) is splinusoid on [0,1], and g(t)=(1/2)f(2(t-1)) on [1,3/2].
// Then h(t) will be f((3/2)t) on [0,2/3] and g((3/2)t) on [2/3,1].
// So replace x with (3/2)x etc in polar forms on [0,2/3], then on [2/3,1] use
// polar forms for g((3/2)t) etc.
// Interval [0,1] is now in thirds, with smallest subdivision size 1/12, so we should
// also have k=12L, so the divisions for polar forms are at j/k with j equal to:
// (1/6)k=2L, (1/2)k=6L, (2/3)k=8L, 9L, 11L, 12L=k.
// Again: t_4=1/k,...,t_{N-4}=t_{k+2}=(k-1)/k, so t_{i+3}=i/k, for 3<i<N-3.

// splinusoid = cubic spline based on p(t) = 6*t - 32*t^3, defined as f(t) =
// p(t) on [0,1/4], -p(t-1/2) on [1/4,3/4], and p(t-1) on [3/4,1].
// f(t) has B-spline coefficients for knot sequence 0,0,0,0,1/4,1/2,3/4,1,1,1,1 given as:
// c_0=0, c_1=1/2, c_2=3/2, c_3=0, c_4=-3/2, c_5=-1/2, c_6=0.
// Below we compute its B-spline coefficients for standard knot sequence 0,0,0,0,1/k,...,(k-1)/k,1,1,1,1
// computed from the polar form for p(t): F[x,y,z] = 2(x+y+z)-32xyz, -p(t-1/2): -F[x-1/2,y-1/2,z-1/2],
// and p(t-1): F[x-1,y-1,z-1].  For simplicity, suppose k=4L, so L/k=1/4, 2L/k=1/2, 3L/k=3/4.
// Also note t_4=1/k,...,t_{N-4}=t_{k+2}=(k-1)/k, so t_{i+3}=i/k, for 3<i<N-3.
// The polar form switches when t_J reaches 1/4=L/k (so J=L+3)
// and then again when t_J reaches 3/4=3L/k (so J=3L+3).
// So for J=0,...,n-1, we compute Bcoeffs c_J as:
// c_J = F[t_{J+1},t_{J+2},t_{J+3}] for J=0,...,L+2,
// c_J = -F[t_{J+1}-1/2,t_{J+2}-1/2,t_{J+3}-1/2] for J=L+3,...,3L+2,
// c_J = F[t_{J+1}-1,t_{J+2}-1,t_{J+3}-1] for J=3L+3,...,n-1=k+2=N-d-1.
// ...
// this one is not implemented yet ...

void MainContentComponent::computeSplinusoid2Bcoeffs()
{
    int k = kVal;
    int d = dVal;
    int n = k + d;
    int N = n + d;
    int L = k / 4;
    float gain = 0.7;
    splinusoid2.k = k;
    splinusoid2.d = d;
    splinusoid2.inputs.clear();
    splinusoid2.targets.clear();
    splinusoid2.knots.clear();
    splinusoid2.bcoeffs.clear();
    float incr = 1 / (float) k;
    
    splinusoid2.knots.add(0);  // knot t_0
    splinusoid2.knots.add(0);  // knot t_1
    splinusoid2.knots.add(0);  // knot t_2
    splinusoid2.knots.add(0);  // knot t_3
    for (int i=4; i<N-3; i++) {
        splinusoid2.knots.add(splinusoid2.knots[i-1]+incr);
    }
    splinusoid2.knots.add(1);  // knot t_{N-3}
    splinusoid2.knots.add(1);  // knot t_{N-2}
    splinusoid2.knots.add(1);  // knot t_{N-1}
    splinusoid2.knots.add(1);  // knot t_N
    
    float val = 0, a = 0, b = 0, c = 0;
    for (int i=0; i<n; i++) {
        a = splinusoid2.knots[i+1];
        b = splinusoid2.knots[i+2];
        c = splinusoid2.knots[i+3];
        if (i < L+3) {
            val = 2*(a+b+c) - 32*a*b*c;
        } else {
            if (i < 3*L+3) {
                a -= 0.5; b -= 0.5; c -= 0.5;
                val = -2*(a+b+c) + 32*a*b*c;
            } else {
                a -= 1; b -= 1; c -= 1;
                val = 2*(a+b+c) - 32*a*b*c;
            }
        }
        if ((i>1) && (i<n-2)) {
            val += (i % 2) * 0.2;
            val -= ((i+1) % 2) * 0.2;
        }
        val *= gain;
        splinusoid2.bcoeffs.set(i, val);
    }
}


// cubic polynomial q(t) = 3*t^2-2*t^3 has B-spline coefficients
// computed from its polar form F[x,y,z] = x*y + x*z + y*z - 2*x*y*z
void MainContentComponent::computeCubicPolyBcoeffs()
{
    int k = kVal;
    int d = dVal;
    int n = k + d;
    int N = n + d;
    cubicDeltaSpline.k = k;
    cubicDeltaSpline.d = d;
    cubicDeltaSpline.inputs.clear();
    cubicDeltaSpline.targets.clear();
    cubicDeltaSpline.knots.clear();
    cubicDeltaSpline.bcoeffs.clear();
    float incr = 1 / (float) k;
    
    cubicDeltaSpline.knots.add(0);  // knot t_0
    cubicDeltaSpline.knots.add(0);  // knot t_1
    cubicDeltaSpline.knots.add(0);  // knot t_2
    cubicDeltaSpline.knots.add(0);  // knot t_3
    for (int i=4; i<N-3; i++) {
        cubicDeltaSpline.knots.add(cubicDeltaSpline.knots[i-1]+incr);
    }
    cubicDeltaSpline.knots.add(1);  // knot t_{N-3}
    cubicDeltaSpline.knots.add(1);  // knot t_{N-2}
    cubicDeltaSpline.knots.add(1);  // knot t_{N-1}
    cubicDeltaSpline.knots.add(1);  // knot t_N
    
    float val = 0, a = 0, b = 0, c = 0;
    for (int i=0; i<n; i++) {
        a = cubicDeltaSpline.knots[i+1];
        b = cubicDeltaSpline.knots[i+2];
        c = cubicDeltaSpline.knots[i+3];
        val += a*b + a*c + b*c - 2*a*b*c;
        cubicDeltaSpline.bcoeffs.set(i, val);
    }
}

void MainContentComponent::addDeltaToCycleBcoeffs(CycleSpline cycle)
{
    const float delta = cycle.y1 - cycle.y0;
    if (delta > 0) {
        const int k = kVal;
        const int d = dVal;
        const int n = k + d;
        float val = 0;
        for (int i=0; i<n; i++) {
            val = delta * cubicDeltaSpline.bcoeffs[i];
            cycle.bcoeffs.set(i, val);
        }
    }
}

// cubic polynomial p(t, delta) = y_0 + delta * (3*t^2-2*t^3) as B-spline sum has B-spline coefficients
// c_0 = y_0 = c_1, c_{n-2} = y_1 = c_{n-1}, where y_1 = y_0 + delta.
// We solve for the bcoeffs of p(t, delta) - y_0 using the polar form for q(t) = 3*t^2-2*t^3,
// which also happens to be just the sum of the last two cubic Bernstein polynomials
// 3*(1-t)*t^2 + t^3 = B^3_2(t) + B^3_3(t) = 1 - B^3_0(t) - B^3_1(t) = 1 - (1-t)^3 - 3*(1-t)^2*t.
// F[x,y,z] = x*y + x*z + y*z - 2*x*y*z.
void MainContentComponent::computeCubicDeltaBcoeffs(float y0, float delta, CycleSpline cycle)
{
    int k = cycle.k;
    int d = cycle.d;
    int n = k + d;
    int N = n + d;
    cycle.inputs.clear();
    cycle.targets.clear();
    cycle.knots.clear();
    cycle.bcoeffs.clear();
    float incr = 1 / (float) k;

    cycle.knots.add(0);  // knot t_0
    cycle.knots.add(0);  // knot t_1
    cycle.knots.add(0);  // knot t_2
    cycle.knots.add(0);  // knot t_3
    for (int i=4; i<N-3; i++) {
        cycle.knots.add(cycle.knots[i-1]+incr);
    }
    cycle.knots.add(1);  // knot t_{N-3}
    cycle.knots.add(1);  // knot t_{N-2}
    cycle.knots.add(1);  // knot t_{N-1}
    cycle.knots.add(1);  // knot t_N

    float val = 0, a = 0, b = 0, c = 0;
    for (int i=0; i<n; i++) {
        val = y0;
        a = cycle.knots[i+1];
        b = cycle.knots[i+2];
        c = cycle.knots[i+3];
        val += delta * (a*b + a*c + b*c - 2*a*b*c);
        cycle.bcoeffs.set(i, val);
    }
}

// cycle1 is cycle j-1, cycle2 is cycle j, cycle(test) is the test for cycle j.
// Specifically, cycle(test) uses the interval [a,b] and y0, y1 from cycle j, but uses cycle j-1 bcoeffs.
// This function computes new b = endpoint of cycle j, based on comparison of slopes using cycle(test).
// Choose the minimum such error within p samples of cycle[j].b, with end of cycle j lifted to y1.
// Previous cycle j-1 will have endpoints (a(j-1),y0(j-1)) and (b(j-1),y1(j-1)), so cycle(test)
// will have endpoints (a(j),y0(j)) and (b(j),y1(j)), with a(j)=b(j-1) and y0(j)=y1(j-1).
// This way cycle(test) will start where cycle j-1 ends.  Next, use bcoeffs(j-1) and b(test)
// on cycle(test), and compute outputs on cycle(test) as y0(j) + values from bcoeffs(j-1).
// The latter are computed by letting delta equal y1(j)-y0(j) and computing the cubic delta bcoeffs,
// and adding to bcoeffs(test).  When min value b(min) for b(test) error is found, set b(j) = b(min).
// Then do regular B-spline fit on cycle j with this b=b(min).
float MainContentComponent::bWithMinError(CycleSpline& cycle1, CycleSpline& cycle2, int p)
{
    float error = 0, min = 100, minIndex = 0;
    float b0 = cycle2.b - (float)p;
    float b = b0;
    float a = cycle2.a;
    float k = cycle2.k;
    CycleSpline tempCycle2 = CycleSpline(k, a, b);
    tempCycle2.y0 = cycle2.y0;
    for (int i=0; i<2*p; i++) {
        b = b0 + (float)i;
        tempCycle2.b = b;
        tempCycle2.y1 = interpFloat(tempCycle2.b, floatBuffer);
//        addDeltaToCycleBcoeffs(tempCycle2);
        computeCycleSplineOutputs(tempCycle2); // this uses y0 and y1
        error = compareSlopeError(cycle1, tempCycle2, floatBuffer);
        if (error < min) {
            min = error;
            minIndex = i;
        }
//        DBG("error for offset " << i-p << " and b = " << b << " :  " << error);
    }
//    DBG("b value offset: " << minIndex - p);
    b = b0 + (float)minIndex;
    return b;
}

// Next compute discrete first derivative and measure error which should detect the anomoly of cycle shape discontinuity.
// Measure the discrete derivative at a sample by the slope of the line through the next sample.  For a cycle spline these
// are just the outputs. This treats cycle1 and cycle2 as read-only.
float MainContentComponent::compareSlopeError(CycleSpline& cycle1,  CycleSpline& cycle2, AudioBuffer<float>& samples)
{
//    float slopeWeight = 100;
    float y1Weight = 100;
    float derivWeight = 100;
    float output = 0, sample1 = 0, sample2 = 0, output1 = 0, output2 = 0, val = 0;
    float sampleSlope = 0, splineSlope = 0;
    float a = cycle2.a;
    float b = cycle2.b;
    float y0 = cycle2.y0;
    float y1 = cycle2.y1;
    int k = cycle2.k;
    int d = cycle2.d, n = k+d;
    int j = (int) a + 1;
    CycleSpline tempCycle2 = CycleSpline(k, a, b);
    tempCycle2.y0 = y0;
    tempCycle2.y1 = y1;
    // now set bcoeffs of tempCycle2 to those of cycle1 and compute outputs
    for (int p=0; p<n; p++) {
        tempCycle2.bcoeffs.set(p, cycle1.bcoeffs[p]);
    }
    computeCycleSplineOutputs(tempCycle2);  // this uses y0 and y1 with delta cubic
    int N = tempCycle2.outputs.size();
    // compute error for tempCycle2 using cycle1 bcoeffs
    for (int i=0; i<N-1; i++) {
        sample1 = samples.getSample(0, i+j);
        sample2 = samples.getSample(0, i+j+1);
        sampleSlope = sample2 - sample1;
        output1 = tempCycle2.outputs[i];
        output2 = tempCycle2.outputs[i+1];
        splineSlope = output2 - output1;
        val = derivWeight * (sampleSlope - splineSlope);
        output += val * val;
    }
    output = output / (float) N;

    // Also trying another error measure:  difference between first and last slopes.  Need to keep this value small
    // in order to avoid the "drift" of left endpoint b past critical points in the cycle shape.  This requires
    // the computation of the new cycle spline fit, or the use of samples of the audio data.  We will try the latter
    // first, to save computation time, and simply estimate the slope at the first and last samples.

    // This one didn't work so well, so leaving it out now.
//    float slopea = samples.getSample(0, j+1) - samples.getSample(0, j);
//    float slopeb = samples.getSample(0, j+N-1) - samples.getSample(0, j+N-2);
//    float diff = abs(slopeb - slopea) * slopeWeight;
//    output += diff;
    
    // Another measure:  the distance of y1 from 0, possibly to pull the cycle back closer to zeros.
    // This one works quite well ...
    output += y1Weight * y1 * y1;
    return output;
}


// takes input r for elementary CA(r), and input array of initial state, outputs final state
Array<int> MainContentComponent::computeCA(int r, Array<int> input) {
    int n = input.size();
    Array<int> output;
    for (int j=0; j<n; j++) {
        output.add(0);
    }
    int MASK = 1;
    MASK <<= 8;
    Array<int> local;   // array for local rule, reads left to right as binary digits of r
    for (int i=0; i<8; i++) {
        local.add(0);
        MASK >>= 1;
        if (r & MASK) {
            local.set(i,1);
        }
    }
    // for each element of input, do local computation for CA(r)
    int a=0, b=0, c=0;
    for (int k=0; k<n; k++) {
        if (k == 0) {
            b = input[0]; c = input[1];
        } else {
            if (k == n-1) {
                a = input[n-2]; b = input[n-1];
            } else {
                a = input[k-1]; b = input[k]; c = input[k+1];
            }
        }
        int sum = a+b+c;
        if (sum == 3) {
            output.set(k, local[0]);
        }
        if (sum == 0) {
            output.set(k, local[7]);
        }
        if (sum == 1) {
            if (a == 1) {
                output.set(k, local[3]);
            }
            if (b == 1) {
                output.set(k, local[5]);
            }
            if (c == 1) {
                output.set(k, local[6]);
            }
        }
        if (sum == 2) {
            if (c == 0) {
                output.set(k, local[1]);
            }
            if (b == 0) {
                output.set(k, local[2]);
            }
            if (a == 0) {
                output.set(k, local[4]);
            }
        }
        a = 0; b = 0; c = 0;
    }
    return output;
}

void MainContentComponent::pushNextSampleIntoFifo (float sample) noexcept
{
    // if the fifo contains enough data, set a flag to say
    // that the next line should now be rendered..
    if (fifoIndex == fftSize)
    {
        if (! nextFFTBlockReady)
        {
            std::fill (fftData.begin(), fftData.end(), 0.0f);
            std::copy (fifo.begin(), fifo.end(), fftData.begin());
            nextFFTBlockReady = true;
        }
        fifoIndex = 0;
    }
    fifo[(size_t) fifoIndex++] = sample;
}

void MainContentComponent::fftError(AudioBuffer<float>& samples)
{
    // fill fifo (float buffer of size fftSize) with some samples and compute fft
    
}

// observations with error weights:  currently derivWeight = y1Weight = 100 (ignoring slopeWeight)
// which gives a good output for A450.wav and piano1.wav, handling the shape continuity quite well.
// Next need to combine these models with cycle interpolation.  The new approach should be to work
// with the fully computed model, with each cycle stored in keyCycleArray, then choose from these
// already computed key cycles and do interpolation quickly for new model.  It is also a good time
// to store, or serialize these models, for fast retrieval.

// compute model with lift delta on each cycle
void MainContentComponent::computeModelWithDelta()
{
    // Model without cycle interp but use delta and cubic bcoeffs for lift parameter.
    // The idea is to avoid shape discontinuity by maintaining similarity from one cycle to the next.
    // The zero crossings are no longer used, but rather recompute end points with the period length guess,
    // followed by adjustment of right endpoint based on minimizing some error functions.
    // We start by using the zero crossings and compute error comparing cycle j-1 to j,
    // then use average error as a reference, and when this is exceeded we introduce delta.
    // In addition to comparison with cycle j-1 values and derivatives, we also use y1*y1.
    // Each cycle will have some new parameters, in particular y0 and y1, which are the signal
    // values at the end points a and b.  Also, use new knot sequence, 0,0,0,0,1/k,...,(k-1)/k,1,1,1,1
    
    computeCubicPolyBcoeffs();  // this is q(t) = 3*t^2-2*t^3
//    numCycles -= 3;
    int n = kVal + 3;
    DBG("computing model with delta");
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
//    DBG("buffers set");
    float a = 0, b = 1, error = 0, avgError = 0;
    CycleSpline cycle  = CycleSpline(kVal, a, b);
    CycleSpline cycle1 = CycleSpline(kVal, a, b);
    CycleSpline cycle2 = CycleSpline(kVal, a, b);
    // compute each cycle and write out the computed samples to writeBuffer
    auto* leftBuffer = writeBuffer.getWritePointer (0, 0);
    int currentSampleIndex = 0;
    cycleBreakPoints.set(0, 0);
    cycleErrors.set(0, 0);
//    DBG("first cycle breakpoint and error set");
    int B = 6;
    int counter = 0;
//    initECA();
//    printECA();
    for (int i=0; i<B; i++)  // compute cycles 0 to B-1 based on zero-crossings
    {
        cycleBreakPoints.set(i+1, cycleZeros[i+1]);  // set cycleBreakPoints for i up to B
        // for cycle[6] we need cycleBreakPoints[7] = cycleBreakPoints[6] + samplesPerCycleGuess
        a = cycleZeros[i];
        b = cycleZeros[i+1];
        cycle = CycleSpline(kVal, a, b);
//        DBG("set cycle spline " << i);
        // these have y0=y1=0, so don't use delta
        computeCycleBcoeffsWithDelta(cycle, floatBuffer);
        computeCycleSplineOutputs(cycle);
        int M = cycle.outputs.size();
//        DBG("M:  " << M);
//        DBG("computed outputs for cycle spline " << i);
        for (int j=0; j<M; ++j) {
            leftBuffer[currentSampleIndex] = cycle.outputs[j];
            currentSampleIndex += 1;
        }
        deltaCycleArray.add(cycle);
        // now cycleZeros equal cycleBreakPoints up to i<B
    }
//    DBG("first " << B-1 << " cycles set ...");
    for (int i=1; i<B; i++)  // set errors on cycles 0 to B-1
    {
        if (i < B-1) {
            cycleErrors.set(i, 0);  // set error for cycles 0 through B-2 to be 0
        } else { // i = B-1
            // compute error on cycle B-1 comparing cycle B-1 to B-2
            cycle1 = keyCycleArray[i-1];
            cycle2 = keyCycleArray[i];
            error = compareSlopeError(cycle1, cycle2, floatBuffer);  // error for i=B-1
            avgError = error;
//            DBG("error using previous bcoeffs for cycle " << i << ": " << error);
            cycleErrors.set(i, error);
//            DBG("setting error for cycle " << B-1);
        }
    }
//    for (int i=B; i<numCycles; i++)  // compute cycles > B-1
//    DBG("numCycles = " << numCycles);
    for (int i=B; i<numCycles-1; i++)  // compute cycles > B-1
    {
        DBG("cycle: " << i);
        counter += 1;
        a = cycleBreakPoints[i];
        b = a + samplesPerCycleGuess;   // for piano1.wav this is 155.8 (for freq 283 Hz)
        cycleBreakPoints.set(i+1, b);   // will also use for a on next cycle, unless error is too high
        cycleZeros.set(i+1, b);         // keep cycleZeros in synch with cycleBreakPoints for cycle interpolation
        cycle2 = CycleSpline(kVal, a, b);  // bcoeffs are set to 0's
        cycle1 = deltaCycleArray[i-1];
        cycle2.y0 = cycle1.y1;
        cycle2.y1 = interpFloat(b, floatBuffer);
//        DBG("cycle " << i << "  y1: " << cycle2.y1);
        error = compareSlopeError(cycle1, cycle2, floatBuffer);  // error for current cycle i
        // Next, if (error > 2 * avgError), or some other bound, need to choose
        // various new b values and compute error using previous cycle bcoeffs and
        // the cubic delta function to accomodate the right end point, and choose
        // the new b value, within some number of samples of old b, with minimum error
        float bMin = b;
        if (error > 0.001 * avgError) {
//            DBG("initial error for cycle: " << i << ":  " <<  error);
//            DBG("data for error: a = " << cycle2.a << " b = " << cycle2.b);
//            DBG("computing new b for cycle: " << i);
            bMin = bWithMinError(cycle1, cycle2, 10);
        }
        cycle2.b = bMin;
        b = bMin;
        cycleBreakPoints.set(i+1, b);
        cycleZeros.set(i+1, b);
        // current cycle interval is now set as [a,b] with new b=bMin, next compute spline model:
        cycle2.y1 = interpFloat(cycle2.b, floatBuffer);
        cycle2.y0 = cycle1.y1;
        computeCycleBcoeffsWithDelta(cycle2, floatBuffer);
        
        // change one bcoeff:
//        cycle2.bcoeffs.set(3, cycle2.bcoeffs[3] + 0.03);
//        ------------------------------------------------------------------------
        // modifying bcoeffs with ECA
//        for (int p=0; p<n; p++) {
//            float mult = 0.8 + vectorECA[p] * 0.3;
//            cycle2.bcoeffs.set(p, cycle2.bcoeffs[p] * mult);
//        }
//        if (counter == 5) {
//            vectorECA = computeCA(rVal, vectorECA);
//            printECA();
//            counter = 0;
//        }
//        ------------------------------------------------------------------------
//        if ((i > 10) && (i<50)) {
//            DBG("bcoeff[11] for cycle " << i << ": " << cycle2.bcoeffs[11]);
//        }
//        int n = cycle2.k + 3;
//        for (int q=0; q<n; q++) {
//            if ((q<5) || (q>17)) {
//                cycle2.bcoeffs.set(q, 0);
//            }
//            if ((q>9) && (q<13)) {
//                cycle2.bcoeffs.set(q, 0);
//            }
//        }
//        ------------------------------------------------------------------------
        computeCycleSplineOutputs(cycle2);
//        cycle2.printData();
//        DBG("outputs computed for cycle " << i);
        error = compareSlopeError(cycle1, cycle2, floatBuffer);
        cycleErrors.set(i, error);
        int M = cycle2.outputs.size();
//        DBG("M:  " << M);
        samplesPerCycle.set(i, M);
        for (int j=0; j<M; ++j) {
            leftBuffer[currentSampleIndex] = cycle2.outputs[j];
            currentSampleIndex += 1;
        }
        avgError = 0;
        deltaCycleArray.add(cycle2);
//        if (i>45 && i<50) {
//            DBG("------------");
//            DBG("cycle " << i);
//            DBG("a: " << cycle2.a << " b: " << cycle2.b << " y0: " << cycle2.y0 << " y1: " << cycle2.y1);
//            DBG("first 2 outputs: " << cycle2.outputs[0] << "  " << cycle2.outputs[1]);
//            int N = cycle2.outputs.size();
//            DBG("last  2 outputs: " << cycle2.outputs[N-2] << "  " << cycle2.outputs[N-1]);
//            DBG("------------");
//        }
        for (int j=5; j<i; j++) {
            avgError += cycleErrors[j];
        }
        avgError = avgError / (float)(i-5);
//        DBG("error " << i << ": " << error << "  avgError: " << avgError);
//        DBG("samples per cycle " << i << ": " << samplesPerCycle[i]);
//        int j = (int) a;
//        float slopea = floatBuffer.getSample(0, j+1) - floatBuffer.getSample(0, j);
//        float slopeb = floatBuffer.getSample(0, j+M-1) - floatBuffer.getSample(0, j+M-2);
//        float diff = abs(slopeb - slopea) * 100;
//        DBG("slope at b minus slope at a:  " << diff);
    }
    graphView.setBreakPointsForGraph(cycleBreakPoints, samplesPerCycle);
    graphView.setDeltaCycleArray(deltaCycleArray);
    DBG("finished computing model with delta");
}

// compute spline model without interpolation, write to buffer, and store in allCycleArray
void MainContentComponent::modelWithoutCycleInterp()
{
    DBG("computing model without cycle interp");
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    float a = 0, b = 1;
//    float error1 = 0, error2 = 0, error3 = 0;
    CycleSpline cycle = CycleSpline(kVal, a, b);
    // compute each cycle and write out the computed samples to writeBuffer
    auto* leftBuffer = writeBuffer.getWritePointer (0, 0);
    int currentSampleIndex = 0;

    for (int i=0; i<numCycles-1; i++)  // for (int i=0; i<10; i++)
    {
        a = cycleZeros[i];
        b = cycleZeros[i+1];
        cycle = CycleSpline(kVal, a, b);
//        computeCycleBcoeffs(cycle, floatBuffer);
        setNewTargets(cycle);
        computeNewBcoeffs(cycle);
        computeCycleSplineOutputs(cycle);
//        if (i > 105 && i < 112) {
//            error1 = computeCycleSplineError(cycle, floatBuffer);
//            DBG("error1 for cycle " << i << ":  " << error1);
//            CycleSpline cycle1 = keyCycleArray.getLast();
//            error2 = compareCycleSplineError(cycle1, cycle, floatBuffer);
//            DBG("error2 for cycle " << i << ":  " << error2);
//            error3 = compareSlopeError(cycle1, cycle, floatBuffer);
//            DBG("error3 for cycle " << i << ":  " << error3);
//        }
        int M = cycle.outputs.size();
        for (int j=0; j<M; ++j) {
            leftBuffer[currentSampleIndex] = cycle.outputs[j];
            currentSampleIndex += 1;
        }
        allCycleArray.add(cycle);
//        DBG("allCycleArray.size: " << allCycleArray.size());
    }
}

void MainContentComponent::setKeyCycleIndices() {
    if (otherCycleInterp) {
            // leave keys to reuse and modify
        } else {
            keys.clear();
        }
    if (expCycleInterp) {
//        DBG("setting key cycles for expCycleInterp");
//        keys.add(0); keys.add(1);
//        int key = 1, j = 1;
//        while (key < numCycles) {
//            key = keys[j] * 2;
//            if (key < numCycles) {
//                keys.add(key);
//                j++;
//            }
//        }
        // resetting for special cases
        // this is first instrument case with 1 sec audio, 18 key cycles:
        // 0,5,10,15,20,25,30,40,50,60,70,80,100,120,150,180,220,last
        // the above works for middle C 262 Hz, but not for lower than 220
        // so if f_0 is less than 220, change to:
        // 0,2,4,6,8,10,15,20,25,30,40,50,60,70,80,100,120,last
        if (freqGuess < 220) {
//        if (freqGuess < 300) {
            keys.add(0); keys.add(2); keys.add(4); keys.add(6); keys.add(8); keys.add(10);
            keys.add(15); keys.add(20); keys.add(25); keys.add(30); keys.add(40); keys.add(50);
            keys.add(60); keys.add(70); keys.add(80); keys.add(100); keys.add(120);
        } else {
            keys.add(0); keys.add(5); keys.add(10); keys.add(15); keys.add(20); keys.add(25);
            keys.add(30); keys.add(40); keys.add(50); keys.add(60); keys.add(70); keys.add(80);
            keys.add(100); keys.add(120); keys.add(150); keys.add(180); keys.add(220);
        }
//        int i = 0;
//        int m = 5;
//        while (m * i < 31) {
//            keys.add(m * i);
//            i++;
//        }
//        i = 4;
//        m = 10;
//        while (m * i < 81) {
//            keys.add(m * i);
//            i++;
//        }
//        i = 5;
//        m = 20;
//        while (m * i < 121) {
//            keys.add(m * i);
//            i++;
//        }
//        keys.add(150);
//        keys.add(180);
//        keys.add(220);
    }
    if (endsOnlyCycleInterp) {
        DBG("setting key cycles for endsOnlyCycleInterp");
        keys.add(0);
    }
    if (fibonacciCycleInterp) {
        DBG("setting key cycles for fibonacciCycleInterp");
        keys.add(0); keys.add(1); keys.add(2);
        int key = 2, j = 2;
        while (key < numCycles) {
            key = keys[j] + keys[j-1];
//            if (key < numCycles) {
            // special case for MCM paper:
//            DBG("Note: Fibonacci Cycle Interp is using special case for MCM paper");
//            if (key < 234) {
//                keys.add(key);
//                j++;
//            }
            if (key < numCycles) {
                keys.add(key);
                j++;
            }
        }
    }
    if (regularCycleInterp) {
        DBG("setting key cycles for regularCycleInterp");
        int i = 0;
        while (mVal * i < numCycles) {
            keys.add(mVal * i);
            i++;
        }
    }
    if (keys[keys.size()-1] < numCycles-1) {
        if (expCycleInterp) {
            keys.add(int(freqGuess));
        } else {
            keys.add(numCycles-1);
        }
    }
    if (keysCA > 0) {   // need to restrict number of keys to be < keysCA
        Array<int> temp;
        for (int i=0; i<keys.size(); i++) {
            temp.add(keys[i]);
        }
        if (modelWithDelta) {
            // leave keys as is, in effect for instruments
        } else { // limit keys
            keys.clear();
            for (int j=0; j<keysCA; j++) {
                keys.add(temp[j]);
            }
        }
    }
    bool modified = false;
    int Add = graphView.keysToAdd.size();
    if (Add > 0) {
        for (int i=0; i<Add; i++) {
            int newKey = graphView.keysToAdd[i];
            if (!isKey(newKey)) {
                keys.add(newKey);
                modified = true;
            }
        }
    }
    int Remove = graphView.keysToRemove.size();
    if (Remove > 0) {
        for (int i=0; i<Remove; i++) {
            int oldKey = graphView.keysToRemove[i];
            if (isKey(oldKey)) {
                int j = keyIndex(oldKey);
                keys.remove(j);
                modified = true;
            }
        }
    }
    keys.sort();
    if (modified) {
        setInterpSelectionsFalse();
        otherCycleInterp = true;
        interpSelector.setSelectedId (6);
    }
    graphView.keysToAdd.clear();
    graphView.keysToRemove.clear();
    for (int k=0; k<keys.size(); k++) {
        DBG("keys[" << k << "]: " << keys[k]);
    }
    DBG("there are " << keys.size() << " key cycles");
}


// not using this function, refactored into code above
void MainContentComponent::nextCAcycle(int r)
{
    int n = kVal + 3;  // dimension of splineCycles and size of vectorECA
    computeCA(r, vectorECA);
    for (int i=0; i<n; i++) {
        cycleCA.bcoeffs.set(i, cycleCA.bcoeffs[i] * vectorECA[i]);
    }
}


// compute key cycle bcoeffs, write to keyBcoeffs array
// but do not compute outputs and do not write to buffer
void MainContentComponent::computeKeyCycles()
{
    // compute key cycles and add to array keyCycleArray
    DBG("computing (only) key cycles");
    keyCycleArray.clear();
    writeBuffer.clear();
//    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    DBG("size of writeBuffer: " << lastSample + 1);
    float a = 0, b = 1;
    CycleSpline cycle = CycleSpline(kVal, a, b);
    int i = 0;
    int k = kVal;
    DBG("kVal = " << kVal);
    int d = dVal;
    int n = k + d;  // dim of cycle splines
    // compute key cycles and write bcoeffs to keyBcoeffs array
    for (i=0; i<keys.size(); i++)
    {
        a = cycleZeros[keys[i]];
        b = cycleZeros[keys[i]+1];
        cycle = CycleSpline(kVal, a, b);
        // compute each key cycle bcoeffs and store in array of floats
        setNewTargets(cycle);
        computeNewBcoeffs(cycle);
        for (int p=0; p<n; p++) {
            keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
        }
        keyCycleArray.add(cycle);
    }
}

// We now use this function to write key cycles which are taken from a delta model.
// In this case, the delta model has already been computed for all cycles, allowing cycle endpoints a and b
// to move away from zeros with y values y0 and y1.  All delta cycles are now stored in the array deltaCycleArray.

// same as above compute bcoeffs of key cycles,
// but also compute outputs and write to buffer
void MainContentComponent::writeKeyCyclesToBuffer()
{
    DBG("computing and writing key cycles");
    keyCycleArray.clear();
//    allCycleArray.clear();
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    float a = 0, b = 1;
    CycleSpline cycle = CycleSpline(kVal, a, b);
    int k = kVal;
    int d = dVal;
    int n = k + d;  // dim of cycle splines
//    if (addSubharmonic) {
//        cycleParabola = CycleSpline(kVal, a, b);
//        setParabolicTargets(0.005);
//        computeParabolicBcoeffs();
//    }
    // compute key cycles and write outputs to buffer
    for (int i=0; i<keys.size(); i++)
    {
        if (modelWithDelta || modelWithoutDelta) {
            // key cycles are now chosen with indices from full set of cycles
            // if using delta model key cycles also have y0, y1
            if (modelWithDelta) {
                cycle = deltaCycleArray[keys[i]];
            } else {
                cycle = allCycleArray[keys[i]];
            }
            if (preserveDelta) {
                // leave key cycles with y0, y1
            } else {
                cycle.y0 = 0;
                cycle.y1 = 0;
            }
            a = cycle.a;
            b = cycle.b;
            keyCycleArray.add(cycle);  // comes from delta model so includes y0, y1
        } else {
            a = cycleZeros[keys[i]];
            b = cycleZeros[keys[i]+1];
            cycle = CycleSpline(kVal, a, b);
            // compute each key cycle bcoeffs and store in array of floats
            setNewTargets(cycle);
            computeNewBcoeffs(cycle);
            computeCycleSplineOutputs(cycle);
            keyCycleArray.add(cycle);
        }
        for (int p=0; p<n; p++) {
            float val = cycle.bcoeffs[p];
            keyBcoeffs.set(i*n+p, val);
//            if (addSubharmonic) {
//                if (i % 2 == 0) {
//                    val += cycleParabola.bcoeffs[p];
//                } else {
//                    val -= cycleParabola.bcoeffs[p];
//                }
//                cycle.bcoeffs.set(p, val);
//            }
        }
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
//        DBG("M:  " << M << "  L:  " << L);
        for (int j=0; j<M; j++) {
            float value = cycle.outputs[j];
            // adding subharmonic
//            float t = (float) (a+j+L);
//            if (keys[i] % 2 == 0) {
//                value += 0.01 * std::sinf(Pi * (t-a)/(b-a));
//            } else {
//                value -= 0.01 * std::sinf(Pi * (t-a)/(b-a));
//            }
            writeBuffer.setSample(0, j+L, value);
        }
    }
}

// this sets cycle to use new knot sequence and appropriate targets
void MainContentComponent::setNewTargets(CycleSpline& cycle)
{
    int k = cycle.k;
    int d = cycle.d;
    int n = k + d;
    int N = n + d;
    cycle.inputs.clear();
    cycle.targets.clear();
    cycle.knots.clear();
    cycle.bcoeffs.clear();
//    DBG("printData after clear:");
//    cycle.printData();
    
    float incr = 1 / (float) k;
    cycle.inputs.add(0.5*incr);
    for (int i=1; i<k; i++) {
        cycle.inputs.add(i*incr);
    }
    cycle.inputs.add(1-0.5*incr);
    
    // New knot sequence: 0,0,0,0,1/k,2/k,...,(k-1)/k,1,1,1,1
    for (int i=0; i<d+1; i++) {
        cycle.knots.add(0);
    }
    for (int i=d+1; i<N-d; i++) {
        cycle.knots.add(cycle.knots[i-1]+incr);
    }
    for (int i=0; i<d+1; i++) {
        cycle.knots.add(1);
    }
    
    // set bcoeffs to 0 then will compute n-2 of them c_1,...,c_{n-1}
    // in function computeNewBcoeffs
    for (int i=0; i<n; i++) {
        cycle.bcoeffs.add(0);
    }
    
    float a = cycle.a;
    float b = cycle.b;
    for (int i=0; i<n-2; i++) {
        float t = cycle.inputs[i];
        float output = interpFloat(a + t * (b-a), floatBuffer);
        cycle.targets.set(i, output);
    }
    
//    DBG("New params set for CycleSpline knot sequence: ");
//    DBG("N = " << N << " n = " << n << " k = " << k);
//    for (int i=0; i<N+1; i++) {
//        std::cout << cycle.knots[i] << " ";
//    }
//    DBG("");
}

// uses new knots, inputs, and targets
void MainContentComponent::computeNewBcoeffs(CycleSpline& cycle)
{
    // assume bcoeffs[0] = bcoeffs[n-1] = 0, and compute the other n-2
    // this means B-spline graph will hit target 0 at the ends
    // also bcoeffs[1] and bcoeffs[n-2] control the derivatives at the ends
    int k = cycle.k, d = cycle.d;
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
            val = newBsplineVal(k, j+1, cycle.inputs[i]);  // A[i,j]
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
    y = multMatCol(n-2, B, cycle.targets);
    for (int i=1; i<n-1; i++) {
        cycle.bcoeffs.set(i, y[i-1]);
    }
}

void MainContentComponent::loadSmallModel()
{
    setInterpSelectionsFalse();
    expCycleInterp = true;
    CAmodelButtonClicked();
}

void MainContentComponent::computeMetaSplines()
{
    int n = kVal + dVal;
    if (modelWithSmall) {
        // don't need Delta Model, but should write something to keyBcoeffs array
    } else {
        auto file = File::getSpecialLocation(File::userHomeDirectory).getChildFile("keyBcoeffs");
        FileOutputStream output (file);
        if (output.openedOk())
        {
            output.setPosition (0);  // default would append
            output.truncate();
        }
        if (modelWithDelta) {
            // write key cycle bcoeffs from allCycleArray to keyBcoeffs array
            DBG("writing keyBcoeffs from deltaCycleArray");
            for (int i=0; i<keys.size(); i++) {
                CycleSpline cycle = deltaCycleArray[keys[i]];
                for (int p=0; p<n; p++) {
                    keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
                    output.writeFloat(keyBcoeffs[i*n+p]);
                }
            }
        }
    }
    // this is doing piecewise linear metasplines for now, so bcoeffs for metasplines are skipped
    DBG("computing metasplines");
    // compute metasplines based on distribution of keycycles
    // targets of metasplines are keycycle bcoeffs
    metaSplineArray.clear();
    MetaSpline spline = MetaSpline(keys.size(), mVal);
//    n = kVal + dVal;  // dim of cycle splines, = # bcoeffs for each cycle spline
//    int size = keys[keys.size()-1]+1;  // size = numCycles when modelWithCA
    for (int i=0; i<n; i++) {  // loop on each bcoeff of cycles with same fixed kVal
        spline = MetaSpline(keys.size(), mVal, numCycles); // numCycles = outputs.size()
        if (modelWithCA) {
            spline = MetaSpline(keys.size(), mVal, numCycles);
        }
        for (int j=0; j<keys.size(); j++) {  // j loop on key cycles
            if (regularCycleInterp) {
                // already set uniform inputs, so nothing to do
            } else {
                // put metaspline inputs at values keys/(numCycles-1)
                if (modelWithCA) {
                    spline.inputs.set(j, (float)keys[j] / (float)(numCycles-1));
                } else {
                    spline.inputs.set(j, (float)keys[j] / (float)(numCycles-1));
                }
            }
            spline.targets.set(j, keyBcoeffs[j*n+i]);
        }
//        computeMetaSplineBcoeffs(spline);
        computeLinearMetaSplineOutputs(spline);
//        spline.printData();
        metaSplineArray.add(spline);
    }  // now we have n=k+d meta-splines, one for each of n bcoeffs
}

// not using this function, did this in computeMetaSplines
void MainContentComponent::writeKeyBcoeffsToFile()
{
    int n = kVal + dVal;
    auto file = File::getSpecialLocation(File::userHomeDirectory).getChildFile("keyBcoeffs");
    FileOutputStream output (file);
    if (output.openedOk())
    {
        output.setPosition (0);  // default would append
        output.truncate();
        for (int i=0; i<n; i++) {
//            output.writeFloat(bcoeffs[i]);
        }
    }

}

void MainContentComponent::modelWithCAplusCycleInterp()
{
    DBG("computing CA model");
    if (modelWithDelta) {
        // use key cycles from deltaCycleArray
    } else {
        generateKeyCyclesCA();
    }
    computeMetaSplines();
    writeNormalizedCyclesToBuffer();
    graphView.graphCAmodel = true;
    graphView.repaint();
}

// Need to fix constant cycle length, but can first generate key cycles with a=0, b=1
// then when filling out all cycles can put the correct a, b for constant cycle length.
// Also, we can decouple this CA model from number of cycles in audio file sample.
// CycleSpline(kVal, a, b) uses new knot sequence 0,0,0,0,1/k,...,1,1,1,1 and inputs.
// Added May '22:  keysCA is now bound on number of keys, and last key is just splinusoid
void MainContentComponent::generateKeyCyclesCA()
{
    keyCycleArray.clear();
//    writeBuffer.clear();
//    writeBuffer.setSize(1, sampleRate * lengthInSecondsCA);  // channels, samples
//    DBG("writeBuffer.size(): " << writeBuffer.getNumSamples());
    float a = 0, b = 1;
    CycleSpline cycle = CycleSpline(kVal, a, b);  // uses new knot sequence
    int i = 0;
    int k = kVal;
    int d = dVal;
    int n = k + d;  // dim of cycle splines
    bool usingSplinusoid = false;
    // add first key cycle as graphView.cycleNew or as cycle 10
    if (audioDataLoaded) {
        if (graphView.highlightCycle > -1) {
            cycleCA = graphView.cycleNew;  // this uses new knot sequence
            DBG("cycleCA:");
            cycleCA.printData();
        } else {
            DBG("no cycle spline selected, using cycle 10");
            // use audio graph key cycle 10
            a = cycleZeros[10];
            b = cycleZeros[11];
            cycleCA = CycleSpline(kVal, a, b);
            setNewTargets(cycleCA);
            computeNewBcoeffs(cycleCA);
        }
    } else {
        computeSplinusoidBcoeffs();
        usingSplinusoid = true;
        initECA();
    }
    for (int p=0; p<n; p++) {
        float mult = 0.3 + 0.9 * vectorECA[p];
        // ignore with mult = 1
        mult = 1;
        if (usingSplinusoid) {
            cycle.bcoeffs.set(p, splinusoid.bcoeffs[p] * mult);
            keyBcoeffs.set(i*n+p, splinusoid.bcoeffs[p]);
        } else {
            cycle.bcoeffs.set(p, cycleCA.bcoeffs[p] * mult);
            keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
        }
    }
    keyCycleArray.add(cycle);  // this is selected cycle now with initial vectorECA applied
    // now apply vectorCA to get new key cycles
    for (i=1; i<keys.size(); i++)
    {
        rVal = 30;
        vectorECA = computeCA(rVal, vectorECA);
        for (int p=0; p<n; p++) {
            float r = juce::Random::getSystemRandom().nextFloat();
            r = 2*r-1;  // random float in [-1,1]
            r *= 0.5;   // random float in [-0.3,0.3]
            float mult = 1 + r * vectorECA[p];
            // ignore r
            mult = 2.5 - vectorECA[p];
            // trying plus/minus
//            mult = 5 * r;
//            if (vectorECA[p] < 1) {
//                mult = -5 * r;
//            }
            // ignore all CA with mult = 1
            mult = 1;
            if (i == keys.size()-1) {
                mult = 1;   // this sets last key cycle to splinusoid or cycleCA
            }
            if (usingSplinusoid) {
                cycle.bcoeffs.set(p, splinusoid.bcoeffs[p] * mult);
            } else {
                cycle.bcoeffs.set(p, cycleCA.bcoeffs[p] * mult);
            }
            keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
        }
        keyCycleArray.add(cycle);  // this is selected cycle now with iterated vectorECA applied
    }
}

void MainContentComponent::initECA()
{
    int k = kVal, d = dVal;
    int n = k + d;
    int mid = (int) (n - 1) / 2;   // if n is even, this is not in the middle
//    int quart = (int) mid / 2;
    // initialize vectorECA for ECA iteration
    vectorECA.clear();
    for (int i=0; i<n; i++) {
        vectorECA.add(0);
    }
    vectorECA.set(31, 1);
//    vectorECA.set(mid, 1);
//    vectorECA.set(quart, 1);
//    vectorECA.set(mid + quart, 1);
}

void MainContentComponent::printECA()
{
    int k = kVal, d = dVal, val = 0;
    int n = k + d;
    for (int j=0; j<n; j++) {
        val = vectorECA[j];
        if (val == 0) {
            cout << " ";
        }
        if (val == 1) {
            cout << "1";
        }
    }
    std::cout << endl;
}


void MainContentComponent::writeNormalizedCyclesToBuffer()
{
    // write all cycles using constant cycle length averageSamplesPerCycle
    DBG("computing and writing normalized cycles");
    DBG("freqGuess: " << freqGuess);
    samplesPerCycleCA = (int)(float(sampleRate) / freqGuess);
    DBG("samplesPerCycleCA: " << samplesPerCycleCA);
    writeBuffer.setSize(1, samplesPerCycleCA * numCycles + 1);  // channels, samples
    // choose number of cycles for Attack part A to be keysCA * mVal
    // for example, 5*5=25 works well for A450
    // Note: we are not using keysCA other than to compute A if modelWithDelta is true
    keysCA = 5;
    float A = (keysCA * mVal) * samplesPerCycleCA;
    const float D = A / 2;
    const float S = (numCycles * samplesPerCycleCA - A - D) / 4 ;
    const float R = 3 * S;
    DBG("A: " << A << "  D: " << D << "  S: " << S << "  R: " << R);
    int i = 0;
    float a = 0, b=1, ratio = 1;
    CycleSpline cycle;
    int n = kVal + dVal;
    float total = (float) numCycles;
//    DBG("numCycles: " << numCycles);
//    DBG("key cycles:");
//    for (int h=0; h<keys.size(); h++) {
//        DBG("keys[" << h << "]: " << keys[h]);
//    }
//    DBG("allCycleArray.size:  " << allCycleArray.size());
    while (i < numCycles) {  // i is cycle number whether key or non-key
        DBG("i: " << i);
        // the next block forces each cycle to have averageSamplesPerCycle samples
        // as computed by computeCycleSplineOutputs(cycle) using only a and b
        if (i == 0) {
            a = 0;
            b = (float)((int)samplesPerCycleCA) - 0.5;
        } else {
            a = b;
            b = a + (int)samplesPerCycleCA;
        }
        if (isKey(i)) {   // use key cycle i
            int j = keyIndex(i);
            if (modelWithSmall) {
                cycle = CycleSpline(kVal, a, b);  // this is key cycle j
                for (int p=0; p<n; p++) {
                    cycle.bcoeffs.set(p, keyBcoeffs[j*n+p]);
                }
            } else {
                if (modelWithDelta) {
                    cycle = deltaCycleArray[i];
                } else {
                    cycle = keyCycleArray[j];
                }
            }
            cycle.a = a; cycle.b = b;   // reset to normalized a and b
            cycle.y0 = 0; cycle.y1 = 0; // reset to avoid Delta floating spline
//            DBG("writing key cycle: " << i);
        } else {
            // set bcoeffs for (non-key) cycle i using meta-spline when i < index of last key cycle
            if (i < keys[keys.size()-1]) {
                cycle = CycleSpline(kVal, a, b);  // this is cycle i
                for (int j=0; j<n; j++) {
//                    DBG("getting bcoeffs from metasplines with j = " << j);
                    cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
                }
            } else {  // will not use this block if modelWithDelta true
                // use last key cycle
                cycle = keyCycleArray[keyCycleArray.size()-1];
//                cycle = CycleSpline(kVal, a, b);
//                for (int j=0; j<n; j++) {
                    // use splinusoid
//                     cycle.bcoeffs.set(j, splinusoid.bcoeffs[j]);
//                }
            }
        }
        // changing each cycle to have bcoeffs set here:
//        for (int i=0; i<n; i++) {
//            float y = (float)i;
//            cycle.bcoeffs.set(i, 0);
//            if (i > 5) {
//                y *= 0.05;
//                cycle.bcoeffs.set(i, y);
//            }
//            if (i > 15) {
//                y *= -1;
//                cycle.bcoeffs.set(i, y);
//            }
//            if (i == 20) {
//                cycle.bcoeffs.set(i, 0);
//            }
//        }
//        cycle.bcoeffs.set(16, splinusoid.bcoeffs[16] - 0.5);
        computeCycleSplineOutputs(cycle);
        allCycleArray.add(cycle);
        // for now both randomizeBcoeffs and useEnvelope are being ignored
        if (randomizeBcoeffs) {
            ratio = 1;
            if (i > 10) {
                ratio = (total + 10 - i) / total;
            }
        } else {
            if (useEnvelope) {
                if (cycle.maxVal > 0) {
                    ratio = maxSampleValues[i] / cycle.maxVal;
                } else {
                    DBG("cycle [" << i << "] has maxVal 0");
                }
            }
        }
        ratio = 1;  // this overrides the previous block
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
        DBG("outputs.size: " << M);
        DBG("first sample: " << L);
        if (i == 0) {
            L = 0;
        }
        for (int j=0; j<M; j++) {
            if (useADSR) {
                ratio = 0.5;  // Sustain value 1/2
                if (j+L < A) {
                    ratio = (1/A) * (j+L);  // Attack value goes 0 to 1
                } else {
                    if (j+L < A+D) {
                        ratio = 1 - (1/(2*D))*(j+L-A);  // Decay value, 1 to 1/2
                    } else {
                        if (j+L > A+D+S) {
                            ratio = -(1/(2*R))*(j+L-(A+D+S+R));  // Release value, 1/2 to 0
                        }
                    }
                }
                // here we restrict envelope only to the release part
                ratio = 1;
                if (j+L > A+D+S) {
                    ratio = -(1/R)*(j+L-(A+D+S+R));  // Release ratio, 1 to 0
                }
            }
            writeBuffer.setSample(0, j+L, ratio * cycle.outputs[j]);
        }
        i += 1;
    }
}

void MainContentComponent::writeNonKeyCyclesToBuffer()
{
    allCycleArray.clear();
    DBG("computing and writing non-key cycles");
    int i = 0;
    float a = 0, b=1, y0 = 0, y1 = 0, ratio = 1;
    CycleSpline cycle;
//    CycleSpline secondLastKeyCycle = keyCycleArray[keyCycleArray.size()-2];
    int n = kVal + dVal;
    if (addSubharmonic) {
        cycleParabola = CycleSpline(kVal, a, b);
        setParabolicTargets(0.005);
        computeParabolicBcoeffs();
    }
    while (i < numCycles) {
        if (isKey(i)) {
            // key cycles already written
            int j = keyIndex(i);
            cycle = keyCycleArray[j];
            allCycleArray.add(cycle);
            i += 1;
        } else {
            // set bcoeffs for cycle i using meta-spline and write outputs
            // cycle i starting sample is (int) cycle.a + 1
            a = cycleZeros[i];
            b = cycleZeros[i+1];
            if (preserveDelta) {
                y0 = deltaCycleArray[i].y0;
                y1 = deltaCycleArray[i].y1;
            } else {
                y0 = 0;
                y1 = 0;
            }
    // here we set the cycle length to be the same as second last key cycle
//            int Tail = keys[keys.size()-2];  // secondLastKeyCycleIndex
//            float bTail = cycleZeros[Tail+1];
//            if (i > Tail) {
//                a = bTail + (i-Tail-1) * 200;
//                b = a + 200;
//            }
    // end setting normalized cycle length for tail
            cycle = CycleSpline(kVal, a, b);  // this is cycle i
            cycle.y0 = y0;
            cycle.y1 = y1;
            // here use only the second last key cycle for the tail:
//            for (int j=0; j<n; j++) {
//                if (i < keys[keys.size()-2]) {
//                    cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
//                } else {
//                    cycle.bcoeffs.set(j, secondLastKeyCycle.bcoeffs[j]);
//                }
//            }
            for (int j=0; j<n; j++) {
                cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
                if (addSubharmonic) {
                    float val = cycle.bcoeffs[j];
                    if (i % 2 == 0) {
                        val += cycleParabola.bcoeffs[j];
                    } else {
                        val -= cycleParabola.bcoeffs[j];
                    }
                    cycle.bcoeffs.set(j, val);
                }
            }
            computeCycleSplineOutputs(cycle);
            if (useEnvelope) {
                if (cycle.maxVal > 0) {
                    ratio = maxSampleValues[i] / cycle.maxVal;
                } else {
                    DBG("cycle [" << i << "] has maxVal 0");
                }
            }
            allCycleArray.set(i, cycle);
            // DBG("cycle number:  " << i);
            ratio = 1;  // to ignore ratio factor for cycles
            int M = cycle.outputs.size();
            int L = ((int)a) + 1;  // first sample index to write
            for (int j=0; j<M; j++) {
                float value = ratio * cycle.outputs[j];
                // adding subharmonic
//                float t = (float) (a+j+L);
//                if (i % 2 == 0) {
//                    value += 0.01 * std::sinf(Pi * (t-a)/(b-a));
//                } else {
//                    value -= 0.01 * std::sinf(Pi * (t-a)/(b-a));
//                }
                writeBuffer.setSample(0, j+L, value);
            }
            i += 1;
        }
    }
}


void MainContentComponent::computeCycleWithEnv()
{
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    int i = 0;
    CycleSpline cycle;
    int n = kVal + dVal;
    float a = 0, b = 1, ratio = 1;
    while (i < numCycles) {
        if (normalizeCycleLength) {
            if (i == 0) {
                a = 0;
                b = (float)((int)samplesPerCycleGuess) + 0.5;
            } else {
                a = b;
                b = a + (int)samplesPerCycleGuess;
            }
        } else {
            a = cycleZeros[i];
            b = cycleZeros[i+1];
        }
        cycle = CycleSpline(kVal, a, b);  // this is cycle i
        for (int j=0; j<n; j++) {
            cycle.bcoeffs.set(j, graphView.cycleToGraph.bcoeffs[j]);
        }
        computeCycleSplineOutputs(cycle);
        if (cycle.maxVal > 0) {
            ratio = maxSampleValues[i] / cycle.maxVal;
        } else {
            DBG("cycle [" << i << "] has maxVal 0");
        }
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
        for (int j=0; j<M; j++) {
            writeBuffer.setSample(0, j+L, ratio * cycle.outputs[j]);
        }
        i += 1;
    }
}

// compute all cycles from model and adjust max values to match audio data
void MainContentComponent::computeAllCycles()
{
    DBG("computing all cycles");
    int i = 0;
    float a = 0, b = 1, ratio = 1;
    CycleSpline cycle;
    int n = kVal + dVal;
    while (i < numCycles) {
        // compute bcoeffs for cycle i from meta-spline and write outputs
        // cycle i starting sample is (int) cycle.a + 1
        a = cycleZeros[keys[i]];
        b = cycleZeros[keys[i]+1];
        cycle = CycleSpline(kVal, a, b);  // this is cycle i
        for (int j=0; j<n; j++) {
            cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
        }
        // next function now includes computing maxVal and maxIndex for cycle outputs
        computeCycleSplineOutputs(cycle);
        if (useEnvelope) {
            if (cycle.maxVal > 0) {
                ratio = maxSampleValues[i] / cycle.maxVal;
            } else {
                DBG("cycle [" << i << "] has maxVal 0");
            }
        }
//        DBG("cycle number:  " << i);
//        cycle.printData();
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
        for (int j=0; j<M; j++) {
            writeBuffer.setSample(0, j+L, ratio * cycle.outputs[j]);
        }
        i += 1;
    }
}

// for rendering audio with multiple cycles for transitioning and feeding buffer in Callback
// for now, set first cycle to be cycleToGraph (double-clicked), then randomize bcoeffs
void MainContentComponent::setSplineArrays()
{
    int d = 3;
    int k = kVal;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
    // compute cubic delta spline here
    computeCubicSinusoidBcoeffs();
    // changing the array of control coeffs to be 3*n rows and 4 cols
    // so c_i,j is now c[4*(i+m*n)+j], for cols j=0,1,2,3 and 3 sets of rows for m=0,1,2
    for (int i=0; i<12*n; i++) {  // (3*n)*4 = rows*cols
        controlCoeffs.set(i, 0);
    }
    
    // New knot sequence: 0,0,0,0,1/k,2/k,...,(k-1)/k,1,1,1,1
    float incr = 1 / (float) k;
    if (graphView.graphNewSplineCycle) {
        for (int i=0; i<d+1; i++) {  // d+1
            knotVals.set(i, 0);
        }
        for (int i=4; i<N-d; i++) {  // N-d-1-3
            knotVals.set(i, knotVals[i-1]+incr);
        }
        for (int i=N-d; i<N+1; i++) {  // d+1
            knotVals.set(i, 1);
        }                            // total = N-d-4+2d+2 = N+1

    // Old knot sequence: -3/k,-2/k,-1/k,0,1/k,2/k,...,(k-1)/k,1,1+1/k,1+2/k,1+3/k
    } else {
        for (int i=0; i<N+1; i++) {
            knotVals.set(i, (i-d) * incr);
    //        DBG("knot vals [" << i << "] = " << knotVals[i]);
        }
    }

    for (int i=0; i<n; i++) {
        if (graphView.graphNewSplineCycle) {
//            controlCoeffs.set(4*i, graphView.cycleNew.bcoeffs[i]);
            // change here to use cubic spline sinusoid
            controlCoeffs.set(4*i, cubicSinusoidSpline.bcoeffs[i]);
            graphView.cycleNew.bcoeffs.set(i, cubicSinusoidSpline.bcoeffs[i]);
            graphView.repaint();
        } else {
            controlCoeffs.set(4*i, graphView.cycleToGraph.bcoeffs[i]);
        }
    }
//    for (int m=0; m<10; m++) {
//        for (int i=0; i<n; i++) {
//            controlCoeffs.set(4*(i+m*n), graphView.cycles[m].bcoeffs[i]);
//        }
//    }
}

// randomize bcoeffs in controlCoeffs array depending on which part is in use
// so that we can change the cycle with UI as it is playing
// for this to work, we need to first specify continuity conditions at the ends
// of the interval. This could be simply that we preserve the zero value at
// the ends, and also the slopes.  It might work to also randomize in such a way
// that the slope can change but still remain consistent at the ends.
void MainContentComponent::randomizeCycleBcoeffs(int cycleInUse)
{
    // randomize bcoeffs for cycle[i], where i is not in use
    // cycleInUse is 0 if controlCoeffs first part is playing
    int i = 0, j = 1;
    if (cycleInUse == 0) {
        i = 1;
        j = 0;
    }
    int n = dVal + kVal;
    for (int p=0; p<n; p++) {
        float r = juce::Random::getSystemRandom().nextFloat();
        r = 2 * r - 1; // random float in [-1,1]
        r *= 0.05;
        controlCoeffs.set(4*(p+i*n), controlCoeffs[4*(p+j*n)] + r);
    }
    if (cycleInUse == 0) {
        cycleInUse = 1;
    } else {
        cycleInUse = 0;
    }
}

// previous version
//void MainContentComponent::randomizeCycleBcoeffs(int cycleRendered)
//{
//    // randomize bcoeffs for cycle[i], i = cycleRendered + 5 (mod 10)
//    int i = cycleRendered;
//    i += 5;
//    i = i % 10;  // i is now index of cycle to randomize in controlCoeffs
////    CycleSpline cycle = graphView.cycles[i];
////    for (int j=0; j<cycle.bcoeffs.size(); j++) {
////        float r = juce::Random::getSystemRandom().nextFloat();
////        r = 2 * r - 1; // random float in [-1,1]
////        r *= 0.001;
////        cycle.bcoeffs.set(j, cycle.bcoeffs[i] + r);
////    }
////    // copy cycle.bcoeffs to controlCoeffs:
//    int n = dVal + kVal;
//    for (int p=0; p<n; p++) {
//        float r = juce::Random::getSystemRandom().nextFloat();
//        r = 2 * r - 1; // random float in [-1,1]
//        r *= 0.05;
//        controlCoeffs.set(4*(p+i*n), controlCoeffs[4*(p+i*n)] + r);
//    }
//}

void MainContentComponent::randBcoeff(int n)
{
    float r = juce::Random::getSystemRandom().nextFloat();
    r = 2 * r - 1; // random float in [-1,1]
    r *= 0.1;
    float val = graphView.cycleNew.bcoeffs[n] + r;
    if (abs(val) < 0.8) {
        graphView.cycleNew.bcoeffs.set(n, val);
    }
}
    
void MainContentComponent::fillBcoeffs()
{
//    if (cycleRendered == 0) {
//        // fill controlCoeffs with bcoeffs in first column
//    }
}

float MainContentComponent::computeSpline(int control, float t)
{
    // assume t is in [0,1] and output is 0 at the ends
    // this function is called in the audio callback getNextAudioBlock()

    int m = control;
//    m = 0;  // set this to use only cycles[0] in graphView setSplineArrays
            // ignore the other cycles in controlCoeffs
    int d = 3;
    int k = kVal;
    float output = 0;
    int J = 0;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
    if ((t > 0) && (t < 1)) {
//        for (int i=1; i<N; i++) {
//            if (t < knotVals[i]) {
//              J = i-1;
//                if (J > n-1) {
//                    J = n-1;
//                }
//              break;
//            }
//        }
        int i = 0;
        while (t > knotVals[i]) {
            i += 1;
        }
        J = i - 1;
        if (J > n-1) {
            J = n-1;
        }
        // controlCoeffs[i,j] is controlCoeffs[4*i+j], row=i,col=j
        // m=0,1,2 shifts controlCoeffs by m*n
        for (int p=1; p<d+1; p++)
        {
            for (int i=J-d+p; i<J+1; i++)
            {
              float denom = (knotVals[i+d-(p-1)]-knotVals[i]);
              float fac1 = (t-knotVals[i]) / denom;
              float fac2 = (knotVals[i+d-(p-1)]-t) / denom;
              controlCoeffs.set(4*(i+m*n)+p, fac1 * controlCoeffs[4*(i+m*n)+(p-1)]
                  + fac2 * controlCoeffs[4*(i+m*n-1)+(p-1)]);
            }
        }
        output = controlCoeffs[4*(J+m*n)+d];
    }
    return output;
}

// returns true if i is the value of some key cycle
bool MainContentComponent::isKey(int i)
{
    for (int j=0; j<keys.size(); j++) {
        if (keys[j] == i) {
            return true;
        }
    }
    return false;
}

// returns the index of key cycle with value i
int MainContentComponent::keyIndex(int i)
{
    for (int j=0; j<keys.size(); j++) {
        if (keys[j] == i) {
            return j;
        }
    }
    return -1;
}


// for envelope, find max values from audio data per cycle
void MainContentComponent::findMaxValuesPerCycle(Array<float>&  maxSampleIndices, Array<float>&  maxSampleValues, Array<float>& cycleZeros, AudioBuffer<float>& samples)
{
    int n = cycleZeros.size();
    for (int i=0; i<n-1; i++) {
        float a = cycleZeros[i];
        float b = cycleZeros[i+1];
        float maxValue = 0;
        int maxIndex = 0;
        int startIndex = 0;
        float temp = 0;
        if (a > 0) {
            startIndex = (int)a + 1;
        }
        int endIndex = (int)b;
        int numSamples = endIndex - startIndex + 1;
        for (int j=0; j<numSamples; j++) {
            temp = abs(samples.getSample(0,startIndex + j));
            if (temp > maxValue) {
                maxValue = temp;
                maxIndex = j;
            }
        }
        maxSampleIndices.set(i, maxIndex);
        maxSampleValues.set(i, maxValue);
    }
}

int MainContentComponent::getCycleNumber(float t)
{
    int i = 0;
    while (cycleZeros[i] < t) {
        i++;
    }
    return i-1;
}


// Write model to buffer with uniform (or other) cycle interpolation
void MainContentComponent::writeCycleInterpModelToBuffer()
{
    if (noCycleInterp) {
        mVal = 1;
        modelWithoutCycleInterp();
        return;
    }
    setKeyCycleIndices();
    if (expCycleInterp) {
        writeKeyCyclesToBuffer();
        computeMetaSplines();
        writeNonKeyCyclesToBuffer();
    }
    if (mVal > 1) {
        // Do cycle interpolation with key cycle indices 0,m,2m,...
        // This is now if function writeKeyCyclesToBuffer()
        writeBuffer.clear();
        writeBuffer.setSize(1, lastSample+1);  // channels, samples
//        Array<float> keyBcoeffs;
        int numKeyCycles = numCycles / mVal;  // dim of meta splines
        float a = 0, b = 1;
        CycleSpline cycle = CycleSpline(kVal, a, b);
        auto* leftBuffer  = writeBuffer.getWritePointer (0, 0);
        for (int t=0; t<lastSample+1; t++) {
//            writeBuffer.setSample(0, t, 0);
            leftBuffer[t] = 0;
        }
        int i = 0;
        int k = kVal;
        int d = dVal;
        int n = k + d;  // dim of cycle splines
        // compute key cycles and write outputs to buffer
        while (i < numKeyCycles)
        {
            a = cycleZeros[mVal*i];
            b = cycleZeros[mVal*i+1];
            cycle = CycleSpline(kVal, a, b);
            // compute each key cycle bcoeffs and store in array of floats
            computeCycleBcoeffs(cycle, floatBuffer);
            for (int p=0; p<n; p++) {
                keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
            }
            computeCycleSplineOutputs(cycle);
            int M = cycle.outputs.size();
            int L = (int) a + 1;  // first sample index to write
            for (int j=0; j<M; j++) {
                writeBuffer.setSample(0, j+L, cycle.outputs[j]);
            }
            i += 1;
        }
        // Now key cycles are written to writeBuffer, and key cycle bcoeffs are in array keyBcoeffs
        // next compute bcoeffs for intermediate cycles (from meta-splines) and write outputs
        // meta-splines use the keyBcoeffs as targets
        
        
        MetaSpline spline = MetaSpline(numKeyCycles, mVal);
        Array<MetaSpline> metaSplineArray;
        for (int i=0; i<n; i++) {  // loop on each bcoeff of cycles with same fixed kVal
            spline = MetaSpline(numKeyCycles, mVal,numCycles);   // n = numKeyCycles (targets)
            for (int j=0; j<numKeyCycles; j++) {
                spline.targets.set(j, keyBcoeffs[j*n+i]);
            }
            computeMetaSplineBcoeffs(spline);
//            spline.printData();
            computeMetaSplineOutputs(spline);
            spline.printData();
            metaSplineArray.add(spline);
//            DBG("size of metaSplineArray:  " << metaSplineArray.size());
//            spline.printData();
//            DBG("mVal :  " << mVal << "  kVal :  " << kVal);
//            DBG("numKeyCycles: " << numKeyCycles);
        }  // now we have n meta-splines, one for each of n bcoeffs
        
        // compute bcoeffs and write outputs for intermediate cycles
        i = 0;  // i is cycle number, now compute all cycles up to the last key cycle
        n = kVal + dVal;
        while (i < mVal * numKeyCycles) {
            if (i % mVal == 0) { // key cycles already written
            } else {
                // compute bcoeffs for cycle i from meta-spline and write outputs
                // cycle i starting sample is (int) cycle.a + 1
                a = cycleZeros[i];
                b = cycleZeros[i+1];
                cycle = CycleSpline(kVal, a, b);  // this is cycle i
                for (int j=0; j<n; j++) {
                    cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
                }
                computeCycleSplineOutputs(cycle);
                DBG("cycle number:  " << i);
                cycle.printData();
//                metaSplineArray[5].printData();
                int M = cycle.outputs.size();
                int L = (int) a + 1;  // first sample index to write
                for (int j=0; j<M; j++) {
                    writeBuffer.setSample(0, j+L, cycle.outputs[j]);
                }
            }
            i += 1;
        }
    }
}

// previous version used iterateCA from graphView
//        a = cycleZeros[keys[i]];
//        b = cycleZeros[keys[i]+1];
//        int r = 30;
//        cycle = CycleSpline(kVal, a, b);
//        // iterate CA to compute next key cycle bcoeffs
//        if (modelWithCA) {
//            nextCAcycle(r);  // uses ECA(r)
//        } else {
//            graphView.iterateCA();
//        }
//        for (int p=0; p<n; p++) {
//            keyBcoeffs.set(i*n+p, graphView.cycleNew.bcoeffs[p]);
//            cycle.bcoeffs.set(p, graphView.cycleNew.bcoeffs[p]);
//        }
//        computeCycleSplineOutputs(cycle);
//        keyCycleArray.add(cycle);
//        int M = cycle.outputs.size();
//        int L = ((int)a) + 1;  // first sample index to write
//        for (int j=0; j<M; j++) {
//            float value = cycle.outputs[j];
//            writeBuffer.setSample(0, j+L, value);
//        }


