/*
  ==============================================================================

    CycleInterp.cpp
    Created: 8 Jan 2022 3:50:44pm
    Author:  Matt Klassen

  ==============================================================================
*/

// Cycle Interpolation ToDo List:

// The cycle interpolation model should have chosen key cycles and meta data which shows how
// to construct the interpolating cycles, including meta-splines used to generate the B-spline
// coefficients of non-key cycles, and lengths of intervals for each cycle (if changing).

// If we create the simplest model first, it should be to use the cycle lengths and samples per cycle
// from a wav file, so the model just follows the original data but creates a spline for each cycle
// and then generates sample data from those splines.  This method is used if m=1 to compute model.

// If m>1 then we are doing cycle interpolation with key cycles j*m, j=0,1,2,...

// Suppose for now that we are going to use constant cycle length, and only the B-spline
// coefficients are changing from cycle to cycle, not the cycle lengths or subinterval sizes.

// To deal with the realtime rendering of spline model cycle by cycle:
// I think I will have three sets of arrays which are initialized with the B-spline coeffs
// for the first three cycles.  When the first cycle is finished rendering, we can start to
// write the coefficients of the 4th cycle over those for the first cycle.  The AudioCallBack
// can then pass an int variable outside to say which cycle has finished rendering.  This
// could be a Value type which has a listener that responds by writing the next B-spline
// coefficients.  So cycles are organized mod 3. When cycle i finishes rendering, overwrite
// its B-spline coefficients with those for cycle i + 3.
// write cycleBcoeffs for i=0,1,2
// render 0, overwrite with 3
// render 1, overwrite with 4
// render 2, overwrite with 5 etc.
// The time it takes to write one cycle's bcoeffs should not be longer than the render/playback
// of one cycle.  This way the new coeffs will be in place before needed.
// If new cycle bcoeffs are being computed from a meta-spline, these should still be computable
// in realtime in the above scheme.

// The selection of key cycles should be done with some preset schemes.
// 1. (Simple regular) Select every m^th cycle as key, and interpolate all others
// - this method will kill oscillating patterns between key cycles, and possibly also subharmonics
// * meta-splines should be used to generate B-spline coefficients of non-key cycles.
// * meta-splines should be applied in groups of cycles of reasonable length.
// ex: 400 cycles in 1 sec, use every 10th cycle as key, so 40 key cycles.
// - suppose sample rate 48 kHz, f_0 guess 400, so 120 samples per cycle
// - suppose each key cycle has n=30 interpolation points (targets), so there are 30 meta-splines
// which describe the coefficients as they evolve through the non-key cycles.
// The number of data points in this example is 1) the B-spline coeffs of key cycles: 40*30=1200
// 2) the meta-spline coeffs for each of 30 meta-splines, with 40 targets (from key cycles): 1200
// so 2400 float values, compared with 48000 original samples, or 1/20 = 5% of the original data.

// 2. (Exponential scheme) Select cycles with indices 0,1,2,4,8,16,32...
// - meta-splines could still be used even though the separation of targets is nonuniform.
// - in fact, one could even still use uniform separation for the splines but then apply them
//   to the nonuniform intervals.
// ex. as above, now with only 8 key cycles: 0,1,2,4,8,16,32,40 (include last one).
// So data is now 1) 8*30=240 2) also 240, so 480 floats, or data reduction of 480/48000=.01
// for a reduction to 1% of original data.

// An important point to consider in the above examples is that one can make modifications to
// the meta-splines without changing the reduction of data, for instance to give some oscillations
// to the pattern of non-key cycles, or to experiment with the degrees of the splines, etc.

// A meta-spline can simply be defined on the interval [0,1] with equally spaced intermediate
// points defining subintervals as before.  If we are working on the spline for bcoeff b_0, and
// we have k+1 key cycles, then we can define the subintervals, knot sequence, and target points
// similarly to how we did for the cycle splines.  The target values are now B-spline coeffs at each of
// of the key cycles.

// Note regarding computation: we can choose a smaller cycle length in order to lower the
// size of linear systems n=k+d but then create a sequence of these smaller cycles in chunks to
// perform the role of key cycles.  For cycle interpolation, however, this would only work if the
// chunk of cycles is still part of a sort of periodic pattern.

#include "MainContentComponent.h"

void MainContentComponent::writeModelToBuffer()
{
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
        modelWithoutCycleInterp();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
        return;
    }
    if (normalizeCycleLength) {
        uint32 time1 = Time::getMillisecondCounter();
    setKeyCycleIndices();
    computeKeyCycles();
//    writeKeyCyclesToBuffer();
    computeMetaSplines();
//    computeAllCycles();
//    writeNonKeyCyclesToBuffer();
    writeNormalizedCyclesToBuffer();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
    } else {
        uint32 time1 = Time::getMillisecondCounter();
    setKeyCycleIndices();
//    computeKeyCycles();
    writeKeyCyclesToBuffer();
    computeMetaSplines();
//    computeAllCycles();
    writeNonKeyCyclesToBuffer();
//    writeNormalizedCyclesToBuffer();
        uint32 time2 = Time::getMillisecondCounter();
        float timeVal = (float(time2)-float(time1)) * 0.001;
        DBG("time to compute model (sec):  " << timeVal);
    }
}

// compute spline model without interpolation and write to buffer
void MainContentComponent::modelWithoutCycleInterp()
{
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    float a = 0, b = 1;
    CycleSpline cycle = CycleSpline(kVal, a, b);
    // compute each cycle and write out the computed samples to writeBuffer
    auto* leftBuffer = writeBuffer.getWritePointer (0, 0);
    int currentSampleIndex = 0;
    for (int i=0; i<numCycles-1; i++)  // for (int i=0; i<10; i++)
    {
        a = cycleZeros[i];
        b = cycleZeros[i+1];
        cycle = CycleSpline(kVal, a, b);
        computeCycleBcoeffs(cycle, floatBuffer);
        computeCycleSplineOutputs(cycle);
        int M = cycle.outputs.size();
        for (int j=0; j<M; ++j) {
            leftBuffer[currentSampleIndex] = cycle.outputs[j];
            currentSampleIndex += 1;
        }
    }
}

void MainContentComponent::setKeyCycleIndices() {
    if (otherCycleInterp) {
            // leave keys to reuse and modify
        } else {
            keys.clear();
        }
    if (expCycleInterp) {
        DBG("setting key cycles for expCycleInterp");
        keys.add(0); keys.add(1);
        int key = 1, j = 1;
        while (key < numCycles) {
            key = keys[j] * 2;
            if (key < numCycles) {
                keys.add(key);
                j++;
            }
        }
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
        keys.add(numCycles-1);
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
        DBG("keys[" << k << "]; " << keys[k]);
    }
    DBG("there are " << keys.size() << " key cycles");
}

// need to decide if we are fixing constant cycle length ???
void MainContentComponent::generateKeyCyclesCA()
{
    keyCycleArray.clear();
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    float a = 0, b = 1;
    CycleSpline cycle = CycleSpline(kVal, a, b);
    int i = 0;
    int k = kVal;
    DBG("kVal = " << kVal);
    int d = dVal;
    int n = k + d;  // dim of cycle splines
    // add first key cycle as cycleNew
    if (graphView.highlightCycle > -1) {
        keyCycleArray.add(graphView.cycleNew);
    } else {
        // use audio graph key cycle 1
        a = cycleZeros[1];
        b = cycleZeros[2];
        cycle = CycleSpline(kVal, a, b);
        computeCycleBcoeffs(cycle, floatBuffer);
        graphView.cycleNew = cycle;
//        graphView.iterateCA();
        cycle = graphView.cycleNew;
        for (int p=0; p<n; p++) {
            keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
        }
        keyCycleArray.add(cycle);
        computeCycleSplineOutputs(cycle);
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
        for (int j=0; j<M; j++) {
            float value = cycle.outputs[j];
            writeBuffer.setSample(0, j+L, value);
        }
    }
    for (i=1; i<keys.size(); i++)
    {
        a = cycleZeros[keys[i]];
        b = cycleZeros[keys[i]+1];
        cycle = CycleSpline(kVal, a, b);
        // iterate CA to compute next key cycle bcoeffs
        graphView.iterateCA();
        for (int p=0; p<n; p++) {
            keyBcoeffs.set(i*n+p, graphView.cycleNew.bcoeffs[p]);
            cycle.bcoeffs.set(p, graphView.cycleNew.bcoeffs[p]);
        }
        keyCycleArray.add(cycle);
        computeCycleSplineOutputs(cycle);
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
        for (int j=0; j<M; j++) {
            float value = cycle.outputs[j];
            writeBuffer.setSample(0, j+L, value);
        }
    }
}

void MainContentComponent::modelWithCAplusCycleInterp()
{
    generateKeyCyclesCA();
    computeMetaSplines();
    writeNormalizedCyclesToBuffer();
    graphView.graphCAmodel = true;
    graphView.repaint();
}

// compute key cycle bcoeffs, write to keyBcoeffs array
// but do not compute outputs and do not write to buffer
void MainContentComponent::computeKeyCycles()
{
    // compute key cycles and add to array keyCycleArray
    keyCycleArray.clear();
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
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
        computeCycleBcoeffs(cycle, floatBuffer);
        for (int p=0; p<n; p++) {
            keyBcoeffs.set(i*n+p, cycle.bcoeffs[p]);
        }
        keyCycleArray.add(cycle);
    }
}

// same as above compute bcoeffs of key cycles,
// but also compute outputs and write to buffer
void MainContentComponent::writeKeyCyclesToBuffer()
{
    DBG("computing and writing key cycles");
//    float Pi = juce::MathConstants<float>::pi;
    keyCycleArray.clear();
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    float a = 0, b = 1;
    CycleSpline cycle = CycleSpline(kVal, a, b);
    int k = kVal;
    int d = dVal;
    int n = k + d;  // dim of cycle splines
    if (addSubharmonic) {
        cycleParabola = CycleSpline(kVal, a, b);
        setParabolicTargets(0.005);
        computeParabolicBcoeffs();
    }
    // compute key cycles and write outputs to buffer
    for (int i=0; i<keys.size(); i++)
    {
        a = cycleZeros[keys[i]];
        b = cycleZeros[keys[i]+1];
        cycle = CycleSpline(kVal, a, b);
        // compute each key cycle bcoeffs and store in array of floats
        computeCycleBcoeffs(cycle, floatBuffer);
        for (int p=0; p<n; p++) {
            float val = cycle.bcoeffs[p];
            keyBcoeffs.set(i*n+p, val);
            if (addSubharmonic) {
                if (i % 2 == 0) {
                    val += cycleParabola.bcoeffs[p];
                } else {
                    val -= cycleParabola.bcoeffs[p];
                }
                cycle.bcoeffs.set(p, val);
            }
        }
        computeCycleSplineOutputs(cycle);
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
        keyCycleArray.add(cycle);
    }
}

void MainContentComponent::computeMetaSplines()
{
    // this is doing piecewise linear metasplines for now, so bcoeffs are skipped
    DBG("computing metasplines");
    // compute metasplines based on distribution of keycycles
    // targets of metasplines are keycycle coeffs for first n-2, 0 for the last two
    // which are for the second derivative conditions at the ends
    metaSplineArray.clear();
    MetaSpline spline = MetaSpline(keys.size(), mVal);
    int n = kVal + dVal;  // dim of cycle splines, = # bcoeffs for each cycle spline
    for (int i=0; i<n; i++) {  // loop on each bcoeff of cycles with same fixed kVal
        spline = MetaSpline(keys.size(), mVal, numCycles); // numCycles = outputs.size()
        for (int j=0; j<keys.size(); j++) {  // j loop on key cycles
            if (regularCycleInterp) {
                // already set uniform inputs
            } else {
                // put metaspline inputs at values keys*(1/numCycles-1)
                spline.inputs.set(j, (float)keys[j] / (float)numCycles-1);
            }
            spline.targets.set(j, keyBcoeffs[j*n+i]);
        }
        // the last two targets are already set to zero, for the end derivatives
//        computeMetaSplineBcoeffs(spline);
        computeLinearMetaSplineOutputs(spline);
//        spline.printData();
        metaSplineArray.add(spline);
    }  // now we have n=k+d meta-splines, one for each of n bcoeffs
}


void MainContentComponent::writeNormalizedCyclesToBuffer()
{
    // write all cycles using constant cycle length samplesPerCycleGuess
    DBG("computing and writing normalized cycles");
    int i = 0;
    float a = 0, b=1, ratio = 1;
    CycleSpline cycle;
    int n = kVal + dVal;
    float total = (float) numCycles;
    while (i < numCycles) {
        // the next block forces each cycle to have samplesPerCycleGuess samples
        // as computed by computeCycleSplineOutputs(cycle) using only a and b
        if (i == 0) {
            a = 0;
            b = (float)((int)samplesPerCycleGuess) + 0.5;
        } else {
            a = b;
            b = a + (int)samplesPerCycleGuess;
        }
        if (isKey(i)) {   // use key cycle i
            int j = keyIndex(i);
            cycle = keyCycleArray[j];
            cycle.a = a; cycle.b = b;   // reset to normalized a and b
            DBG("writing key cycle: " << i);
        } else {
            // set bcoeffs for (non-key) cycle i using meta-spline
            DBG("writing cycle: " << i);
            cycle = CycleSpline(kVal, a, b);  // this is cycle i
            for (int j=0; j<n; j++) {
                cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
            }
            cycle.printData();
        }
        computeCycleSplineOutputs(cycle);

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
        int M = cycle.outputs.size();
        int L = ((int)a) + 1;  // first sample index to write
        for (int j=0; j<M; j++) {
            writeBuffer.setSample(0, j+L, ratio * cycle.outputs[j]);
        }
        i += 1;
    }
}

void MainContentComponent::writeNonKeyCyclesToBuffer()
{
    DBG("computing and writing non-key cycles");
//    float Pi = juce::MathConstants<float>::pi;
    int i = 0;
    float a = 0, b=1, ratio = 1;
    CycleSpline cycle;
    CycleSpline secondLastKeyCycle = keyCycleArray[keyCycleArray.size()-2];
    int n = kVal + dVal;
    if (addSubharmonic) {
        cycleParabola = CycleSpline(kVal, a, b);
        setParabolicTargets(0.005);
        computeParabolicBcoeffs();
    }
    while (i < numCycles) {
        if (isKey(i)) {
            // key cycles already written
            i += 1;
        } else {
            // set bcoeffs for cycle i using meta-spline and write outputs
            // cycle i starting sample is (int) cycle.a + 1
            a = cycleZeros[i];
            b = cycleZeros[i+1];
    // here we set the cycle length to be the same as second last key cycle
//            int Tail = keys[keys.size()-2];  // secondLastKeyCycleIndex
//            float bTail = cycleZeros[Tail+1];
//            if (i > Tail) {
//                a = bTail + (i-Tail-1) * 200;
//                b = a + 200;
//            }
    // end setting normalized cycle length for tail
            cycle = CycleSpline(kVal, a, b);  // this is cycle i
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
//    int k = graphView.cycleToGraph.k;
    int k = kVal;
    int n = k + d;  // dimension of splines, also number of inputs
    int N = n + d;  // last index of knot sequence t_0,...,t_N
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
//        for (int i=0; i<N+1; i++) {
//            DBG("knotVals [" << i << "] = " << knotVals[i]);
//        }
    // Old knot sequence: -3/k,-2/k,-1/k,0,1/k,2/k,...,(k-1)/k,1,1+1/k,1+2/k,1+3/k
    } else {
        for (int i=0; i<N+1; i++) {
            knotVals.set(i, (i-d) * incr);
    //        DBG("knot vals [" << i << "] = " << knotVals[i]);
        }
    }

    for (int i=0; i<n; i++) {
        if (graphView.graphNewSplineCycle) {
            controlCoeffs.set(4*i, graphView.cycleNew.bcoeffs[i]);
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
