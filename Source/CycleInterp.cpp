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
        uint32 time1 = Time::getMillisecondCounter();
        setKeyCycleIndices();
        // when computing from deltaModel the next function changes to use already computed key cycles:
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
// F[x,y,z] = x*y + x*z + y*z - x*y*z.
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
    // The first B-1 cycles are modeled without delta (ie. use zero-crossings), then the remaining cycles use delta.
    // The delta model tries to avoid shape discontinuity by maintaining similarity from one cycle to the next.
    // Starting with cycle B, we recompute end point point guess with the period length guess,
    // followed by adjustment of right endpoint based on minimizing some error functions comparing cycle j-1 to j,
    // then use average error as a reference, and when this is exceeded we introduce delta.
    // In addition to comparison with cycle j-1 values and derivatives, we also use y1*y1.
    // Each cycle will have some new parameters, in particular y0 and y1, which are the signal
    // values at the end points a and b.  Also, use new knot sequence, 0,0,0,0,1/k,...,(k-1)/k,1,1,1,1
    
    computeCubicPolyBcoeffs();  // this is q(t) = 3*t^2-2*t^3
//    numCycles -= 3;
    
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
    int B = 8;
    for (int i=0; i<B; i++)  // compute cycles 0 to B-1 based on zeros
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
        allCycleArray.add(cycle);
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
    for (int i=B; i<numCycles; i++)  // compute cycles > B-1
    {
        a = cycleBreakPoints[i];
        b = a + samplesPerCycleGuess;   // for piano1.wav this is 155.8 (for freq 283 Hz)
        cycleBreakPoints.set(i+1, b);   // will also use for a on next cycle, unless error is too high
        cycleZeros.set(i+1, b);         // keep cycleZeros in synch with cycleBreakPoints for cycle interpolation
        cycle2 = CycleSpline(kVal, a, b);  // bcoeffs are set to 0's
        cycle1 = allCycleArray[i-1];
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
//        addDeltaToCycleBcoeffs(cycle2);   // this was redundant
        computeCycleSplineOutputs(cycle2);
//        cycle2.printData();
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
        allCycleArray.add(cycle2);
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
    DBG("finished computing model with delta");
}

// compute spline model without interpolation and write to buffer
void MainContentComponent::modelWithoutCycleInterp()
{
    DBG("computing model without cycle interp");
    writeBuffer.clear();
    writeBuffer.setSize(1, lastSample+1);  // channels, samples
    float a = 0, b = 1, error1 = 0, error2 = 0, error3 = 0;
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
        keyCycleArray.add(cycle);
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
        DBG("keys[" << k << "]: " << keys[k]);
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
    DBG("computing (only) key cycles");
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
// to move away from zeros with y values y0 and y1.  All cycles are now stored in the array allCycleArray.

// same as above compute bcoeffs of key cycles,
// but also compute outputs and write to buffer
void MainContentComponent::writeKeyCyclesToBuffer()
{
    DBG("computing and writing key cycles");
    keyCycleArray.clear();
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
        if (modelWithDelta) {  // in initial delta model all cycles are key and have y0, y1
            // key cycles are now chosen with indices from full set of cycles
            cycle = allCycleArray[keys[i]];
            a = cycle.a;
            b = cycle.b;
//            if (keys[i] == 360) {
//                DBG("cycle 360:");
//                cycle.printData();
//            }
        } else {
            a = cycleZeros[keys[i]];
            b = cycleZeros[keys[i]+1];
            cycle = CycleSpline(kVal, a, b);
            // compute each key cycle bcoeffs and store in array of floats
            setNewTargets(cycle);
            computeNewBcoeffs(cycle);
            computeCycleSplineOutputs(cycle);
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

void MainContentComponent::computeMetaSplines()
{
    // this is doing piecewise linear metasplines for now, so bcoeffs for metasplines are skipped
    DBG("computing metasplines");
    // compute metasplines based on distribution of keycycles
    // targets of metasplines are keycycle bcoeffs
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
                spline.inputs.set(j, (float)keys[j] / (float)(numCycles-1));
            }
            spline.targets.set(j, keyBcoeffs[j*n+i]);
        }
//        computeMetaSplineBcoeffs(spline);
        computeLinearMetaSplineOutputs(spline);
//        spline.printData();
        metaSplineArray.add(spline);
    }  // now we have n=k+d meta-splines, one for each of n bcoeffs
}


void MainContentComponent::writeNormalizedCyclesToBuffer()
{
    // write all cycles using constant cycle length averageSamplesPerCycle
    DBG("computing and writing normalized cycles");
    DBG("freqGuess: " << freqGuess);
    DBG("samplesPerCycleGuess: " << samplesPerCycleGuess);
    DBG("averageSamplesPerCycle: " << averageSamplesPerCycle);
    int i = 0;
    float a = 0, b=1, ratio = 1;
    CycleSpline cycle;
    int n = kVal + dVal;
    float total = (float) numCycles;
    while (i < numCycles-20) {
        DBG("i: " << i);
        // the next block forces each cycle to have averageSamplesPerCycle samples
        // as computed by computeCycleSplineOutputs(cycle) using only a and b
        if (i == 0) {
            a = 0;
            b = (float)((int)averageSamplesPerCycle) + 0.5;
        } else {
            a = b;
            b = a + (int)averageSamplesPerCycle;
        }
        if (isKey(i)) {   // use key cycle i
            int j = keyIndex(i);
            cycle = keyCycleArray[j];
            cycle.a = a; cycle.b = b;   // reset to normalized a and b
//            DBG("writing key cycle: " << i);
        } else {
            // set bcoeffs for (non-key) cycle i using meta-spline
//            DBG("writing cycle: " << i);
            cycle = CycleSpline(kVal, a, b);  // this is cycle i
            for (int j=0; j<n; j++) {
                cycle.bcoeffs.set(j, metaSplineArray[j].outputs[i]);
            }
//            cycle.printData();
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
            i += 1;
        } else {
            // set bcoeffs for cycle i using meta-spline and write outputs
            // cycle i starting sample is (int) cycle.a + 1
            a = cycleZeros[i];
            b = cycleZeros[i+1];
            y0 = allCycleArray[i].y0;
            y1 = allCycleArray[i].y1;
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
