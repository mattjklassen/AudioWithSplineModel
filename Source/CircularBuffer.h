/*
  ==============================================================================

    CircularBuffer.h
    Created: 21 Feb 2022 5:36:45pm
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once
#include <stddef.h>
#include <math.h>
#include <vector>
#include "../JuceLibraryCode/JuceHeader.h"

class CircularBuffer
{
public:
    // ctor
    // pre-allocates the internal buffer_ to be inMaxNumElements in size
    CircularBuffer(unsigned inMaxNumElements);

    // interface for pushing data to the buffer (data in)
    // updates currWriteIndex_
    // asserts if currWriteIndex_ overtakes currReadIndex_ (inMaxNumElements was too small!)
    //     this would mean we are overwriting data that didn't get a chance to be read
    void pushSample(const float inSample); // hacky, slow (unless this is for a delay line)
    void pushBuffer(const float* inSamples, unsigned inNumSamples); // optimized
    
    // interface for popping data off the buffer (data out)
    // updates currReadIndex_
    // asserts if currReadIndex_ overtakes currWriteIndex_ (inMaxNumElements was too small)
    //     this would mean we are reading stale/garbage data that hasn't been (re)written yet
    float popSample(); // hacky, slow (unless this is for a delay line)
    //  (assumes inBufferToFill is pre-allocated)
    void popBuffer(float* inBufferToFill, unsigned inNumSampesToPop); // optimized

private:
    // underlying buffer
    std::vector<float> buffer_;

    // cached size of the buffer
    unsigned bufferSize_;

    // the first index we will read from next
    unsigned currReadIndex_ = 0;

    // the first index we will write to next
    unsigned currWriteIndex_ = 0;

}; // class CircularBuffer


