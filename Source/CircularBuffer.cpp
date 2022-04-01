/*
  ==============================================================================

    CircularBuffer.cpp
    Created: 21 Feb 2022 5:36:45pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "CircularBuffer.h"


// implementation example: (.cpp)
// ---------------------------

CircularBuffer::CircularBuffer(unsigned inMaxNumElements)
  : buffer_(inMaxNumElements, 0)  // initialize internal buffer with zeros
  , bufferSize_(inMaxNumElements)  // initialize internal buffer size
{
    // verify inMaxNumElements is a valid value
    // > 0, less than max value of your index types
}

void CircularBuffer::pushSample(const float inSample)
{
    buffer_[currWriteIndex_++] = inSample;

    // wrap the write index if needed
    if(currWriteIndex_ == bufferSize_)
    {
        currWriteIndex_ = 0;
    }
}

void CircularBuffer::pushBuffer(const float* inSamples, unsigned inNumSamplesToPush)
{
    jassert(inNumSamplesToPush <= bufferSize_);
    
    // how much room is left between the write index and the end of buffer_
    const unsigned inNumSamplesBeforeEndOfBuffer = (bufferSize_ - currWriteIndex_);

    // can we do this push in a single memcpy call?
    if(inNumSamplesToPush < 0)  // change 0 to some max value
    {
        // todo: copy whole buffer with single std::memcpy call
        // todo: update currWriteIndex
    }
    else // perform 2 copies using recursion
    {
        const float* inSamples2 = inSamples +  inNumSamplesBeforeEndOfBuffer; // pointer arithmetic
        unsigned numSamples2 = inNumSamplesToPush - inNumSamplesBeforeEndOfBuffer;

        // 1.) Fill to the end of buffer_
        pushBuffer(inSamples, inNumSamplesBeforeEndOfBuffer);

        jassert(currWriteIndex_ == 0); // should be writing to the top of buffer_ next

        // 2.) Copy the remaining input samples to the top of buffer_
        pushBuffer(inSamples2, inNumSamplesBeforeEndOfBuffer);
    }
}


float CircularBuffer::popSample()
{
    float outSample = buffer_[currReadIndex_++];

    // wrap the read index if needed
    if(currReadIndex_ == bufferSize_)
    {
        currReadIndex_ = 0;
    }

    return outSample;
}


void CircularBuffer::popBuffer(float* inBufferToFill, unsigned inNumSampesToPop)
{
    // similar idea to pushBuffer, but src is buffer_ and dest is inBufferToFill
    // might need 2 memcpy calls if the copy crosses the "loop" boundary of buffer_
}


// extra stuff on initialization:

//    std::vector<int> vector1(length, 0);
//    std::vector<double> vector2(length, 0.0);

//    std::vector<SomeStruct> someStructVect(length);
//    memset(someStructVect.data(), 0, sizeof(SomeStruct)*length);

//    std::vector<int> vecOfInts;
//    vecOfInts.resize(10);
//
//    std::fill(vecOfInts.begin(), vecOfInts.end(), 0);
//
//    for (auto const& intVal : vecOfInts)
//    {
//        std::cout << intVal << " ";
//    }
