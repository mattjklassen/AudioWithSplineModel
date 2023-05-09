/*
  ==============================================================================

    zeros.h
    Created: 9 Nov 2021 3:57:35pm
    Author:  Matt Klassen

  ==============================================================================
*/

#pragma once

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <vector>
#include "../JuceLibraryCode/JuceHeader.h"

float zerocross(int I, short *samples);
float zerocrossFloat(int I, AudioBuffer<float>& floatBuffer);
float interp(float input, short *samples);
float interpFloat(float input, AudioBuffer<float>& floatBuffer);
//void FindAllZerosFloat(int sampleRate, AudioBuffer<float>& floatBuffer,  Array<float>& allZeros);
//void FindZerosClosestToPeriods(int sampleRate, int Periods, float Freq, Array<float>& cycleZeros, Array<float>& allZeros, Array<int>& samplesPerCycle);
int FindAllZerosFloat(int sampleRate, int Periods, float Freq, AudioBuffer<float>& floatBuffer,  Array<float>& allZeros);
int FindZerosClosestToPeriods(int sampleRate, int Periods, float Freq, Array<float>& cycleZeros, Array<float>& allZeros, Array<int>& samplesPerCycle, float LastZero);
