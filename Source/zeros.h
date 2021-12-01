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

using namespace std;

float zerocross(int I, short *samples);
float interp(float input, short *samples);
int FindAllZeros(int Periods, float Freq, char *data, float *allzeros);
int FindZerosClosestToPeriods(int Periods, float Freq, float *zeros, float *allzeros, float *SPP, float LastZero);
