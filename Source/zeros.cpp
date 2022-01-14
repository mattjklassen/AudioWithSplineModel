/*
  ==============================================================================

    zeros.cpp
    Created: 9 Nov 2021 4:02:40pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "zeros.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>


//float zerocross(int I, short *samples);
//float interp(float input, short *samples);
//int FindAllZeros(int Periods, float Freq, char *data, float *allzeros);
//int FindZerosClosestToPeriods(int Periods, float Freq, float *zeros, float *allzeros, float *SPP, float LastZero);
//

// computes the zero crossing between two consecutive samples of opposite sign
float zerocross(int I, short *samples)
{
    // find zero between samples I and I+1:
    // let x1 = 0, y1 = samples[I] < 0 or > 0,
    // and x2 = 1, y2 = samples[I+1] > 0 or < 0 resp.
    // m = (y2-y1)/(x2-x1) = y2-y1, y-y1 = m*(x-x1), 0-y1 = m*(x-x1)
    // zero crossing: x = x1 - y1/m = -y1/(y2-y1) = y1/(y1-y2)
    // = samples[I]/(samples[I]-samples[I+1])
    if ((float)samples[I]*(float)samples[I+1] > 0)
    {
    // cout << "zerocross called with samples of same sign" << endl;
    return 0.0f;
    }
    return (float)(samples[I]) / ((float)samples[I]-(float)samples[I+1]);
}

float zerocrossFloat(int I, AudioBuffer<float>& floatBuffer)
{
    // find zero between samples I and I+1:
    // let x1 = 0, y1 = samples[I] < 0 or > 0,
    // and x2 = 1, y2 = samples[I+1] > 0 or < 0 resp.
    // m = (y2-y1)/(x2-x1) = y2-y1, y-y1 = m*(x-x1), 0-y1 = m*(x-x1)
    // zero crossing: x = x1 - y1/m = -y1/(y2-y1) = y1/(y1-y2)
    // = samples[I]/(samples[I]-samples[I+1])
//    if (samples[I]*samples[I+1] > 0)
    if (floatBuffer.getSample(0, I) * floatBuffer.getSample(0, I+1) > 0)
    {
    // cout << "zerocross called with samples of same sign" << endl;
    return 0.0f;
    }
    return (floatBuffer.getSample(0, I)) / (floatBuffer.getSample(0, I)-floatBuffer.getSample(0, I+1));
}

// interpolate linearly between samples at input
// output is float in interval [-1,1]
float interp(float input, short *samples)
{
    // find interpolated output value bewteen samples I and I+1
    // where I=(int)floor(input) is the greatest integer less or equal to input
    int I = (int)floor(input);
    float x = input - (float)I;
    float output = (1.0f - x)*((float)samples[I]) + x*((float)samples[I+1]);
    return output / 32768.0f ;
}

float interpFloat(float input, AudioBuffer<float>& floatBuffer)
{
    // find interpolated output value bewteen samples I and I+1
    // where I=(int)floor(input) is the greatest integer less or equal to input
    int I = (int) floor(input);
    float x = input - (float) I;
    float output = (1.0f - x) * floatBuffer.getSample(0, I) + x * floatBuffer.getSample(0, I+1);
    return output;
}

// Find all zero crossings in wave file data

int FindAllZerosFloat(int sampleRate, int Periods, float Freq, AudioBuffer<float>& floatBuffer,  Array<float>& allZeros)
{
    float Spp = (float)sampleRate / (float)Freq;  // samples per period guess
    int S = (int)Spp;  // generic guess at (int) samples per period
//    DBG("samples per period guess: " << S);
    int J1 = S*Periods; // max number of samples to scan
    int lastSample = floatBuffer.getNumSamples();
  
    // float Guess;  // for next zero close to period
    int J = 0;  // for big loop on samples
    // float J2 = 0.0f;  // total count of samples
  
//    short *samples = reinterpret_cast<short*>(data);
    // unsigned n = size/2; // number of data samples

    int NumAllZeros = 0;
    // float LastZero = 0.0f;
  
    float zero = 0.0f;
    int i = 0;
    int Toggle = 0;
  
//    std::cout << "All Zero Crossings:" << std::endl;

    float zero0 = 0.0f;  // previous zero
    float zero1 = 0.0f;  // current zero
    int pos0=0;  // previous sample was positive if pos0=1
    int pos1=0;  // current sample was positive if pos1=1
    int neg0=0;  // previous sample was negative if neg0=1
    int neg1=0;  // current sample was negative if neg1=1
    int foundzero=0;
    
    if (floatBuffer.getSample(0, 0) == 0)
    {
        allZeros.add(0);  // set first zero to 0
    }
    if (floatBuffer.getSample(0, 0)>0)
    {
        pos0=1;
    }
    if (floatBuffer.getSample(0, 0)<0)
    {
        neg0=1;
    }
    
    // set max for J approximately by number of cycles or periods
    // to be scanned times the samples per period, or just
    // the number of samples
    
    // i is index for allzeros[]
    
//    for (J=1; J<J1; J++)  // loop on samples
    for (J=1; J<lastSample; J++)  // loop on samples
    {
       foundzero = 0;
       if (floatBuffer.getSample(0, J) == 0)
       {
           foundzero = 1;
           zero1 = (float)(J);
       }
    
       if (floatBuffer.getSample(0, J) > 0) pos1 = 1;  // current sample pos
       if (floatBuffer.getSample(0, J) < 0) neg1 = 1;  // current sample neg
       
       if ((pos1 && neg0) || (neg1 && pos0)) // found zero between J-1 and J
       {
           foundzero = 1;
           zero = zerocrossFloat(J-1, floatBuffer);  // zero crossing between samples
           // added Dec 1, 2021:
           // keep zeros away from endpoints of interval between samples:
           if (zero < 0.01) {
//               std::cout << "zero:  " << zero << std::endl;
               zero = 0.01;
           }
           if (zero > 0.99) {
//               std::cout << "zero:  " << zero << std::endl;
               zero = 0.99;
           }
           // so zeros are not exactly at sample values
           zero1 = (float)(J-1) + zero;  // set current zero sample value
       }
       if (foundzero)
       {
           allZeros.add(zero1);  // set zero number i
       // Toggle = 1;
           if (Toggle == 0)
           {
           //    cout << setprecision(8) << setw(10) << zero1 << "  "
           //    << setprecision(4) << zero1 - zero0 << endl;
           }
           zero0 = zero1;  // prior zero is set to current zero
           i++;            // advance index of allzeros[]
       }
    
       neg0 = neg1;  // set prior sample value neg to current
       pos0 = pos1;  // set prior sample value pos to current
       neg1 = 0;     // reset for current sample value
       pos1 = 0;     // reset for current sample value
    
    }  // end loop on samples to find all zeros
    
    NumAllZeros = i-1;

    return NumAllZeros;
}   // end function FindAllZeros()


// Find zero crossings closest to periods, return the number of zeros
// which will also equal the number of values filled in samplesPerCycle

int FindZerosClosestToPeriods(int sampleRate, int Periods, float Freq, Array<float>& cycleZeros, Array<float>& allZeros, Array<int>& samplesPerCycle, float LastZero)
{
    //  cout << endl << "Zeros Closest to Periods:" << endl << endl;
    //  cout << "Index    Sample   Difference" << endl << endl;
    
    // Now find zeros[] which are the closest to periods
    // zeros[0]=0, Guess for zeros[1] is zeros[0]+SPP
    // Step through allzeros[], find one closest to Guess
    // To do this just check the conditional (Guess > allzeros[])
    // when this fails set z1 = allzeros[i] and z0 = allzeros[i-1]
    // then pick the one that is closest to Guess for zeros[1].
    // Next set Guess to zeros[1]+SPP.  Step through allzeros[]
    // starting with allzeros[i].  It's possible that Guess < allzeros[i]
    // already, which case again z1 = allzeros[i] and z0 = allzeros[i-1]
    // and we will choose closest to Guess which will be z1 and the
    // previous would be z0.  In most cases we will have more zeros.
    
    int i = 0;   // index of allzeros[]
    // zeros[0] = allzeros[0];
    cycleZeros.add(0.0f);
//    SPP[0] = 0.0f;
    float Spp = (float)sampleRate / (float)Freq;  // samples per period guess
    // S = (int)Spp;  // generic guess at (int) samples per period
    // J1 = S*Periods; // max number of samples to scan
    float Guess = cycleZeros[0] + Spp; // for next zero close to period
    int J = 0;  // for loop on Periods
    float zero0, zero1;
    
    for (J=1; J<Periods+10; J++)
    {
        while (Guess > allZeros[i])
        {
            i++;
        }
        zero1 = allZeros[i];
        zero0 = allZeros[i-1];
        if ( abs(Guess - zero1) < abs(Guess - zero0) )
        {
            cycleZeros.add(zero1);
        } else {
            cycleZeros.add(zero0);
        }
        samplesPerCycle.add(int(cycleZeros[J]) - int(cycleZeros[J-1]));
        // cout << setw(4) << J << ":  "
          //   << fixed
          //   << setw(8) << setprecision(2) << zeros[J] << "  "
          //   << fixed
          //   << setw(6) << setprecision(2) << SPP[J] << endl;
        Guess = cycleZeros[J] + Spp;
        if (Guess > LastZero) {
            DBG("Guess " << Guess << " is > LastZero " << LastZero);
//            cycleZeros.add(LastZero);
            DBG("J: " << J << "  cycZeros.size(): " << cycleZeros.size());
            DBG("last cycleZero: " << cycleZeros[cycleZeros.size()-1]);
            break;
        }
    } // end for J
    
        return J;   // This is number of zeros found by checking close to periods
                    // and not exceeding the last found in allzeros[]
    
}  // end function FindZerosClosestToPeriods

