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

using namespace std;

float zerocross(int I, short *samples);
float interp(float input, short *samples);
int FindAllZeros(int Periods, float Freq, char *data, float *allzeros);
int FindZerosClosestToPeriods(int Periods, float Freq, float *zeros, float *allzeros, float *SPP, float LastZero);


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

// Find all zero crossings in wave file data

int FindAllZeros(int Periods, float Freq, char *data, float *allzeros)
{
    float Spp = 44100.0f / (float)Freq;  // samples per period guess
    int S = (int)Spp;  // generic guess at (int) samples per period
    int J1 = S*Periods; // max number of samples to scan
  
    // float Guess;  // for next zero close to period
    int J = 0;  // for big loop on samples
    // float J2 = 0.0f;  // total count of samples
  
    short *samples = reinterpret_cast<short*>(data);
    // unsigned n = size/2; // number of data samples

    int NumAllZeros = 0;
    // float LastZero = 0.0f;
  
    float zero = 0.0f;
    int i = 0;
    int Toggle = 0;
  
    // cout << "All Zero Crossings:" << endl << endl;

    float zero0 = 0.0f;  // previous zero
    float zero1 = 0.0f;  // current zero
    int pos0=0;  // previous sample was positive if pos0=1
    int pos1=0;  // current sample was positive if pos1=1
    int neg0=0;  // previous sample was negative if neg0=1
    int neg1=0;  // current sample was negative if neg1=1
    int foundzero=0;
    
    if (samples[0] == 0)
    {
        allzeros[0] = 0;  // set first zero to 0
    }
    if (samples[0]>0)
    {
        pos0=1;
    }
    if (samples[0]<0)
    {
        neg0=1;
    }
    
    // set max for J approximately by number of cycles or periods
    // to be scanned times the samples per period, or just
    // the number of samples
    
    // i is index for allzeros[]
    
    for (J=1; J<J1; J++)  // loop on samples
    {
       foundzero = 0;
       if (samples[J] == 0)
       {
           foundzero = 1;
           zero1 = (float)(J);
       }
    
       if (samples[J] > 0) pos1 = 1;  // current sample pos
       if (samples[J] < 0) neg1 = 1;  // current sample neg
       
       if ((pos1 && neg0) || (neg1 && pos0)) // zero is between J-1 and J
       {
           foundzero = 1;
           zero = zerocross(J-1, samples);  // zero crossing between samples
           zero1 = (float)(J-1) + zero;  // set current zero sample value
       }
       if (foundzero)
       {
           allzeros[i] = zero1;  // set zero number i
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
// which will also equal the number of values filled in SPP[] = samples per period

int FindZerosClosestToPeriods(int Periods, float Freq, float *zeros, float *allzeros, float *SPP, float LastZero)
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
    zeros[0] = 0.0f;
    SPP[0] = 0.0f;
    float Spp = 44100.0f / (float)Freq;  // samples per period guess
    // S = (int)Spp;  // generic guess at (int) samples per period
    // J1 = S*Periods; // max number of samples to scan
    float Guess = zeros[0] + Spp; // for next zero close to period
    int J = 0;  // for loop on Periods
    float zero0, zero1;
    
    for (J=1; J<Periods; J++)
    {
        while (Guess > allzeros[i])
        {
            i++;
        }
        zero1 = allzeros[i];
        zero0 = allzeros[i-1];
        zeros[J] = zero0;
        if ( abs(Guess - zero1) < abs(Guess - zero0) )
        {
            zeros[J] = zero1;
        }
        SPP[J] = zeros[J] - zeros[J-1];
        // cout << setw(4) << J << ":  "
          //   << fixed
          //   << setw(8) << setprecision(2) << zeros[J] << "  "
          //   << fixed
          //   << setw(6) << setprecision(2) << SPP[J] << endl;
        Guess = zeros[J] + Spp;
        if (Guess > LastZero) break;
    } // end for J
    
        return J;   // This is number of zeros found by checking close to periods
                  // and not exceeding the last found in allzeros[]
    
}  // end function FindZerosClosestToPeriods
