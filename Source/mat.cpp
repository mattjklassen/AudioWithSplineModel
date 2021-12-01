/*
  ==============================================================================

    mat.cpp
    Created: 13 Nov 2021 1:48:59pm
    Author:  Matt Klassen

  ==============================================================================
*/

#include "mat.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

void fillbmat(int k, int d, float *bmat, float *knots);
void fillmat(int k, int d, float *mat, float *inputs, float *knots);
void gausselim(int n, float *mat, float *matinv);
int  getpivot(int n, int j, float *mat);
void getzero(int n, int i, int j, float *mat, float *matinv);
void initbmat(int k, int d, float *bmat);
void matcopy(int n, float *mat, float *temp);
void printmat(int n, float *mat);
void printmult(int n, float *mat, float *matinv);
void rowmult(int n, int i, float C, float *mat, float *matinv);
void rowswap(int n, int m, int r, float *mat, float *matinv);
void setidentity(int n, float *mat);


void fillbmat(int k, int d, float *bmat, float *knots)
{
// matrix bmat describes B-splines as column vectors or coordinate
// vectors with respect to the truncated power function basis.
// Each coordinate vector comes from a quotient of Vandermonde
// determinants which interprets the divided difference formulation
// of B-splines through Cramer's Rule.
  int sign = 1;
  int i, j, m;

  for (j=0; j<k+d; j++)   // columns j for B-splines B^d_j(t)
  {
      sign = 1;
      for (i=j; i<j+d+2; i++)  // rows i for truncated powers (t-ti)^d_+
      {
	  if (i>k+d-1)
	  {
	      break;
	  }
	  bmat[(k+d)*i+j] = knots[j+d+1]-knots[j];
	  for (m=j; m<i; m++)
	  {
	      bmat[(k+d)*i+j] /= knots[i]-knots[m];
	  }
	  for (m=i+1; m<j+d+2; m++)
	  {
	      bmat[(k+d)*i+j] /= knots[m]-knots[i];
	  }
	  bmat[(k+d)*i+j] *= sign;
	  sign *= -1;
      }
  }
}  // end function fillbmat

void fillmat(int k, int d, float *mat, float *inputs, float *knots)
{
  int i, j;
  for (i=0; i<k+d; i++)  // rows i for inputs[i]
  {
      for (j=0; j<k+d; j++)   // columns j for function (t-tj)^d_+
      {
	  if (inputs[i] < knots[j])
	  {
	      mat[(k+d)*i+j] = 0;
	  } else {
	      mat[(k+d)*i+j] = pow(inputs[i]-knots[j],d);
	  }
	  // if (i == j)
	  // {
	  //  cout << "i,j:  " << i << "," << j << endl;
	  //  cout << "mat[(k+d)*i+j] : " << mat[(k+d)*i+j] << endl;
	  // }
      }
  }
}   // end function fillmat


void gausselim(int n, float *mat, float *matinv)
{

  int i, j, m=0;

  for (j=0; j<n; j++)
  {
      // get pivot for column j, choosing max from columns
      m = getpivot(n, j, mat);
      // cout << "pivot row is m =  " << m << endl;
      // swap row m (with max abs) for row j 
      if (m > j)
      {
         rowswap(n, m, j, mat, matinv);
      }
      float pivot = mat[n*j+j];
      // cout << "pivot value is =  " << pivot << endl;
      if (pivot > 0 || pivot < 0)
      {
         // get 1 on diagonal
         rowmult(n, j, 1/pivot, mat, matinv);
      } else 
      {
	//  cout << "No pivot for column " << j << endl;
	  return;
      }
      for (i=0; i<n; i++)
      {
	  // get zeros in rest of column j
	  if (i < j || i > j)
	  {
	     getzero(n, i, j, mat, matinv);  // skips case i=j
	  }
      }
  }

}  // end function gausselim



// find pivot row for col j with max abs value
// assume nxn matrix mat
int getpivot(int n, int j, float *mat)
{
    int i;
    int m = j;
    float max = abs(mat[n*j+j]);
    for (i=j+1; i<n; i++)
    {
	if (abs(mat[n*i+j]) > max)
	{
	    max = abs(mat[n*i+j]);
	    m = i;
	}
    }
    return m;
}


// get zero in col j for row i, for mat nxn
// assumes pivot is 1
void getzero(int n, int i, int j, float *mat, float *matinv)
{
    int m;
    float C = mat[n*i+j];
    // replace row i with itself plus (-C)*(pivot row)
    for (m=0; m<n; m++)
    {
	mat[n*i+m] += (-C)*mat[n*j+m];
	matinv[n*i+m] += (-C)*matinv[n*j+m];
    }
}

void initbmat(int k, int d, float *bmat)
{
  int i, j;
  for (j=0; j<k+d; j++)   // columns j for B-splines B^d_j(t)
  {
      for (i=0; i<k+d; i++)  // rows i for truncated powers (t-ti)^d_+
      {
	  bmat[(k+d)*i+j] = 0;
      }
  }
}  // end function initbmat


void matcopy(int n, float *mat, float *temp)
{
    int i, j;

    for (i=0; i<n; i++)
    {
       for (j=0; j<n; j++)
       {
           temp[n*i+j] = mat[n*i+j];
       }
    }
}  // end function matcopy

void printmat(int n, float *mat)
{
    int i, j;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
    	    cout << i << "," << j << " :  " << mat[n*i+j] << "  ";
        }
	cout << endl;
    }
}

void printmult(int n, float *mat, float *matinv)
{
   float prod;
   int i, j, m;

   for (i=0; i<n; i++)
   {
      for (j=0; j<n; j++)
      {
         prod = 0;
	 for (m=0; m<n; m++)
	 {
	     prod += mat[n*i+m]*matinv[n*m+j];
	 }
         cout << "i,j : " << i << "," << j << " :  " << prod << endl;
      }
    }
}

void rowmult(int n, int i, float C, float *mat, float *matinv)
{
    int j;
    for (j=0; j<n; j++)
    {
	mat[n*i+j] = C*mat[n*i+j];
	matinv[n*i+j] = C*matinv[n*i+j];
    }
}

// swap row m with row r
void rowswap(int n, int m, int r, float *mat, float *matinv)
{
    float temp;
    int j;
    for (j=0; j<n; j++)
    {
	temp = mat[n*m+j];
	mat[n*m+j] = mat[n*r+j];
	mat[n*r+j] = temp;
	temp = matinv[n*m+j];
	matinv[n*m+j] = matinv[n*r+j];
	matinv[n*r+j] = temp;
    }
}

void setidentity(int n, float *mat)
{
  int i, j;
  for (i=0; i<n; i++) 
  {
      for (j=0; j<n; j++) 
      {
	  if (i == j)
	  {
	      mat[(n)*i+j] = 1;
	  } 
	  else
	  {
	      mat[(n)*i+j] = 0;
	  }
      }
  }
}   // end function setidentity

