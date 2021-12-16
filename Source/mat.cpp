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
#include "../JuceLibraryCode/JuceHeader.h"

using namespace std;

//void fillbmat(int k, int d, double *bmat, double *knots);
//void fillmat(int k, int d, double *mat, double *inputs, double *knots);
//void gausselim(int n, double *mat, double *matinv);
//int  getpivot(int n, int j, double *mat);
//void getzero(int n, int i, int j, double *mat, double *matinv);
//void initbmat(int k, int d, double *bmat);
//void matcopy(int n, double *mat, double *temp);
//void printmat(int n, double *mat);
//void printmult(int n, double *mat, double *matinv);
//void rowmult(int n, int i, double C, double *mat, double *matinv);
//void rowswap(int n, int m, int r, double *mat, double *matinv);
//void setidentity(int n, double *mat);


void fillbmat(int k, int d, double *bmat, double *knots)
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

void fillmat(int k, int d, double *mat, double *inputs, double *knots)
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

Array<float> identityMatrix(int n) {
    Array<float> B;
    for (int i=0; i<n*n; i++) {
        B.add(i, 0);
    }
    for (int i=0; i<n; i++) {
        B.set(i*n+i, 1.0);
    }
    return B;
}

void matrixPrint(int n, Array<float>& A) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            std::cout << A[i*n+j];
        }
        std::cout << std::endl;
    }
}

void gaussElim(int n, Array<float>& A, Array<float>& B)
{
    // B will be the inverse of A both assumed n x n with A invertible
    
    int i, j, m=0;

    for (j=0; j<n; j++)
    {
        // get pivot for column j, choosing max from columns
        m = getPivot(n, j, A);
        // cout << "pivot row is m =  " << m << endl;
        // swap row m (with max abs) for row j
        if (m > j) {
           rowSwap(n, m, j, A, B);
        }
        double pivot = A[n*j+j];
            // cout << "pivot value is =  " << pivot << endl;
            if (pivot > 0 || pivot < 0) {
            // get 1 on diagonal
            rowMult(n, j, 1/pivot, A, B);
        } else {
            //  cout << "No pivot for column " << j << endl;
            return;
        }
        for (i=0; i<n; i++) {
            // get zeros in rest of column j
            if (i < j || i > j) {
                   getZero(n, i, j, A, B);  // skips case i=j
                }
        }
    }
    
}

void gausselim(int n, double *mat, double *matinv)
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
      double pivot = mat[n*j+j];
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
int getPivot(int n, int j, Array<float>& A)
{
    int i;
    int m = j;
    float max = abs(A[n*j+j]);
    for (i=j+1; i<n; i++) {
        if (abs(A[n*i+j]) > max) {
            max = abs(A[n*i+j]);
            m = i;
        }
    }
    return m;
}


// find pivot row for col j with max abs value
// assume nxn matrix mat
int getpivot(int n, int j, double *mat)
{
    int i;
    int m = j;
    double max = abs(mat[n*j+j]);
    for (i=j+1; i<n; i++) {
        if (abs(mat[n*i+j]) > max) {
            max = abs(mat[n*i+j]);
            m = i;
        }
    }
    return m;
}

// get zero in col j for row i, for mat nxn
// assumes pivot is 1
void getZero(int n, int i, int j, Array<float>& A, Array<float>& B)
{
    int m;
    float C = A[n*i+j];
    // replace row i with itself plus (-C)*(pivot row)
    for (m=0; m<n; m++) {
        A.set(n*i+m, A[n*i+m] + (-C)*A[n*j+m]);
        B.set(n*i+m, B[n*i+m] + (-C)*B[n*j+m]);
    }
}

// get zero in col j for row i, for mat nxn
// assumes pivot is 1
void getzero(int n, int i, int j, double *mat, double *matinv)
{
    int m;
    double C = mat[n*i+j];
    // replace row i with itself plus (-C)*(pivot row)
    for (m=0; m<n; m++)
    {
	mat[n*i+m] += (-C)*mat[n*j+m];
	matinv[n*i+m] += (-C)*matinv[n*j+m];
    }
}

void initbmat(int k, int d, double *bmat)
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


void matcopy(int n, double *mat, double *temp)
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

void printMatrix(int n, Array<float>& A)
{
    int i, j;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            cout << i << "," << j << " :  " << A[n*i+j] << "  ";
        }
    cout << endl;
    }
}

Array<float> multMatCol(int n, Array<float>& A, Array<float>& x)
{
    Array<float> y;
    for (int i=0; i<n; i++) {
        y.add(0);
    }
    if (A.size() != n * n) {
        DBG("matrix is not nxn");
        return y;
    }
    if (x.size() != n) {
        DBG("vector is not nxn");
        return y;
    }
    int i, j;
    for (i=0; i<n; i++)
    {
        float val = 0;
        for (j=0; j<n; j++)
        {
            val += A[i * n + j] * x[j];
        }
        y.set(i, val);
    }
    return y;
}

Array<float> matMult(int n, Array<float>& A, Array<float>& B)
{
    // multiply nxn matrices A * B = C
    Array<float> C;
    for (int i=0; i<n*n; i++) {
        C.add(0);
    }
    if (A.size() != n*n) {
        DBG("matrix A is not nxn");
        return C;
    }
    if (B.size() != n*n) {
        DBG("vector is not nxn");
        return C;
    }
    int i, j, k;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            float val = 0;
            for (k=0; k<n; k++) {
                val += A[i * n + k] * B[k * n + j];
            }
            C.set(i * n + j, val);
        }
    }
    return C;
}

void printmat(int n, double *mat)
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

void printMult(int n, Array<float>& A, Array<float>& B)
{
   double prod;
   int i, j, m;

   for (i=0; i<n; i++)
   {
      for (j=0; j<n; j++)
      {
         prod = 0;
         for (m=0; m<n; m++)
         {
             prod += A[n*i+m] * B[n*m+j];
         }
         cout << "i,j : " << i << "," << j << " :  " << prod << endl;
      }
    }
}

void printmult(int n, double *mat, double *matinv)
{
   double prod;
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

void rowMult(int n, int i, float C, Array<float>& A, Array<float>& B)
{
    int j;
    for (j=0; j<n; j++) {
        A.set(n*i+j, C*A[n*i+j]);
        B.set(n*i+j, C*B[n*i+j]);
    }
}

void rowmult(int n, int i, double C, double *mat, double *matinv)
{
    int j;
    for (j=0; j<n; j++)
    {
	mat[n*i+j] = C*mat[n*i+j];
	matinv[n*i+j] = C*matinv[n*i+j];
    }
}

// swap row m with row r
void rowSwap(int n, int m, int r, Array<float>& A, Array<float>& B)
{
    float temp;
    int j;
    for (j=0; j<n; j++) {
        temp = A[n*m+j];
        A.set(n*m+j, A[n*r+j]);
        A.set(n*r+j, temp);
        temp = B[n*m+j];
        B.set(n*m+j, B[n*r+j]);
        B.set(n*r+j,  temp);
    }
}

// swap row m with row r
void rowswap(int n, int m, int r, double *mat, double *matinv)
{
    double temp;
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

void setidentity(int n, double *mat)
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

