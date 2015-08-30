/**
 
 Copyright (c) 2009 Christopher Pugh
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 **/

#include "sh.h"

#include <iostream>

double  doubleFactorial(unsigned int in)
{
    double result = 1;
    for(int i = in; i > 1; i-= 2)
    {
        result *= i;
    }
    return result;
}

double factorial(unsigned int in)
{
    double result = 1;
    for(int i = in; i > 0; i--)
    {
        result *= i;
    }
    
    return result;
}


void SHCoeffCountToBandIndex(unsigned int i, int x[2])
{
    int& l = x[0];
    int& m = x[1];

    l = floor(sqrt(i));
    m = (i - (l*l)) - l;
}

///\todo we might be able to use template metaprogramming to automate concatenating the polynomial this way
double sphericalHarmonicsLowOrderFastEvaluate(unsigned int l, int m, double theta, double phi)
{
    double x = sin(phi) * cos(theta);
    double y = sin(phi) * sin(theta);
    double z = cos(phi);
    
    if(l == 1)
    {
        switch(m)
        {
            case -1 : return .488603 * y;
            case 0 : return (.488603 * z);
            case 1 : return (.488603 * x);
        }
        
    }
    if(l == 2)
    {
        switch(m)
        {
                
            case -2 : return (1.092548 * x * y);
            case -1: return (1.092548 * y * z);
            case 0 : return (.315392 * (3 * z * z - 1.0));
            case 1 : return (1.092548 * x * z);
            case 2 : return (0.546274 * (x * x - y * y));
        }
    }
    if( l ==0)
    {
        return 0.282095;
    }
    
    return 0.0;
}

double sphericalHarmonicsEvaluate(unsigned int l, int m, double theta, double phi)
{
    
    
#ifndef SH_NO_FAST_LOW_ORDER_EVALUATE
    if(l <= 2)
    {
        return sphericalHarmonicsLowOrderFastEvaluate(l,m,theta,phi);
    }
#endif
    
    const double root2 =sqrt(2);

	double result = 0.0;

    if(m >0)
    {
        result = root2 * shNormConstant(l,m) * cos(double(m)*theta) * legendreEvaluate(l,m,cos(phi));
    }
    else if(m < 0)
    {
        result = root2 * shNormConstant(l,-m) * sin(double(-m)*theta) * legendreEvaluate(l,-m,cos(phi));
    }
    else
    {

        result= shNormConstant(l, 0) * legendreEvaluate(l, 0, cos(phi));
    }
    
    return result;
}


//evaluates the real legendre polynomial of band l, index m, at x
double legendreEvaluate(unsigned int l, unsigned int m, double x)
{
    double r = 1.0;

	if(m == 0) r = 1; //base case legendre polynomial(0,0)

	else
    {
	   r =  pow(-1.0, m) * doubleFactorial(2*m - 1) * pow(1.0 - x*x, .5 * double(m));
	}

    if(l==m) return r;

    //calculate the legendre polynomial at the point up one band
    double rmPlusOne = x* double(2*m + 1) * r;

    if(l == m+1) return rmPlusOne; //lf l is m+1 we're done

    for(unsigned int i = m+2; i <= l; i++)
    {
       double temp = rmPlusOne;
       rmPlusOne =  x*double(2*i-1)*rmPlusOne - double((i+m-1))*r;
       rmPlusOne /= double(i-m);
       r = temp;
    }

    return rmPlusOne;
}


double cosine_lobe_coefficient_premultiplied(int l)
{
    if(!l)
    {
        return M_PI;
    }
    
    if(l == 1)
    {   const double A1 = (2.0 / 3.0) * M_PI;
        return A1;
    }
    
    if(l & 1) // l > 1, odd
    {
        return 0.0;
    }
    
    double num_b = ((l >>1) - 1) & 1 ? -1.0 : 1.0;
    double den_b = static_cast<double>((l + 2) * (l-1));
    double b = num_b / den_b;
    
    unsigned int fac = factorial(l>>1);
    double c = static_cast<double>(factorial(l)) / ((1<<l) * fac* fac);
    
    
    return (2.0 * M_PI) * b * c;
}

/// unnormalized real sh convolution; useful if your kernel coefficients are already scaled
void convolutionRealSHUnnormalized(double* out, double* f, double* kernel, int bands, int channels_f)
{
    int ctr = 0;
    for(int l = 0; l < bands; l++)
    {
        int zonalHarmonicIndex = indexOfZonalHarmonic(l);
        double kernelVal = kernel[zonalHarmonicIndex];
        
        for(int m = -l; m <= l; m++, ctr++)
        {
            for(int ch = 0; ch < channels_f; ch++)
            {
                out[ctr * channels_f + ch] = kernelVal * f[ctr * channels_f + ch];
            }
        }
    }
}

///evaluate the 1D gaussian curve at a particular point
inline double gaussian1D(double distanceX, double stdev)
{
    double exponent = (-distanceX * distanceX ) / (2.0 * stdev * stdev);
    return  exp(exponent);
}


void gaussianWindowing(double* out, double* f, int bands, int channels_f)
{
    double targetSTDEV = (bands/3.0);
    
    int ctr = 0;
    for(int l = 0; l < bands; l++)
    {
        double kernelVal = gaussian1D(l, targetSTDEV);
        
        for(int m = -l; m <= l; m++, ctr++)
        {
            for(int ch = 0; ch < channels_f; ch++)
            {
                out[ctr * channels_f + ch] = kernelVal * f[ctr * channels_f + ch];
            }
        }
    }
}

void hanningWindowing(double* out, double* f, int bands, int channels_f)
{
    int ctr = 0;
    for(int l = 0; l < bands; l++)
    {
        double kernelVal = .5 *(1.0 + (cos((M_PI * l) / bands)));
        
        for(int m = -l; m <= l; m++, ctr++)
        {
            for(int ch = 0; ch < channels_f; ch++)
            {
                out[ctr * channels_f + ch] = kernelVal * f[ctr * channels_f + ch];
            }
        }
    }
}


void lanczosWindowing(double* out, double* f, int bands, int channels_f)
{
    int ctr = 0;
    for(int l = 0; l < bands; l++)
    {
        double x = std::max<double>(l, .001) * (M_PI / (bands));
    
        double kernelVal = sin(x) / x;

        for(int m = -l; m <= l; m++, ctr++)
        {
            for(int ch = 0; ch < channels_f; ch++)
            {
                out[ctr * channels_f + ch] = kernelVal * f[ctr * channels_f + ch];
            }
        }
    }
}


/// unnormalized real sh convolution; useful if your kernel coefficients are already scaled
void convolutionRealSH(double* out, double* f, double* kernel, int bands, int channels_f)
{
    int ctr = 0;
    for(int l = 0; l < bands; l++)
    {
        int zonalHarmonicIndex = indexOfZonalHarmonic(l);
        double kernelVal = shNormConstantInv(l) * kernel[zonalHarmonicIndex];
        
        for(int m = -l; m <= l; m++, ctr++)
        {
            for(int ch = 0; ch < channels_f; ch++)
            {
                out[ctr * channels_f + ch] = kernelVal * f[ctr * channels_f + ch];
            }
        }
    }
}




