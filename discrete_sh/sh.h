/**
 
 Copyright (c) 2009 Christopher Pugh
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 **/


#ifndef VIRTUOSO_SPHERICAL_HARMONICS_H_INCLUDED
#define VIRTUOSO_SPHERICAL_HARMONICS_H_INCLUDED


#include <math.h>

static const double fourPi = 4* M_PI;

//evaluates the real spherical harmonic of band l, index m, at spherical coordinates theta, phi
double sphericalHarmonicsEvaluate(unsigned int l, int m, double theta,  double phi);

//evaluates the real legendre polynomial of band l, index m, at x
double legendreEvaluate(unsigned int l, unsigned int m, double x);

///converts a coefficient index to an l,m pair
void SHCoeffCountToBandIndex(unsigned int i, int x[2]);

double doubleFactorial(unsigned int in);

//integer factorial function
double factorial( unsigned int in);


double shNormConstant(int l, int m);

inline double shNormConstant(int l)
{
    return sqrt((double)(2 * l + 1) / fourPi);
}

/// gives the SH coefficient for the projection of cos(phi) where phi is the polar angle.  Takes only "l" because this function is radially symmetric and therefore only the zonal harmonics are nonzero
/// for purposes of convolution or band limiting, note that coefficients beyond l = 2 quickly approach 0
double cosine_lobe_coefficient_premultiplied(int l);

/// normalized, eg, for if we want a lambertian BRDF
inline double cosine_lobe_normalized_coefficient_premultiplied(int l)
{
    return cosine_lobe_coefficient_premultiplied(l) / M_PI;
}

inline double cosine_lobe_coefficient(int l)
{
    return shNormConstant(l) * cosine_lobe_coefficient_premultiplied(l);
}

inline double cosine_lobe_normalized_coefficient(int l)
{
    return shNormConstant(l) * cosine_lobe_normalized_coefficient_premultiplied(l);
}

/// gives the SH coefficient for the approximate projection of cos(phi)^pow where phi is the polar angle.  Takes only "l" because this function is radially symmetric and therefore only the zonal harmonics are nonzero
/// the "premultiplied" version is already scaled by a constant necessary for normalized spherical convolution
inline double cosine_lobe_normalized_coefficient_premultiplied(int l, double pow)
{
    return exp((double)(-l * l) / (double)(2 * pow));
}

inline double shNormConstantInv(int l)
{
    return sqrt(fourPi / (2 * l + 1));
}


/// gives the SH coefficient for the approximate projection of cos(phi)^pow where phi is the polar angle.  Takes only "l" because this function is radially symmetric and therefore only the zonal harmonics are nonzero
inline double cosine_lobe_normalized_coefficient(int l, double pow)
{
    return shNormConstantInv(l) * cosine_lobe_normalized_coefficient_premultiplied(l, pow);
}

/// how many bands are desired for accurate SH representation of cos(phi)^pow where phi is the polar angle
/// This is a conservative function that covers all exponents from 1-256 fairly accurately
inline int cosine_lobe_target_bands(double pow, double err = 1e-5)
{
    return sqrt(- 2.0 * pow * std::log(err));
}

inline double shNormConstant(int l, int m)
{
    return sqrt(double ((2*l+1)*factorial(l - m)) / (fourPi * factorial(l+m)));
}

inline double shNormConstantInv(int l, int m)
{
    return sqrt(double ( (fourPi * factorial(l+m)) / ((2*l+1)*factorial(l - m))));
}

inline int totalCoeffs(int l)
{
    int rval = (l + 1);
    return rval * rval;
}

/// zero indexed
inline int indexOfLM(int l, int m)
{
    return (l + 1) * l + m;
}

/// zero indexed
inline int indexOfZonalHarmonic(int l)
{
    return (l + 1) * l;
}

/// unnormalized real sh convolution; useful if your kernel coefficients are already scaled
void convolutionRealSHUnnormalized(double* out, double* f, double* kernel, int bands, int channels_f);

/// unnormalized real sh convolution; useful if your kernel coefficients are already scaled
void convolutionRealSH(double* out, double* f, double* kernel, int bands, int channels_f);

void gaussianWindowing(double* out, double* f, int bands, int channels_f);
void hanningWindowing(double* out, double* f, int bands, int channels_f);
void lanczosWindowing(double* out, double* f, int bands, int channels_f);

#endif
