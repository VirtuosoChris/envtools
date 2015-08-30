/**
 
 Copyright (c) 2009 Christopher Pugh
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 **/


#ifndef VIRTUOSO_SH_PROJECTION_H_INCLUDED
#define VIRTUOSO_SH_PROJECTION_H_INCLUDED

#include <image.h>

///it's important to know that if you're using rgb, each "coefficient" is an rgb triple, so you need to pass in as coefficients an array that holds numCoefficients*channels values

void projectLatLongMap(const float* image, unsigned int width, unsigned int height, unsigned int channels, double* coefficients, unsigned int numCoefficients);

void createLatLongMap(float* pixels, unsigned int width, unsigned int height, unsigned int channels, double* coefficients, unsigned int numCoefficients);


void projectLatLongMap(const FloatImage& latLong, double* coefficientsBegin, unsigned int numCoefficients);

FloatImage createLatLongMap(const FloatImage::index_type&, double* coefficientsBegin, unsigned int numCoefficients);


#endif
