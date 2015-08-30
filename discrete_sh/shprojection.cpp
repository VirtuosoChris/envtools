/**
 
Copyright (c) 2009 Christopher Pugh

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
**/

#include "shprojection.h"
#include "sh.h"
#include <image.h>

void createLatLongMap(float* pixels, unsigned int width, unsigned int height, unsigned int channels, double* coefficients, unsigned int numCoefficients)
{
    const unsigned int numPixels = width*height;

    for(unsigned int pixel = 0; pixel < numPixels; pixel++)
    {
        //clear the pixel prior to summation
        for(unsigned int ch = 0; ch < channels; ch++)
        {
            pixels[pixel * channels + ch] =0.0;
        }

        double u = double(pixel%width) / width;
        double v = double(pixel /width) / height;

        double th =  2*3.14159*u;
        double phi = 3.14159*v;

        int band_index[2]; //the band and index array to store for the SH coefficients
        
        for( int coeff = 0; coeff < numCoefficients; coeff++)
        {
            SHCoeffCountToBandIndex(coeff,band_index);

            double sh = sphericalHarmonicsEvaluate(band_index[0],band_index[1], th, phi);

            for(unsigned int ch = 0; ch < channels; ch++)
            {
                pixels[pixel * channels + ch] += coefficients[coeff*channels +ch] * sh;
            }

        }


    }
}

FloatImage createLatLongMap(const FloatImage::index_type& dimensions, double* coefficients, unsigned int numCoefficients)
{
    FloatImage latLong(dimensions);

    const FloatImage::index_type dims = latLong.getDimensions();

    createLatLongMap(latLong.dataPtr(),dims[1],dims[2],dims[0],coefficients,numCoefficients);

    return latLong;

}

void projectLatLongMap(const FloatImage& latLong, double* coefficientsBegin, unsigned int coefficientsPerChannel)
{

    const unsigned int& channels = latLong.getDimensions()[0];
    const unsigned int& width = latLong.getDimensions()[1];
    const unsigned int& height = latLong.getDimensions()[2];

    projectLatLongMap(latLong.dataPtr(), width, height, channels, coefficientsBegin, coefficientsPerChannel);

}

void projectLatLongMap(const float* image, unsigned int width, unsigned int height, unsigned int channels, double* coefficients, unsigned int numCoefficients)
{

    const unsigned int numPixels = width*height;

    //zero the coefficients prior to summation
    for(unsigned int i = 0; i < numCoefficients * channels; i++)
    {
        coefficients[i] = 0.0;
    }

    int band_index[2]; //the band and index array to store for the SH coefficients

    //iterate through the pixels and compute the coefficients
    for(unsigned int pixel = 0; pixel < numPixels; pixel++)
    {
        //get normalized image coordinates

        double u = double(pixel % width)  / width;
        double v = double(pixel / width ) / height;

        //convert normalized image coordinates to angular coordinates
        double th =  2*3.14159*u;
        double phi = 3.14159*v;
        
        double dTh = 1.0 / (2 * M_PI);
        double dPhi = 1.0 / M_PI;

        for( int coeff = 0; coeff < numCoefficients; coeff++)
        {
            SHCoeffCountToBandIndex(coeff,band_index);

            double add =  (1.0 / numPixels) * sphericalHarmonicsEvaluate(band_index[0],band_index[1],th,phi);
        
            
            add *= sin(phi);
            add /= dTh * dPhi;

            for(unsigned int ch = 0; ch < channels; ch++)
            {
                float imageRead = image[pixel*channels + ch];

                coefficients[ coeff*channels + ch] += add * imageRead;
            }

        }
    }
}
