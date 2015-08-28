// Copyright (c) 2013 Robert Kooima
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include <getopt.h>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../sh.hpp"
#include "../image/inout.h"


//------------------------------------------------------------------------------
#if 0
template <typename real> void grace_cathedral(sht<real>& T)
{
    T.F(0,  0, 0) =  0.79; T.F(0,  0, 1) =  0.44; T.F(0,  0, 2) =  0.54;
    T.F(1, -1, 0) =  0.39; T.F(1, -1, 1) =  0.35; T.F(1, -1, 2) =  0.60;
    T.F(1,  0, 0) = -0.34; T.F(1,  0, 1) = -0.18; T.F(1,  0, 2) = -0.27;
    T.F(1,  1, 0) = -0.29; T.F(1,  1, 1) = -0.06; T.F(1,  1, 2) =  0.01;
    T.F(2, -2, 0) = -0.11; T.F(2, -2, 1) = -0.05; T.F(2, -2, 2) = -0.12;
    T.F(2, -1, 0) = -0.26; T.F(2, -1, 1) = -0.22; T.F(2, -1, 2) = -0.47;
    T.F(2,  0, 0) = -0.16; T.F(2,  0, 1) = -0.09; T.F(2,  0, 2) = -0.15;
    T.F(2,  1, 0) =  0.56; T.F(2,  1, 1) =  0.21; T.F(2,  1, 2) =  0.14;
    T.F(2,  2, 0) =  0.21; T.F(2,  2, 1) = -0.05; T.F(2,  2, 2) = -0.30;
}

template <typename real> void eucalyptus_grove(sht<real>& T)
{
    T.F(0,  0, 0) =  0.38; T.F(0,  0, 1) =  0.43; T.F(0,  0, 2) =  0.45;
    T.F(1, -1, 0) =  0.29; T.F(1, -1, 1) =  0.36; T.F(1, -1, 2) =  0.41;
    T.F(1,  0, 0) =  0.04; T.F(1,  0, 1) =  0.03; T.F(1,  0, 2) =  0.01;
    T.F(1,  1, 0) = -0.10; T.F(1,  1, 1) = -0.10; T.F(1,  1, 2) = -0.09;
    T.F(2, -2, 0) = -0.06; T.F(2, -2, 1) = -0.06; T.F(2, -2, 2) = -0.04;
    T.F(2, -1, 0) =  0.01; T.F(2, -1, 1) = -0.01; T.F(2, -1, 2) = -0.05;
    T.F(2,  0, 0) = -0.09; T.F(2,  0, 1) = -0.13; T.F(2,  0, 2) = -0.15;
    T.F(2,  1, 0) = -0.06; T.F(2,  1, 1) = -0.05; T.F(2,  1, 2) = -0.04;
    T.F(2,  2, 0) =  0.02; T.F(2,  2, 1) = -0.00; T.F(2,  2, 2) = -0.05;
}

template <typename real> void st_peters_basilica(sht<real>& T)
{
    T.F(0,  0, 0) =  0.36; T.F(0,  0, 1) =  0.26; T.F(0,  0, 2) =  0.23;
    T.F(1, -1, 0) =  0.18; T.F(1, -1, 1) =  0.14; T.F(1, -1, 2) =  0.13;
    T.F(1,  0, 0) = -0.02; T.F(1,  0, 1) = -0.01; T.F(1,  0, 2) = -0.00;
    T.F(1,  1, 0) =  0.03; T.F(1,  1, 1) =  0.02; T.F(1,  1, 2) =  0.01;
    T.F(2, -2, 0) =  0.02; T.F(2, -2, 1) =  0.01; T.F(2, -2, 2) =  0.00;
    T.F(2, -1, 0) = -0.05; T.F(2, -1, 1) = -0.03; T.F(2, -1, 2) = -0.01;
    T.F(2,  0, 0) = -0.09; T.F(2,  0, 1) = -0.08; T.F(2,  0, 2) = -0.07;
    T.F(2,  1, 0) =  0.01; T.F(2,  1, 1) =  0.00; T.F(2,  1, 2) =  0.00;
    T.F(2,  2, 0) = -0.08; T.F(2,  2, 1) = -0.06; T.F(2,  2, 2) =  0.00;
}
#endif

//------------------------------------------------------------------------------

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s [-o output] [-b b] [-l l] [-m m] [-n n]\n"
            "\t-o ... Output file name (out.tif)\n"
            "\t-b ... Output depth     (4)\n"
            "\t-c ... Output channels  (3)\n"
            "\t-l ... Harmonic degree  (0)\n"
            "\t-m ... Harmonic order   (0)\n"
            "\t-n ... Synthesis degree (256)\n\n", exe);

    return -1;
}

int main(int argc, char **argv)
{
    // Set default options.
    try{
    const char *out = "out.tif";
    int         n   = 256;
    int         b   = 1;
    int         c   = 3;
    int         l   = 0;
    int         m   = 0;

    // Parse the options.

    int o;

    while ((o = getopt(argc, argv, "b:c:l:m:n:o:")) != -1)
        switch (o)
        {
            case 'b': b = strtol(optarg, 0, 0); break;
            case 'c': c = strtol(optarg, 0, 0); break;
            case 'l': l = strtol(optarg, 0, 0); break;
            case 'm': m = strtol(optarg, 0, 0); break;
            case 'n': n = strtol(optarg, 0, 0); break;
            case 'o': out = optarg;             break;

            default: return usage(argv[0]);
        }

    // Confirm a reasonable output request.

    if (n > 0 && b > 0 && m >= -l && m <= l && l < n)
    {
        // Allocate source and destination buffers.

        if (float *dst = (float *) calloc(4 * n * n * c, sizeof (float)))
        {
            // Instance the transformer and construct the input.

            sht<long double> T(n, c);

            if (c == 1)
                T.F(l, m, 0) =  1;
            else
            {
                T.F(l, m, 0) = -1;
                T.F(l, m, 1) =  1;
            }

            // Synthesize and store the output.

            T.syn();
            T.S.get(dst, 2 * n);

            
            ///(out, 2 * n, 2 * n, c, b, dst);
            image img;
            img.w = img.h = 2*n;
            img.c = c;
            img.p = dst;
           
            image_writer(out, &img, 1);

            free(dst);
        }
    }
    else usage(argv[0]);
    }catch(const std::runtime_error& e)
    {
        std::clog<<e.what()<<std::endl;
    }
    return 0;
}

//------------------------------------------------------------------------------
