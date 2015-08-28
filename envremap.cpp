/* Copyright (c) 2010-2013 Robert Kooima                                      */
/*                                                                            */
/* Permission is hereby granted, free of charge, to any person obtaining a    */
/* copy of this software and associated documentation files (the "Software"), */
/* to deal in the Software without restriction, including without limitation  */
/* the rights to use, copy, modify, merge, publish, distribute, sublicense,   */
/* and/or sell copies of the Software, and to permit persons to whom the      */
/* Software is furnished to do so, subject to the following conditions:       */
/*                                                                            */
/* The above copyright notice and this permission notice shall be included in */
/* all copies or substantial portions of the Software.                        */
/*                                                                            */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    */
/* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    */
/* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        */
/* DEALINGS IN THE SOFTWARE.                                                  */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <math.h>

#include "envremap.h"

#define SAMPLEFORMAT_UINT 1
#define SAMPLEFORMAT_INT 2
#define SAMPLEFORMAT_IEEEFP 3
#define SAMPLEFORMAT_VOID 4
#define SAMPLEFORMAT_COMPLEXINT 5
#define SAMPLEFORMAT_COMPLEXIEEEFP 6

/*----------------------------------------------------------------------------*/
/* A small set of single precision mathematical utilities.                    */

#define PI2 1.5707963f
#define PI  3.1415927f
#define TAU 6.2831853f

static inline float lerp(float a, float b, float k)
{
    return a * (1.f - k) + b * k;
}

static inline float clamp(float f, float a, float z)
{
    if      (f < a) return a;
    else if (f > z) return z;
    else            return f;
}

static inline float length(float i, float j)
{
    return sqrtf(i * i + j * j);
}

static inline void normalize(float *v)
{
    const float k = 1.f / sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] *= k;
    v[1] *= k;
    v[2] *= k;
}

inline void add(float *a, const float *b, const float *c)
{
    a[0] = b[0] + c[0];
    a[1] = b[1] + c[1];
    a[2] = b[2] + c[2];
}

/*----------------------------------------------------------------------------*/
/* Pixel block and line copying procedures used for cube map bordering.       */

static void blit(void *dst, int dw, int dx, int dy,
                 void *src, int sw, int sx, int sy,
                        int w, int h, int c, int b)
{
    int i;
    for (i = 0; i < h; i++)
        memcpy((char *) dst + ((i + dy) * dw + dx) * c * b,
               (char *) src + ((i + sy) * sw + sx) * c * b,
                                                 w * c * b);
}

#if 0
static void line(void *dst, int dw, int dx, int dy, int ddx, int ddy,
                 void *src, int sw, int sx, int sy, int sdx, int sdy,
                                                 int n, int c, int b)
{
    int i;
    for (i = 0; i < n; i++)
        memcpy((char *) dst + ((dy + ddy * i) * dw + dx + ddx * i) * c * b,
               (char *) src + ((sy + sdy * i) * sw + sx + sdx * i) * c * b,
                                                                     c * b);
}
#endif

/*----------------------------------------------------------------------------*/

#define SAMP(img, i, j, k) img.p[img.c * (img.w * i + j) + k]

static float *rotN(image *img, int i, int j)
{
    const int ii = i;
    const int jj = j;
    return img->p + img->c * (img->w * ii + jj);
}

static float *rotL(image *img, int i, int j)
{
    const int ii = j;
    const int jj = img->h - i - 1;
    return img->p + img->c * (img->w * ii + jj);
}

static float *rotR(image *img, int i, int j)
{
    const int ii = img->w - j - 1;
    const int jj = i;
    return img->p + img->c * (img->w * ii + jj);
}

typedef float *(*rot)(image *, int, int);

static void border(image *a, rot rota, image *b, rot rotb, int d)
{
    const size_t s = b->c * sizeof (float);
    const int    n = b->h;

    for     (int i = d; i < n - d; i++)
        for (int j = 0; j <     d; j++)
        {
            memcpy(rota(a, i, n - d + j), rotb(b, i,         d + j), s);
            memcpy(rotb(b, i,         j), rota(a, i, n - d - d + j), s);
        }
}

image *image_alloc(int n, int h, int w, int c)
{
    image *img;
    int f;
    
    if ((img = (image *) calloc(n, sizeof (image))))
        for (f = 0; f < n; f++)
        {
            img[f].p = (float *) calloc(w * h * c, sizeof (float));
            img[f].w = w;
            img[f].h = h;
            img[f].c = c;            
        }
    
    return img;
}


/* Add borders to a cubemap image. Assume the given image pointer is an array */
/* of six images. Copy each to a new sef of six images, each two pixels wider */
/* and higher. Also copy the borders. This is necessary for correct cubemap   */
/* sampling.                                                                  */
image *image_border(image *src)
{
    image *dst = 0;
    
    const int d = 1;
    
    if ((src) && (dst = image_alloc(6, src[0].w + 2 * d,
                                    src[0].w + 2 * d,
                                    src[0].c
                                    )))
    {
        const int n = src[0].w;
        const int c = src[0].c;
        const int b = 4;
        
        const int N = n + 2 * d;
        
        /* Copy all page data. */
        
        blit(dst[0].p, N, d, d, src[0].p, n, 0, 0, n, n, c, b);
        blit(dst[1].p, N, d, d, src[1].p, n, 0, 0, n, n, c, b);
        blit(dst[2].p, N, d, d, src[2].p, n, 0, 0, n, n, c, b);
        blit(dst[3].p, N, d, d, src[3].p, n, 0, 0, n, n, c, b);
        blit(dst[4].p, N, d, d, src[4].p, n, 0, 0, n, n, c, b);
        blit(dst[5].p, N, d, d, src[5].p, n, 0, 0, n, n, c, b);
        
        border(dst + 0, rotN, dst + 5, rotN, d);
        border(dst + 5, rotN, dst + 1, rotN, d);
        border(dst + 1, rotN, dst + 4, rotN, d);
        border(dst + 4, rotN, dst + 0, rotN, d);
        
        border(dst + 1, rotR, dst + 2, rotN, d);
        border(dst + 1, rotL, dst + 3, rotN, d);
        
        border(dst + 2, rotN, dst + 0, rotL, d);
        border(dst + 3, rotN, dst + 0, rotR, d);
        
        border(dst + 2, rotL, dst + 4, rotL, d);
        border(dst + 2, rotR, dst + 5, rotL, d);
        border(dst + 3, rotL, dst + 5, rotR, d);
        border(dst + 3, rotR, dst + 4, rotR, d);
        
#if 0
        /* Corner patch hack. */
        
        for     (f = 0; f < 6; f++)
            for (k = 0; k < c; k++)
            {
                SAMP(dst[f], 0, 0, k) = (SAMP(dst[f], 1, 0, k) +
                                         SAMP(dst[f], 0, 1, k) +
                                         SAMP(dst[f], 1, 1, k)) / 3.0f;
                SAMP(dst[f], 0, M, k) = (SAMP(dst[f], 1, M, k) +
                                         SAMP(dst[f], 0, L, k) +
                                         SAMP(dst[f], 1, L, k)) / 3.0f;
                SAMP(dst[f], M, 0, k) = (SAMP(dst[f], L, 0, k) +
                                         SAMP(dst[f], M, 1, k) +
                                         SAMP(dst[f], L, 1, k)) / 3.0f;
                SAMP(dst[f], M, M, k) = (SAMP(dst[f], L, M, k) +
                                         SAMP(dst[f], M, L, k) +
                                         SAMP(dst[f], L, L, k)) / 3.0f;
            }
#endif
    }
    return dst;
}

/*----------------------------------------------------------------------------*/

/* Sample an image at row i column j using linear interpolation.              */

void filter_linear(const image *img, float i, float j, float *p)
{
    const float ii = clamp(i - 0.5f, 0.0f, img->h - 1.0f);
    const float jj = clamp(j - 0.5f, 0.0f, img->w - 1.0f);

    const long  i0 = lrintf(floorf(ii)), i1 = lrintf(ceilf(ii));
    const long  j0 = lrintf(floorf(jj)), j1 = lrintf(ceilf(jj));

    const float di = ii - i0;
    const float dj = jj - j0;

    int k;

    for (k = 0; k < img->c; k++)
        p[k] += lerp(lerp(img->p[(img->w * i0 + j0) * img->c + k],
                          img->p[(img->w * i0 + j1) * img->c + k], dj),
                     lerp(img->p[(img->w * i1 + j0) * img->c + k],
                          img->p[(img->w * i1 + j1) * img->c + k], dj), di);
}

/* Sample an image at row i column j using nearest neighbor.                  */

void filter_nearest(const image *img, float i, float j, float *p)
{
    const float ii = clamp(i - 0.5f, 0.0f, img->h - 1.0f);
    const float jj = clamp(j - 0.5f, 0.0f, img->w - 1.0f);

    const long  i0 = lrintf(ii);
    const long  j0 = lrintf(jj);

    int k;

    for (k = 0; k < img->c; k++)
        p[k] += img->p[(img->w * i0 + j0) * img->c + k];
}

/*----------------------------------------------------------------------------*/

int cube_to_img(int *f, float *i, float *j, int h, int w, const float *v)
{
    const float X = fabsf(v[0]);
    const float Y = fabsf(v[1]);
    const float Z = fabsf(v[2]);

    float x;
    float y;

    if      (v[0] > 0 && X >= Y && X >= Z) { *f = 0; x = -v[2] / X; y = -v[1] / X; }
    else if (v[0] < 0 && X >= Y && X >= Z) { *f = 1; x =  v[2] / X; y = -v[1] / X; }
    else if (v[1] > 0 && Y >= X && Y >= Z) { *f = 2; x =  v[0] / Y; y =  v[2] / Y; }
    else if (v[1] < 0 && Y >= X && Y >= Z) { *f = 3; x =  v[0] / Y; y = -v[2] / Y; }
    else if (v[2] > 0 && Z >= X && Z >= Y) { *f = 4; x =  v[0] / Z; y = -v[1] / Z; }
    else if (v[2] < 0 && Z >= X && Z >= Y) { *f = 5; x = -v[0] / Z; y = -v[1] / Z; }
    else return 0;

    *i = 1.0f + (h - 2) * (y + 1.0f) / 2.0f;
    *j = 1.0f + (w - 2) * (x + 1.0f) / 2.0f;

    return 1;
}

int dome_to_img(int *f, float *i, float *j, int h, int w, const float *v)
{
    if (v[1] >= 0)
    {
        const float d = sqrtf(v[0] * v[0] + v[2] * v[2]);
        const float r = acosf(v[1]) / PI2;

        *f = 0;
        *i = h * (1.0f - r * v[2] / d) / 2.0f;
        *j = w * (1.0f + r * v[0] / d) / 2.0f;

        return 1;
    }
    return 0;
}

int hemi_to_img(int *f, float *i, float *j, int h, int w, const float *v)
{
    if (v[2] <= 0)
    {
        const float d = sqrtf(v[0] * v[0] + v[1] * v[1]);
        const float r = acosf(-v[2]) / PI2;

        *f = 0;
        *i = h * (1.0f - r * v[1] / d) / 2.0f;
        *j = w * (1.0f + r * v[0] / d) / 2.0f;

        return 1;
    }
    return 0;
}

int ball_to_img(int *f, float *i, float *j, int h, int w, const float *v)
{
    const float d = sqrtf(v[0] * v[0] + v[1] * v[1]);
    const float r = sinf(acosf(v[2]) * 0.5f);

    *f = 0;
    *i = h * (1.0f - r * v[1] / d) / 2.0f;
    *j = w * (1.0f + r * v[0] / d) / 2.0f;

    return 1;
}

int rect_to_img(int *f, float *i, float *j, int h, int w, const float *v)
{
    *f = 0;
    *i = h * (       acosf (v[1])        / PI);
    *j = w * (0.5f + atan2f(v[0], -v[2]) / TAU);

    return 1;
}

/*----------------------------------------------------------------------------*/

int cube_to_env(int f, float i, float j, int h, int w, float *v)
{
    const int p[6][3][3] = {
        {{  0,  0, -1 }, {  0, -1,  0 }, {  1,  0,  0 }},
        {{  0,  0,  1 }, {  0, -1,  0 }, { -1,  0,  0 }},
        {{  1,  0,  0 }, {  0,  0,  1 }, {  0,  1,  0 }},
        {{  1,  0,  0 }, {  0,  0, -1 }, {  0, -1,  0 }},
        {{  1,  0,  0 }, {  0, -1,  0 }, {  0,  0,  1 }},
        {{ -1,  0,  0 }, {  0, -1,  0 }, {  0,  0, -1 }},
    };

    const float y = 2.0f * i / h - 1.0f;
    const float x = 2.0f * j / w - 1.0f;

    v[0] = p[f][0][0] * x + p[f][1][0] * y + p[f][2][0];
    v[1] = p[f][0][1] * x + p[f][1][1] * y + p[f][2][1];
    v[2] = p[f][0][2] * x + p[f][1][2] * y + p[f][2][2];

    normalize(v);
    return 1;
}

int dome_to_env(int f, float i, float j, int h, int w, float *v)
{
    const float y = 2.0f * i / h - 1.0f;
    const float x = 2.0f * j / w - 1.0f;

    if (length(x, y) <= 1.0f)
    {
        const float lat = PI2 - PI2 * length(x, y);
        const float lon =             atan2f(x, y);

        v[0] =  sinf(lon) * cosf(lat);
        v[1] =              sinf(lat);
        v[2] = -cosf(lon) * cosf(lat);

        return 1;
    }
    return 0;
}

int hemi_to_env(int f, float i, float j, int h, int w, float *v)
{
    const float y = 2.0f * i / h  - 1.0f;
    const float x = 2.0f * j / w  - 1.0f;

    if (length(x, y) <= 1.0f)
    {
        const float lat = PI2 - PI2 * length(x, y);
        const float lon =             atan2f(x, y);

        v[0] =  sinf(lon) * cosf(lat);
        v[1] = -cosf(lon) * cosf(lat);
        v[2] =             -sinf(lat);

        return 1;
    }
    return 0;
}

int ball_to_env(int f, float i, float j, int h, int w, float *v)
{
    const float y = 2.0f * i / h  - 1.0f;
    const float x = 2.0f * j / w  - 1.0f;

    if (length(x, y) <= 1.0f)
    {
        const float lat = 2.0f * asin(length(x, y));
        const float lon =             atan2f(x, y);

        v[0] =  sinf(lon) * sinf(lat);
        v[1] = -cosf(lon) * sinf(lat);
        v[2] =              cosf(lat);

        return 1;
    }
    return 0;
}

int rect_to_env(int f, float i, float j, int h, int w, float *v)
{
    const float lat = PI2 - PI * i / h;
    const float lon = TAU * j / w - PI;

    v[0] =  sinf(lon) * cosf(lat);
    v[1] =              sinf(lat);
    v[2] = -cosf(lon) * cosf(lat);

    return 1;
}

/*----------------------------------------------------------------------------*/

static int xfm(const float *rot, float *v)
{
    if (rot[0])
    {
        float s = sinf(rot[0] * PI / 180.0f);
        float c = cosf(rot[0] * PI / 180.0f);
        float y = v[1] * c - v[2] * s;
        float z = v[1] * s + v[2] * c;

        v[1] = y; v[2] = z;
    }
    if (rot[1])
    {
        float s = sinf(rot[1] * PI / 180.0f);
        float c = cosf(rot[1] * PI / 180.0f);
        float z = v[2] * c - v[0] * s;
        float x = v[2] * s + v[0] * c;

        v[0] = x; v[2] = z;
    }
    if (rot[2])
    {
        float s = sinf(rot[2] * PI / 180.0f);
        float c = cosf(rot[2] * PI / 180.0f);
        float x = v[0] * c - v[1] * s;
        float y = v[0] * s + v[1] * c;

        v[0] = x; v[1] = y;
    }
    return 1;
}

void supersample(const image   *src,
                 const image   *dst,
                 const pattern *pat,
                 const float   *rot,
                 filter fil, to_img img, to_env env, int f, int i, int j)
{
    int    F;
    float  I;
    float  J;
    int    k;
    int    c = 0;
    float *p = dst[f].p + dst[f].c * (dst[f].w * i + j);

    /* For each sample of the supersampling pattern... */

    for (k = 0; k < pat->n; k++)
    {
        const float ii = pat->p[k].i + i;
        const float jj = pat->p[k].j + j;

        /* Project and unproject giving the source location. Sample there. */

        float v[3];

        if (env( f, ii, jj, dst->h, dst->w, v) && xfm(rot, v) &&
            img(&F, &I, &J, src->h, src->w, v))
        {
#if 1
            fil(src + F, I, J, p);
#else
            p[0] = (v[0] + 1.0f) / 2.0f;
            p[1] = (v[1] + 1.0f) / 2.0f;
            p[2] = (v[2] + 1.0f) / 2.0f;
#endif
            c++;
        }
    }

    /* Normalize the sample. */

    for (k = 0; k < dst->c; k++)
        p[k] /= c;
}

void process(const image   *src,
             const image   *dst,
             const pattern *pat,
             const float   *rot,
             filter fil, to_img img, to_env env, int n)
{
    int i;
    int j;
    int f;

    /* Sample all destination rows, columns, and pages. */

    #pragma omp parallel for private(j, f)
    for         (i = 0; i < dst->h; i++)
        for     (j = 0; j < dst->w; j++)
            for (f = 0; f <      n; f++)
                supersample(src, dst, pat, rot, fil, img, env, f, i, j);
}

/*----------------------------------------------------------------------------*/

point cent_points[] = {
    { 0.5f, 0.5f },
};

point rgss_points[] = {
    { 0.125f, 0.625f },
    { 0.375f, 0.125f },
    { 0.625f, 0.875f },
    { 0.875f, 0.375f },
};

point box2_points[] = {
    { 0.25f, 0.25f },
    { 0.25f, 0.75f },
    { 0.75f, 0.25f },
    { 0.75f, 0.75f },
};

point box3_points[] = {
    { 0.1666667f, 0.1666667f },
    { 0.1666667f, 0.5000000f },
    { 0.1666667f, 0.8333333f },
    { 0.5000000f, 0.1666667f },
    { 0.5000000f, 0.5000000f },
    { 0.5000000f, 0.8333333f },
    { 0.8333333f, 0.1666667f },
    { 0.8333333f, 0.5000000f },
    { 0.8333333f, 0.8333333f },
};

point box4_points[] = {
    { 0.125f, 0.125f },
    { 0.125f, 0.375f },
    { 0.125f, 0.625f },
    { 0.125f, 0.875f },
    { 0.375f, 0.125f },
    { 0.375f, 0.375f },
    { 0.375f, 0.625f },
    { 0.375f, 0.875f },
    { 0.625f, 0.125f },
    { 0.625f, 0.375f },
    { 0.625f, 0.625f },
    { 0.625f, 0.875f },
    { 0.875f, 0.125f },
    { 0.875f, 0.375f },
    { 0.875f, 0.625f },
    { 0.875f, 0.875f },
};

const pattern cent_pattern = {  1, cent_points };
const pattern rgss_pattern = {  4, rgss_points };
const pattern box2_pattern = {  4, box2_points };
const pattern box3_pattern = {  9, box3_points };
const pattern box4_pattern = { 16, box4_points };



