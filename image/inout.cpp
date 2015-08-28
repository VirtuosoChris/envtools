
#include "inout.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <tiffio.h>


#include <stdlib.h>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#include "gray.h"
#include "sRGB.h"

#ifdef TIFF_SUPPORT

static inline float clamp(float f, float a, float z)
{
    if      (f < a) return a;
    else if (f > z) return z;
    else            return f;
}


/* Read one scanline from the given TIFF file, converting it from the format */
/* of the TIFF to float. The file must have contiguous planar configuration.  */
/* This is a drop-in replacement for TIFFReadScanline.                        */

int TIFFReadFloatScanline(TIFF *T, float *dst, uint32 r)
{
    static tdata_t *src = NULL;
    static tsize_t  len = 0;
    
    uint32 w = 0;
    uint16 c = 0;
    uint16 b = 0;
    uint16 s = 0;
    uint16 p = 0;
    
    TIFFGetField(T, TIFFTAG_IMAGEWIDTH,      &w);
    TIFFGetField(T, TIFFTAG_SAMPLESPERPIXEL, &c);
    TIFFGetField(T, TIFFTAG_BITSPERSAMPLE,   &b);
    TIFFGetField(T, TIFFTAG_SAMPLEFORMAT,    &s);
    TIFFGetField(T, TIFFTAG_PLANARCONFIG,    &p);
    
    if (p == PLANARCONFIG_CONTIG)
    {
        if (len != TIFFScanlineSize(T))
        {
            len  = TIFFScanlineSize(T);
            src  = (tdata_t*)realloc(src, len);
        }
        
        if (src && TIFFReadScanline(T, src, r, 0) > 0)
        {
            if      ((b ==  8) && (s == SAMPLEFORMAT_UINT || s == 0))
                for (uint32 i = 0; i < w * c; i++)
                    dst[i] = ((uint8  *) src)[i] / 255.0f;
            
            else if ((b ==  8) && (s == SAMPLEFORMAT_INT))
                for (uint32 i = 0; i < w * c; i++)
                    dst[i] = ((int8   *) src)[i] / 127.0f;
            
            else if ((b == 16) && (s == SAMPLEFORMAT_UINT || s == 0))
                for (uint32 i = 0; i < w * c; i++)
                    dst[i] = ((uint16 *) src)[i] / 65535.0f;
            
            else if ((b == 16) && (s == SAMPLEFORMAT_INT))
                for (uint32 i = 0; i < w * c; i++)
                    dst[i] = ((int16  *) src)[i] / 32767.0f;
            
            else if ((b == 32) && (s == SAMPLEFORMAT_IEEEFP))
                for (uint32 i = 0; i < w * c; i++)
                    dst[i] = ((float  *) src)[i];
            
            else return -1;
        }
    }
    return +1;
}

/* Write one scanline to the given TIFF file, converting it from float to the */
/* format of the TIFF. The file must have contiguous planar configuration.    */
/* This is a drop-in replacement for TIFFWriteScanline.                       */

int TIFFWriteFloatScanline(TIFF *T, float *src, uint32 r)
{
    static tdata_t *dst = NULL;
    static tsize_t  len = 0;
    
    uint32 w = 0;
    uint16 c = 0;
    uint16 b = 0;
    uint16 s = 0;
    uint16 p = 0;
    
    TIFFGetField(T, TIFFTAG_IMAGEWIDTH,      &w);
    TIFFGetField(T, TIFFTAG_SAMPLESPERPIXEL, &c);
    TIFFGetField(T, TIFFTAG_BITSPERSAMPLE,   &b);
    TIFFGetField(T, TIFFTAG_SAMPLEFORMAT,    &s);
    TIFFGetField(T, TIFFTAG_PLANARCONFIG,    &p);
    
    if (p == PLANARCONFIG_CONTIG)
    {
        if (len != TIFFScanlineSize(T))
        {
            len  = TIFFScanlineSize(T);
            dst  = (tdata_t*)realloc(dst, len);
        }
        
        if (dst)
        {
            if      ((b ==  8) && (s == SAMPLEFORMAT_UINT || s == 0))
                for (uint32 i = 0; i < w * c; i++)
                    ((uint8  *) dst)[i] = clamp(src[i], 0.0f, 1.0f) * 255.0f;
            
            else if ((b ==  8) && (s == SAMPLEFORMAT_INT))
                for (uint32 i = 0; i < w * c; i++)
                    ((int8   *) dst)[i] = clamp(src[i], 0.0f, 1.0f) * 127.0f;
            
            else if ((b == 16) && (s == SAMPLEFORMAT_UINT || s == 0))
                for (uint32 i = 0; i < w * c; i++)
                    ((uint16 *) dst)[i] = clamp(src[i], 0.0f, 1.0f) * 65535.0f;
            
            else if ((b == 16) && (s == SAMPLEFORMAT_INT))
                for (uint32 i = 0; i < w * c; i++)
                    ((int16  *) dst)[i] = clamp(src[i], 0.0f, 1.0f) * 32767.0f;
            
            else if ((b == 32) && (s == SAMPLEFORMAT_IEEEFP))
                for (uint32 i = 0; i < w * c; i++)
                    ((float  *) dst)[i] = src[i];
            
            else return -1;
            
            if (TIFFWriteScanline(T, dst, r, 0) == -1)
                return -1;
        }
    }
    return +1;
}

static image* tiff_reader(const char *name, int n)
{
    image *in = 0;
    TIFF  *T  = 0;
    int    f;
    
    if ((T = TIFFOpen(name, "r")))
    {
        if ((in = (image *) calloc(n, sizeof (image))))
        {
            for (f = 0; f < n; f++)
            {
                uint16 b, c, s = 0;
                uint32 w, h, r;
                float *p;
                
                TIFFSetDirectory(T, f);
                
                TIFFGetField(T, TIFFTAG_IMAGEWIDTH,      &w);
                TIFFGetField(T, TIFFTAG_IMAGELENGTH,     &h);
                TIFFGetField(T, TIFFTAG_SAMPLESPERPIXEL, &c);
                TIFFGetField(T, TIFFTAG_BITSPERSAMPLE,   &b);
                TIFFGetField(T, TIFFTAG_SAMPLEFORMAT,    &s);
                
                if ((p = (float *) malloc(h * w * c * sizeof (float))))
                {
                    for (r = 0; r < h; ++r)
                        TIFFReadFloatScanline(T, p + w * c * r, r);
                    
                    in[f].p = (float *) p;
                    in[f].w = (int)     w;
                    in[f].h = (int)     h;
                    in[f].c = (int)     c;
                }
            }
        }
        TIFFClose(T);
    }
    return in;
}

/* Write n pages to the named TIFF image file.                                */

static void tiff_writer(const char *name, image *out, int n)
{
    TIFF  *T = 0;
    uint32 f, r;
    
    if ((T = TIFFOpen(name, "w")))
    {
        for (f = 0; f < n; ++f)
        {
            TIFFSetField(T, TIFFTAG_IMAGEWIDTH,      out[f].w);
            TIFFSetField(T, TIFFTAG_IMAGELENGTH,     out[f].h);
            TIFFSetField(T, TIFFTAG_SAMPLESPERPIXEL, out[f].c);
            TIFFSetField(T, TIFFTAG_BITSPERSAMPLE,   32);
            TIFFSetField(T, TIFFTAG_SAMPLEFORMAT,    SAMPLEFORMAT_IEEEFP);
            TIFFSetField(T, TIFFTAG_ORIENTATION,     ORIENTATION_TOPLEFT);
            TIFFSetField(T, TIFFTAG_PLANARCONFIG,    PLANARCONFIG_CONTIG);
            
            if (out[f].c == 1)
            {
                TIFFSetField(T, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_MINISBLACK);
                TIFFSetField(T, TIFFTAG_ICCPROFILE, sizeof (grayicc), grayicc);
            }
            else
            {
                TIFFSetField(T, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_RGB);
                TIFFSetField(T, TIFFTAG_ICCPROFILE, sizeof (sRGBicc), sRGBicc);
            }
            
            for (r = 0; r < out[f].h; ++r)
                TIFFWriteFloatScanline(T, out[f].p + out[f].w * out[f].c * r, r);
            
            TIFFWriteDirectory(T);
        }
        TIFFClose(T);
    }
}

#endif

#ifdef EXR_SUPPORT

#include <OpenEXR/ImfCRgbaFile.h>

void *image_read_exr(const char *name, int *w, int *h, int *c, int *b)
{
    ImfInputFile    *file;
    ImfRgba         *data;
    const ImfHeader *head;
    
    float *p = 0;
    
    if ((file = ImfOpenInputFile(name)))
    {
        if ((head = ImfInputHeader(file)))
        {
            int x0;
            int x1;
            int y0;
            int y1;
            int i;
            
            /* Read and extract header info. */
            
            ImfHeaderDataWindow(head, &x0, &y0, &x1, &y1);
            
            *w = x1 - x0 + 1;
            *h = y1 - y0 + 1;
            *c = 4;
            *b = sizeof (float);
            
            /* Allocate temporary storage and read pixel data to it. */
            
            if ((data = (ImfRgba *) malloc((*w) * (*h) * sizeof (ImfRgba))))
            {
                ImfInputSetFrameBuffer(file, data - x0 - y0 * (*w), 1, (*w));
                ImfInputReadPixels    (file, y0, y1);
                
                /* Allocate final storage and copy pixel data to it. */
                
                if ((p = (float *) malloc((*w) * (*h) * (*c) * (*b))))
                {
                    for (i = 0; i < (*w) * (*h); ++i)
                    {
                        p[i * 4 + 0] = ImfHalfToFloat(data[i].r);
                        p[i * 4 + 1] = ImfHalfToFloat(data[i].g);
                        p[i * 4 + 2] = ImfHalfToFloat(data[i].b);
                        p[i * 4 + 3] = ImfHalfToFloat(data[i].a);
                    }
                }
                free(data);
            }
        }
        ImfCloseInputFile(file);
    }
    return p;
}

/// writes half float image in exr format
void image_write_exr(const char *name, int w, int h, int c, void *p)
{
    ImfOutputFile *file;
    ImfRgba       *data;
    ImfHeader     *head;
    
    const float *q = (const float *) p;
    
    /* Allocation and intialize a new header. */
    
    if ((head = ImfNewHeader()))
    {
        ImfHeaderSetDataWindow   (head, 0, 0, w - 1, h - 1);
        ImfHeaderSetDisplayWindow(head, 0, 0, w - 1, h - 1);
        ImfHeaderSetCompression  (head, IMF_ZIP_COMPRESSION);
        ImfHeaderSetLineOrder    (head, IMF_INCREASING_Y);
        
        if ((file = ImfOpenOutputFile(name, head, IMF_WRITE_RGBA)))
        {
            /* Allocate temporary storage and copy pixel data to it. */
            
            if ((data = (ImfRgba *) calloc(w * h, sizeof (ImfRgba))))
            {
                int i;
                
                for (i = 0; i < w * h; ++i)
                {
                    float R = (c > 0) ? q[i * c + 0] : 0.0f;
                    float G = (c > 1) ? q[i * c + 1] : 0.0f;
                    float B = (c > 2) ? q[i * c + 2] : 0.0f;
                    float A = (c > 3) ? q[i * c + 3] : 1.0f;
                    
                    ImfFloatToHalf(R, &data[i].r);
                    ImfFloatToHalf(G, &data[i].g);
                    ImfFloatToHalf(B, &data[i].b);
                    ImfFloatToHalf(A, &data[i].a);
                }
                
                /* Write the file. */
                
                ImfOutputSetFrameBuffer(file, data, 1, w);
                ImfOutputWritePixels   (file, h);
                
                free(data);
            }
            ImfCloseOutputFile(file);
        }
        ImfDeleteHeader(head);
    }
}

#endif

static bool extcmp(const char *name, const char *ext)
{
    int off = strlen(name) - strlen(ext);
    if(off<0)return false;
    
    int rval = strcmp(name + off, ext);
    
    return rval;
    
}

std::string splice(const std::string& str, const std::string& ext, const std::string& splice)
{
    std::string str2(str);
    std::string newsuff = splice + ext;
    std::size_t loc = str.find(ext, str.length() - ext.length() -1 );
    return str2.replace(loc, ext.length(), newsuff);
}


image *image_reader(const std::vector<const char *>& inputFiles, int n)
{
    if(!(inputFiles.size() ==1 || inputFiles.size() == 6))
    {
        throw std::runtime_error("Error: Supported input configuration is one input file (cubes supported in multilayer tif) or 6 input files for each cube face.");
    }
    
    //input 6 layer tif for cube map
    if((inputFiles.size() ==1) && (n ==6))
    {
        if(!extcmp(inputFiles[0], ".tif") || !extcmp(inputFiles[0], ".tiff"))
        {
            return tiff_reader(inputFiles[0], n);
        }
    }
    
    image* rval = (image*)calloc(n, sizeof(image));
    
    for(int page =0; page < n; page++)
    {
        const char* name = inputFiles[page];
        image* im = &rval[page];
        
        if(!extcmp(name, ".hdr"))
        {
            im->p = stbi_loadf(name, &im->w, &im->h, &im->c, 0);
        }
#ifdef EXR_SUPPORT
        else if(!extcmp(name,".exr"))
        {
            int b;
            im->p = (float*)image_read_exr(name, &im->w, &im->h, &im->c, &b);
        }
#endif
#ifdef TIFF_SUPPORT
        else if(!extcmp(name, ".tif") || !extcmp(name, ".tiff"))
        {
            image* readVal = tiff_reader(name, n);
            *im = *readVal;
            free(readVal);
        }
#endif
        else if (!extcmp(name, ".bin"))
        {
            std::ifstream file(name);
            
            file.read((char*)(&(im->w)), sizeof(int));
            file.read((char*)(&(im->h)), sizeof(int));
            file.read((char*)(&(im->c)), sizeof(int));
            
            std::size_t bytesToRead = sizeof(float) * im->w * im->h * im->c;
            im->p = (float*)malloc(bytesToRead);
            
            file.read((char*)(im->p), bytesToRead);
        }
    }
    
    return rval;
}

image *image_reader(const char *name, int n)
{
    if(!extcmp(name, ".hdr"))
    {
        image* rval = (image*)calloc(1, sizeof(image));
        rval->p = stbi_loadf(name, &rval->w, &rval->h, &rval->c, 0);
        return rval;
    }
#ifdef EXR_SUPPORT
    else if(!extcmp(name,".exr"))
    {
        image* rval = (image*) calloc(1, sizeof(image));
        
        int b;
        rval->p = (float*)image_read_exr(name, &rval->w, &rval->h, &rval->c, &b);
        return rval;
        
    }
#endif
#ifdef TIFF_SUPPORT
    else if(!extcmp(name, ".tif") || !extcmp(name, ".tiff"))
    {
        return tiff_reader(name, n);
    }
#endif
    else if (!extcmp(name, ".bin"))
    {
        image* rval = (image*)calloc(1, sizeof(image));
        
        std::ifstream file(name);
        
        file.read((char*)(&(rval->w)), sizeof(int));
        file.read((char*)(&(rval->h)), sizeof(int));
        file.read((char*)(&(rval->c)), sizeof(int));
        
        std::size_t bytesToRead = sizeof(float) * rval->w * rval->h * rval->c;
        rval->p = (float*)malloc(bytesToRead);
        
        file.read((char*)(rval->p), bytesToRead);
        return rval;
    }
    
    return NULL;
}


void image_writer(const char *name, image *out, int n)
{
    const char* CUBE_FACES[] = {"_RIGHT", "_LEFT", "_TOP", "_BOTTOM", "_FRONT", "_BACK"};
    
    /// absolute value each image face
    for(int k =0; k < n; k++)
    {
        for(int i =0; i < out->w * out->h * out->c; i++)
        {
            out[k].p[i] = std::max<float>(0.0,out[k].p[i]);
        }
    }
    
    if(!extcmp(name, ".hdr"))
    {
        if(n != 1)
        {
            for(int i =0; i < n; i++)
            {
                
                std::string filename = splice(name, ".hdr", CUBE_FACES[i]);
                stbi_write_hdr(filename.c_str(), out[i].w, out[i].h, out[i].c, out[i].p);
            }
        }
        else
        {
            stbi_write_hdr(name, out[0].w, out[0].h, out[0].c, out[0].p);
        }
    }
    
#ifdef EXR_SUPPORT
    else if(!extcmp(name,".exr"))
    {
        if(n != 1)
        {
            for(int i =0; i < n; i++)
            {
                
                std::string filename = splice(name, ".exr", CUBE_FACES[i]);
                image_write_exr(filename.c_str(), out->w, out->h, out->c, out->p);
            }
        }
        else
        {
            image_write_exr(name, out->w, out->h, out->c, out->p);
        }
    }
#endif
#ifdef TIFF_SUPPORT
    else if(!extcmp(name, ".tif") || !extcmp(name, ".tiff"))
    {
        tiff_writer(name, out, n);
    }
#endif
    else if (!extcmp(name, ".bin"))
    {
        std::ofstream file(name);
        
        file.write((char*)(&(out->w)), sizeof(int));
        file.write((char*)(&(out->h)), sizeof(int));
        file.write((char*)(&(out->c)), sizeof(int));
        
        std::size_t bytesToRead = sizeof(float) * out->w * out->h * out->c;
        
        file.write((char*)(out->p), bytesToRead);
    }
    else throw std::runtime_error("Unable to write output image: unknown extension");
}