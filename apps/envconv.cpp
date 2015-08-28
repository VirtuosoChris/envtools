#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <math.h>

#include <sstream>
#include <iostream>

#include <fstream>

#include "../envremap.h"
#include "../image/inout.h"

#include <vector>

#include "../sh.hpp"

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s [-i format] [-p pattern] [-f filter] [-n n] [[-o dst] [-d] [-s exponent] [-g width] [-h width] [-l width]]+ src\n"
            "\t-i ... Input (and output) file type: cube, dome, hemi, ball, latlong, llsquare [latlong]\n"
            "\t-p ... Sample pattern: cent, rgss, box2, box3, box4    [rgss]\n"
            "\t-f ... Filter type: nearest, linear                  [linear]\n"
            "\t-n ... Output size [1024]\n"
            "\t-d ... diffuse convolution.\n"
            "\t-s ... specular (phong) convolution.  parameter is the exponent of the phong lobe.\n"
            "\t-g ... gaussian window function.  parameter is the width of the window.\n"
            "\t-h ... hanning window function.  parameter is the width of the window.\n"
            "\t-l ... lanczos window function. parameter is the width of the window.\n"
            "\n\nHints:\nThis program allows chaining multiple convolutions together.\nFor RGBE (.hdr) images, it is advisable to set the filter type to nearest, or to convert to a different format for intermediate data processing.\n\n"
            "\n\n Cube face order is right, left, top, bottom, front back ({+-x,+-y,+-z})\n"
            ,
            exe);
    return 0;
}


enum SupportedConvolutionTypes
{
    DIFFUSE_CONVOLUTION,
    PHONG_CONVOLUTION,
    GAUSSIAN_CONVOLUTION,
    HANNING_CONVOLUTION,
    LANCZOS_CONVOLUTION
};

typedef std::pair<SupportedConvolutionTypes, double> ConvolutionParams;
typedef std::vector <ConvolutionParams> Convolutions;
typedef std::tuple<std::string, Convolutions, int> ConvolutionsForTarget;

bool addConvolution(SupportedConvolutionTypes convType,  std::vector<ConvolutionsForTarget>& convolutionTargets)
{
    try
    {
        double convparameter = std::stod(optarg, 0);
        std::get<1>(convolutionTargets.back()).push_back({convType, convparameter});
        return true;
    }
    catch (const std::exception& e)
    {
        return false;
    }
}


///\todo optimize reprojections for resolution targets smaller than original
bool makeConvolutionTargets(sht<double>& shTransform, std::vector<ConvolutionsForTarget>& convolutionTargets, image* reprojected, const char* o,  int intermediateSize, const pattern* patt, const float* rot, filter fil, to_img img, to_env env, int num, bool outputIsLLSquare=false)
{
    image* originalSH = image_alloc(1, intermediateSize>>1, intermediateSize>>1, reprojected->c);
    shTransform.F.get(originalSH->p, intermediateSize>>1);
    
    for(int i =0; i < convolutionTargets.size(); i++)
    {
        
        image* dst = NULL;
        
        int dim = std::get<2>(convolutionTargets[i]);

        if (!strcmp(o, "latlong"))
        {
            ///\todo release mem
            dst = image_alloc(num, dim, 2*dim, reprojected->c);
        }
        else
        {
            ///\todo release mem
            dst = image_alloc(num, dim, dim, reprojected->c);
        }
        
        Convolutions& convolutions = std::get<1>(convolutionTargets[i]);

        std::clog<<"working on "<<std::get<0>(convolutionTargets[i])<<std::endl;
        
        std::function< double (int)> kernel = [convolutions](int l) -> double
        {
            double kernelVal = 1.0;
            
            for(int c =0; c < convolutions.size(); c++)
            {
                switch(convolutions[c].first)
                {
                    case DIFFUSE_CONVOLUTION:
                        kernelVal *= diffuse_kernel(l);
                        break;
                    case PHONG_CONVOLUTION:
                        kernelVal *= phong_kernel(l, convolutions[c].second);
                        break;
                    case LANCZOS_CONVOLUTION:
                        kernelVal *= lanczos_kernel(l, convolutions[c].second);
                        break;
                    case HANNING_CONVOLUTION:
                        kernelVal *= hanning_kernel(l, convolutions[c].second);
                        break;
                    case GAUSSIAN_CONVOLUTION:
                        kernelVal *= gaussian_kernel(l, convolutions[c].second);
                        break;
                }
            }
            return kernelVal;
        };
        
        // put back the original, unfiltered, projected coefficients for each new output target
        if(i)
        {
            shTransform.F.set(originalSH->p, intermediateSize>>1);
        }
        
        shTransform.F.convolution(kernel);
        
        std::clog<<"creating spatial output for "<<std::get<0>(convolutionTargets[i])<<std::endl;
        
        // project to spat
        shTransform.syn(); /// to spat
        
        shTransform.S.get(reprojected->p, intermediateSize);
        
        // image_writer("REPROJECTED.hdr", reprojected, 1);
        
        // project the convolved image back to the original format and write the output file
        if (reprojected && dst && (reprojected->w == dst->w) && outputIsLLSquare)
        {
            for(int op = 0; op < (reprojected->w * reprojected->h * reprojected->c); op++)
            {
                dst->p[op] = std::max(0.0f, dst->p[op]);
            }
            image_writer(std::get<0>(convolutionTargets[i]).c_str(), reprojected, num);
        }
        else
        {
            process(reprojected, dst, patt, rot, fil, img, env, num);
            
            for(int op = 0; op < (dst->w * dst->h * dst->c); op++)
            {
                dst->p[op] = std::max(0.0f, dst->p[op]);
            }

            image_writer(std::get<0>(convolutionTargets[i]).c_str(), dst, num);
        }
        
        std::clog<<"wrote "<<std::get<0>(convolutionTargets[i])<<std::endl;
        
    } /// end for - targets
    
    return true;
}


int main(int argc, char **argv)
{
    
    try{
        
    
    std::vector<ConvolutionsForTarget> convolutionTargets;
    
    /* Set some default behaviors. */
    
    const char *i = "latlong";
    const char *o = "latlong";
    const char *p = "rgss";
    const char *f = "linear";
    
    float rot[3] = { 0.f, 0.f, 0.f };
    
    int n = 1024;
    int c;
    
    /* Parse the command line options. */
    
    while ((c = getopt(argc, argv, "i:o:p:n:f:d:s:g:h:l:x:y:z:")) != -1)
    {
        switch (c)
    {
            
        case 'i': // for convolution program we want input and output formats to match
        {
            i = optarg;
            o = optarg;
            
            break;
        }
        case 'o': //allow multiple output files generated from different convolutions of the same sh projection
        {
            convolutionTargets.push_back({optarg, Convolutions(), n});
            //currentOutput = optarg;
            break;
        }
        case 'p': p      = optarg;               break;
        case 'f': f      = optarg;               break;
        case 'x': rot[0] = strtod(optarg, 0);    break;
        case 'y': rot[1] = strtod(optarg, 0);    break;
        case 'z': rot[2] = strtod(optarg, 0);    break;
        case 'n':
        {
            n  = strtol(optarg, 0, 0);
            std::get<2>(convolutionTargets.back()) = n;
            break;
        }
        
        //parse convolutions and parameters.  program allows us to chain together multiple convolutions
            
        case 'd': //parameter is the width of the kernel
            if(!addConvolution(DIFFUSE_CONVOLUTION,convolutionTargets))return usage(argv[0]); break;
            /*if(!addConvolution(DIFFUSE_CONVOLUTION,   convolutions))return usage(argv[0]); break;*/
        case 's': if(!addConvolution(PHONG_CONVOLUTION,convolutionTargets))return usage(argv[0]); break;
        case 'g': if(!addConvolution(GAUSSIAN_CONVOLUTION,  convolutionTargets))return usage(argv[0]); break;
        case 'h': if(!addConvolution(HANNING_CONVOLUTION,   convolutionTargets))return usage(argv[0]); break;
        case 'l': if(!addConvolution(LANCZOS_CONVOLUTION,   convolutionTargets))return usage(argv[0]); break;
        
        default: return usage(argv[0]);
    }
    }
    int      num = 1;
    image   *src = 0;
    //image   *dst = 0;
    image   *tmp = 0;
    to_img   img;
    to_env   env;
    filter   fil;
    
    /* Select the sampler. */
    
    if      (!strcmp(f, "linear"))  fil = filter_linear;
    else if (!strcmp(f, "nearest")) fil = filter_nearest;
    else return usage(argv[0]);
    
    /* Read the input image. */
    
    std::vector<const char*> inputFiles;
    inputFiles.insert(inputFiles.begin(), &argv[optind], &argv[argc]);
        
    if (optind + 1 <= argc)
    {
        if (!strcmp(i, "cube"))
        {
            tmp = image_reader(inputFiles, 6);
            src = image_border(tmp);
            img = cube_to_img;
        }
        else if (!strcmp(i, "dome"))
        {
            src = image_reader(inputFiles, 1);
            img = dome_to_img;
        }
        else if (!strcmp(i, "hemi"))
        {
            src = image_reader(inputFiles, 1);
            img = hemi_to_img;
        }
        else if (!strcmp(i, "ball"))
        {
            src = image_reader(inputFiles, 1);
            img = ball_to_img;
        }
        else if (!strcmp(i, "latlong") || !strcmp(i, "llsquare"))
        {
            src = image_reader(inputFiles, 1);
            img = rect_to_img;
        }
        else return usage(argv[0]);
    }
    else return usage(argv[0]);
    
    /* Prepare the output image. */
    
    if (src)
    {
        if (!strcmp(o, "cube"))
        {
            num=6;
           /// dst = image_alloc((num = 6), n, n, src->c);
            env = cube_to_env;
        }
        else if (!strcmp(o, "dome"))
        {
           /// dst = image_alloc((num = 1), n, n, src->c);
            env = dome_to_env;
        }
        else if (!strcmp(o, "hemi"))
        {
            ///dst = image_alloc((num = 1), n, n, src->c);
            env = hemi_to_env;
        }
        else if (!strcmp(o, "ball"))
        {
           /// dst = image_alloc((num = 1), n, n, src->c);
            env = ball_to_env;
        }
        else if (!strcmp(o, "latlong"))
        {
            ///dst = image_alloc((num = 1), n, 2 * n, src->c);
            env = rect_to_env;
        }
        else if (!strcmp(o, "llsquare"))
        {
          ///  dst = image_alloc((num = 1), n, n, src->c);
            env = rect_to_env;
        }
        else return usage(argv[0]);
    }
    else
    {
        throw std::runtime_error("Unable to read input image");
    }
        
    const pattern* patt;
    
    if (!strcmp(p, "cent"))
        patt = &cent_pattern;
    
    else if (!strcmp(p, "rgss"))
        patt = &rgss_pattern;
    
    else if (!strcmp(p, "box2"))
        patt = &box2_pattern;
    
    else if (!strcmp(p, "box3"))
        patt = &box3_pattern;
    
    else if (!strcmp(p, "box4"))
        patt = &box4_pattern;
    
    else return usage(argv[0]);
    
    /**** perform SH projection and convolution ***/
    
    // to do convolution we will first convert to llsquare.
    // When we are done we will project the convolved LLSquare back to the original format
    unsigned int intermediateSize = std::max(src->w, src->h);
        
    if(!strcmp(o, "cube")) // cube face is NxN.  Horizontally, we have 4N pixels to cover full 2PI radians, and vertically we have 2N pixels to cover.  To preserve info we have 4Nx4N for the square lat long map
    {
        intermediateSize -=2; //account for padding
        intermediateSize *= 4;
    }
    
    image* intermediate=NULL;
    
    if (!strcmp(i, "llsquare"))
    {
        intermediate = src;
    }
    else
    {
        intermediate = image_alloc(1, intermediateSize, intermediateSize, src->c);
        process(src, intermediate, patt, rot, fil, img, rect_to_env, 1);
    }
        
    /***/
     //   image_writer("intermediate.hdr", intermediate, 1);
    /***/
        
    
    /// create intermediate, llsquare image
    
    image* reprojected = image_alloc(1, intermediateSize, intermediateSize, src->c);
    
    sht<double> shTransform(intermediateSize>>1, src->c);
    
    shTransform.S.set(intermediate->p);
    
    std::clog<<"Projecting input to SH."<<std::endl;
    
    shTransform.ana();  // to sh
    
    std::clog<<"SH Projection done"<<std::endl;;
    
    if(!makeConvolutionTargets(shTransform, convolutionTargets,reprojected, /*dst,*/ o, intermediateSize, patt, rot, fil, rect_to_img, env, num, !strcmp(o, "llsquare")))
    {
        return usage(argv[0]);
    }
        
    }catch(const std::runtime_error& e)
    {
        std::clog<<e.what()<<std::endl;
    }
    
    return 0;
}