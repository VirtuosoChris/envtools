#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <math.h>

#include <sstream>
#include <iostream>

#include <fstream>

#include <vector>

#include "../envremap.h"
#include "../image/inout.h"

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s [-i input] [-o output] [-p pattern] [-f filter] [-n n] src dst\n"
            "\t-i ... Input  file type: cube, dome, hemi, ball, latlong, llsquare  [latlong]\n"
            "\t-o ... Output file type: cube, dome, hemi, ball, latlong, llsquare  [latlong]\n"
            "\t-p ... Sample pattern: cent, rgss, box2, box3, box4    [rgss]\n"
            "\t-f ... Filter type: nearest, linear                  [linear]\n"
            "\t-n ... Output size                                     [1024]\n"
            "\n\n Cube face order is right, left, top, bottom, front back ({+-x,+-y,+-z})\n"
            ,
            exe);
    return 0;
}

int main(int argc, char **argv)
{
try{
    /* Set some default behaviors. */
    
    const char *i = "latlong";
    const char *o = "latlong";
    const char *p = "rgss";
    const char *f = "linear";
    
    float rot[3] = { 0.f, 0.f, 0.f };
    
    int n = 1024;
    int c;
    
    /* Parse the command line options. */
    
    const char* convtype=NULL;
    double convparameter = 0.0;
    
    
    while ((c = getopt(argc, argv, "i:o:p:n:f:c:x:y:z:")) != -1)
        switch (c)
    {
        case 'i': i      = optarg;               break;
        case 'o': o      = optarg;               break;
        case 'p': p      = optarg;               break;
        case 'f': f      = optarg;               break;
        case 'x': rot[0] = strtod(optarg, 0);    break;
        case 'y': rot[1] = strtod(optarg, 0);    break;
        case 'z': rot[2] = strtod(optarg, 0);    break;
        case 'n': n      = strtol(optarg, 0, 0); break;
        case 'c':
        {
            convtype = optarg;
            
            if(!strcmp(convtype, "phong") || !strcmp(convtype, "gauss") || !strcmp(convtype, "hanning") || !strcmp(convtype, "lanczos"))
            {
                try {
                    convparameter = std::stod(argv[optind], 0);
                }
                catch (const std::exception& e)
                {
                    return usage(argv[0]);
                }

            }
            else
            {
                return usage(argv[0]);
            }
            
            break;
        }
            
        default: return usage(argv[0]);
    }
    
    int      num = 1;
    image   *src = 0;
    image   *dst = 0;
    image   *tmp = 0;
    to_img   img;
    to_env   env;
    filter   fil;
    
    /* Select the sampler. */
    
    if      (!strcmp(f, "linear"))  fil = filter_linear;
    else if (!strcmp(f, "nearest")) fil = filter_nearest;
    else return usage(argv[0]);
    
    /* Read the input image. */
    
    const int fileArgCt = 2;
    
    std::vector<const char*> inputFiles;
    inputFiles.insert(inputFiles.begin(), &argv[optind], &argv[argc-1]);
    
    if (optind + fileArgCt <= argc)
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
            dst = image_alloc((num = 6), n, n, src->c);
            env = cube_to_env;
        }
        else if (!strcmp(o, "dome"))
        {
            dst = image_alloc((num = 1), n, n, src->c);
            env = dome_to_env;
        }
        else if (!strcmp(o, "hemi"))
        {
            dst = image_alloc((num = 1), n, n, src->c);
            env = hemi_to_env;
        }
        else if (!strcmp(o, "ball"))
        {
            dst = image_alloc((num = 1), n, n, src->c);
            env = ball_to_env;
        }
        else if (!strcmp(o, "latlong"))
        {
            dst = image_alloc((num = 1), n, 2 * n, src->c);
            env = rect_to_env;
        }
        else if (!strcmp(o, "llsquare"))
        {
            dst = image_alloc((num = 1), n, n, src->c);
            env = rect_to_env;
        }
        else return usage(argv[0]);
    }
    else
    {
        throw std::runtime_error("Failed to load file");
    }
    
    /* Perform the remapping using the selected pattern. */
    
    if (src && dst)
    {
        if (!strcmp(p, "cent"))
            process(src, dst, &cent_pattern, rot, fil, img, env, num);
        
        else if (!strcmp(p, "rgss"))
            process(src, dst, &rgss_pattern, rot, fil, img, env, num);
        
        else if (!strcmp(p, "box2"))
            process(src, dst, &box2_pattern, rot, fil, img, env, num);
        
        else if (!strcmp(p, "box3"))
            process(src, dst, &box3_pattern, rot, fil, img, env, num);
        
        else if (!strcmp(p, "box4"))
            process(src, dst, &box4_pattern, rot, fil, img, env, num);
        
        else return usage(argv[0]);
    
        /* Write the output. */
        
        //image_writer(argv[optind + 1], dst, num);
        image_writer(argv[argc - 1], dst, num);
    }
    
}catch (const std::runtime_error& e)
{
    std::clog<<e.what()<<std::endl;
}
    return 0;
}