#ifndef __envtools__inout__
#define __envtools__inout__

#include "../envremap.h"
#include <vector>

image *image_reader(const std::vector<const char *>& inputFiles, int n);
image *image_reader(const char *name, int n);

void image_writer(const char *outputFiles, image *out, int n);


#endif /* defined(__envtools__inout__) */
