#ifndef envtools_envremap_h
#define envtools_envremap_h

/* In image structure represents an input or output raster.                   */

struct image
{
    float *p;  // data
    int    h;  // height
    int    w;  // width
    int    c;  // sample count
    /// int    b;  // sample depth
    /// int    s;  // sample format
};

typedef struct image image;


/* A pattern structure represents a supersampling pattern.                    */

struct point
{
    float i;
    float j;
};

typedef struct point point;

struct pattern
{
    int    n;
    point *p;
};

typedef struct pattern pattern;

/*----------------------------------------------------------------------------*/

typedef void (*filter)(const image *, float, float, float *);

typedef int (*to_img)(int *, float *, float *, int, int, const float *);
typedef int (*to_env)(int,   float,   float,   int, int,       float *);

/// Sample an image at row i column j using linear interpolation
void filter_linear(const image *img, float i, float j, float *p);

/// Sample an image at row i column j using nearest neighbor
void filter_nearest(const image *img, float i, float j, float *p);

int ball_to_env(int f, float i, float j, int h, int w, float *v);

int hemi_to_env(int f, float i, float j, int h, int w, float *v);

int dome_to_env(int f, float i, float j, int h, int w, float *v);

int cube_to_env(int f, float i, float j, int h, int w, float *v);

int rect_to_env(int f, float i, float j, int h, int w, float *v);

int ball_to_img(int *f, float *i, float *j, int h, int w, const float *v);

int rect_to_img(int *f, float *i, float *j, int h, int w, const float *v);

int hemi_to_img(int *f, float *i, float *j, int h, int w, const float *v);

int cube_to_img(int *f, float *i, float *j, int h, int w, const float *v);

int dome_to_img(int *f, float *i, float *j, int h, int w, const float *v);

int dome_to_img(int *f, float *i, float *j, int h, int w, const float *v);

image *image_alloc(int n, int h, int w, int c);

void process(const image   *src,
             const image   *dst,
             const pattern *pat,
             const float   *rot,
             filter fil, to_img img, to_env env, int n);

/// Add borders to a cubemap image. Assume the given image pointer is an array
/// of six images. Copy each to a new sef of six images, each two pixels wider
/// and higher. Also copy the borders. This is necessary for correct cubemap sampling.
image* image_border(image *src);

extern point cent_points[] ;
extern point rgss_points[];
extern point box2_points[];
extern point box3_points[];
extern point box4_points[];

extern const pattern cent_pattern;
extern const pattern rgss_pattern;
extern const pattern box2_pattern;
extern const pattern box3_pattern;
extern const pattern box4_pattern;

#endif
