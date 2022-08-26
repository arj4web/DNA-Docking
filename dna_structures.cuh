#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <cufft.h>
#include <cuda.h>

static const dim3 threadperblock2D(32,32);
static const dim3 threadperblock3D(10,10,10);

struct Atom
{
    int serial;
    char atom_name[5];
    float coord[4];
    float occupancy;
    float temp_factor;
    float charge;
    char molecule[2];
};

struct Nucleic_Acid
{
    char nucleicAcid_name[4];
    char chainID[2];
    char nuCode[6];
    int nc;
    int size;
    struct Atom *atom;
};

struct DNA_structure{
    char name[256];
    int length;
    struct Nucleic_Acid *nucleotide;
};

//Macros

#define gaddress(x,y,z,grid_size) ( (z) + ( 2 * ( (grid_size) / 2 + 1 ) ) * ( (y) + (grid_size) * (x) ) )

#define gcentre(ordinate,grid_span,grid_size) ( (float)(ordinate) + .5 ) * (float)( (grid_span) / (float)(grid_size) ) - (float)( (grid_span) / 2 )

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
extern __device__ __host__ int gord( float position , float grid_span , int grid_size ) ;
extern __device__ float pythagoras( float x1 , float y1 , float z1 , float x2 , float y2 , float z2 ) ;
