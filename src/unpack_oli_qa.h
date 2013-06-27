#ifndef _UNPACK_OLI_QA_H_
#define _UNPACK_OLI_QA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xtiffio.h"
#include "geotiffio.h"
#include "bool.h"
#include "error_handler.h"

#define STR_SIZE 1024

/* Set up the enumerated types for the quality bands */
typedef enum
{
    FILL = 0,
    DROPPED_FRAME,
    TERRAIN_OCCL,
    WATER,
    CLOUD_SHADOW,
    VEG,
    SNOW_ICE,
    CIRRUS,
    CLOUD,
    NQUALITY_TYPES
} Quality_t;

/* Set up the enumerated types for low, medium, or high confidence */
typedef enum
{
    UNDEFINED = 0,
    LOW,
    MED,
    HIGH,
    NCONF_TYPES
} Confidence_t;

/* Set up local defines for the UTM and PS projections */
#define UNDEFINED_PROJ -99
#define UTM_PROJ 1
#define PS_PROJ 2

/* Prototypes */
void usage ();

short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    bool *combine_bits,   /* O: should the QA bits be combined? */
    char **infile,        /* O: address of input filename */
    char **outfile,       /* O: address of output filename */
    bool qa_specd[NQUALITY_TYPES],  /* O: array to specify which of the QA
                                          bands was specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES]
                          /* O: array to specify the confidence level for
                                each of the quality fields */
);

TIFF *create_tiff
(
    char *tiffile,         /* I: geotiff filename */
    int proj,              /* I: projection type */
    uint32 nlines,         /* I: number of lines in tiff image */
    uint32 nsamps,         /* I: number of samples in tiff image */
    double *tie_point,     /* I: corner tie points for projection [3] is ULx
                                 and [4] is ULy (pass address of this array, as
                                 memory is allocated) */
    double *pixel_size,    /* I: pixel size array (x, y) (pass address of this
                                 array, as memory is allocated) */
    uint16 coord_sys,      /* I: coordinate system used (PixelIsArea or
                                 PixelIsPoint) */
    uint16 model_type,     /* I: geokey for the model type */
    uint16 linear_units,   /* I: geokey for the linear units */
    uint16 angular_units,  /* I: geokey for the angular units */
    uint16 projected_type, /* I: geokey for the angular units */
    uint16 proj_linear_units, /* O: geokey for proj linear units (PS proj) */
    double proj_parms[15],  /* O: projection parameters (PS proj) */
    char *citation         /* I: citation string */
);

short read_attributes
(
    char *infile,      /* I: input geotiff filename */
    int *proj,         /* O: projection type */
    uint32 *nlines,    /* O: number of lines in tiff image */
    uint32 *nsamps,    /* O: number of samples in tiff image */
    uint16 *bitspersample,  /* O: bits per sample in tiff image */
    uint16 *sampleformat,   /* O: data type of tiff image */
    double tie_point[6],    /* O: corner tie points for projection [3] is ULx
                                  and [4] is ULy */
    double pixel_size[3],   /* O: pixel size array (x, y, -) */
    uint16 *coord_sys,      /* O: coordinate system used (PixelIsArea or
                                  PixelIsPoint) */
    uint16 *model_type,     /* O: geokey for the model type */
    uint16 *linear_units,   /* O: geokey for the linear units */
    uint16 *angular_units,  /* O: geokey for the angular units */
    uint16 *projected_type, /* O: geokey for the angular units */
    uint16 *proj_linear_units, /* O: geokey for proj linear units (PS proj) */
    double proj_parms[15],  /* O: projection parameters (PS proj) */
    char *citation          /* O: citation string */
);

short unpack_bits
(
    char *qa_infile,      /* I: input QA filename */
    char *qa_outfile,     /* I: output QA base filename */
    bool qa_specd[NQUALITY_TYPES],  /* I: array to specify which QA bands
                                          was specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES]
                          /* I: array to specify the confidence level for
                                each of the quality fields */
);

short unpack_combine_bits
(
    char *qa_infile,      /* I: input QA filename */
    char *qa_outfile,     /* I: output QA filename */
    bool qa_specd[NQUALITY_TYPES],  /* I: array to specify which QA bands
                                          was specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES]
                          /* I: array to specify the confidence level for
                                each of the quality fields */
);

#endif
