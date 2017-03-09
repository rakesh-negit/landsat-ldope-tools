#include "geotiffio.h"
#include "xtiffio.h"
#include "unpack_collection_qa.h"

/* Define the constants used for shifting bits and ANDing with the bits to
   get to the desire quality bits */
const int SINGLE_BIT = 0x01; /* 00000001 */
const int DOUBLE_BIT = 0x03; /* 00000011 */
const int SHIFT[NQUALITY_TYPES] = {0, 1, 2, 4, 5, 7, 9, 11};

/******************************************************************************
MODULE:  read_attributes

PURPOSE:  Read the file attributes and geokey information from the GeoTIFF file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
6/19/2013     Gail Schmidt     Original Development
7/12/2016     Gail Schmidt     Updated to read ProjStraightVertPoleLongGeoKey
                               if ProjStraightVertPoleLongGeoKey does not exist
                               for Polar Stereographic scenes

NOTES:
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "read_attributes"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variable */
    int cit_length;          /* length of citation key */
    int size;                /* size of individual element in citation var */
    int proj_type;           /* projection type geokey value */
    tagtype_t type;          /* tag type of the citation variable */
    double *tiePoint=NULL;   /* pointer for reading the tiepoints */
    double *pixelScale=NULL; /* pointer for reading the pixel size */
    uint16 count;            /* count of attributes to be read from tiff file */
    TIFF *fp_tiff=NULL;      /* tiff file pointer for input file */
    GTIF *fp_gtif=NULL;      /* geotiff key parser for input file */

    /* Open the input tiff file */
    if ((fp_tiff = XTIFFOpen (infile, "r")) == NULL)
    {
        sprintf (errmsg, "Error opening base TIFF file %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    } 
  
    /* Get metadata from tiff file */
    if (TIFFGetField (fp_tiff, TIFFTAG_IMAGELENGTH, nlines) == 0)
    {
        sprintf (errmsg, "Error reading number of lines from base TIFF file "
            "%s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (TIFFGetField (fp_tiff, TIFFTAG_IMAGEWIDTH, nsamps) == 0)
    {
        sprintf (errmsg, "Error reading number of samples from base TIFF file "
            "%s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (TIFFGetField (fp_tiff, TIFFTAG_BITSPERSAMPLE, bitspersample) == 0)
    {
        sprintf (errmsg, "Error reading bitspersample from base TIFF file "
            "%s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    if (TIFFGetField (fp_tiff, TIFFTAG_SAMPLEFORMAT, sampleformat) == 0)
    {
        sprintf (errmsg, "Error reading sampleformat from base TIFF file "
            "%s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    count = 6;
    if (TIFFGetField (fp_tiff, TIFFTAG_GEOTIEPOINTS, &count, &tiePoint) == 0)
    {
        sprintf (errmsg, "Error reading tiepoints from base TIFF file %s",
            infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    for (i = 0; i < count; i++)
        tie_point[i] = tiePoint[i];
  
    count = 3;
    if (TIFFGetField (fp_tiff, TIFFTAG_GEOPIXELSCALE, &count, &pixelScale) == 0)
    {
        sprintf (errmsg, "Error reading pixel size from base TIFF file %s",
            infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    for (i = 0; i < count; i++)
        pixel_size[i] = pixelScale[i];

    /* Open the key parser for the geotiff file */
    fp_gtif = GTIFNew (fp_tiff);

    /* Try to read the coordinate transform geokey.  If it is set and set to
       polar stereographic, then process the PS geokeys.  If it isn't set, then
       assume geokeys for UTM.  If it's set and not PS, then flag it as an
       error. */
    *proj = UNDEFINED_PROJ;
    if (!GTIFKeyGet (fp_gtif, ProjCoordTransGeoKey, &proj_type, 0, 1))
    {  /* assume UTM projection */
        *proj = UTM_PROJ;
    }
    else if (proj_type == CT_PolarStereographic)
    {  /* assume PS projection */
        *proj = PS_PROJ;
    }
    else
    {
        sprintf (errmsg, "Unsupported projection type in the GeoTIFF file.  "
            "If the ProjCoordTransGeoKey is not set, then UTM is assumed.  If "
            "this key is set, then it is expected to be "
            "CT_PolarStereographic.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the various GeoKeys */
    if (!GTIFKeyGet (fp_gtif, GTModelTypeGeoKey, model_type, 0, 1))
    {
        sprintf (errmsg, "Error reading the GTModelTypeGeoKey from the "
            "GeoTIFF file %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    cit_length = GTIFKeyInfo (fp_gtif, GTCitationGeoKey, &size, &type);
    if (GTIFKeyGet (fp_gtif, GTCitationGeoKey, citation, 0, cit_length) !=
        cit_length)
    {
        sprintf (errmsg, "Error reading the GTCitationGeoKey from the "
            "GeoTIFF file %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (!GTIFKeyGet (fp_gtif, GeogLinearUnitsGeoKey, linear_units, 0, 1))
    {
        sprintf (errmsg, "Error reading the GeogLinearUnitsGeoKey from the "
            "GeoTIFF file %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (!GTIFKeyGet (fp_gtif, GeogAngularUnitsGeoKey, angular_units, 0, 1))
    {
        sprintf (errmsg, "Error reading the GeogAngularUnitsGeoKey from the "
            "GeoTIFF file %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (!GTIFKeyGet (fp_gtif, ProjectedCSTypeGeoKey, projected_type, 0, 1))
    {
        sprintf (errmsg, "Error reading the ProjectedCSTypeGeoKey from the "
            "GeoTIFF file %s", infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* GTRasterTypeGeoKey dictates whether the reference coordinate is the UL
       (*RasterPixelIsArea*, code 1) or center (*RasterPixelIsPoint*, code 2)
       of the UL pixel. If this key is missing, the default (as defined by the
       specification) is to be *RasterPixelIsArea*, which is the UL of the UL
       pixel. */
    if (!GTIFKeyGet (fp_gtif, GTRasterTypeGeoKey, coord_sys, 0, 1))
    {
        /* use a flag to specify that it wasn't in the current file */
        *coord_sys = -99;
    }

    /* Read additional geokeys for polar stereographic projection */
    if (*proj == PS_PROJ)
    {
        if (!GTIFKeyGet (fp_gtif, ProjLinearUnitsGeoKey, proj_linear_units,
            0, 1))
        {
            sprintf (errmsg, "Error reading ProjLinearUnitsGeoKey from "
                "GeoTIFF file %s", infile);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (!GTIFKeyGet (fp_gtif, ProjNatOriginLongGeoKey, &proj_parms[4],
            0, 1))
        {
            if (!GTIFKeyGet (fp_gtif, ProjStraightVertPoleLongGeoKey,
                &proj_parms[4], 0, 1))
            {
                sprintf (errmsg, "Error reading ProjNatOriginLongGeoKey or "
                    "ProjStraightVertPoleLongGeoKey from GeoTIFF file %s",
                    infile);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }

        if (!GTIFKeyGet (fp_gtif, ProjNatOriginLatGeoKey, &proj_parms[5], 0, 1))
        {
            sprintf (errmsg, "Error reading ProjNatOrigintLatGeoKey from "
                "GeoTIFF file %s", infile);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (!GTIFKeyGet (fp_gtif, ProjFalseEastingGeoKey, &proj_parms[6], 0, 1))
        {
            sprintf (errmsg, "Error reading ProjFalseEastingGeoKey from "
                "GeoTIFF file %s", infile);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (!GTIFKeyGet (fp_gtif, ProjFalseNorthingGeoKey, &proj_parms[7],
            0, 1))
        {
            sprintf (errmsg, "Error reading ProjFalseNorthingGeoKey from "
                "GeoTIFF file %s", infile);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the geotiff key parser */
    GTIFFree (fp_gtif);

    /* Close the input tiff file */
    XTIFFClose (fp_tiff);

    return (SUCCESS);
}


/******************************************************************************
MODULE:  create_tiff

PURPOSE:  Create the tiff file and set the attributes.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
NULL            An error occurred during processing
TIFF *          Successful creation of the new GeoTIFF file

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
6/19/2013     Gail Schmidt     Original Development

NOTES:
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "create_tiff"; /* function name */
    char errmsg[STR_SIZE];   /* error message */

    uint16 count;            /* count of attributes written to tiff file */
    TIFF *fp_tiff=NULL;      /* tiff file pointer */
    GTIF *fp_gtif=NULL;      /* geotiff key parser */

    /* Create the tiff file */
    if ((fp_tiff = XTIFFOpen (tiffile, "w")) == NULL)
    {
        sprintf (errmsg, "Error creating base TIFF file %s", tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    } 
  
    /* Set metadata for tiff file */
    if (TIFFSetField (fp_tiff, TIFFTAG_IMAGELENGTH, nlines) == 0)
    {
        sprintf (errmsg, "Error setting number of lines to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (TIFFSetField (fp_tiff, TIFFTAG_IMAGEWIDTH, nsamps) == 0)
    {
        sprintf (errmsg, "Error setting number of samples to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (TIFFSetField (fp_tiff, TIFFTAG_BITSPERSAMPLE, 8) == 0)
    {
        sprintf (errmsg, "Error setting bitspersample to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (TIFFSetField (fp_tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT) == 0)
    {
        sprintf (errmsg, "Error setting sampleformat to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (TIFFSetField (fp_tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE) == 0)
    {
        sprintf (errmsg, "Error setting compression to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (TIFFSetField (fp_tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK)
        == 0)
    {
        sprintf (errmsg, "Error setting photometric to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    if (TIFFSetField (fp_tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG) == 0)
    {
        sprintf (errmsg, "Error setting planarconfig to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    count = 6;
    if (TIFFSetField (fp_tiff, TIFFTAG_GEOTIEPOINTS, count, tie_point) == 0)
    {
        sprintf (errmsg, "Error setting tiepoints to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }
  
    count = 3;
    if (TIFFSetField (fp_tiff, TIFFTAG_GEOPIXELSCALE, count, pixel_size) == 0)
    {
        sprintf (errmsg, "Error setting pixel size to base TIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Open the geotiff key parser */
    fp_gtif = GTIFNew (fp_tiff);

    /* Set the GTRasterTypeGeoKey, if it was set in the original product */
    if (coord_sys != -99)
    {
        if (!GTIFKeySet (fp_gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, coord_sys))
        {
            sprintf (errmsg, "Error setting GTRasterTypeGeokey to GeoTIFF file "
                "%s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
    }

    /* Set the rest of the geokeys */
    if (!GTIFKeySet (fp_gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, model_type))
    {
        sprintf (errmsg, "Error setting GTModelTypeGeoKey to GeoTIFF file "
            "%s", tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    if (!GTIFKeySet (fp_gtif, GTCitationGeoKey, TYPE_ASCII, 0, citation))
    {
        sprintf (errmsg, "Error setting GTCitationGeoKey to GeoTIFF file "
            "%s", tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    if (!GTIFKeySet (fp_gtif, GeogLinearUnitsGeoKey, TYPE_SHORT, 1,
        linear_units))
    {
        sprintf (errmsg, "Error setting GeogLinearUnitsGeoKey to GeoTIFF file "
            "%s", tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    if (!GTIFKeySet (fp_gtif, GeogAngularUnitsGeoKey, TYPE_SHORT, 1,
        angular_units))
    {
        sprintf (errmsg, "Error setting GeogAngularUnitsGeoKey to GeoTIFF file "
            "%s", tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    if (!GTIFKeySet (fp_gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1,
        projected_type))
    {
        sprintf (errmsg, "Error setting ProjectedCSTypeGeoKey to GeoTIFF file "
            "%s", tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Set additional geokeys for polar stereographic projection */
    if (proj == PS_PROJ)
    {
        if (!GTIFKeySet (fp_gtif, ProjCoordTransGeoKey, TYPE_SHORT, 1,
            CT_PolarStereographic))
        {
            sprintf (errmsg, "Error setting ProjCoordTransGeoKey for Polar "
                "Stereographic in the GeoTIFF file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        if (!GTIFKeySet (fp_gtif, GeographicTypeGeoKey, TYPE_SHORT, 1,
            GCS_WGS_84))
        {
            sprintf (errmsg, "Error setting GeographicTypeGeoKey for Polar "
                "Stereographic in the GeoTIFF file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        if (!GTIFKeySet (fp_gtif, ProjLinearUnitsGeoKey, TYPE_SHORT, 1,
            proj_linear_units))
        {
            sprintf (errmsg, "Error setting ProjLinearUnitsGeoKey to GeoTIFF "
                "file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        if (!GTIFKeySet (fp_gtif, ProjNatOriginLongGeoKey, TYPE_DOUBLE, 1,
            proj_parms[4]))
        {
            sprintf (errmsg, "Error setting ProjNatOrigintLongGeoKey to "
                "GeoTIFF file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        if (!GTIFKeySet (fp_gtif, ProjNatOriginLatGeoKey, TYPE_DOUBLE, 1,
            proj_parms[5]))
        {
            sprintf (errmsg, "Error setting ProjNatOrigintLatGeoKey to "
                "GeoTIFF file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        if (!GTIFKeySet (fp_gtif, ProjFalseEastingGeoKey, TYPE_DOUBLE, 1,
            proj_parms[6]))
        {
            sprintf (errmsg, "Error setting ProjFalseEastingGeoKey to "
                "GeoTIFF file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }

        if (!GTIFKeySet (fp_gtif, ProjFalseNorthingGeoKey, TYPE_DOUBLE, 1,
            proj_parms[7]))
        {
            sprintf (errmsg, "Error setting ProjFalseNorthingGeoKey to "
                "GeoTIFF file %s", tiffile);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
    }

    /* Write the geotiff keys */
    if (!GTIFWriteKeys (fp_gtif))
    {
        sprintf (errmsg, "Error writing the geokeys to the GeoTIFF file %s",
            tiffile);
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Close the geotiff key parser */
    GTIFFree (fp_gtif);

    /* Pass back the input tiff file pointer */
    return (fp_tiff);
}


/******************************************************************************
MODULE:  unpack_bits

PURPOSE:  Unpack the QA band for the specified quality bits.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
8/31/2016     Ray Dittmeier    Original Development

NOTES:
******************************************************************************/
short unpack_bits
(
    char *qa_infile,      /* I: input QA filename */
    char *qa_outfile,     /* I: output QA base filename */
    bool qa_specd[NQUALITY_TYPES],  /* I: array to specify which QA bands
                                          were specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES],
                          /* I: array to specify the confidence level for
                                each of the quality fields */
    int satellite_number  /* I: number of the satellite, e.g.: 8 */
)
{
    char FUNC_NAME[] = "unpack_bits"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char tmpstr[STR_SIZE];   /* temporary pointer string message */
    char citation[STR_SIZE]; /* geokey for citation string */
    char outfile[NQUALITY_TYPES][STR_SIZE];   /* array of output filenames */
    uint32 nlines, nsamps;   /* number of lines and samples */
    int i;                   /* looping variable */
    int line, samp;          /* current line and sample to be processed */
    uint16 bitspersample;    /* bits per sample in input tiff image */
    uint16 sampleformat;     /* data type of input tiff image */
    int proj_type;           /* projection type */
    uint8 *unpack_buf=NULL;  /* unpacked QA band from the file */
    uint16 coord_sys;        /* geokey for coordinate system */
    uint16 model_type;       /* geokey for the model type */
    uint16 linear_units;     /* geokey for the linear units */
    uint16 angular_units;    /* geokey for the angular units */
    uint16 projected_type;   /* geokey for the angular units */
    uint16 unpack_val;       /* unpacked bit value for current pixel */
    uint16 *qa_buf=NULL;     /* QA band from the file */
    uint16 proj_linear_units; /* geokey for the proj linear units (PS proj) */
    double proj_parms[15];   /* projection parameters (PS proj) */
    double tie_points[6];    /* corner point information */
    double pixel_size[3];    /* pixel size (x, y, -) */
    TIFF *in_fp_tiff=NULL;   /* tiff file pointer for input file */
    TIFF *out_fp_tiff[NQUALITY_TYPES];  /* array of tiff file pointers for each
                                           output file */

    /* Init the output file pointer */
    for (i = 0; i < NQUALITY_TYPES; i++)
        out_fp_tiff[i] = NULL;

    /* Read the file attributes from the input file */
    if (read_attributes (qa_infile, &proj_type, &nlines, &nsamps,
        &bitspersample, &sampleformat, tie_points, pixel_size, &coord_sys,
        &model_type, &linear_units, &angular_units, &projected_type,
        &proj_linear_units, proj_parms, citation) != SUCCESS)
    {
        sprintf (errmsg, "Error reading attributes from geoTIFF file %s",
            qa_infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Check to make sure the product is a 16-bit unsigned integer */
    if (bitspersample != 16)
    {
        sprintf (errmsg, "Input GeoTIFF QA band is expected to be a 16-bit "
            "integer but instead it is a %d-bit product", bitspersample);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (sampleformat != SAMPLEFORMAT_UINT)
    {
        if (sampleformat == SAMPLEFORMAT_INT)
            sprintf (tmpstr, "signed integer");
        else if (sampleformat == SAMPLEFORMAT_IEEEFP)
            sprintf (tmpstr, "float");
        else
            sprintf (tmpstr, "unknown");
        sprintf (errmsg, "Error: input GeoTIFF QA band is expected to be "
            "an unsigned integer but instead it is a %s product", tmpstr);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate memory for the QA band and an unpacked band (one scanline) */
    qa_buf = (uint16 *) calloc (nsamps, sizeof (uint16));
    if (qa_buf == NULL)
    {
        sprintf (errmsg, "Error allocating memory (full scene) for the input "
            "QA band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    unpack_buf = (uint8 *) calloc (nsamps, sizeof (uint8));
    if (unpack_buf == NULL)
    {
        sprintf (errmsg, "Error allocating memory (full scene) for the "
            "unpacked QA band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the input tiff file */
    if ((in_fp_tiff = XTIFFOpen (qa_infile, "r")) == NULL)
    {
        sprintf (errmsg, "Error opening base TIFF file %s", qa_infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    } 
  
    /* Create and open the output tiff files, depending on which QA bits were
       specified to be unpacked */
    if (qa_specd[FILL])
    {
        sprintf (outfile[FILL], "%s_fill.tif", qa_outfile);
        out_fp_tiff[FILL] = create_tiff (outfile[FILL], proj_type, nlines,
            nsamps, tie_points, pixel_size, coord_sys, model_type, linear_units,
            angular_units, projected_type, proj_linear_units, proj_parms,
            citation);
        if (!out_fp_tiff[FILL])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s", outfile[FILL]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* For L4-7 this is dropped pixel.  For L8 it is terrain occlusion. */
    if (qa_specd[OCCLUSION_OR_DROPPED])
    {
        if (satellite_number == 8)
        {
            sprintf (outfile[OCCLUSION_OR_DROPPED], "%s_terrain_occl.tif", 
                qa_outfile);
        }
        else
        {
            sprintf (outfile[OCCLUSION_OR_DROPPED], "%s_dropped_pixel.tif", 
                qa_outfile);
        }
        out_fp_tiff[OCCLUSION_OR_DROPPED] = create_tiff (
            outfile[OCCLUSION_OR_DROPPED], proj_type, nlines, nsamps, 
            tie_points, pixel_size, coord_sys, model_type, linear_units, 
            angular_units, projected_type, proj_linear_units, proj_parms, 
            citation);
        if (!out_fp_tiff[OCCLUSION_OR_DROPPED])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s",
                outfile[OCCLUSION_OR_DROPPED]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (qa_specd[RADIOMETRIC_SAT])
    {
        sprintf (outfile[RADIOMETRIC_SAT], "%s_radiometric_sat.tif", 
            qa_outfile);
        out_fp_tiff[RADIOMETRIC_SAT] = create_tiff (outfile[RADIOMETRIC_SAT],
            proj_type, nlines, nsamps, tie_points, pixel_size, coord_sys,
            model_type, linear_units, angular_units, projected_type,
            proj_linear_units, proj_parms, citation);
        if (!out_fp_tiff[RADIOMETRIC_SAT])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s",
                outfile[RADIOMETRIC_SAT]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (qa_specd[CLOUD])
    {
        sprintf (outfile[CLOUD], "%s_cloud.tif", qa_outfile);
        out_fp_tiff[CLOUD] = create_tiff (outfile[CLOUD], proj_type, nlines,
            nsamps, tie_points, pixel_size, coord_sys, model_type, linear_units,
            angular_units, projected_type, proj_linear_units, proj_parms,
            citation);
        if (!out_fp_tiff[CLOUD])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s", outfile[CLOUD]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (qa_specd[CLOUD_CONFIDENCE])
    {
        sprintf (outfile[CLOUD_CONFIDENCE], "%s_cloud_confidence.tif", 
            qa_outfile);
        out_fp_tiff[CLOUD_CONFIDENCE] = create_tiff (outfile[CLOUD_CONFIDENCE], 
            proj_type, nlines, nsamps, tie_points, pixel_size, coord_sys, 
            model_type, linear_units, angular_units, projected_type, 
            proj_linear_units, proj_parms, citation);
        if (!out_fp_tiff[CLOUD_CONFIDENCE])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s", 
                outfile[CLOUD_CONFIDENCE]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (qa_specd[CLOUD_SHADOW])
    {
        sprintf (outfile[CLOUD_SHADOW], "%s_cloud_shadow.tif", qa_outfile);
        out_fp_tiff[CLOUD_SHADOW] = create_tiff (outfile[CLOUD_SHADOW],
            proj_type, nlines, nsamps, tie_points, pixel_size, coord_sys,
            model_type, linear_units, angular_units, projected_type,
            proj_linear_units, proj_parms, citation);
        if (!out_fp_tiff[CLOUD_SHADOW])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s",
                outfile[CLOUD_SHADOW]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (qa_specd[SNOW_ICE])
    {
        sprintf (outfile[SNOW_ICE], "%s_snow_ice.tif", qa_outfile);
        out_fp_tiff[SNOW_ICE] = create_tiff (outfile[SNOW_ICE], proj_type,
            nlines, nsamps, tie_points, pixel_size, coord_sys, model_type,
            linear_units, angular_units, projected_type, proj_linear_units,
            proj_parms, citation);
        if (!out_fp_tiff[SNOW_ICE])
        {
            sprintf (errmsg, "Error creating geoTIFF file %s",
                outfile[SNOW_ICE]);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    if (satellite_number == 8)
    {
        if (qa_specd[CIRRUS])
        {
            sprintf (outfile[CIRRUS], "%s_cirrus.tif", qa_outfile);
            out_fp_tiff[CIRRUS] = create_tiff (outfile[CIRRUS], proj_type, 
                nlines, nsamps, tie_points, pixel_size, coord_sys, model_type, 
                linear_units, angular_units, projected_type, proj_linear_units,
                proj_parms, citation);
            if (!out_fp_tiff[CIRRUS])
            {
                sprintf (errmsg, "Error creating geoTIFF file %s", 
                    outfile[CIRRUS]);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }
    }

    /* Loop through the QA band and process one scanline at a time */
    for (line = 0; line < nlines; line++)
    {
        if (TIFFReadScanline (in_fp_tiff, qa_buf, line, 0) == -1)
        {
            sprintf (errmsg, "Error reading line %d from the input file", line);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        } 

        /* Unpack each line based on the user input for which quality bits
           should be output */
        if (qa_specd[FILL])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[FILL];
                unpack_buf[samp] = (uint8) (unpack_val & SINGLE_BIT);
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[FILL], unpack_buf, line, 0)
                == -1)
            {
                sprintf (errmsg, "Error writing line %d to the fill file",
                    line);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        /* For L4-7 this is dropped pixel.  For L8 it is terrain occlusion.
           They occupy the same single bit. */
        if (qa_specd[OCCLUSION_OR_DROPPED])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[OCCLUSION_OR_DROPPED];
                unpack_buf[samp] = (uint8) (unpack_val & SINGLE_BIT);
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[OCCLUSION_OR_DROPPED], 
                unpack_buf, line, 0) == -1)
            {
                if (satellite_number == 8)
                {
                    sprintf (errmsg, "Error writing line %d to the terrain "
                        "occlusion file", line);
                }
                else
                {
                    sprintf (errmsg, "Error writing line %d to the dropped "
                        "pixel file", line);
                }
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        if (qa_specd[RADIOMETRIC_SAT])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[RADIOMETRIC_SAT];
                unpack_buf[samp] = (uint8) (unpack_val & DOUBLE_BIT);
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[RADIOMETRIC_SAT], unpack_buf, 
                line, 0) == -1)
            {
                sprintf (errmsg, "Error writing line %d to the radiometric "
                    "saturation file", line);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        if (qa_specd[CLOUD])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[CLOUD];
                unpack_buf[samp] = (uint8) (unpack_val & SINGLE_BIT);
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[CLOUD], unpack_buf, line,
                0) == -1)
            {
                sprintf (errmsg, "Error writing line %d to the cloud file",
                    line);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        if (qa_specd[CLOUD_CONFIDENCE])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[CLOUD_CONFIDENCE];
                if ((uint8) (unpack_val & DOUBLE_BIT) >= 
                    qa_conf[CLOUD_CONFIDENCE])
                    unpack_buf[samp] = 1;
                else
                    unpack_buf[samp] = 0;
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[CLOUD_CONFIDENCE], unpack_buf, 
                line, 0) == -1)
            {
                sprintf (errmsg, "Error writing line %d to the cloud "
                    "confidence file", line);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        if (qa_specd[CLOUD_SHADOW])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[CLOUD_SHADOW];
                if ((uint8) (unpack_val & DOUBLE_BIT) >= qa_conf[CLOUD_SHADOW])
                    unpack_buf[samp] = 1;
                else
                    unpack_buf[samp] = 0;
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[CLOUD_SHADOW], unpack_buf, line,
                0) == -1)
            {
                sprintf (errmsg, "Error writing line %d to the cloud shadow "
                    "file", line);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        if (qa_specd[SNOW_ICE])
        {
            /* Loop through each sample and unpack the pixel */
            for (samp = 0; samp < nsamps; samp++)
            {
                unpack_val = qa_buf[samp] >> SHIFT[SNOW_ICE];
                if ((uint8) (unpack_val & DOUBLE_BIT) >= qa_conf[SNOW_ICE])
                    unpack_buf[samp] = 1;
                else
                    unpack_buf[samp] = 0;
            }

            /* Write line to the output file */
            if (TIFFWriteScanline (out_fp_tiff[SNOW_ICE], unpack_buf, line, 0)
                == -1)
            {
                sprintf (errmsg, "Error writing line %d to the snow/ice file",
                    line);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            } 
        }

        if (satellite_number == 8)
        {
            if (qa_specd[CIRRUS])
            {
                /* Loop through each sample and unpack the pixel */
                for (samp = 0; samp < nsamps; samp++)
                {
                    unpack_val = qa_buf[samp] >> SHIFT[CIRRUS];
                    if ((uint8) (unpack_val & DOUBLE_BIT) >= qa_conf[CIRRUS])
                        unpack_buf[samp] = 1;
                    else
                        unpack_buf[samp] = 0;
                }

                /* Write line to the output file */
                if (TIFFWriteScanline (out_fp_tiff[CIRRUS], unpack_buf, line, 0)
                    == -1)
                {
                    sprintf (errmsg, "Error writing line %d to the cirrus file",
                        line);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                } 
            }
        }
    }

    /* Close the input and output tiff files */
    XTIFFClose (in_fp_tiff);
    if (qa_specd[FILL])
        XTIFFClose (out_fp_tiff[FILL]);
    if (qa_specd[OCCLUSION_OR_DROPPED])
        XTIFFClose (out_fp_tiff[OCCLUSION_OR_DROPPED]);
    if (qa_specd[RADIOMETRIC_SAT])
        XTIFFClose (out_fp_tiff[RADIOMETRIC_SAT]);
    if (qa_specd[CLOUD])
        XTIFFClose (out_fp_tiff[CLOUD]);
    if (qa_specd[CLOUD_CONFIDENCE])
        XTIFFClose (out_fp_tiff[CLOUD_CONFIDENCE]);
    if (qa_specd[CLOUD_SHADOW])
        XTIFFClose (out_fp_tiff[CLOUD_SHADOW]);
    if (qa_specd[SNOW_ICE])
        XTIFFClose (out_fp_tiff[SNOW_ICE]);
    if (satellite_number == 8)
    {
        if (qa_specd[CIRRUS])
            XTIFFClose (out_fp_tiff[CIRRUS]);
    }

    /* Free the buffer pointers */
    if (qa_buf != NULL)
        free (qa_buf);
    if (unpack_buf != NULL)
        free (unpack_buf);

    return (SUCCESS);
}


/******************************************************************************
MODULE:  unpack_combine_bits

PURPOSE:  Unpack the OLI QA band for the specified quality bits and combine
them into one output mask.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
8/31/2016     Ray Dittmeier    Original Development

NOTES:
******************************************************************************/
short unpack_combine_bits
(
    char *qa_infile,      /* I: input QA filename */
    char *qa_outfile,     /* I: output QA filename */
    bool qa_specd[NQUALITY_TYPES],  /* I: array to specify which QA bands
                                          were specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES],
                          /* I: array to specify the confidence level for
                                each of the quality fields */
    int satellite_number  /* I: number of the satellite, e.g.: 8 */
)
{
    char FUNC_NAME[] = "unpack_combine_bits"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char tmpstr[STR_SIZE];   /* temporary pointer string message */
    char citation[STR_SIZE]; /* geokey for citation string */
    uint32 nlines, nsamps;   /* number of lines and samples */
    int line, samp;          /* current line and sample to be processed */
    uint16 bitspersample;    /* bits per sample in input tiff image */
    uint16 sampleformat;     /* data type of input tiff image */
    int proj_type;           /* projection type */
    uint8 *unpack_buf=NULL;  /* unpacked QA band from the OLI file */
    uint16 coord_sys;        /* geokey for coordinate system */
    uint16 model_type;       /* geokey for the model type */
    uint16 linear_units;     /* geokey for the linear units */
    uint16 angular_units;    /* geokey for the angular units */
    uint16 projected_type;   /* geokey for the angular units */
    uint16 unpack_val;       /* unpacked bit value for current pixel */
    uint16 *qa_buf=NULL;     /* QA band from the OLI file */
    uint16 proj_linear_units; /* geokey for the proj linear units (PS proj) */
    double proj_parms[15];   /* projection parameters (PS proj) */
    double tie_points[6];    /* corner point information */
    double pixel_size[3];    /* pixel size (x, y, -) */
    TIFF *in_fp_tiff=NULL;   /* tiff file pointer for input file */
    TIFF *out_fp_tiff=NULL;  /* tiff file pointer for output file */

    /* Read the file attributes from the input file */
    if (read_attributes (qa_infile, &proj_type, &nlines, &nsamps,
        &bitspersample, &sampleformat, tie_points, pixel_size, &coord_sys,
        &model_type, &linear_units, &angular_units, &projected_type,
        &proj_linear_units, proj_parms, citation) != SUCCESS)
    {
        sprintf (errmsg, "Error reading attributes from geoTIFF file %s",
            qa_infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Check to make sure the product is a 16-bit unsigned integer */
    if (bitspersample != 16)
    {
        sprintf (errmsg, "Input GeoTIFF QA band is expected to be a 16-bit "
            "integer but instead it is a %d-bit product", bitspersample);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (sampleformat != SAMPLEFORMAT_UINT)
    {
        if (sampleformat == SAMPLEFORMAT_INT)
            sprintf (tmpstr, "signed integer");
        else if (sampleformat == SAMPLEFORMAT_IEEEFP)
            sprintf (tmpstr, "float");
        else
            sprintf (tmpstr, "unknown");
        sprintf (errmsg, "Error: input GeoTIFF QA band is expected to be "
            "an unsigned integer but instead it is a %s product", tmpstr);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate memory for the QA band and an unpacked band (one scanline) */
    qa_buf = (uint16 *) calloc (nsamps, sizeof (uint16));
    if (qa_buf == NULL)
    {
        sprintf (errmsg, "Error allocating memory (full scene) for the input "
            "QA band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    unpack_buf = (uint8 *) calloc (nsamps, sizeof (uint8));
    if (unpack_buf == NULL)
    {
        sprintf (errmsg, "Error allocating memory (full scene) for the "
            "unpacked QA band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the input tiff file */
    if ((in_fp_tiff = XTIFFOpen (qa_infile, "r")) == NULL)
    {
        sprintf (errmsg, "Error opening base TIFF file %s", qa_infile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    } 
  
    /* Create and open the output tiff file */
    out_fp_tiff = create_tiff (qa_outfile, proj_type, nlines, nsamps,
        tie_points, pixel_size, coord_sys, model_type, linear_units,
        angular_units, projected_type, proj_linear_units, proj_parms,
        citation);
    if (!out_fp_tiff)
    {
        sprintf (errmsg, "Error creating geoTIFF file %s", qa_outfile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Loop through the QA band and process one scanline at a time */
    for (line = 0; line < nlines; line++)
    {
        if (TIFFReadScanline (in_fp_tiff, qa_buf, line, 0) == -1)
        {
            sprintf (errmsg, "Error reading line %d from the input file", line);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        } 

        /* Unpack each pixel based on the user input for which quality bits
           should be output.  Once a bit is turned on, then move on to the
           next pixel. */
        for (samp = 0; samp < nsamps; samp++)
        {
            unpack_buf[samp] = 0;
            if (qa_specd[FILL])
            {
                unpack_val = qa_buf[samp] >> SHIFT[FILL];
                if ((uint8) (unpack_val & SINGLE_BIT) == 1)
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            /* For L4-7 this is dropped pixel.  For L8 it's terrain occlusion.
               They occupy the same single bit. */
            if (qa_specd[OCCLUSION_OR_DROPPED])
            {
                unpack_val = qa_buf[samp] >> SHIFT[OCCLUSION_OR_DROPPED];
                if ((uint8) (unpack_val & SINGLE_BIT) == 1)
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            if (qa_specd[RADIOMETRIC_SAT])
            {
                /* Radiometric saturation is a 2-bit input so any input value 
                   over 0 will result in a 1 output. It's handled like the
                   confidence values, but without a threshold. */ 
                unpack_val = qa_buf[samp] >> SHIFT[RADIOMETRIC_SAT];
                if ((uint8) (unpack_val & DOUBLE_BIT) >= 1)
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            if (qa_specd[CLOUD])
            {
                unpack_val = qa_buf[samp] >> SHIFT[CLOUD];
                if ((uint8) (unpack_val & SINGLE_BIT) == 1)
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            if (qa_specd[CLOUD_CONFIDENCE])
            {
                unpack_val = qa_buf[samp] >> SHIFT[CLOUD_CONFIDENCE];
                if ((uint8) (unpack_val & DOUBLE_BIT) 
                    >= qa_conf[CLOUD_CONFIDENCE])
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            if (qa_specd[CLOUD_SHADOW])
            {
                unpack_val = qa_buf[samp] >> SHIFT[CLOUD_SHADOW];
                if ((uint8) (unpack_val & DOUBLE_BIT) >= qa_conf[CLOUD_SHADOW])
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            if (qa_specd[SNOW_ICE])
            {
                unpack_val = qa_buf[samp] >> SHIFT[SNOW_ICE];
                if ((uint8) (unpack_val & DOUBLE_BIT) >= qa_conf[SNOW_ICE])
                {
                    unpack_buf[samp] = 1;
                    continue;
                }
            }

            if (satellite_number == 8)
            {
                if (qa_specd[CIRRUS])
                {
                    unpack_val = qa_buf[samp] >> SHIFT[CIRRUS];
                    if ((uint8) (unpack_val & DOUBLE_BIT) >= qa_conf[CIRRUS])
                    {
                        unpack_buf[samp] = 1;
                        continue;
                    }
                }
            }
        }  /* end for samp */

        /* Write combined line to the output file */
        if (TIFFWriteScanline (out_fp_tiff, unpack_buf, line, 0) == -1)
        {
            sprintf (errmsg, "Error writing line %d to the cloud file", line);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        } 
    }  /* end for line */

    /* Close the input and output tiff files */
    XTIFFClose (in_fp_tiff);
    XTIFFClose (out_fp_tiff);

    /* Free the buffer pointers */
    if (qa_buf != NULL)
        free (qa_buf);
    if (unpack_buf != NULL)
        free (unpack_buf);

    return (SUCCESS);
}
