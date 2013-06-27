#include "unpack_oli_qa.h"

static const char conf_vals[NCONF_TYPES] = {'u', 'l', 'm', 'h'};

/******************************************************************************
MODULE:  unpack_oli_qa

PURPOSE:  Unpack the specified OLI QA band and write the quality bands to
individual GeoTIFF files using the specified base output filename.

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
6/21/2013     Gail Schmidt     Modified to allow the user to specify the
                               confidence level for processing the 2-bit
                               quality fields

NOTES:
******************************************************************************/
int main (int argc, char *argv[])
{
    bool combine_bits;       /* should the QA bits be combined? */
    bool qa_specd[NQUALITY_TYPES];  /* array to specify which of the QA bands
                                       was specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES]; /* array to specify the confidence
                                       level for each of the quality fields */
    char *qa_infile=NULL;    /* input QA filename */
    char *qa_outfile=NULL;   /* output QA filename or basename */
    char tmp_char;           /* temporary character for each QA band */
    int retval;              /* return status */

    printf ("Unpack of OLI QA band started ...\n");

    /* Read the command-line arguments to determine which file needs to be
       processed and which quality bands will be dumped */
    retval = get_args (argc, argv, &combine_bits, &qa_infile, &qa_outfile,
        qa_specd, qa_conf);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Tell the user which bands will be unpacked */
    printf ("OLI QA input file: %s\n", qa_infile);
    if (!combine_bits)
        printf ("Unpacked QA output file basename: %s\n", qa_outfile);
    else
        printf ("Unpacked and combined QA output filename: %s\n", qa_outfile);
    printf ("Process    Description\n"
            "-------    -----------\n");
    if (qa_specd[FILL])
        tmp_char = 'x';
    else
        tmp_char = ' ';
    printf ("   %c       Fill\n", tmp_char);

    if (qa_specd[DROPPED_FRAME])
        tmp_char = 'x';
    else
        tmp_char = ' ';
    printf ("   %c       Dropped frame\n", tmp_char);

    if (qa_specd[TERRAIN_OCCL])
        tmp_char = 'x';
    else
        tmp_char = ' ';
    printf ("   %c       Terrain occlusion\n", tmp_char);

    if (qa_specd[WATER])
        tmp_char = conf_vals[qa_conf[WATER]];
    else
        tmp_char = ' ';
    printf ("   %c       Water confidence\n", tmp_char);

    if (qa_specd[CLOUD_SHADOW])
        tmp_char = conf_vals[qa_conf[CLOUD_SHADOW]];
    else
        tmp_char = ' ';
    printf ("   %c       Cloud shadow\n", tmp_char);

    if (qa_specd[VEG])
        tmp_char = conf_vals[qa_conf[VEG]];
    else
        tmp_char = ' ';
    printf ("   %c       Vegetation confidence\n", tmp_char);

    if (qa_specd[SNOW_ICE])
        tmp_char = conf_vals[qa_conf[SNOW_ICE]];
    else
        tmp_char = ' ';
    printf ("   %c       Snow/ice confidence\n", tmp_char);

    if (qa_specd[CIRRUS])
        tmp_char = conf_vals[qa_conf[CIRRUS]];
    else
        tmp_char = ' ';
    printf ("   %c       Cirrus confidence\n", tmp_char);

    if (qa_specd[CLOUD])
        tmp_char = conf_vals[qa_conf[CLOUD]];
    else
        tmp_char = ' ';
    printf ("   %c       Cloud confidence\n\n", tmp_char);

    /* Read the input QA band, unpack the bits, combine if specified, and
       write out the desired band(s) */
    if (!combine_bits)
    {
        /* Unpack the bits into individual bands */
        retval = unpack_bits (qa_infile, qa_outfile, qa_specd, qa_conf);
        if (retval != SUCCESS)
        {   /* unpack_bits already printed the error message */
            exit (ERROR);
        }
    }
    else
    {
        /* Unpack the bits and combine into one band */
        retval = unpack_combine_bits (qa_infile, qa_outfile, qa_specd, qa_conf);
        if (retval != SUCCESS)
        {   /* unpack_bits already printed the error message */
            exit (ERROR);
        }
    }

    /* Free the filename pointers */
    if (qa_infile != NULL)
        free (qa_infile);
    if (qa_outfile != NULL)
        free (qa_outfile);

    /* Indicate successful completion of processing */
    printf ("Unpack of OLI QA band complete!\n");
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
6/19/2013   Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("unpack_oli_qa will read the QA band from the Landsat L8 OLI "
            "input file, then unpack this band into individual QA bands "
            "stored in multiple files using the user-specified base output "
            "filename.  The output bands will refer to the QA bits (from "
            "right to left), representing the QA information which is stored "
            "in the QA band.  In some cases a single bit is used to represent "
            "quality data in the OLI QA band and in other cases two bits are "
            "used for the quality info.\n\n"
            "For quality data represented by a single bit, the output values "
            "are as follows:\n"
            "  0 = No, this condition does not exist\n"
            "  1 = Yes, this condition exists\n\n"
            "For quality data represented by two bits, the user will be "
            "allowed to specify the confidence levels included in the mask.  "
            "The current confidence levels in the OLI QA band are as follows:\n"
            " 00 = Algorithm did not determine the status of this condition\n"
            " 01 = Algorithm has low confidence that this condition exists "
            "(0-33 percent confidence)\n"
            " 10 = Algorithm has medium confidence that this condition exists "
            "(34-66 percent confidence)\n"
            " 11 = Algorithm has high confidence that this condition exists "
            "(67-100 percent confidence)\n"
            "If the user specifies a confidence level of 'low' for the 2-bit "
            "confidence field, then the output mask will be flagged as 1 / yes "
            "if the 2-bit confidence value is low, medium, or high.  If the "
            "user specifies a confidence level of 'med' for the confidence "
            "field, then the output mask will be flagged if the confidence "
            "value is medium or high.  And, if the user specifies a confidence "
            "level of 'high', then the output mask will be flagged if the "
            "confidence value is high.\n\n"
            "The following table identifies the output quality band and how it "
            "correlates to the bits in the individual QA bands, when not using "
            "combine bits.  The user may elect to combine the specified QA "
            "bits into one single output file.  In that case, if any of the "
            "QA bits are turned on or meet the specified confidence level, "
            "then the output mask for that pixel will be flagged as "
            "1 / yes.\n\n"
            "QA Band               QA Bit(s)    Description\n"
            "------------------    ---------    -----------\n"
            "_fill.tif                0         Fill\n"
            "_dropped_frame.tif       1         Dropped frame\n"
            "_terrain_occl.tif        2         Terrain occlusion\n"
            "   N/A                   3         Reserved\n"
            "_water.tif               4,5       Water confidence\n"
            "_cloud_shadow.tif        6,7       Cloud shadow\n"
            "_vegetation.tif          8,9       Vegetation confidence\n"
            "_snow_ice.tif            10,11     Snow/ice confidence\n"
            "_cirrus.tif              12,13     Cirrus confidence\n"
            "_cloud.tif               14,15     Cloud confidence\n\n"
            "For more information about the Landsat OLI QA Band file, please "
            "refer to http://landsat.usgs.gov/L8QualityAssessmentBand.php\n\n");
    printf ("unpack_oli_qa --help will print the usage information\n\n");
    printf ("usage: unpack_oli_qa "
            "--ifile=input_QA_filename "
            "--ofile=output_unpacked_QA_filename "
            "[--all=conf_level][--fill=conf_level][--drop_frame=conf_level] "
            "[--terrain_occl=conf_level][--water=conf_level] "
            "[--cloud_shadow=conf_level][--veg=conf_level] "
            "[--snow_ice=conf_level][--cirrus=conf_level] "
            "[--cloud=conf_level] [--combine]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -ifile: name of the input QA file (GeoTIFF product with "
            "uint16 bands)\n");
    printf ("    -ofile: basename of the output unpacked QA file if not "
            "combining the QA bits, otherwise the full filename of the output "
            "file if combining the QA bits (GeoTIFF products with uint8 bands "
            "to match the user-specified quality bits)\n");
    printf ("\nwhere the following is optional:\n");
    printf ("    -combine: indicates the specified QA bits will be combined "
            "into one single output band (default is false)\n");
    printf ("\nwhere the following QA bit parameters are optional:\n");
    printf ("    -all: specifies all the quality bits should be output "
            "(default is true), using the specified confidence level for "
            "2-bit QA fields\n");
    printf ("    -fill: specifies the fill bit should be output\n");
    printf ("    -drop_frame: specifies the dropped frame bit should be "
            "output\n");
    printf ("    -terrain_occl: specifies the terrain occlusion bit should "
            "be output\n");
    printf ("    -water: specifies the water confidence should be output "
            "using the specified confidence level\n");
    printf ("    -cloud_shadow: specifies the cloud shadow confidence should "
            "be output using the specified confidence level\n");
    printf ("    -veg: specifies the vegetation confidence should be output "
            "using the specified confidence level\n");
    printf ("    -snow_ice: specifies the snow/ice confidence should be output "
            "using the specified confidence level\n");
    printf ("    -cirrus: specifies the cirrus confidence should be output "
            "using the specified confidence level\n");
    printf ("    -cloud: specifies the cloud confidence should be output "
            "using the specified confidence level\n");
    printf ("\nwhere the conf_level can be 'low', 'med', or 'high' and the "
            "default, if not specified, is medium confidence.\n");
    printf ("\nunpack_oli_qa --help will print the usage statement\n");
    printf ("\nThe following example will unpack all the QA bits into their "
            "own single-band GeoTIFF files.  This will use the default of "
            "medium confidence (and above) for the 2-bit quality fields.\n");
    printf ("unpack_oli_qa "
            "--ifile=LC80340322013132LGN01_BQA.tif "
            "--ofile=LC80340322013132LGN01 "
            "--all\n");
    printf ("\nThe following example will unpack the fill, water, vegetation, "
            "snow/ice, and cloud quality fields each into their own GeoTIFF "
            "file.  The fill field is a single bit field and does not require "
            "a confidence level.  The water pixels will be masked if their "
            "confidence level is high.  The vegetation pixels will be masked "
            "if their confidence level is low, medium, or high.  The snow/ice "
            "pixels will be masked if their confidence is medium (by default) "
            "or high.  And, the cloud pixels will also be masked if their "
            "confidence level is medium or high.\n");
    printf ("unpack_oli_qa "
            "--ifile=LC80340322013132LGN01_BQA.tif "
            "--ofile=LC80340322013132LGN01 "
            "--fill --water=high --veg=low --snow_ice --cloud=med\n");
    printf ("\nThe following example will unpack the fill, cloud, and cirrus "
            "quality fields each into one combined file.  The fill field is a "
            "single bit field and does not require a confidence level.  The "
            "cloud pixels will be masked if their confidence level is high.  "
            "The cirrus pixels will be masked if their confidence level is "
            "low, medium, or high.\n");
    printf ("unpack_oli_qa "
            "--ifile=LC80340322013132LGN01_BQA.tif "
            "--ofile=LC80340322013132LGN01_mask.tif "
            "--fill --cloud=high --cirrus=low --combine\n");
}
