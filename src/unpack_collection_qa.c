#include "unpack_collection_qa.h"

static const char conf_vals[NCONF_TYPES] = {'u', 'l', 'm', 'h'};

/******************************************************************************
MODULE:  unpack_collection_qa

PURPOSE:  Unpack the specified collection QA band and write the quality bands 
to individual GeoTIFF files using the specified base output filename.

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
    int satellite_number = 0; /* number of the satellite, e.g.: 8 */

    printf ("Unpack of QA band started ...\n");

    /* Read the command-line arguments to determine which file needs to be
       processed and which quality bands will be dumped */
    retval = get_args (argc, argv, &combine_bits, &satellite_number, 
        &qa_infile, &qa_outfile, qa_specd, qa_conf);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Tell the user which bands will be unpacked */
    printf ("QA input file: %s\n", qa_infile);
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

    /* For L4-7 this is dropped pixel.  For L8 it's terrain occlusion */
    if (qa_specd[OCCLUSION_OR_DROPPED])
        tmp_char = 'x';
    else
        tmp_char = ' ';
    if (satellite_number == 8)
    {
        printf ("   %c       Terrain occlusion\n", tmp_char);
    }
    else
    {
        printf ("   %c       Dropped pixel\n", tmp_char);
    }

    if (qa_specd[RADIOMETRIC_SAT])
        tmp_char = 'x';
    else
        tmp_char = ' ';
    printf ("   %c       Radiometric saturation\n", tmp_char);

    if (qa_specd[CLOUD])
        tmp_char = 'x';
    else
        tmp_char = ' ';
    printf ("   %c       Cloud\n", tmp_char);

    if (qa_specd[CLOUD_CONFIDENCE])
        tmp_char = conf_vals[qa_conf[CLOUD_CONFIDENCE]];
    else
        tmp_char = ' ';
    printf ("   %c       Cloud confidence\n", tmp_char);

    if (qa_specd[CLOUD_SHADOW])
        tmp_char = conf_vals[qa_conf[CLOUD_SHADOW]];
    else
        tmp_char = ' ';
    printf ("   %c       Cloud shadow confidence\n", tmp_char);

    if (qa_specd[SNOW_ICE])
        tmp_char = conf_vals[qa_conf[SNOW_ICE]];
    else
        tmp_char = ' ';
    printf ("   %c       Snow/ice confidence\n", tmp_char);

    if (satellite_number == 8)
    {
        if (qa_specd[CIRRUS])
            tmp_char = conf_vals[qa_conf[CIRRUS]];
        else
            tmp_char = ' ';
        printf ("   %c       Cirrus confidence\n", tmp_char);
    }

    /* Read the input QA band, unpack the bits, combine if specified, and
       write out the desired band(s) */
    if (!combine_bits)
    {
        /* Unpack the bits into individual bands */
        retval = unpack_bits (qa_infile, qa_outfile, qa_specd, qa_conf,
            satellite_number);
        if (retval != SUCCESS)
        {   /* unpack_bits already printed the error message */
            exit (ERROR);
        }
    }
    else
    {
        /* Unpack the bits and combine into one band */
        retval = unpack_combine_bits (qa_infile, qa_outfile, qa_specd, qa_conf,
            satellite_number);
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
    printf ("Unpack of QA band complete!\n");
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
8/31/2016     Ray Dittmeier    Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("unpack_collection_qa will read the QA band from the Landsat L4-8 "
            "input file, then unpack this band into individual QA bands "
            "stored in multiple files using the user-specified base output "
            "filename.  The output bands will refer to the QA bits (from "
            "right to left), representing the QA information which is stored "
            "in the QA band.  In some cases a single bit is used to represent "
            "quality data in the QA band and in other cases two bits are used "
            "for the quality info.\n\n"
            "For quality data represented by a single bit, the output values "
            "are as follows:\n"
            "  0 = No, this condition does not exist\n"
            "  1 = Yes, this condition exists\n"
            "\n"
            "For radiometric saturation, the interpretation of the bits is:\n"
            " 00 = No bands contain saturation\n"
            " 01 = 1-2 bands contain saturation (low)\n"
            " 10 = 3-4 bands contain saturation (medium)\n"
            " 11 = 5 or more bands contain saturation (high)\n"
            "\n"
            "For other quality data represented by two bits, the user will be "
            "allowed to specify the confidence levels included in the mask.  "
            "The current confidence levels in the QA band are as follows:\n"
            " 00 = Algorithm did not determine the status of this condition\n"
            " 01 = Algorithm has low to no confidence that this condition "
            "exists (0-33 percent confidence)\n"
            " 10 = Algorithm has medium confidence that this condition exists "
            "(34-66 percent confidence)\n"
            " 11 = Algorithm has high confidence that this condition exists "
            "(67-100 percent confidence)\n"
            "\n"
            "If the user specifies a confidence level of 'low' for the 2-bit "
            "confidence field, then the output mask will be flagged as 1 / yes "
            "if the 2-bit confidence value is low, medium, or high.  If the "
            "user specifies a confidence level of 'med' for the confidence "
            "field, then the output mask will be flagged if the confidence "
            "value is medium or high.  And, if the user specifies a confidence "
            "level of 'high', then the output mask will be flagged if the "
            "confidence value is high.\n"
            "\n"
            "The following tables identify the output quality band and how it "
            "correlates to the bits in the individual QA bands, when not using "
            "combined bits.  The user may elect to combine the specified QA "
            "bits into one single output file.  In that case, if any of the "
            "QA bits are turned on or meet the specified confidence level, "
            "then the output mask for that pixel will be flagged as "
            "1 / yes.\n\n"
            "Landsat 4, 5 and 7\n"
            "\n"
            "QA Band               QA Bit(s)    Description\n"
            "------------------    ---------    -----------\n"
            "_fill.tif                0         Fill\n"
            "_dropped_pixel.tif       1         Dropped pixel\n"
            "_radiometric_sat.tif     2,3       Radiometric saturation\n"
            "_cloud.tif               4         Cloud\n"
            "_cloud_confidence.tif    5,6       Cloud confidence\n"
            "_cloud_shadow.tif        7,8       Cloud shadow\n"
            "_snow_ice.tif            9,10      Snow/ice confidence\n"
            "   N/A                   11 - 15   Reserved\n"
            "\n"
            "Landsat 8\n"
            "\n"
            "QA Band               QA Bit(s)    Description\n"
            "------------------    ---------    -----------\n"
            "_fill.tif                0         Fill\n"
            "_terrain_occl.tif        1         Terrain occlusion\n"
            "_radiometric_sat.tif     2,3       Radiometric saturation\n"
            "_cloud.tif               4         Cloud\n"
            "_cloud_confidence.tif    5,6       Cloud confidence\n"
            "_cloud_shadow.tif        7,8       Cloud shadow\n"
            "_snow_ice.tif            9,10      Snow/ice confidence\n"
            "_cirrus.tif              11,12     Cirrus confidence\n"
            "   N/A                   13,14,15  Reserved\n"
            "\n"
            "For more information about the Landsat QA Band file, please "
            "refer to http://landsat.usgs.gov/collectionqualityband.php\n\n");
    printf ("unpack_collection_qa --help will print the usage information\n\n");
    printf ("usage: unpack_collection_qa "
            "--ifile=input_QA_filename "
            "--ofile=output_unpacked_QA_filename "
            "[--all=conf_level][--fill][--drop_pixel][--terrain_occl] "
            "[--radiometric_sat][--cloud][--cloud_confidence=conf_level] "
            "[--cloud_shadow=conf_level][--snow_ice=conf_level] "
            "[--cirrus=conf_level] [--combine]\n");
    printf ("\nwhere --drop_pixel is only available for Landsat 4-7 files\n"
            "and --terrain_occl and --cirrus are only available for Landsat 8\n"
            "files\n");
    printf ("\nwhere the following parameters are required:\n");
    printf ("    -ifile: name of the input QA file (GeoTIFF product with "
            "uint16 bands).  The name should follow the Landsat collection "
            "filename format\n");
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
    printf ("    -drop_pixel: specifies the dropped pixel bit should be "
            "output (L4-7 scenes only)\n");
    printf ("    -terrain_occl: specifies the terrain occlusion bit should "
            "be output (L8 scenes only)\n");
    printf ("    -radiometric_sat: specifies the radiometric saturation bits "
            "should be output\n");
    printf ("    -cloud: specifies the cloud bit should be output\n");
    printf ("    -cloud_confidence: specifies the cloud confidence should be "
            "output using the specified confidence level\n");
    printf ("    -cloud_shadow: specifies the cloud shadow confidence should "
            "be output using the specified confidence level\n");
    printf ("    -snow_ice: specifies the snow/ice confidence should be output "
            "using the specified confidence level\n");
    printf ("    -cirrus: specifies the cirrus confidence should be output "
            "using the specified confidence level (L8 scenes only)\n");
    printf ("\nwhere the conf_level can be 'low', 'med', or 'high' and the "
            "default, if not specified, is medium confidence.\n");
    printf ("\nunpack_collection_qa --help will print the usage statement\n");
    printf ("\nThe following example will unpack all the QA bits into their "
            "own single-band GeoTIFF files.  This will use the default of "
            "medium confidence (and above) for the 2-bit quality fields, "
            "excluding radiometric saturation which doesn't use a confidence "
            "level.\n");
    printf ("unpack_collection_qa "
            "--ifile=LE07_L1TP_042027_20050927_20160512_01_T1_BQA.TIF "
            "--ofile=LE07_L1TP_042027_20050927_20160512_01_T1_BQA "
            "--all\n");
    printf ("\nThe following example will unpack the fill, radiometric "
            "saturation, cloud shadow, snow/ice, and cloud confidence fields "
            "each into their own GeoTIFF file.  The fill field is a single bit "
            "field and does not require a confidence level.  The radiometric "
            "saturation pixels also do not require a confidence level. The "
            "cloud shadow pixels will be masked if their confidence level is "
            "low, medium, or high.  The snow/ice pixels will be masked if "
            "their confidence is medium (by default) or high.  The cloud "
            "confidence pixels will also be masked if their confidence level "
            "is medium or high.\n");
    printf ("unpack_collection_qa "
            "--ifile=LE07_L1TP_042027_20050927_20160512_01_T1_BQA.TIF "
            "--ofile=LE07_L1TP_042027_20050927_20160512_01_T1_BQA "
            "--fill --radiometric_sat --cloud_shadow=low --snow_ice "
            "--cloud_confidence=med\n");
    printf ("\nThe following example will unpack the fill, cloud shadow, and "
            "cirrus quality fields each into one combined file.  The fill "
            "field is a single bit field and does not require a confidence "
            "level. The cloud shadow pixels will be masked if their confidence "
            "level is high. The cirrus pixels will be masked if their "
            "confidence level is low, medium, or high.\n");
    printf ("unpack_collection_qa "
            "--ifile=LC08_L1GT_029030_20151209_2015013_01_T1_BQA.tif "
            "--ofile=LC08_L1GT_029030_20151209_2015013_01_T1_BQA "
            "--fill --cloud_shadow=high --cirrus=low --combine\n");
}
