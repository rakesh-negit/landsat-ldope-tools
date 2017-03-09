#include <getopt.h>
#include "unpack_collection_qa.h"

/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
8/31/2016   Ray Dittmeier    Original Development based on unpack_oli_get_args.c

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
******************************************************************************/
short get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    bool *combine_bits,   /* O: should the QA bits be combined? */
    int *satellite_number, /* O: number of the satellite, e.g.: 8 */
    char **infile,        /* O: address of input filename */
    char **outfile,       /* O: address of output filename or base filename */
    bool qa_specd[NQUALITY_TYPES],  /* O: array to specify which of the QA
                                          bands was specified for processing */
    Confidence_t qa_conf[NQUALITY_TYPES]
                          /* O: array to specify the confidence level for
                                each of the quality fields */
)
{
    int c;                               /* current argument index */
    int option_index;                    /* index for command-line option */
    int i;                               /* looping variable */
    static int combine_flag=false;       /* all quality bands flag */
    static int all_flag=false;           /* all quality bands flag */
    static int fill_flag=false;          /* fill band flag */
    static int drop_pixel_flag=false;    /* L4-7 dropped pixel band flag */
    static int terrain_occl_flag=false;  /* L8 terrain occlusion band flag */
    static int radiometric_sat_flag=false; /* radiometric saturation band 
                                            flag */
    static int cloud_flag=false;         /* cloud band flag */
    static int cloud_confidence_flag=false; /* cloud confidence band flag */
    static int cloud_shadow_flag=false;  /* cloud shadow band flag */
    static int snow_ice_flag=false;      /* snow/ice confidence band flag */
    static int cirrus_flag=false;        /* cirrus confidence band flag */
    char errmsg[STR_SIZE];               /* error message */
    char satellite_number_str[3];        /* String for the satellite number */

    char FUNC_NAME[] = "get_args";       /* function name */

    static struct option long_options[] =
    {
        {"combine", no_argument, &combine_flag, true},
        {"fill", no_argument, &fill_flag, true},
        {"drop_pixel", no_argument, &drop_pixel_flag, true},
        {"terrain_occl", no_argument, &terrain_occl_flag, true},
        {"cloud", no_argument, &cloud_flag, true},
        {"all", optional_argument, 0, 'a'},
        {"radiometric_sat", no_argument, &radiometric_sat_flag, true},
        {"cloud_confidence", optional_argument, 0, 'n'},
        {"cloud_shadow", optional_argument, 0, 'd'},
        {"snow_ice", optional_argument, 0, 's'},
        {"cirrus", optional_argument, 0, 'r'},
        {"ifile", required_argument, 0, 'i'},
        {"ofile", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Initialize the confidence levels for the QA fields.  Single bit QA
       fields will be undefined.  The two-bit fields will be medium. */
    qa_conf[FILL] = UNDEFINED;
    qa_conf[OCCLUSION_OR_DROPPED] = UNDEFINED;
    qa_conf[RADIOMETRIC_SAT] = UNDEFINED;
    qa_conf[CLOUD] = UNDEFINED;
    qa_conf[CLOUD_CONFIDENCE] = MED;
    qa_conf[CLOUD_SHADOW] = MED;
    qa_conf[SNOW_ICE] = MED;
    qa_conf[CIRRUS] = MED;

    /* Loop through all the cmd-line options */
    opterr = 0;   /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {   /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
     
            case 'h':  /* help */
                usage ();
                return (ERROR);
                break;

            case 'i':  /* ifile */
                *infile = strdup (optarg);
                break;
     
            case 'o':  /* ofile */
                *outfile = strdup (optarg);
                break;
     
            case 'a':  /* all */
                all_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                    {
                        qa_conf[CLOUD_CONFIDENCE] = LOW;
                        qa_conf[CLOUD_SHADOW] = LOW;
                        qa_conf[SNOW_ICE] = LOW;
                        qa_conf[CIRRUS] = LOW;
                    }
                    else if (!strcmp (optarg, "high"))
                    {
                        qa_conf[CLOUD_CONFIDENCE] = HIGH;
                        qa_conf[CLOUD_SHADOW] = HIGH;
                        qa_conf[SNOW_ICE] = HIGH;
                        qa_conf[CIRRUS] = HIGH;
                    }
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "water band", optarg);
                        error_handler (true, FUNC_NAME, errmsg);
                        usage ();
                        return (ERROR);
                    }
                }
                break;

            case 'n':  /* cloud_confidence */
                cloud_confidence_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[CLOUD_CONFIDENCE] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[CLOUD_CONFIDENCE] = HIGH;
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "cloud confidence band", optarg);
                        error_handler (true, FUNC_NAME, errmsg);
                        usage ();
                        return (ERROR);
                    }
                }
                break;

            case 'd':  /* cloud_shadow */
                cloud_shadow_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[CLOUD_SHADOW] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[CLOUD_SHADOW] = HIGH;
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "cloud shadow band", optarg);
                        error_handler (true, FUNC_NAME, errmsg);
                        usage ();
                        return (ERROR);
                    }
                }
                break;

            case 's':  /* snow_ice */
                snow_ice_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[SNOW_ICE] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[SNOW_ICE] = HIGH;
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "snow/ice band", optarg);
                        error_handler (true, FUNC_NAME, errmsg);
                        usage ();
                        return (ERROR);
                    }
                }
                break;

            case 'r':  /* cirrus */
                cirrus_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[CIRRUS] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[CIRRUS] = HIGH;
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "cirrus band", optarg);
                        error_handler (true, FUNC_NAME, errmsg);
                        usage ();
                        return (ERROR);
                    }
                }
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind-1]);
                error_handler (true, FUNC_NAME, errmsg);
                usage ();
                return (ERROR);
                break;
        }
    }

    /* Make sure the infiles and outfiles were specified */
    if (*infile == NULL)
    {
        sprintf (errmsg, "Input QA file is a required argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    if (*outfile == NULL)
    {
        sprintf (errmsg, "Unpacked bits output base QA file is a required "
            "argument");
        error_handler (true, FUNC_NAME, errmsg);
        usage ();
        return (ERROR);
    }

    /* Assume the input file follows the Landsat product identifier format.
       In the collection era, that starts with LXSS, where SS is the satellite
       number.  Extract the satellite number */
    strncpy(satellite_number_str, *infile + 2,
        sizeof(satellite_number_str) - 1);
    satellite_number_str[2] = '\0';
    *satellite_number = atoi(satellite_number_str);
    if (*satellite_number != 4 && *satellite_number != 5 &&
        *satellite_number != 7 && *satellite_number != 8)
    {
        sprintf (errmsg, "Error with filename format: the filename should "
            "adhere to the Landsat collection filename format with satellite "
            "number in positions 3 and 4.  This tool supports satellites 4, "
            "5, 7, and 8 with format 04, 05, 07, and 08.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (*satellite_number == 8)
    {
        if (drop_pixel_flag)
        {
            sprintf (errmsg, "Dropped pixel is not supported for this "
                "satellite.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    } 
    else
    {
        if (cirrus_flag)
        {
            sprintf (errmsg, "Cirrus is not supported for this satellite.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        if (terrain_occl_flag)
        {
            sprintf (errmsg, "Terrain occlusion is not supported for this "
                "satellite.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Check if the bits are to be combined */
    *combine_bits = false;
    if (combine_flag)
        *combine_bits = true;

    /* If none of the quality band flags were specified, then set the all_flag
       to true as the default */
    if (*satellite_number == 8)
    {
        if (!all_flag && !fill_flag && !terrain_occl_flag &&
            !radiometric_sat_flag && !cloud_flag && !cloud_confidence_flag &&
            !cloud_shadow_flag && !snow_ice_flag && !cirrus_flag)
            all_flag = true;
    }
    else
    {
        if (!all_flag && !fill_flag && !drop_pixel_flag &&
            !radiometric_sat_flag && !cloud_flag && !cloud_confidence_flag &&
            !cloud_shadow_flag && !snow_ice_flag)
            all_flag = true;
    }

    /* Initialize the QA array to false */
    for (i = 0; i < NQUALITY_TYPES; i++)
        qa_specd[i] = false;

    /* If the all flag was specified, then turn all the quality bands on for
       processing */
    if (all_flag)
    {
        qa_specd[FILL] = true;
        qa_specd[OCCLUSION_OR_DROPPED] = true;
        qa_specd[RADIOMETRIC_SAT] = true;
        qa_specd[CLOUD] = true;
        qa_specd[CLOUD_CONFIDENCE] = true;
        qa_specd[CLOUD_SHADOW] = true;
        qa_specd[SNOW_ICE] = true;
        if (*satellite_number == 8)
        {
            qa_specd[CIRRUS] = true;
        }
        return (SUCCESS);
    }

    /* Set up the array to depict which quality bands need to be processed */
    if (fill_flag)
        qa_specd[FILL] = true;
    if (*satellite_number != 8)
    {
        if (drop_pixel_flag)
            qa_specd[OCCLUSION_OR_DROPPED] = true;
    }
    if (*satellite_number == 8)
    {
        if (terrain_occl_flag)
            qa_specd[OCCLUSION_OR_DROPPED] = true;
    }
    if (radiometric_sat_flag)
        qa_specd[RADIOMETRIC_SAT] = true;
    if (cloud_flag)
        qa_specd[CLOUD] = true;
    if (cloud_confidence_flag)
        qa_specd[CLOUD_CONFIDENCE] = true;
    if (cloud_shadow_flag)
        qa_specd[CLOUD_SHADOW] = true;
    if (snow_ice_flag)
        qa_specd[SNOW_ICE] = true;
    if (*satellite_number == 8)
    {
        if (cirrus_flag)
            qa_specd[CIRRUS] = true;
    }

    return (SUCCESS);
}
