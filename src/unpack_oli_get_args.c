#include <getopt.h>
#include "unpack_oli_qa.h"

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
6/19/2013   Gail Schmidt     Original Development
6/21/2013   Gail Schmidt     Added support for confidence levels to be
                             specified for the 2-bit quality fields

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
    static int drop_frame_flag=false;    /* dropped frame band flag */
    static int terrain_occl_flag=false;  /* terrain occlusion band flag */
    static int water_flag=false;         /* water confidence band flag */
    static int cloud_shadow_flag=false;  /* cloud shadow band flag */
    static int veg_flag=false;           /* vegetation confidence band flag */
    static int snow_ice_flag=false;      /* snow/ice confidence band flag */
    static int cirrus_flag=false;        /* cirrus confidence band flag */
    static int cloud_flag=false;         /* cloud confidence band flag */
    char errmsg[STR_SIZE];               /* error message */
    char FUNC_NAME[] = "get_args";       /* function name */

    static struct option long_options[] =
    {
        {"combine", no_argument, &combine_flag, true},
        {"fill", no_argument, &fill_flag, true},
        {"drop_frame", no_argument, &drop_frame_flag, true},
        {"terrain_occl", no_argument, &terrain_occl_flag, true},
        {"all", optional_argument, 0, 'a'},
        {"water", optional_argument, 0, 'w'},
        {"cloud_shadow", optional_argument, 0, 'd'},
        {"veg", optional_argument, 0, 'v'},
        {"snow_ice", optional_argument, 0, 's'},
        {"cirrus", optional_argument, 0, 'r'},
        {"cloud", optional_argument, 0, 'c'},
        {"ifile", required_argument, 0, 'i'},
        {"ofile", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Initialize the confidence levels for the QA fields.  Single bit QA
       fields will be undefined.  The two-bit fields will be medium. */
    qa_conf[FILL] = UNDEFINED;
    qa_conf[DROPPED_FRAME] = UNDEFINED;
    qa_conf[TERRAIN_OCCL] = UNDEFINED;
    qa_conf[WATER] = MED;
    qa_conf[CLOUD_SHADOW] = MED;
    qa_conf[VEG] = MED;
    qa_conf[SNOW_ICE] = MED;
    qa_conf[CIRRUS] = MED;
    qa_conf[CLOUD] = MED;

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
                        qa_conf[WATER] = LOW;
                        qa_conf[CLOUD_SHADOW] = LOW;
                        qa_conf[VEG] = LOW;
                        qa_conf[SNOW_ICE] = LOW;
                        qa_conf[CIRRUS] = LOW;
                        qa_conf[CLOUD] = LOW;
                    }
                    else if (!strcmp (optarg, "high"))
                    {
                        qa_conf[WATER] = HIGH;
                        qa_conf[CLOUD_SHADOW] = HIGH;
                        qa_conf[VEG] = HIGH;
                        qa_conf[SNOW_ICE] = HIGH;
                        qa_conf[CIRRUS] = HIGH;
                        qa_conf[CLOUD] = HIGH;
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
            case 'w':  /* water */
                water_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[WATER] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[WATER] = HIGH;
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

            case 'v':  /* vegetation */
                veg_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[VEG] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[VEG] = HIGH;
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "vegetation band", optarg);
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

            case 'c':  /* cloud */
                cloud_flag = true;
                /* Validate the confidence argument.  Default of medium is
                   already set */
                if (optarg != NULL)
                {
                    if (!strcmp (optarg, "low"))
                        qa_conf[CLOUD] = LOW;
                    else if (!strcmp (optarg, "high"))
                        qa_conf[CLOUD] = HIGH;
                    else if (strcmp (optarg, "med"))
                    {
                        sprintf (errmsg, "Unknown confidence level of %s for "
                            "cloud band", optarg);
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
        sprintf (errmsg, "Input OLI QA file is a required argument");
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

    /* Check if the bits are to be combined */
    *combine_bits = false;
    if (combine_flag)
        *combine_bits = true;

    /* If none of the quality band flags were specified, then set the all_flag
       to true as the default */
    if (!all_flag && !fill_flag && !drop_frame_flag && !terrain_occl_flag &&
        !water_flag && !cloud_shadow_flag && !veg_flag && !snow_ice_flag &&
        !cirrus_flag && !cloud_flag)
        all_flag = true;

    /* If the all flag was specified, then turn all the quality bands on for
       processing */
    if (all_flag)
    {
        for (i = 0; i < NQUALITY_TYPES; i++)
            qa_specd[i] = true;
        return (SUCCESS);
    }

    /* Initialize the QA array to false */
    for (i = 0; i < NQUALITY_TYPES; i++)
        qa_specd[i] = false;

    /* Set up the array to depict which quality bands need to be processed */
    if (fill_flag)
        qa_specd[FILL] = true;
    if (drop_frame_flag)
        qa_specd[DROPPED_FRAME] = true;
    if (terrain_occl_flag)
        qa_specd[TERRAIN_OCCL] = true;
    if (water_flag)
        qa_specd[WATER] = true;
    if (cloud_shadow_flag)
        qa_specd[CLOUD_SHADOW] = true;
    if (veg_flag)
        qa_specd[VEG] = true;
    if (snow_ice_flag)
        qa_specd[SNOW_ICE] = true;
    if (cirrus_flag)
        qa_specd[CIRRUS] = true;
    if (cloud_flag)
        qa_specd[CLOUD] = true;

    return (SUCCESS);
}
