/****************************************************************************
!C

!File: reduce_sds.c

!Description:
  Contains routines for sds reduction 

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  This software was developed by:
    Land Data Operational Product Evaluation (LDOPE) Team for the
    National Aeronautics and Space Administration, Goddard Space Flight
    Center

!References and Credits:

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

    Yi Zhang 
    LDOPE                             Science Systems and Applications Inc.
    sophyi@ltpmail.gsfc.nasa.gov      NASA/GSFC Code 922
    phone: 301-614-5497               Greenbelt, MD 20771

    David Roy 
    LDOPE                             University of Maryland 
    droy@kratmos.gsfc.nasa.gov	      NASA/GSFC Code 922
    phone: 301-614-5571               Greenbelt, MD 20771

!Design Notes: (none)

!END
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mfhdf.h"
#include "qa_tool.h"
#include "str_op.h"
#include "sds_rw.h"
#include "meta.h"
#include "main_util.h"
#include "alloc_mem.h"

#define HELP \
"NAME \n" \
"    reduce_sds - reduce the spatial dimensions of one or more SDS of a\n" \
"    MODIS Land HDF-EOS data product.\n" \
" \n" \
"SYNOPSIS \n" \
"    reduce_sds -help [filename] \n" \
"    reduce_sds -of=<output_filename> -rf=<reduction factor> -sub|avg|cnt|cl\n"\
"           [-sds=<SDSname1>[<,SDSname2>...]] \n" \
"           [-bit=[<bit range>]<opr><value>[,[<bit range>]<opr><value>]...]\n"\
"           [-min] [-max] [-std] [-num] [-meta] [-float] filename\n" \
" \n" \
"DESCRIPTION \n" \
"    The spatial dimensions of input SDS(s) may be reduced using one of four\n"\
"    different methods. The reduction factor (rf) must be a non-zero\n"\
"    positive integer. The output SDS x dimension will be ((x div rf) + \n"\
"    (x mod rf)) and similarly the y dimension will be ((y div rf) +\n"\
"    (y mod rf)). If the input SDS list contains SDSs with different spatial\n"\
"    dimensions the reduction factor will be applied to the SDS with the\n"\
"    smallest spatial dimension and the other SDS(s) will be reduced to have\n"\
"    the same output dimension. \n" \
" \n" \
"    All SDS fill values are ignored.\n" \
" \n" \
"    This tool may be used to reduce data volumes, and to enable analysis of\n"\
"    the different MODIS Land data product spatial resolutions\n"\
"    (250m, 500m, 1km), or to enable quick comparison with other coarser\n"\
"    spatial resolution data sets. \n" \
" \n" \
"    This tool compliments the tool enlarge_sds. \n" \
"  \n" \
"    This tool supports 2D/3D/4D SDS(s).\n" \
" \n" \
"    The tool command arguments can be specified in any order. \n" \
" \n" \
"ARGUMENTS \n" \
"    -help  [filename]        Display this help message. If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDS in the file are\n" \
"                             displayed. \n" \
"    -of=<filename>           Output filename \n" \
"    -rf=<reduction factor>   Reduction factor (a non-zero positive integer)\n"\
"    -sub                     Reduce by sub-sampling. Pixel value at (i, j)\n"\
"                             in the output SDS is copied from the pixel at\n"\
"                             (rf*i + rf/2, rf*j + rf/2) in the input SDS.\n" \
"    -avg                     Reduce by averaging. Pixel value at (i, j) in\n"\
"                             the output SDS is the average of pixel values\n"\
"                             in a sub window defined from \n" \
"                             (rf*i, rf*j) to (rf*i + rf - 1, rf*j + rf -1)\n"\
"                             in the input SDS.  Optional minimum, maximum,\n"\
"                             standard deviation, and number of averaged\n"\
"                             pixels may be output as separate SDS.\n" \
"    -cnt                     Reduce by pixel counting. Pixel value at\n"\
"                             (i, j) in the output SDS is the number of\n"\
"                             pixels in a sub window defined from \n" \
"                             (rf*i, rf*j) to (rf*i + rf -1, rf*j + rf -1)\n"\
"                             in the input SDS with bit value equal to the\n"\
"                             user specified value. Relational operators\n" \
"                             may be used. The output SDS value is in the \n" \
"                             range {0 : rf x rf}.\n" \
"    -cl                      Reduce by majority class value. Pixel value\n"\
"                             at (i, j) in the output SDS is set to the\n"\
"                             majority value of the pixels in a sub window\n" \
"                             defined from (rf*i, rf*j) to\n"\
"                             (rf*i + rf -1, rf*j + rf -1).  This option is\n"\
"                             best used for data with a small number of\n" \
"                             nominal pixel values. \n" \
"    -sds=<SDS list>          List of SDS to reduce. SDS names are separated\n"\
"                             by commas with no space. By default all SDSs\n"\
"                             are processed maintaining the input SDS\n"\
"                             interleaving. \n" \
"                             To process a specific layer of a 3D SDS\n"\
"                             specify the element number of the third\n"\
"                             dimension as a dot extension of the SDS name:\n"\
"                             sds_name.n (e.g., sur_refl_b02.1 = the layer\n" \
"                             defined by the 1st element of the 3rd\n"\
"                             dimension of the 3D SDS sur_refl_b02). \n" \
"                             To process a specific layer of a 4D SDS,\n"\
"                             specify the higher dimension element number(s)\n"\
"                             as a dot extension of the SDS name:\n"\
"                             sds_name.n.m (e.g., Surface_Refl.1.2 = the\n"\
"                             layer defined by the 1st element of the 3rd\n"\
"                             dimension and the 2nd element of the 4th\n"\
"                             dimension of the 4D SDS Surface_Refl). \n" \
"                             Note that wildcards and ranges of element\n"\
"                             values may be specified as sds_name.* and as\n"\
"                             sds_name.n1-n2.m respectively. \n" \
"    -bit=[<bit range>]<operator><value>,[<bit range>]<operator><value>, ...\n"\
"                             This option is applicable only with the -cnt\n"\
"                             option. The SDS bit range and corresponding\n"\
"                             value are specified in decimal separated by a\n"\
"                             relational operator. If the bit range is not\n"\
"                             specified the tool considers all the bits in\n"\
"                             the SDS.  Multiple range, operator and value\n"\
"                             combinations separated by commas will result\n"\
"                             in separate output SDS for each of such\n"\
"                             combination. Valid relational operators are:\n"\
"                             ==, <, >, <=, >=, != \n" \
"    -std                     Compute the standard deviation in each sub\n"\
"                             window  (for -avg option only).\n" \
"    -min                     Compute the minimum value in each sub window\n"\
"                             (for -avg option only). \n" \
"    -max                     Compute the maximum value in each sub window\n"\
"                             (for -avg option only). \n" \
"    -num                     Compute the number of averaged pixels in each\n"\
"                             sub window (for -avg option only). \n" \
"    -meta                    Copy metadata from input file to output file.\n" \
"    -float                   Output average value SDS in float data type,\n"\
"                             default is the input data type (for -avg\n"\
"                             option only). \n" \
"    filename                 Input filename. \n" \
"\n" \
"EXAMPLES \n" \
"    reduce_sds -sds=Fpar_1km -sub -rf=10 -of=sub_fpar.hdf \n" \
"               MOD15A1.A2001193.h09v05.004.2002198025239.hdf \n" \
"    {Note: This example reduces the spatial dimensions of SDS Fpar_1km by\n"\
"           10 using the sub-sampling method. The output file will have\n"\
"           spatial dimensions approximately 10 times smaller.} \n" \
" \n" \
"    reduce_sds -sds=Lai_1km -avg -min -max -std -num -rf=10 \n" \
"               -of=avg_lai.hdf MOD15A1.A2001193.h09v05.004.2002198025239.hdf\n"\
" \n" \
"    reduce_sds -sds=Cloud_Mask -sub -rf=5 -of=cloud_sub.hdf \n" \
"               MYD35_L2.A2002189.2040.003.2002191125354.hdf \n" \
" \n" \
"    reduce_sds -sds=Cloud_Mask.1-2 -avg -rf=5 -of=cloud_avg.hdf \n" \
"               MYD35_L2.A2002189.2040.003.2002191125354.hdf \n" \
" \n" \
"    reduce_sds -cnt -bit=\"0-3<=2,0-3==3,0-3==4,0-3==5,0-3==6,0-3==7\"\n"\
"               -meta -sds=\"most confident detected fire\"\n"\
"               -rf=5 -of=fire_class.hdf \n" \
"               MOD14A1.A2002185.h30v11.003.2002204204451.hdf  \n" \
" \n" \
"    reduce_sds -rf=2 -cl -of=land_cover_class.hdf \n" \
"               MOD12Q1.A2000289.h01v11.003.2002171021653.hdf \n" \
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
"\n" \
"Version 1.0, 08/08/2002 \n"

#define USAGE \
"usage: \n" \
"    reduce_sds -help [filename] \n" \
"    reduce_sds -of=<output_filename> -rf=<reduction factor> -sub|avg|cnt|cl\n"\
"           [-sds=<SDSname1>[<,SDSname2>...]] \n" \
"           [-bit=[<bit range>]<opr><value>[,[<bit range>]<opr><value>]...]\n"\
"[-min] [-max] [-std] [-num] [-meta] [-float] filename \n" \
" \n" \
"OPTIONS \n" \
"    -help  [filename]        Display this help message. If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDS in the file are\n" \
"                             displayed. \n" \
"    -of=<filename>           Output filename \n" \
"    -rf=<reduction factor>   Reduction factor (a non-zero positive integer)\n"\
"    -sub                     Reduce by sub-sampling. Pixel value at (i, j)\n"\
"                             in the output SDS is copied from the pixel at\n"\
"                             (rf*i + rf/2, rf*j + rf/2) in the input SDS.\n" \
"    -avg                     Reduce by averaging. Pixel value at (i, j) in\n"\
"                             the output SDS is the average of pixel values\n"\
"                             in a sub window defined from \n" \
"                             (rf*i, rf*j) to (rf*i + rf - 1, rf*j + rf -1)\n"\
"                             in the input SDS.  Optional minimum, maximum,\n"\
"                             standard deviation, and number of averaged\n"\
"                             pixels may be output as separate SDS.\n" \
"    -cnt                     Reduce by pixel counting. Pixel value at\n"\
"                             (i, j) in the output SDS is the number of\n"\
"                             pixels in a sub window defined from \n" \
"                             (rf*i, rf*j) to (rf*i + rf -1, rf*j + rf -1)\n"\
"                             in the input SDS with bit value equal to the\n"\
"                             user specified value. Relational operators\n" \
"                             may be used. The output SDS value is in the \n" \
"                             range {0 : rf x rf}.\n" \
"    -cl                      Reduce by majority class value. Pixel value\n"\
"                             at (i, j) in the output SDS is set to the\n"\
"                             majority value of the pixels in a sub window\n" \
"                             defined from (rf*i, rf*j) to\n"\
"                             (rf*i + rf -1, rf*j + rf -1).  This option is\n"\
"                             best used for data with a small number of\n" \
"                             nominal pixel values. \n" \
"    -sds=<SDS list>          List of SDS to reduce. SDS names are separated\n"\
"                             by commas with no space. By default all SDSs\n"\
"                             are processed maintaining the input SDS\n"\
"                             interleaving. \n" \
"                             To process a specific layer of a 3D SDS\n"\
"                             specify the element number of the third\n"\
"                             dimension as a dot extension of the SDS name:\n"\
"                             sds_name.n (e.g., sur_refl_b02.1 = the layer\n" \
"                             defined by the 1st element of the 3rd\n"\
"                             dimension of the 3D SDS sur_refl_b02). \n" \
"                             To process a specific layer of a 4D SDS,\n"\
"                             specify the higher dimension element number(s)\n"\
"                             as a dot extension of the SDS name:\n"\
"                             sds_name.n.m (e.g., Surface_Refl.1.2 = the\n"\
"                             layer defined by the 1st element of the 3rd\n"\
"                             dimension and the 2nd element of the 4th\n"\
"                             dimension of the 4D SDS Surface_Refl). \n" \
"                             Note that wildcards and ranges of element\n"\
"                             values may be specified as sds_name.* and as\n"\
"                             sds_name.n1-n2.m respectively. \n" \
"    -bit=[<bit range>]<operator><value>,[<bit range>]<operator><value>, ...\n"\
"                             This option is applicable only with the -cnt\n"\
"                             option. The SDS bit range and corresponding\n"\
"                             value are specified in decimal separated by a\n"\
"                             relational operator. If the bit range is not\n"\
"                             specified the tool considers all the bits in\n"\
"                             the SDS.  Multiple range, operator and value\n"\
"                             combinations separated by commas will result\n"\
"                             in separate output SDS for each of such\n"\
"                             combination. Valid relational operators are:\n"\
"                             ==, <, >, <=, >=, != \n" \
"    -std                     Compute the standard deviation in each sub\n"\
"                             window  (for -avg option only).\n" \
"    -min                     Compute the minimum value in each sub window\n"\
"                             (for -avg option only). \n" \
"    -max                     Compute the maximum value in each sub window\n"\
"                             (for -avg option only). \n" \
"    -num                     Compute the number of averaged pixels in each\n"\
"                             sub window (for -avg option only). \n" \
"    -meta                    Copy metadata from input file to output file.\n" \
"    -float                   Output average value SDS in float data type,\n"\
"                             default is the input data type (for -avg\n"\
"                             option only). \n" \
"    filename                 Input filename. \n" \
"\n"


static unsigned int BIT[] = {
  0x1,        0x3,        0x7,        0xf,
  0x1f,       0x3f,       0x7f,       0xff,
  0x1ff,      0x3ff,      0x7ff,      0xfff,
  0x1fff,     0x3fff,     0x7fff,     0xffff,
  0x1ffff,    0x3ffff,    0x7ffff,    0xfffff,
  0x1fffff,   0x3fffff,   0x7fffff,   0xffffff,
  0x1ffffff,  0x3ffffff,  0x7ffffff,  0xfffffff,
  0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff
};

int parse_cmd_reduce_sds(int argc, char **argv, char **sds_names, int *sds_cnt, char **arg_list);
void reduce_mult_sds(char **sds_names, int sds_cnt, int32 out_sd_id, char **arg_list);
int get_bit_opt(char *cmp_str, int **bn_arr, int *bn_val, int *bn_opr, int *bn_cnt);
void reduce_an_sds_by_sub(sds_t *in_sds_info, int32 out_sd_id, int res);
void reduce_an_sds_by_avg(sds_t *in_sds_info, int32 out_sd_id, int res, int *out_flag, char avg_dt);
void reduce_an_sds_by_cnt(sds_t *in_sds_info, int32 out_sd_id, int res, char *cmp_str);
void reduce_an_sds_by_class(sds_t *in_sds_info, int32 out_sd_id, int res);

int main(int argc, char *argv[])
/****************************************************************************
!C

!Description:
  Function main() for sds reduction 

!Input Parameters:
  Command line arguments (See help for details)

!Input/Output Parameters: (none)

!Output Parameters:
  Exit status 0 on completion, and a completion message.

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)

!END
*****************************************************************************/
{
  int i, k, status;
  int nsds, isds;
  int32 out_sd_id;
  int32 sd_id, sds_id, dt;
  int32 msds, nattr, rank, dim_size[4];
  char **arg_list, **sds_names;
  char dim_str[MAX_STR_LEN], in_fname[MAX_PATH_LENGTH];
  char name[MAX_SDS_NAME_LEN], sds_name[MAX_SDS_NAME_LEN];

  if ((argc==2) && (strcmp(argv[1], "-help")==0))
  {
    /* Print command line help when -help option is input */
    fprintf(stderr, "%s\n", HELP);
    exit(0);
  }
  else
  {
    /* Print all SDS names when -help option is input with the input filename */
    for (i=0, status=0; i<argc; i++)
      if (strcmp(argv[i], "-help") == 0) status = i;
    if (status != 0)
    {
      /* Print all the SDS names in the input file */
      msds = 0;
      for (i=1; i<argc; i++)
      {
        if (argv[i][0] != '-')
        {
          if ((sd_id = SDstart(argv[i], DFACC_READ)) == FAIL)
            fprintf(stderr, "Cannot open the HDF file %s\n", argv[i]);
          else
          {
            if (SDfileinfo(sd_id, &msds, &nattr) == FAIL)
              fprintf(stderr, "Error reading information for HDF file %s\n", argv[i]);
            else
            {
              fprintf(stdout, "Valid SDS names, dimension and data type in file: %s\n", argv[i]);
              for (isds=0; isds<msds; isds++) {
                if ((sds_id = SDselect(sd_id, isds)) != FAIL) {
                  if (SDgetinfo(sds_id, name, &rank, dim_size, &dt, &nattr) != FAIL) {
                    sprintf(sds_name, "%s (" LONG_INT_FMT, name, dim_size[0]);
                    for (k=1; k<rank; k++) {
                      sprintf(dim_str, " x " LONG_INT_FMT, dim_size[k]);
                      strcat(sds_name, dim_str);
                    }
                    switch(dt) {
                      case  5: strcat(sds_name, ") FLOAT32"); break;
                      case  6: strcat(sds_name, ") FLOAT64"); break;
                      case 20: strcat(sds_name, ") INT8"); break;
                      case 21: strcat(sds_name, ") UINT8"); break;
                      case 22: strcat(sds_name, ") INT16"); break;
                      case 23: strcat(sds_name, ") UINT16"); break;
                      case 24: strcat(sds_name, ") INT32"); break;
                      case 25: strcat(sds_name, ") UINT32"); break;
                      default: strcat(sds_name, ") Unknown"); break;
                    }
                    fprintf(stdout, "\t%s\n", sds_name);
                    SDendaccess(sds_id);
                  }
                }
                else fprintf(stdout, "Error opening SDS with index %d\n", isds);
              }
            }
            SDend(sd_id);
          }
          if (msds == 0)
            fprintf(stdout, "Input files do not contain any valid SDS\n");
          break;
        }
      } /* for (i=2; . .  ) */
      exit(0);
    }
  }

  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_names in read_l2g(main)\n");
  if ((arg_list = (char **)Calloc2D(11, MAX_PARAM_LENGTH, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for arg_list in read_l2g(main)\n");
  if ((sds_names != NULL) && (arg_list != NULL))
  {
    fprintf(stdout, "Started process reduce_sds . . . . \n");

    /* Parse command arguments */
    status = parse_cmd_reduce_sds(argc, argv, sds_names, &nsds, arg_list);
    if (status == -1)
      {
	fprintf(stderr, "%s\n", USAGE);
	exit(0);
      }
    else if (status != 0)
    {
	/* open output HDF file */
      if ((out_sd_id =SDstart(arg_list[0], DFACC_CREATE)) == FAIL)
        fprintf(stderr, "Cannot create the output file %s\n", arg_list[0]);
      else
      {
	strcpy(in_fname, arg_list[10]);
	fprintf(stderr, "Processing input HDF file %s\n", in_fname);

	/* reformat SDS names with proper name extension for layer numbers */
	update_nd_sdsnames(sds_names, &nsds, in_fname);

	/* reduce all the input SDS */
        reduce_mult_sds(sds_names, nsds, out_sd_id, arg_list); 

	/* close the output HDF file */
	SDend(out_sd_id);
      }
    }
    Free2D((void **)arg_list);
    Free2D((void **)sds_names);
    fprintf(stdout, "Processing done !\n");
  }
  return 0;
}

int parse_cmd_reduce_sds(int argc, char **argv, char **sds_names, int *sds_cnt, char **arg_list)
/****************************************************************************
!C

!Description:
  Parse the command line arguments.

!Input Parameters:
  argc:	number of arguments in the command line.
  argv: string array containing the input arguments.

!Input/Output Parameters: (none)

!Output Parameters:
  sds_names: string array containing the input SDS names.
  sds_cnt: number of SDS input
  arg_list: return argument list.
  returns 1 on success and -1 on error.

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)

!END
***************************************************************************/
{
  int i;
  int if_cnt;
  int status;

  status = 1;
  if_cnt = *sds_cnt = 0;
  for (i=0; i<11; i++)
    arg_list[i][0] = '\0';

  /* parse arguments and argument values */
  for (i=1; i<argc; i++)
  {
    if (is_arg_id(argv[i], "-sds=") == 0) 
      get_arg_val_arr(argv[i], sds_names, sds_cnt);
    else if (is_arg_id(argv[i], "-sds") == 0) ;
    else if (is_arg_id(argv[i], "-of=") == 0)
      get_arg_val(argv[i], arg_list[0]);
    else if (is_arg_id(argv[i], "-rf=") == 0)
      get_arg_val(argv[i], arg_list[1]);
    else if ((strcmp(argv[i], "-sub") == 0) || (strcmp(argv[i], "-avg") == 0) || 
		(strcmp(argv[i], "-cnt") == 0) || (strcmp(argv[i], "-cl") == 0))
      strcpy(arg_list[2], argv[i]);
    else if (strcmp(argv[i], "-min") == 0) strcpy(arg_list[3], "y");
    else if (strcmp(argv[i], "-max") == 0) strcpy(arg_list[4], "y");
    else if (strcmp(argv[i], "-std") == 0) strcpy(arg_list[5], "y");
    else if (strcmp(argv[i], "-num") == 0) strcpy(arg_list[6], "y");
    else if (is_arg_id(argv[i], "-bit=") == 0) get_arg_val(argv[i], arg_list[7]);
    else if (strcmp(argv[i], "-float") == 0) strcpy(arg_list[8], "y");
    else if (strcmp(argv[i], "-meta") == 0) strcpy(arg_list[9], "y");
    else if (argv[i][0] == '-') fprintf(stderr, "Unknonwn option %s\n", argv[i]);
    else {
      if (if_cnt == 1)
	fprintf(stdout, "Only one input file is accepted, input file %s is ignored.\n", argv[i]);
      else {
	strcpy(arg_list[10], argv[i]);
        ++if_cnt;
      }
    }
  }

  /* check for required arguments */
  if (if_cnt == 0)
  {
    status = -1;
    fprintf(stderr, "No input filename specified\n");
  }
  for (i=0; i<3; i++)
    if (arg_list[i][0] == '\0') 
    {
      status = -1;
      if (i == 0) fprintf(stderr, "No output filename specified\n");
      else if (i == 1) fprintf(stderr, "No reduction factor specified\n");
      else if (i == 2) fprintf(stderr, "No reduction option sub|avg|cnt specified\n");
    }

  /* if required arguments are missing update the status to failure */
  if (status != -1)
  {
    if (strcmp(arg_list[2], "-cnt") == 0)
    {
      if (strlen(arg_list[7]) <= 0) 
      {
	fprintf(stderr, "Missing -bn option for sds reduce by cnt\n");
        status = -1;
      }
    }
    else if (strlen(arg_list[7]) != 0)
      fprintf(stderr, "Ignoring -val option for reduce option %s\n", arg_list[2]);
  }
  if ((status != -1) && (atoi(arg_list[1]) <= 0))
  {
    fprintf(stderr, "Invalid value for reduction factor \n");
    status = 0;
  }

  /* if SDS names not input set to default all SDS */
  if ((status > 0) && (*sds_cnt == 0))
  {
    fprintf(stderr, "No SDS name input. Reading all SDS . . \n");
    *sds_cnt = 1;
    strcpy(sds_names[0], "all");
  }
  return status;
}

void reduce_mult_sds(char **sds_names, int sds_cnt, int32 out_sd_id, char **arg_list) 
/****************************************************************************
!C

!Description:
  Routine to initiate reduce sds based on user options for all the input SDSs.

!Input Parameters:
  sds_names: string array containing the input SDS names.
  sds_cnt: number of SDS input
  out_sd_id: ouput SDS id.
  arg_list: argument list containing paramters.

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)

!END
***************************************************************************/
{
  int rank, out_flag[5];
  int i, isds, in_res, res;
  int xdim_min, rf[MAX_NUM_SDS];
  sds_t in_sds_info;
  char infile[MAX_PATH_LENGTH];

  /* set input SDS names */
  strcpy(infile, arg_list[10]);
  if ((sds_cnt == 1) && (strcmp(sds_names[0], "all") == 0))
  {
    sds_cnt = get_sds_names(infile, sds_names);
    if (sds_cnt == 0)
    {
      fprintf(stderr, "No SDS found in %s\n", infile);
      return;
    }
  }

  /* open the input HDF file */
  if ((in_sds_info.sd_id = SDstart(infile, DFACC_READ)) == FAIL)
  {
    fprintf(stderr, "Cannot open input HDF file in reduce_mult_sds\n");
    return;
  }

  /* open all input SDS and obtain the dimension information */
  for (isds=0; isds<sds_cnt; isds++)
  {
    strcpy(in_sds_info.name, sds_names[isds]);
    in_sds_info.sds_id = -1;
    if (get_sds_info((char *)NULL, &in_sds_info) != -1)
    {
      rank = in_sds_info.rank;
      rf[isds] = (in_sds_info.dim_size[rank-1] < in_sds_info.dim_size[0]) ?
			in_sds_info.dim_size[0] : in_sds_info.dim_size[rank-2];
      SDendaccess(in_sds_info.sds_id);
    }
    else rf[isds] = -1;
  }
  xdim_min = 9999;
  in_res = (int)atoi(arg_list[1]);
  for (isds=0; isds<sds_cnt; isds++)
    if ((rf[isds] != -1) && (rf[isds] < xdim_min))
      xdim_min = rf[isds];
  for (isds=0; isds<sds_cnt; isds++)
    if (rf[isds] != -1)
      rf[isds] = rf[isds]/xdim_min;

  /* process each SDS according to input option */
  for (isds=0; isds<sds_cnt; isds++)
  {
    res = in_res*rf[isds];
    strcpy(in_sds_info.name, sds_names[isds]);
    fprintf(stderr, "\tProcessing SDS %s\n", in_sds_info.name);
    in_sds_info.sds_id = -1;
    if (get_sds_info((char *)NULL, &in_sds_info) == -1)
    {
      fprintf(stderr, "Ignoring invalid input sds %s\n", in_sds_info.name);
      if (in_sds_info.sds_id == FAIL)
        SDendaccess(in_sds_info.sds_id);
      continue;
    }
    if (strcmp(arg_list[2], "-sub") == 0)
      reduce_an_sds_by_sub(&in_sds_info, out_sd_id, res);
    else if (strcmp(arg_list[2], "-cnt") == 0)
      reduce_an_sds_by_cnt(&in_sds_info, out_sd_id, res, arg_list[7]);
    else if (strcmp(arg_list[2], "-cl") == 0)
      reduce_an_sds_by_class(&in_sds_info, out_sd_id, res);
    else if (strcmp(arg_list[2], "-avg") == 0)  
    {
      out_flag[0] = 1;
      for (i=3; i<7; i++)
        out_flag[i-2] = (arg_list[i][0] == 'y') ? 1 : 0;
      reduce_an_sds_by_avg(&in_sds_info, out_sd_id, res, out_flag, arg_list[8][0]);
    }
    SDendaccess(in_sds_info.sds_id);
  }

  /* copy metadata from input to output */
  if (arg_list[9][0] == 'y') 
    copy_metadata(in_sds_info.sd_id, out_sd_id);

  /* close the output HDF file */
  SDend(in_sds_info.sd_id);
}

void reduce_an_sds_by_sub(sds_t *in_sds_info, int32 out_sd_id, int res)
/****************************************************************************
!C

!Description:
  Routine to reduce an sds by sub-sampling.

!Input Parameters:
  in_sds_info: structure containing input SDS information.
  out_sd_id: ouput SDS id.
  res: reduction factor.

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)

!END
****************************************************************************/
{
  int	i, j, k, n, m; 
  int   st_c, st1_c, offset;
  int   ic, ic1, ic2, rank;
  int   bsq, res2, res_by_2; 
  int   ndata_in = 0, ndata_out = 0;
  int   ndata_sm_in = 0, ndata_sm_out;
  int	irow, nrow = 0, ncol = 0, out_ncol = 0;
  int32 attr_type, attr_cnt;
  int32 in_edge[4] = {0, 0, 0, 0};
  int32 out_edge[4] = {0, 0, 0, 0};
  int32 in_start[4] = {0, 0, 0, 0};
  int32 out_start[4] = {0, 0, 0, 0};
  sds_t out_sds_info;
  char sds_name[MAX_SDS_NAME_LEN];
  void  *attr_val, *data_in, *data_out;

  /* create output SDS with proper dimension and other attributes */
  if ((attr_val = (void *)get_sds_attr(in_sds_info->sds_id, "_FillValue", &attr_type, &attr_cnt)) == NULL) 
    fprintf(stderr, "Attribute "_FillValue" not defined for output\n");
  rank = in_sds_info->rank;
  out_sds_info.sds_id = -1;
  out_sds_info.sd_id = out_sd_id;
  out_sds_info.data_size = in_sds_info->data_size;
  out_sds_info.data_type = in_sds_info->data_type;
  bsq = ((rank == 2) || (in_sds_info->dim_size[0] < in_sds_info->dim_size[rank-1])) ? 1 : 0; 
  if (sd_charpos(in_sds_info->name, '.', 0) != -1)
  {
    out_sds_info.rank = 2;
    if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
    {
      out_sds_info.dim_size[0] = (int)ceil((float)in_sds_info->dim_size[0]/res);
      out_sds_info.dim_size[1] = (int)ceil((float)in_sds_info->dim_size[1]/res);
    }
    else
    {
      out_sds_info.dim_size[0] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
      out_sds_info.dim_size[1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
    }
  }                                           
  else
  {
    out_sds_info.rank = in_sds_info->rank;
    for (i=0; i<rank; i++)
      out_sds_info.dim_size[i] = in_sds_info->dim_size[i];
    if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
    {
      out_sds_info.dim_size[0] = (int)ceil((float)in_sds_info->dim_size[0]/res);
      out_sds_info.dim_size[1] = (int)ceil((float)in_sds_info->dim_size[1]/res);
    }
    else
    {
      out_sds_info.dim_size[rank-2] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
      out_sds_info.dim_size[rank-1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
    }
  }
  sprintf(out_sds_info.name, "%s_sub", in_sds_info->name);
  if (open_sds((char *)NULL, &out_sds_info, 'W') == -1)
  {
    if (attr_val != NULL) free(attr_val);
    return;
  }
  else if (attr_val != NULL)
  {
    if (SDsetattr(out_sds_info.sds_id, "_FillValue", attr_type, 1, (VOIDP)attr_val) == FAIL)
      fprintf(stderr, "Cannot write sds attribute _FillValue in reduce_sds\n");
    free(attr_val);
  }

  /* get input SDS name and layer options */
  /* get the number of observations in the input and output */

  get_sdsname_dim(in_sds_info->name, sds_name, &n, &m);
  if (rank == 2)
  {
    nrow = in_sds_info->dim_size[0];
    ncol = ndata_in = in_sds_info->dim_size[1];
    out_ncol = ndata_out = out_sds_info.dim_size[1];
  }                                                                            
  else if (rank > 2)
  {
    if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1])
    {
      nrow = in_sds_info->dim_size[0];
      ncol = ndata_in = in_sds_info->dim_size[1];
      out_ncol = out_sds_info.dim_size[1];
      for (i=2; i<rank; i++)
        ndata_in *= in_sds_info->dim_size[i];
      ndata_out = (int)ceil((float)in_sds_info->dim_size[1]/res);
      if ((n == -1) && (m == -1))
        for (i=2; i<rank; i++)
          ndata_out *= in_sds_info->dim_size[i];
    }
    else
    {
      nrow = in_sds_info->dim_size[rank-2];
      ncol = ndata_in = in_sds_info->dim_size[rank-1];
      out_ncol = ((n == -1)&&(m == -1)) ? out_sds_info.dim_size[rank-1] : out_sds_info.dim_size[1];
      for (i=0; i<rank-2; i++)
        ndata_in *= in_sds_info->dim_size[i];
      ndata_out = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
      if ((n == -1) && (m == -1))
        for (i=0; i<rank-2; i++)
          ndata_out *= in_sds_info->dim_size[i];
    }
  }

  if ((data_in = (void *)calloc(ndata_in, in_sds_info->data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_in in reduce_an_sds_by_sub\n");
  if ((data_out = (void *)calloc(ndata_out, out_sds_info.data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_out in reduce_an_sds_by_sub\n");
  if ((data_in != NULL) && (data_out != NULL))
  {
    if (rank == 2)
    {
      in_edge[0] = out_edge[0] = 1;
      in_edge[1] = in_sds_info->dim_size[1];
      out_edge[1] = out_sds_info.dim_size[1];
    }                       
    else if (rank > 2)
    {
      for (i=0; i<rank; i++)
        in_edge[i] = out_edge[i] = in_sds_info->dim_size[i];
      if ((n == -1) && (m == -1))
      {
        if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1])
        {
          in_edge[0] = out_edge[0] = 1;
          out_edge[1] = out_sds_info.dim_size[1];
        }
        else
        {
          in_edge[rank-2] = 1; out_edge[rank-2] = 1;
          out_edge[rank-1] = out_sds_info.dim_size[rank-1];
        }
      }
      else
      {
        out_edge[0] = 1;
        out_edge[1] = out_sds_info.dim_size[1];
        if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1])
          in_edge[0] = 1;
        else
          in_edge[rank-2] = 1;
      }
    }                                

    /* compute start of the requested SDS layer and offset of each
       observations */

    compute_sds_start_offset(in_sds_info, n, m, &st_c, &offset);   

    /* compute the number of output data samples per line */

    if (bsq == 1)
    {
      ndata_sm_out = 1;
      if ((n == -1) && (m == -1))
      {
        for (i=0; i<rank-2; i++)
	  ndata_sm_out *= out_sds_info.dim_size[i];
      }
    }
    else
    {
      ndata_sm_in = 1;
      for (i=2; i<rank; i++)
	ndata_sm_in *= in_sds_info->dim_size[i]; 
      ndata_sm_out = ((n == -1) && (m == -1)) ? ndata_sm_in : 1;
    }
   
    /* process the requested layer of the input SDS */
 
    res_by_2 = res/2;
    if ((rank>2) && (n == -1) && (m == -1) && (in_sds_info->dim_size[rank-1]>in_sds_info->dim_size[0]))
      out_start[rank-2] = 0;
    else
      out_start[0] = 0;

    /* process one line per res lines of input */
    /* process only one sample per res samples of input samples */

    for (irow=0; irow<nrow; irow += res) 
    {
      res2 = irow + res_by_2;
      if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank- 1]))
        in_start[0] = (res2 < nrow) ? res2 : irow+(nrow-irow)/2;
      else in_start[rank-2] = (res2 < nrow) ? res2 : irow+(nrow-irow)/2;

      if (SDreaddata(in_sds_info->sds_id, in_start, NULL, in_edge, data_in) == FAIL)
      {
	fprintf(stderr, "Error reading input HDF file in reduce_sds\n");
	break;
      }
      if (bsq == 1)
      {
	st1_c = st_c;
	for (k=0, ic=0; k<ndata_sm_out; k++, st1_c += ncol)
        {
          for (i=0, ic1=0; i<out_ncol; i++, ic1 += res, ic++)
          {
	    res2 = ic1 + res_by_2;
	    j = (res2 < ncol) ? st1_c + res2 : st1_c + ic1 + (ncol - ic1)/2;
	    switch(in_sds_info->data_type)
	    {
	      case  5: ((float32 *)data_out)[ic] = ((float32 *)data_in)[j]; break;
	      case  6: ((float64 *)data_out)[ic] = ((float64 *)data_in)[j]; break;
	      case 20: ((int8 *)data_out)[ic] = ((int8 *)data_in)[j]; break;
	      case 21: ((uint8 *)data_out)[ic] = ((uint8 *)data_in)[j]; break;
	      case 22: ((int16 *)data_out)[ic] = ((int16 *)data_in)[j]; break;
	      case 23: ((uint16 *)data_out)[ic] = ((uint16 *)data_in)[j]; break;
	      case 24: ((int32 *)data_out)[ic] = ((int32 *)data_in)[j]; break;
	      case 25: ((uint32 *)data_out)[ic] = ((uint32 *)data_in)[j]; break;
	    }
	  } /* for (i=0; . . .) */
	} /* for (k=0; . . .) */
      }
      else
      {
        for (i=0, ic1=0, ic=0; i<out_ncol; i++, ic1 += res)
        {
	  res2 = ic1 + res_by_2;
	  ic2 = (res2 < ncol) ? res2 : ic1 + (ncol - ic1)/2;
	  j = st_c + ic2*ndata_sm_in; 
	  for (k=0; k<ndata_sm_out; k++, j++, ic++)
	  {
            switch(in_sds_info->data_type)
            {
	      case  5: ((float32 *)data_out)[ic] = ((float32 *)data_in)[j]; break;
	      case  6: ((float64 *)data_out)[ic] = ((float64 *)data_in)[j]; break;
              case 20: ((int8 *)data_out)[ic] = ((int8 *)data_in)[j]; break;
              case 21: ((uint8 *)data_out)[ic] = ((uint8 *)data_in)[j]; break;
              case 22: ((int16 *)data_out)[ic] = ((int16 *)data_in)[j]; break;
              case 23: ((uint16 *)data_out)[ic] = ((uint16 *)data_in)[j]; break;
              case 24: ((int32 *)data_out)[ic] = ((int32 *)data_in)[j]; break;
              case 25: ((uint32 *)data_out)[ic] = ((uint32 *)data_in)[j]; break;
            }
	  } /* for (k=0; . . ) */
        } /* for (i=0; . . .) */
      }
      if (SDwritedata(out_sds_info.sds_id, out_start, NULL, out_edge, data_out) == FAIL)
      {
  	fprintf(stderr, "Error writing output to the HDF file in reduce_sds\n");
	break;
      }
      if ((rank>2) && (n == -1) && (m == -1) && (in_sds_info->dim_size[rank-1]>in_sds_info->dim_size[0]))
        ++out_start[rank-2];
      else ++out_start[0];
    } /* for (irow=0; . . )	*/
  }
  SDendaccess(out_sds_info.sds_id);
  if (data_in != NULL) free(data_in);
  if (data_out != NULL) free(data_out);
}

void reduce_an_sds_by_avg(sds_t *in_sds_info, int32 out_sd_id, int res, int *out_flag, char avg_dt)
/****************************************************************************
!C

!Description:
  Routine to reduce an sds by local averaging.

!Input Parameters:
  in_sds_info: structure containing input SDS information.
  out_sd_id: ouput SDS id.
  res: reduction factor.
  out_flag: optional output flag (min, max, avg, std_dev, num)
	    (0: no, 1: yes for each output)
  avg_dt: output data type for average data out.
	    (0: integer data type, 1: float data type)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)

!END
****************************************************************************/
{ 
  int status = 1;
  int n, m, k, bsq;
  int n_res, m_res;
  int ic, ic1, ic2;
  int rank, out_rank;
  int out_dim_size[4];
  int c = 0, i, j, ii, jj;
  int ndata_in, ndata_out;
  int fill_val, npix_in_avg;
  int ndata_sm_in = 0, ndata_sm_out;
  int st_c, st_c1, sh_sm, offset;
  int irow, icol, jcol, ncol, out_nrow, out_ncol;
  int32 attr_type1 = 0;
  int32 in_edge[4] = {0, 0, 0, 0};
  int32 out_edge[4] = {0, 0, 0, 0};
  int32 in_start[4] = {0, 0, 0, 0};
  int32 out_start[4] = {0, 0, 0, 0};
  sds_t out_sds_info[5];
  double sum2;
  float *sig_buf = NULL, avg, value = 0.0, sum, min, max;
  void *avg_buf = NULL, *num_buf = NULL, *min_buf = NULL, *max_buf = NULL;
  void *data_in;
  char sds_name[MAX_SDS_NAME_LEN];

  /* set the output SDS attributes */
  /* set the number of input and output samples per line */

  rank = in_sds_info->rank;
  fill_val = in_sds_info->fill_val;
  get_sdsname_dim(in_sds_info->name, sds_name, &n, &m);
  out_rank = ((n == -1) && (m == -1)) ? rank : 2;
  if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
  {
    out_nrow = out_dim_size[0] = (int)ceil((float)in_sds_info->dim_size[0]/res);
    out_ncol = out_dim_size[1] = (int)ceil((float)in_sds_info->dim_size[1]/res);
    for (i=2; i<out_rank; i++)
      out_dim_size[i] = in_sds_info->dim_size[i];
    ncol = ndata_in = in_sds_info->dim_size[1];
    for (i=2; i<rank; i++) 
      ndata_in *= in_sds_info->dim_size[i];
    ndata_out = out_dim_size[1];
    for (i=2; i<out_rank; i++) 
      ndata_out *= out_dim_size[i];
  }
  else
  {
    if (out_rank == 2)
    {
      out_nrow = out_dim_size[0] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
      out_ncol = out_dim_size[1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
      ndata_out = out_dim_size[1];
    }
    else
    {
      out_nrow = out_dim_size[rank-2] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
      out_ncol = out_dim_size[rank-1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
      for (i=0; i<rank-2; i++) 
        out_dim_size[i] = in_sds_info->dim_size[i];
      ndata_out = out_dim_size[rank-1];
      for (i=0; i<out_rank-2; i++) 
        ndata_out *= out_dim_size[i];
    }
    ncol = ndata_in = in_sds_info->dim_size[rank-1];
    for (i=0; i<rank-2; i++) 
      ndata_in *= in_sds_info->dim_size[i];
  }
  ndata_in *= res;

  /* set the SDS data types for optional output SDS */

  for (i=0; i<5; i++)
    if (out_flag[i] == 1)
    {
      out_sds_info[i].sds_id = -1;
      out_sds_info[i].sd_id = out_sd_id;
      out_sds_info[i].rank = out_rank;
      for (j=0; j<out_rank; j++)
        out_sds_info[i].dim_size[j] = out_dim_size[j];
      switch(i)
      {
	case 0: sprintf(out_sds_info[i].name, "%s_avg", in_sds_info->name);
	        out_sds_info[i].data_type = (avg_dt == 'y') ? DFNT_FLOAT32 : in_sds_info->data_type;
	        out_sds_info[i].data_size = DFKNTsize(out_sds_info[i].data_type);
                if ((avg_buf = (void *)calloc(ndata_out, out_sds_info[i].data_size)) == NULL)
		{
      		  fprintf(stderr, "Cannot allocate memory for avg_buf in reduce_an_sds_by_avg\n");
		  status = -1;
		}
	    	break;
	case 1: sprintf(out_sds_info[i].name, "%s_min", in_sds_info->name);
	        out_sds_info[i].data_type = in_sds_info->data_type;
	        out_sds_info[i].data_size = DFKNTsize(out_sds_info[i].data_type);
                if ((min_buf = (void *)calloc(ndata_out, out_sds_info[i].data_size)) == NULL)
		{
      		  fprintf(stderr, "Cannot allocate memory for min_buf in reduce_an_sds_by_avg\n");
		  status = -1;
		}
		break;
	case 2: sprintf(out_sds_info[i].name, "%s_max", in_sds_info->name);
	        out_sds_info[i].data_type = in_sds_info->data_type;
	        out_sds_info[i].data_size = DFKNTsize(out_sds_info[i].data_type);
                if ((max_buf = (void *)calloc(ndata_out, out_sds_info[i].data_size)) == NULL)
		{
      		  fprintf(stderr, "Cannot allocate memory for max_buf in reduce_an_sds_by_avg\n");
		  status = -1;
		}
		break;
	case 3: sprintf(out_sds_info[i].name, "%s_sig", in_sds_info->name);
	        out_sds_info[i].data_type = DFNT_FLOAT32;
	        out_sds_info[i].data_size = DFKNTsize(out_sds_info[i].data_type);
                if ((sig_buf = (void *)calloc(ndata_out, out_sds_info[i].data_size)) == NULL)
		{
      		  fprintf(stderr, "Cannot allocate memory for sig_buf in reduce_an_sds_by_avg\n");
		  status = -1;
		}
		break;
	case 4: sprintf(out_sds_info[i].name, "%s_num", in_sds_info->name);
	        out_sds_info[i].data_type = (res*res < 256) ? DFNT_INT8 : DFNT_INT16;
	        out_sds_info[i].data_size = DFKNTsize(out_sds_info[i].data_type);
                if ((num_buf = (void *)calloc(ndata_out, out_sds_info[i].data_size)) == NULL)
		{
      		  fprintf(stderr, "Cannot allocate memory for num_buf in reduce_an_sds_by_avg\n");
		  status = -1;
		}
		break;
      }
      if (open_sds((char *)NULL, &out_sds_info[i], 'W') == -1) status = -1;
    } /* if (out_flag . . . 	*/

  /* compute start and offset for the selected layer of the input SDS */

  compute_sds_start_offset(in_sds_info, n, m, &st_c, &offset);   

  /* compute the number of output samples per line */

  bsq = ((rank == 2) || (in_sds_info->dim_size[0] < in_sds_info->dim_size[rank-1])) ? 1 : 0;
  if (bsq == 1)
  {
    ndata_sm_out = 1;
    if ((n == -1) && (m == -1))
    {
      for (i=0; i<rank-2; i++)
        ndata_sm_out *= out_sds_info[0].dim_size[i];
    }
    else
    {
      if (m == -1)
	st_c = n*res*ncol;
      else
        st_c = res*ncol*in_sds_info->dim_size[0]*in_sds_info->dim_size[1];
    }
  }
  else
  {
    ndata_sm_in = 1;
    for (i=2; i<rank; i++)
      ndata_sm_in *= in_sds_info->dim_size[i]; 
    ndata_sm_out = ((n == -1) && (m == -1)) ? ndata_sm_in : 1;
  }

  if ((data_in = (void *)calloc(ndata_in, in_sds_info->data_size)) == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for data_in in reduce_an_sds_by_avg\n");
    status = -1;
  }

  if ((data_in != NULL) && (status != -1))
  {
    if (rank == 2)
    {
      in_edge[0] = res; in_edge[1] = in_sds_info->dim_size[1];
      out_edge[0] = 1; out_edge[1] = out_sds_info[0].dim_size[1];
    }
    else if (rank > 2)
    {
      for (i=0; i<rank; i++) in_edge[i] = in_sds_info->dim_size[i];
      for (i=0; i<out_rank; i++) out_edge[i] = out_sds_info[0].dim_size[i];
      if ((n == -1) && (m == -1))
      {
        if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1])
	{
	  in_edge[0] = res; out_edge[0] = 1;
	}
	else 
	{
	  in_edge[rank-2] = res; out_edge[rank-2] = 1;
	}
      }
      else
      {
        out_edge[0] = 1;
        if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]) in_edge[0] = res;
        else in_edge[rank-2] = res;
      }
    }

    /* reads res input lines and average res x res samples */
    /* output one line per res input lines. */
    /* if required create min, max, std, and npix SDS */

    for (irow=0; irow<out_nrow; irow++)
    {
      if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
      {
        in_start[0] = irow*res;
	if ((irow+1)*res <= in_sds_info->dim_size[0]) in_edge[0] = res;
	else in_edge[0] = in_sds_info->dim_size[0] - in_start[0];
	n_res = in_edge[0];
      }
      else
      {
        in_start[rank-2] = irow*res;
        if ((irow+1)*res <= in_sds_info->dim_size[rank-2]) in_edge[rank-2] = res;
	else in_edge[rank-2] = in_sds_info->dim_size[rank-2] - in_start[rank-2];
	n_res = in_edge[rank-2];
      }
      if ((rank>2) && (n == -1) && (m == -1) && (in_sds_info->dim_size[rank-1]>in_sds_info->dim_size[0]))
        out_start[rank-2] = irow;
      else out_start[0] = irow;

      if (SDreaddata(in_sds_info->sds_id, in_start, NULL, in_edge, data_in) == FAIL)
      {
  	fprintf(stderr, "Error reading data line from SDS %s in reduce_an_sds_by_avg\n", in_sds_info->name);
        break;
      }
      else 
      {
	if (bsq == 1)
	{
          if ((n == -1) && (m == -1)) offset = ncol;
          else
          {
	    if (m == -1) offset = ncol*n_res;
	    else offset = ncol*(in_sds_info->dim_size[0]*in_sds_info->dim_size[1])*res;
          }
	  sh_sm = n_res*ncol;
	  st_c1 = st_c;
	  for (k=0, ic=0; k<ndata_sm_out; k++, st_c1 += sh_sm)
	  {
	    for (ic1=0, ic2=0; ic1<out_ncol; ic1++, ic2 += res, ic++)
	    {
	      sum = 0.0;
	      sum2 = 0.0;
	      min = 120000;
	      max = -120000;
	      npix_in_avg = 0;
	      m_res = (ic2 < ncol) ? res : (ncol - ic2 + res);
	      for (ii=0, icol=ic2; ii<n_res; ii++, icol += ncol)
	      {
	        jcol = st_c1 + icol;
	        for (jj=0; jj<m_res; jj++, jcol++)
	        {
                  switch (in_sds_info->data_type)
                  {
                    case  5: value = (float)(((float32 *)data_in)[jcol]); break;
                    case  6: value = (float)(((float64 *)data_in)[jcol]); break;
                    case 20: value = ((int8 *)data_in)[jcol]; break;
                    case 21: value = ((uint8 *)data_in)[jcol]; break;
                    case 22: value = ((int16 *)data_in)[jcol]; break;
                    case 23: value = ((uint16 *)data_in)[jcol]; break;
                    case 24: value = (float)(((int32 *)data_in)[jcol]); break;
                    case 25: value = (float)(((uint32 *)data_in)[jcol]); break;
                  }      
                  if (value != fill_val)
                  {
	            npix_in_avg++;
		    sum += value;
                    if ((out_flag[1] == 1) && (min > value)) min = value;
                    if ((out_flag[2] == 1) && (max < value)) max = value;
                    if (out_flag[3] == 1) sum2 += value*value;
                  }
	        }  /* for (jj=0; . . . ) */
	      }  /* for (ii=0; . . .) */

	      avg = (npix_in_avg == 0) ? (float32)fill_val : sum/npix_in_avg;
	      if (avg_dt == 'y')
		((float32 *)avg_buf)[ic] = (float32)avg;
	      else
	      {
                switch (in_sds_info->data_type)
                {
                  case  5: ((float32 *)avg_buf)[ic] = (float32)avg; break;
                  case  6: ((float64 *)avg_buf)[ic] = (float64)avg; break;
                  case 20: ((int8 *)avg_buf)[ic] = (int8)(avg + 0.5); break;
                  case 21: ((uint8 *)avg_buf)[ic] = (uint8)(avg + 0.5); break;
                  case 22: ((int16 *)avg_buf)[ic] = (int16)(avg + 0.5); break;
                  case 23: ((uint16 *)avg_buf)[ic] = (uint16)(avg + 0.5); break;
                  case 24: ((int32 *)avg_buf)[ic] = (int32)(avg + 0.5); break;
                  case 25: ((uint32 *)avg_buf)[ic] = (uint32)(avg + 0.5); break;
                }
	      }
              if (out_flag[1] == 1)
              {
	        if (npix_in_avg == 0) min = (float)fill_val;
                switch (in_sds_info->data_type)
                {
                  case  5: ((float32 *)min_buf)[ic] = (float32)min; break;
                  case  6: ((float64 *)min_buf)[ic] = (float64)min; break;
                  case 20: ((int8 *)min_buf)[ic] = (int8)min; break;
                  case 21: ((uint8 *)min_buf)[ic] = (uint8)min; break;
                  case 22: ((int16 *)min_buf)[ic] = (int16)min; break;
                  case 23: ((uint16 *)min_buf)[ic] = (uint16)min; break;
                  case 24: ((int32 *)min_buf)[ic] = (int32)min; break;
                  case 25: ((uint32 *)min_buf)[ic] = (uint32)min; break;
                }
              }
              if (out_flag[2] == 1)
              {
	        if (npix_in_avg == 0) max = (float)fill_val;
                switch (in_sds_info->data_type)
                {
                  case  5: ((float32 *)max_buf)[ic] = (float32)max; break;
                  case  6: ((float64 *)max_buf)[ic] = (float64)max; break;
                  case 20: ((int8 *)max_buf)[ic] = (int8)max; break;
                  case 21: ((uint8 *)max_buf)[ic] = (uint8)max; break;
                  case 22: ((int16 *)max_buf)[ic] = (int16)max; break;
                  case 23: ((uint16 *)max_buf)[ic] = (uint16)max; break;
                  case 24: ((int32 *)max_buf)[ic] = (int32)max; break;
                  case 25: ((uint32 *)max_buf)[ic] = (uint32)max; break;
                }
              }
              if (out_flag[3] == 1)
		sig_buf[ic] = (npix_in_avg == 0) ? (float32)fill_val : (float32)sqrt(sum2/npix_in_avg - avg*avg);
              if (out_flag[4] == 1)
              {
                if (out_sds_info[4].data_type == DFNT_INT8)
                  ((int8 *)num_buf)[ic] = (int8)npix_in_avg;
                else if (out_sds_info[4].data_type == DFNT_INT16)
                  ((int16 *)num_buf)[ic] = (int16)npix_in_avg;
              }
	    }  /* for (ic1=0; . . .) */
	  }  /* for (k=0; . . .) */
	}  /* if (bsq == 1)   */
	else
	{
	  sh_sm = ncol*ndata_sm_in;
	  for (ic1=0, ic2=0, ic=0; ic1<out_ncol; ic1++, ic2 += res)
	  {
	    st_c1 = st_c + ic1*res*ndata_sm_in;
	    m_res = (ic2 < ncol) ? res : (ncol - ic2 + res);
	    for (k=0; k<ndata_sm_out; k++, ic++, st_c1++)
	    {
	      sum = 0.0;
	      sum2 = 0.0;
	      min = 120000;
	      max = -120000;
	      npix_in_avg = 0;
	      for (ii=0; ii<n_res; ii++)
	      {
		icol = st_c1 + ii*sh_sm;
	        for (jj=0, jcol=icol; jj<m_res; jj++, jcol += ndata_sm_in)
	        {
                  switch (in_sds_info->data_type)
                  {
                    case  5: value = ((float32 *)data_in)[jcol]; break;
                    case  6: value = (float)(((float64 *)data_in)[jcol]); break;
                    case 20: value = ((int8 *)data_in)[jcol]; break;
                    case 21: value = ((uint8 *)data_in)[jcol]; break;
                    case 22: value = ((int16 *)data_in)[jcol]; break;
                    case 23: value = ((uint16 *)data_in)[jcol]; break;
                    case 24: value = (float)(((int32 *)data_in)[jcol]); break;
                    case 25: value = (float)(((uint32 *)data_in)[jcol]); break;
                  }      
                  if (value != fill_val)
                  {
	            npix_in_avg++;
                    sum += value;
                    if ((out_flag[1] == 1) && (min > value)) min = value;
                    if ((out_flag[2] == 1) && (max < value)) max = value;
	            if (out_flag[3] == 1) sum2 += value*value;
                  }
	        }  /* for (jj=0; . . . ) */
	      }  /* for (ii=0; . . .) */

	      avg = (npix_in_avg == 0) ? (float32)fill_val : sum/npix_in_avg;
	      if (avg_dt == 'y')
		((float32 *)avg_buf)[ic] = (float32)avg;
	      else
	      {
                switch (in_sds_info->data_type)
                {
                  case  5: ((float32 *)avg_buf)[ic] = (float32)avg; break;
                  case  6: ((float64 *)avg_buf)[ic] = (float64)avg; break;
                  case 20: ((int8 *)avg_buf)[ic] = (int8)(avg + 0.5); break;
                  case 21: ((uint8 *)avg_buf)[ic] = (uint8)(avg + 0.5); break;
                  case 22: ((int16 *)avg_buf)[ic] = (int16)(avg + 0.5); break;
                  case 23: ((uint16 *)avg_buf)[ic] = (uint16)(avg + 0.5); break;
                  case 24: ((int32 *)avg_buf)[ic] = (int32)(avg + 0.5); break;
                  case 25: ((uint32 *)avg_buf)[ic] = (uint32)(avg + 0.5); break;
                }
	      }
              if (out_flag[1] == 1)
              {
	        if (npix_in_avg == 0) min = (float)fill_val;
                switch (in_sds_info->data_type)
                {
                  case  5: ((float32 *)min_buf)[ic] = (float32)min; break;
                  case  6: ((float64 *)min_buf)[ic] = (float64)min; break;
                  case 20: ((int8 *)min_buf)[ic] = (int8)min; break;
                  case 21: ((uint8 *)min_buf)[ic] = (uint8)min; break;
                  case 22: ((int16 *)min_buf)[ic] = (int16)min; break;
                  case 23: ((uint16 *)min_buf)[ic] = (uint16)min; break;
                  case 24: ((int32 *)min_buf)[ic] = (int32)min; break;
                  case 25: ((uint32 *)min_buf)[ic] = (uint32)min; break;
                }
              }
              if (out_flag[2] == 1)
              {
	        if (npix_in_avg == 0) max = (float)fill_val;
                switch (in_sds_info->data_type)
                {
                  case  5: ((float32 *)max_buf)[ic] = (float32)max; break;
                  case  6: ((float64 *)max_buf)[ic] = (float64)max; break;
                  case 20: ((int8 *)max_buf)[ic] = (int8)max; break;
                  case 21: ((uint8 *)max_buf)[ic] = (uint8)max; break;
                  case 22: ((int16 *)max_buf)[ic] = (int16)max; break;
                  case 23: ((uint16 *)max_buf)[ic] = (uint16)max; break;
                  case 24: ((int32 *)max_buf)[ic] = (int32)max; break;
                  case 25: ((uint32 *)max_buf)[ic] = (uint32)max; break;
                }
              }
              if (out_flag[3] == 1)
		sig_buf[ic] = (npix_in_avg == 0) ? (float32)fill_val : (float32)sqrt(sum2/npix_in_avg - avg*avg);
              if (out_flag[4] == 1)
              {
                if (out_sds_info[4].data_type == DFNT_INT8)
                  ((int8 *)num_buf)[ic] = (int8)npix_in_avg;
                else if (out_sds_info[4].data_type == DFNT_INT16)
                  ((int16 *)num_buf)[ic] = (int16)npix_in_avg;
              }
	    }  /* for (k=0; . . .) */
	  }  /* for (ic1=0; . . .) */
	} /* else  . . .  */
        if (SDwritedata(out_sds_info[0].sds_id, out_start, NULL, out_edge, (VOIDP)avg_buf) == FAIL)
	  fprintf(stderr, "Error writing line of SDS to HDF file in reduce_an_sds_by_avg\n");
        if (out_flag[1] == 1)
        {
          if (SDwritedata(out_sds_info[1].sds_id, out_start, NULL, out_edge, (VOIDP)min_buf) == FAIL)
	    fprintf(stderr, "Error writing line of SDS to HDF file in reduce_an_sds_by_avg\n");
        }
        if (out_flag[2] == 1)
        {
          if (SDwritedata(out_sds_info[2].sds_id, out_start, NULL, out_edge, (VOIDP)max_buf) == FAIL)
	    fprintf(stderr, "Error writing line of SDS to HDF file in reduce_an_sds_by_avg\n");
        }
        if (out_flag[3] == 1)
        {
          if (SDwritedata(out_sds_info[3].sds_id, out_start, NULL, out_edge, (VOIDP)sig_buf) == FAIL)
	    fprintf(stderr, "Error writing line of SDS to HDF file in reduce_an_sds_by_avg\n");
        }
        if (out_flag[4] == 1)
        {
          if (SDwritedata(out_sds_info[4].sds_id, out_start, NULL, out_edge, (VOIDP)num_buf) == FAIL)
	    fprintf(stderr, "Error writing line of SDS to HDF file in reduce_an_sds_by_avg\n");
        }
      }  /* else		*/
    }  /* for (ir=0; . . . 	*/ 
    for (i=0; i<5; i++)
      if (out_flag[i] == 1)
      {
	switch(i)
	{
      	  case 0: attr_type1 = (avg_dt == 'y') ? DFNT_FLOAT32 : in_sds_info->data_type; 
		  c = 1; break;
	  case 1: attr_type1 = in_sds_info->data_type; c = 1; break;
	  case 2: attr_type1 = in_sds_info->data_type; c = 1; break;
	  case 3: attr_type1 = DFNT_FLOAT32; c = 1; break;
	  case 4: attr_type1 = DFNT_INT8; c = 0; break;
	}
	write_attr_fval(out_sds_info[i].sds_id, attr_type1, c, fill_val, ATTR_FILL_NAME);
      }
  }  /* if ((data_in != NULL) . . . 	*/
  SDendaccess(in_sds_info->sds_id);
  SDendaccess(out_sds_info[0].sds_id);
  free(data_in);
  free(avg_buf);
  if (out_flag[1] == 1) { free(min_buf); SDendaccess(out_sds_info[1].sds_id); }
  if (out_flag[2] == 1) { free(max_buf); SDendaccess(out_sds_info[2].sds_id); }
  if (out_flag[3] == 1) { free(sig_buf); SDendaccess(out_sds_info[3].sds_id); }
  if (out_flag[4] == 1) { free(num_buf); SDendaccess(out_sds_info[4].sds_id); }
}

void reduce_an_sds_by_cnt(sds_t *in_sds_info, int32 out_sd_id, int res, char *cmp_str)
/****************************************************************************
!C

!Description:
  Routine to reduce an sds by outputing the count of the number of pixels
  that satisfy input data value constraint.

!Input Parameters:
  in_sds_info: structure containing input SDS information.
  out_sd_id: ouput SDS id.
  res: reduction factor.
  cmp_str: string containing the data value constraint.

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Note(s): (none)

!End
******************************************************************************/
{ 
  int **bn_arr;
  int n, m, bsq;
  int status = 1;
  int st_c, st_c1;
  int ic, ic1, ic2;
  int n_res, m_res;
  int sh_sm, offset;
  int rank, out_rank;
  int out_dim_size[4];
  int i, j, k, ii, jj;
  int cnt[MAX_NUM_SDS];
  int isds, nsds, opr_id;
  int out_nrow, out_ncol;
  int ndata_in, ndata_out;
  int bn_val[MAX_NUM_BITS];
  int bn_cnt[MAX_NUM_BITS];
  int bn_opr[MAX_NUM_BITS];
  int ncol, irow, icol, jcol;
  int value = 0, fill_val = 0, sds_val;
  int ndata_sm_in = 0, ndata_sm_out;
  int32 attr_type, attr_cnt;
  int32 in_edge[4] = {0, 0, 0, 0};
  int32 out_edge[4] = {0, 0, 0, 0};
  int32 in_start[4] = {0, 0, 0, 0};
  int32 out_start[4] = {0, 0, 0, 0};
  sds_t *out_sds_info = NULL;
  char sds_name[MAX_SDS_NAME_LEN];
  void *attr_val, *data_in, *data_out[MAX_NUM_SDS];
  char *opr[] = {"<=", ">=", "==", "!=", "<", ">"};

  if ((bn_arr = (int **)Calloc2D(MAX_NUM_BITS, 16, sizeof(int))) == NULL)
    fprintf(stderr, "Cannot allocate memory for bn_arr in reduce_sds_by_cnt\n");

  /* get the number of ouptut SDS based on the bit constraint */
  /* also get the bit number array and associated value */

  if ((nsds = get_bit_opt(cmp_str, bn_arr, bn_val, bn_opr, bn_cnt)) == 0)
    fprintf(stderr, "-bn option value %s is incorrect or empty\n", cmp_str);
  if (nsds != 0)
    if ((out_sds_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
      fprintf(stderr, "Cannot allocate memory for out_sds_info in reduce_sds_by_cnt\n");

  /* set output SDS attributes */
  /* compute number of input and output samples per line */

  if ((out_sds_info != NULL) && (bn_arr != NULL) && (nsds != 0))
  {
    if ((attr_val = get_sds_attr(in_sds_info->sds_id, "_FillValue", &attr_type, &attr_cnt)) == NULL)
      fprintf(stderr, "Attribute fill value not defined for output\n");
    else
      switch(attr_type)
      {
        case 20: fill_val = ((int8 *)attr_val)[0]; break;
        case 21: fill_val = ((uint8 *)attr_val)[0]; break;
        case 22: fill_val = ((int16 *)attr_val)[0]; break;
        case 23: fill_val = ((uint16 *)attr_val)[0]; break;
        case 24: fill_val = ((int32 *)attr_val)[0]; break;
        case 25: fill_val = ((uint32 *)attr_val)[0]; break;
      }

    rank = in_sds_info->rank;
    get_sdsname_dim(in_sds_info->name, sds_name, &n, &m);
    out_rank = ((n == -1) && (m == -1)) ? rank : 2;
    if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
    {
      out_nrow = out_dim_size[0] = (int)ceil((float)in_sds_info->dim_size[0]/res);
      out_ncol = out_dim_size[1] = (int)ceil((float)in_sds_info->dim_size[1]/res);
      for (i=2; i<out_rank; i++)
        out_dim_size[i] = in_sds_info->dim_size[i];
      ncol = ndata_in = in_sds_info->dim_size[1];
      for (i=2; i<rank; i++)
        ndata_in *= in_sds_info->dim_size[i];
      ndata_out = out_dim_size[1];
      for (i=2; i<out_rank; i++)
        ndata_out *= out_dim_size[i];
    }
    else
    {
      if (out_rank == 2)
      {
        out_nrow = out_dim_size[0] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
        out_ncol = out_dim_size[1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
        ndata_out = out_dim_size[1];
      }
      else
      {
        out_nrow = out_dim_size[rank-2] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
        out_ncol = out_dim_size[rank-1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
        for (i=0; i<rank-2; i++)
          out_dim_size[i] = in_sds_info->dim_size[i];
        ndata_out = out_dim_size[rank-1];
        for (i=0; i<out_rank-2; i++)
          ndata_out *= out_dim_size[i];
      }
      ncol = ndata_in = in_sds_info->dim_size[rank-1];
      for (i=0; i<rank-2; i++)
        ndata_in *= in_sds_info->dim_size[i];
    }
    ndata_in *= res;

    /* set the output SDS name and open the output SDS */

    for (isds=0; isds<nsds; isds++)
    {
      out_sds_info[isds].sds_id = -1;
      out_sds_info[isds].sd_id = out_sd_id;
      out_sds_info[isds].rank = out_rank;
      for (j=0; j<out_rank; j++)
        out_sds_info[isds].dim_size[j] = out_dim_size[j];
      opr_id = bn_opr[isds];
      if (bn_cnt[isds] == 0)
	sprintf(out_sds_info[isds].name, "%s_cnt_sds%s%d", in_sds_info->name, opr[opr_id], 
		bn_val[isds]);
      else if (bn_cnt[isds] == 1) 
	sprintf(out_sds_info[isds].name, "%s_cnt_bit%d%s%d", in_sds_info->name, 
		bn_arr[isds][0], opr[opr_id], bn_val[isds]);
      else
	sprintf(out_sds_info[isds].name, "%s_cnt_bits%d-%d%s%d", in_sds_info->name, 
		bn_arr[isds][0], bn_arr[isds][0]+bn_cnt[isds]-1, opr[opr_id], bn_val[isds]);
      out_sds_info[isds].data_type = (res <= 16) ? DFNT_UINT8 : DFNT_UINT16;
      out_sds_info[isds].data_size = DFKNTsize(out_sds_info[isds].data_type);
      if (open_sds((char *)NULL, &out_sds_info[isds], 'W') != -1)
      {
        if (res <= 16) write_attr_fval(out_sds_info[isds].sds_id, DFNT_UINT8, 0, fill_val, ATTR_FILL_NAME);
        else write_attr_fval(out_sds_info[isds].sds_id, DFNT_UINT16, 0, fill_val, ATTR_FILL_NAME);
      }
    } /* for (isds=0; . . . ) */

    /* compute the start and offset of requested layer */

    compute_sds_start_offset(in_sds_info, n, m, &st_c, &offset);
    bsq = ((rank == 2) || (in_sds_info->dim_size[0] < in_sds_info->dim_size[rank-1])) ? 1 : 0;
    if (bsq == 1)
    {
      ndata_sm_out = 1;
      if ((n == -1) && (m == -1))
      {
        for (i=0; i<rank-2; i++)
          ndata_sm_out *= out_sds_info[0].dim_size[i];
      }
      else
      {
        if (m == -1)
          st_c = n*res*ncol;
        else
          st_c = res*ncol*in_sds_info->dim_size[0]*in_sds_info->dim_size[1];
      }
    }
    else
    {
      ndata_sm_in = 1;
      for (i=2; i<rank; i++)
        ndata_sm_in *= in_sds_info->dim_size[i];
      ndata_sm_out = ((n == -1) && (m == -1)) ? ndata_sm_in : 1;
    }

    if ((data_in = (void *)calloc(ndata_in, in_sds_info->data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_in in reduce_an_sds_by_cnt\n");

    for (isds=0; isds<nsds; isds++)
      if ((data_out[isds] = (void *)calloc(ndata_out, out_sds_info[isds].data_size)) == NULL)
      {
        fprintf(stderr, "Cannot allocate memory for data_out in reduce_an_sds_by_cnt\n");
	status = -1;
      }

    if ((data_in != NULL) && (status != -1))
    {
      if (rank == 2)
      {
        in_edge[0] = res; in_edge[1] = in_sds_info->dim_size[1];
        out_edge[0] = 1; out_edge[1] = out_sds_info[0].dim_size[1];
      }
      else if (rank > 2)
      {
        for (i=0; i<rank; i++) in_edge[i] = in_sds_info->dim_size[i];
        for (i=0; i<out_rank; i++) out_edge[i] = out_sds_info[0].dim_size[i];
        if ((n == -1) && (m == -1))
        {
          if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1])
          {
            in_edge[0] = res; out_edge[0] = 1;
          }
          else
          {
            in_edge[rank-2] = res; out_edge[rank-2] = 1;
          }
        }
        else
        {
          out_edge[0] = 1;
          if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]) in_edge[0] = res;
          else in_edge[rank-2] = res;
        }
      }

      /* read res lines from input and output one sample per res x res input */
      /* output one line per res line of input data */

      for (irow=0; irow<out_nrow; irow++)
      {
        if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
        {
          in_start[0] = irow*res;
          if ((irow+1)*res <= in_sds_info->dim_size[0]) in_edge[0] = res;
          else in_edge[0] = in_sds_info->dim_size[0] - in_start[0];
          n_res = in_edge[0];
        }
        else
        {
          in_start[rank-2] = irow*res;
          if ((irow+1)*res <= in_sds_info->dim_size[rank-2]) in_edge[rank-2] = res;
          else in_edge[rank-2] = in_sds_info->dim_size[rank-2] - in_start[rank-2];
          n_res = in_edge[rank-2];
        }
        if ((rank>2) && (n == -1) && (m == -1) && (in_sds_info->dim_size[rank-1]>in_sds_info->dim_size[0]))
          out_start[rank-2] = irow;
        else out_start[0] = irow;
  
        if (SDreaddata(in_sds_info->sds_id, in_start, NULL, in_edge, data_in) == FAIL)
        {
          fprintf(stderr, "Error reading data line from SDS %s in reduce_an_sds_by_avg\n", in_sds_info->name);
          break;
        }
	else
	{
          if (bsq == 1)
          {
            if ((n == -1) && (m == -1))
              offset = ncol;
            else
            {
              if (m == -1)
                offset = ncol*n_res;
              else
                offset = ncol*(in_sds_info->dim_size[0]*in_sds_info->dim_size[1])*res;
            }
            sh_sm = n_res*ncol;
            st_c1 = st_c;
            for (k=0, ic=0; k<ndata_sm_out; k++, st_c1 += sh_sm)
            {
              for (ic1=0, ic2=0; ic1<out_ncol; ic1++, ic2 += res, ic++)
              {
		for (isds=0; isds<nsds; isds++)
		  cnt[isds] = 0;
                m_res = (ic2 < ncol) ? res : (ncol - ic2 + res);
                for (ii=0, icol=ic2; ii<n_res; ii++, icol += ncol)
	        {
                  jcol = st_c1 + icol;
                  for (jj=0; jj<m_res; jj++, jcol++)
		  {
		    switch (in_sds_info->data_type)
                    {
                      case 20: value = ((int8 *)data_in)[jcol]; break;
                      case 21: value = ((uint8 *)data_in)[jcol]; break;
                      case 22: value = ((int16 *)data_in)[jcol]; break;
                      case 23: value = ((uint16 *)data_in)[jcol]; break;
                      case 24: value = ((int32 *)data_in)[jcol]; break;
                      case 25: value = ((uint32 *)data_in)[jcol]; break;
		    }	
		    if (value != fill_val)
		      for (isds=0; isds<nsds; isds++)
		      {
			sds_val = value;
		        if (bn_cnt[isds] != 0)
		        {
			  if (bn_arr[isds][0] > 0) 
			    sds_val = sds_val >> bn_arr[isds][0];
		          sds_val &= BIT[bn_cnt[isds]-1];
		        }
		        switch(bn_opr[isds])
		        {
		          case 0: if (sds_val <= bn_val[isds]) cnt[isds]++; break;
		          case 1: if (sds_val >= bn_val[isds]) cnt[isds]++; break;
		          case 2: if (sds_val == bn_val[isds]) cnt[isds]++; break;
		          case 3: if (sds_val != bn_val[isds]) cnt[isds]++; break;
		          case 4: if (sds_val < bn_val[isds]) cnt[isds]++; break;
		          case 5: if (sds_val > bn_val[isds]) cnt[isds]++; break;
		        }
		      }  /* for (isds=0; . . . ) */
		  }  /* for (jj=0; . . .) */
		}  /* for (ii=0; . . . ) */
		for (isds=0; isds<nsds; isds++)
		{
	          if (res <= 16) ((uint8 *)data_out[isds])[ic] = (uint8)cnt[isds];
	          else ((uint16 *)data_out[isds])[ic] = (uint16)cnt[isds];
		}
	      }  /* for (ic1=0; . . .)  */
	    }  /* for (k=0; . . . . ) 	*/
	  } 
	  else
	  {
            sh_sm = ncol*ndata_sm_in;
            for (ic1=0, ic2=0, ic=0; ic1<out_ncol; ic1++, ic2 += res)
            {
              st_c1 = st_c + ic1*res*ndata_sm_in;
              m_res = (ic2 < ncol) ? res : (ncol - ic2 + res);
              for (k=0; k<ndata_sm_out; k++, ic++, st_c1++)
              {
		for (isds=0; isds<nsds; isds++)
		  cnt[isds] = 0;
                for (ii=0; ii<n_res; ii++)
                {
                  icol = st_c1 + ii*sh_sm;
                  for (jj=0, jcol=icol; jj<m_res; jj++, jcol += ndata_sm_in)
                  {
                    switch (in_sds_info->data_type)
                    {
                      case 20: value = ((int8 *)data_in)[jcol]; break;
                      case 21: value = ((uint8 *)data_in)[jcol]; break;
                      case 22: value = ((int16 *)data_in)[jcol]; break;
                      case 23: value = ((uint16 *)data_in)[jcol]; break;
                      case 24: value = ((int32 *)data_in)[jcol]; break;
                      case 25: value = ((uint32 *)data_in)[jcol]; break;
                    }
                    if (value != fill_val)
                      for (isds=0; isds<nsds; isds++)
                      {
                        sds_val = value;
                        if (bn_cnt[isds] != 0)
                        {
                          if (bn_arr[isds][0] > 0)
                            sds_val = sds_val >> bn_arr[isds][0];
                          sds_val &= BIT[bn_cnt[isds]-1];
                        }
                        switch(bn_opr[isds])
                        {
                          case 0: if (sds_val <= bn_val[isds]) cnt[isds]++; break;
                          case 1: if (sds_val >= bn_val[isds]) cnt[isds]++; break;
                          case 2: if (sds_val == bn_val[isds]) cnt[isds]++; break;
                          case 3: if (sds_val != bn_val[isds]) cnt[isds]++; break;
                          case 4: if (sds_val < bn_val[isds]) cnt[isds]++; break;
                          case 5: if (sds_val > bn_val[isds]) cnt[isds]++; break;
                        }
                      } /* for (isds=0; . . . ) */
                  }  /* for (jj=0; . . . ) */
                }  /* for (ii=0; . . .) */
		for (isds=0; isds<nsds; isds++)
		{
	          if (res <= 16) ((uint8 *)data_out[isds])[ic] = (uint8)cnt[isds];
	          else ((uint16 *)data_out[isds])[ic] = (uint16)cnt[isds];
		}
	      }  /* for (k=0; . . . ) */
	    }  /* for (ic1=0; . .  ) */
	  } /* else . . . .  */
	  for (isds=0; isds<nsds; isds++)
            if (SDwritedata(out_sds_info[isds].sds_id, out_start, NULL, out_edge, data_out[isds]) == FAIL)
	      fprintf(stderr, "Error writing line of SDS to file in reduce_an_sds_by_cnt\n");
        }  /* else		*/
      }  /* for (ir=0; . . . 	*/ 
      free(data_in);
      for (isds=0; isds<nsds; isds++)
        free(data_out[isds]);
    }  /* if ((data_in != NULL) . . . 	*/
    free(attr_val);
    SDendaccess(in_sds_info->sds_id);
    for (isds=0; isds<nsds; isds++)
      if (out_sds_info[isds].sds_id != -1) SDendaccess(out_sds_info[isds].sds_id);
  }  
  free(out_sds_info);
  return;
}

int get_bit_opt(char *cmp_str, int **bn_arr, int *bn_val, int *bn_opr, int *bn_cnt)
/****************************************************************************
!C

!Description:
  Obtain the bit array and the corresponding bit data value for use in 
  reduce an sds by cnt option.

!Input Parameters:
  cmp_str: input bit number and bit value option string 

!Input/Output Parameters: (none)

!Output Parameters: (none)
  bn_arr: 2-d array containing the start and end bit numbers of
	  each bit range array.
  bn_val: bit range value for each of the range in bn_arr
  bn_opr: bit relational operation for each bn_arr and bn_val
  bn_cnt: number of bit relational constraints 

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Note(s): (none)

!End
******************************************************************************/
{
  int p1, p2, p3 = 0;
  int i, j, k, j1; 
  int len, len1, cnt;
  int ibit, st_bit, end_bit;
  char val_str[10], num_str[10], bn_str[20];
  char *opr[] = {"<=", ">=", "==", "!=", "<", ">"};

  len = (int)strlen(cmp_str);
  p1 = p2 = k = 0;
  while (p2 != -1)
  {
    p2 = sd_charpos(cmp_str, ',', p1);
    if (p2 != -1) sd_strmid(cmp_str, p1, p2-p1, bn_str);
    else sd_strmid(cmp_str, p1, len-p1, bn_str);
    for (i=0; i<6; i++)
      if ((p3 = sd_strpos(bn_str, opr[i], 0)) != -1) break;
    if (i >= 6)
      fprintf(stderr, "Ignoring invalid value %s in -bit option \n", bn_str);
    else
    {
      bn_opr[k] = i;
      cnt = 0;
      if (p3 != 0)
      {
        for (j=0; j<p3; j++) 
          if (bn_str[j] == '-') break; 
	  else num_str[j] = bn_str[j];
        num_str[j] = '\0';
	st_bit = (int)atoi(num_str);
        if (j >= p3) end_bit = st_bit;
	else
        {
	  j++;
	  for (j1=j; j1<p3; j1++)
	    num_str[j1-j] = bn_str[j1];
	  num_str[j1-j] = '\0';
	  end_bit = (int)atoi(num_str);
        }
	for (ibit=st_bit; ibit<=end_bit; ibit++, cnt++)
	  bn_arr[k][cnt] = ibit;
      }
      bn_cnt[k] = cnt;
      len1 = (int)strlen(bn_str);
      p3 = (i > 3) ? (p3 + 1) : (p3 + 2);
      sd_strmid(bn_str, p3, len1-p3, val_str);
      bn_val[k] = (int)atoi(val_str);
      k++;
    }
    p1 = p2 + 1;
  }
  return k;
}

void reduce_an_sds_by_class(sds_t *in_sds_info, int32 out_sd_id, int res)
/****************************************************************************
!C

!Description:
  Routine to reduce an sds by classifying the local area to one of the set
  of classes in the input.

!Input Parameters:
  in_sds_info: structure containing input SDS information.
  out_sd_id: ouput SDS id.
  res: reduction factor.

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original October 2002. Version 1.0

!Team-unique Header:
  See file prologue.

!Developer(s):

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Note(s): (none)

!End
******************************************************************************/
{ 
  int n, m, bsq;
  int st_c, st_c1;
  int ic, ic1, ic2;
  int n_res, m_res;
  int sh_sm, offset;
  int rank, out_rank;
  int out_dim_size[4];
  int i, j, k, ii, jj;
  int value = 0, fill_val = 0;
  int out_nrow, out_ncol;
  int ndata_in, ndata_out;
  int ncol, irow, icol, jcol;
  int ndata_sm_in = 0, ndata_sm_out;
  int32 attr_type, attr_cnt;
  int32 in_edge[4] = {0, 0, 0, 0};
  int32 out_edge[4] = {0, 0, 0, 0};
  int32 in_start[4] = {0, 0, 0, 0};
  int32 out_start[4] = {0, 0, 0, 0};
  sds_t out_sds_info;
  char sds_name[MAX_SDS_NAME_LEN];
  void *attr_val, *data_in, *data_out;
  int cl_num, class_cnt[MAX_NUM_CLASS], max_cl_id;

  if ((attr_val = get_sds_attr(in_sds_info->sds_id, "_FillValue", &attr_type, &attr_cnt)) == NULL)
    fprintf(stderr, "Attribute fill value not defined for output\n");
  else
    switch(attr_type)
    {
      case 20: fill_val = ((int8 *)attr_val)[0]; break;
      case 21: fill_val = ((uint8 *)attr_val)[0]; break;
      case 22: fill_val = ((int16 *)attr_val)[0]; break;
      case 23: fill_val = ((uint16 *)attr_val)[0]; break;
      case 24: fill_val = ((int32 *)attr_val)[0]; break;
      case 25: fill_val = ((uint32 *)attr_val)[0]; break;
    }

  /* set the output SDS attributes */
  /* compute number of input and output samples per output line */

  rank = in_sds_info->rank;
  get_sdsname_dim(in_sds_info->name, sds_name, &n, &m);
  out_rank = ((n == -1) && (m == -1)) ? rank : 2;
  if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
  {
    out_nrow = out_dim_size[0] = (int)ceil((float)in_sds_info->dim_size[0]/res);
    out_ncol = out_dim_size[1] = (int)ceil((float)in_sds_info->dim_size[1]/res);
    for (i=2; i<out_rank; i++)
      out_dim_size[i] = in_sds_info->dim_size[i];
    ncol = ndata_in = in_sds_info->dim_size[1];
    for (i=2; i<rank; i++)
      ndata_in *= in_sds_info->dim_size[i];
    ndata_out = out_dim_size[1];
    for (i=2; i<out_rank; i++)
      ndata_out *= out_dim_size[i];
  }
  else
  {
    if (out_rank == 2)
    {
      out_nrow = out_dim_size[0] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
      out_ncol = out_dim_size[1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
      ndata_out = out_dim_size[1];
    }
    else
    {
      out_nrow = out_dim_size[rank-2] = (int)ceil((float)in_sds_info->dim_size[rank-2]/res);
      out_ncol = out_dim_size[rank-1] = (int)ceil((float)in_sds_info->dim_size[rank-1]/res);
      for (i=0; i<rank-2; i++)
        out_dim_size[i] = in_sds_info->dim_size[i];
      ndata_out = out_dim_size[rank-1];
      for (i=0; i<out_rank-2; i++)
        ndata_out *= out_dim_size[i];
    }
    ncol = ndata_in = in_sds_info->dim_size[rank-1];
    for (i=0; i<rank-2; i++)
      ndata_in *= in_sds_info->dim_size[i];
  }
  ndata_in *= res;

  /* open output SDS */

  out_sds_info.sds_id = -1;
  out_sds_info.sd_id = out_sd_id;
  out_sds_info.rank = out_rank;
  for (j=0; j<out_rank; j++)
    out_sds_info.dim_size[j] = out_dim_size[j];
  sprintf(out_sds_info.name, "%s", in_sds_info->name);
  out_sds_info.data_type = in_sds_info->data_type; 
  out_sds_info.data_size = DFKNTsize(out_sds_info.data_type);
  if (open_sds((char *)NULL, &out_sds_info, 'W') != -1)
    write_attr_fval(out_sds_info.sds_id, in_sds_info->data_type, 1, fill_val, ATTR_FILL_NAME);

  /* compute start and offset of requested layer */

  compute_sds_start_offset(in_sds_info, n, m, &st_c, &offset);
  bsq = ((rank == 2) || (in_sds_info->dim_size[0] < in_sds_info->dim_size[rank-1])) ? 1 : 0;
  if (bsq == 1)
  {
    ndata_sm_out = 1;
    if ((n == -1) && (m == -1))
    {
      for (i=0; i<rank-2; i++)
        ndata_sm_out *= out_sds_info.dim_size[i];
    }
    else
    {
      if (m == -1) st_c = n*res*ncol;
      else st_c = res*ncol*in_sds_info->dim_size[0]*in_sds_info->dim_size[1];
    }
  }
  else
  {
    ndata_sm_in = 1;
    for (i=2; i<rank; i++)
      ndata_sm_in *= in_sds_info->dim_size[i];
    ndata_sm_out = ((n == -1) && (m == -1)) ? ndata_sm_in : 1;
  }

  if ((data_in = (void *)calloc(ndata_in, in_sds_info->data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_in in reduce_an_sds_by_cnt\n");

  if ((data_out = (void *)calloc(ndata_out, out_sds_info.data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_out in reduce_an_sds_by_cnt\n");

  if ((data_in != NULL) && (data_out != NULL))
  {
    if (rank == 2)
    {
      in_edge[0] = res; in_edge[1] = in_sds_info->dim_size[1];
      out_edge[0] = 1; out_edge[1] = out_sds_info.dim_size[1];
    }
    else if (rank > 2)
    {
      for (i=0; i<rank; i++) in_edge[i] = in_sds_info->dim_size[i];
      for (i=0; i<out_rank; i++) out_edge[i] = out_sds_info.dim_size[i];
      if ((n == -1) && (m == -1))
      {
        if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1])
        {
          in_edge[0] = res; out_edge[0] = 1;
        }
        else
        {
          in_edge[rank-2] = res; out_edge[rank-2] = 1;
        }
      }
      else
      {
        out_edge[0] = 1;
        if (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]) in_edge[0] = res;
        else in_edge[rank-2] = res;
      }
    }

    /* read res lines of input at a time. Process res x res */
    /* block and ouput one line per res line of input */

    for (irow=0; irow<out_nrow; irow++)
    {
      if ((rank == 2) || (in_sds_info->dim_size[0] > in_sds_info->dim_size[rank-1]))
      {
        in_start[0] = irow*res;
        if ((irow+1)*res <= in_sds_info->dim_size[0]) in_edge[0] = res;
        else in_edge[0] = in_sds_info->dim_size[0] - in_start[0];
        n_res = in_edge[0];
      }
      else
      {
        in_start[rank-2] = irow*res;
        if ((irow+1)*res <= in_sds_info->dim_size[rank-2]) in_edge[rank-2] = res;
        else in_edge[rank-2] = in_sds_info->dim_size[rank-2] - in_start[rank-2];
        n_res = in_edge[rank-2];
      }
      if ((rank>2) && (n == -1) && (m == -1) && (in_sds_info->dim_size[rank-1]>in_sds_info->dim_size[0]))
        out_start[rank-2] = irow;
      else out_start[0] = irow;
  
      if (SDreaddata(in_sds_info->sds_id, in_start, NULL, in_edge, data_in) == FAIL)
      {
        fprintf(stderr, "Error reading data line from SDS %s in reduce_an_sds_by_avg\n", in_sds_info->name);
        break;
      }
      else
      {
        if (bsq == 1)
        {
          if ((n == -1) && (m == -1)) offset = ncol;
          else
          {
            if (m == -1) offset = ncol*n_res;
            else offset = ncol*(in_sds_info->dim_size[0]*in_sds_info->dim_size[1])*res;
          }
          sh_sm = n_res*ncol;
          st_c1 = st_c;
          for (k=0, ic=0; k<ndata_sm_out; k++, st_c1 += sh_sm)
          {
            for (ic1=0, ic2=0; ic1<out_ncol; ic1++, ic2 += res, ic++)
            {
	      for (cl_num=0; cl_num<MAX_NUM_CLASS; cl_num++)
		class_cnt[cl_num] = 0;
              m_res = (ic2 < ncol) ? res : (ncol - ic2 + res);
              for (ii=0, icol=ic2; ii<n_res; ii++, icol += ncol)
	      {
                jcol = st_c1 + icol;
                for (jj=0; jj<m_res; jj++, jcol++)
		{
		  switch (in_sds_info->data_type)
                  {
                    case 20: value = ((int8 *)data_in)[jcol]; break;
                    case 21: value = ((uint8 *)data_in)[jcol]; break;
                    case 22: value = ((int16 *)data_in)[jcol]; break;
                    case 23: value = ((uint16 *)data_in)[jcol]; break;
                    case 24: value = ((int32 *)data_in)[jcol]; break;
                    case 25: value = ((uint32 *)data_in)[jcol]; break;
		  }	
		  if ((value != fill_val) && (value < MAX_NUM_CLASS))
		    ++class_cnt[value];
		}  /* for (jj=0; . . . ) */
	      }  /* for (ii=0; . . . ) */
	      for (cl_num=0, max_cl_id=0; cl_num<MAX_NUM_CLASS; cl_num++)
	        if (class_cnt[cl_num] > class_cnt[max_cl_id]) 
		  max_cl_id = cl_num;
	      if (class_cnt[max_cl_id] == 0) max_cl_id = fill_val;
	      switch (out_sds_info.data_type)
              {
                case 20: ((int8 *)data_out)[ic] = (int8)max_cl_id; break;
                case 21: ((uint8 *)data_out)[ic] = (uint8)max_cl_id; break;
                case 22: ((int16 *)data_out)[ic] = (int16)max_cl_id; break;
                case 23: ((uint16 *)data_out)[ic] = (uint16)max_cl_id; break;
                case 24: ((int32 *)data_out)[ic] = max_cl_id; break;
                case 25: ((uint32 *)data_out)[ic] = max_cl_id; break;
	      }	
	    }  /* for (ic1=0; . . .)  */
	  }  /* for (k=0; . . . . ) 	*/
	} 
	else
	{
          sh_sm = ncol*ndata_sm_in;
          for (ic1=0, ic2=0, ic=0; ic1<out_ncol; ic1++, ic2 += res)
          {
            st_c1 = st_c + ic1*res*ndata_sm_in;
            m_res = (ic2 < ncol) ? res : (ncol - ic2 + res);
            for (k=0; k<ndata_sm_out; k++, ic++, st_c1++)
            {
	      for (cl_num=0; cl_num<MAX_NUM_CLASS; cl_num++)
		class_cnt[cl_num] = 0;
              for (ii=0; ii<n_res; ii++)
              {
                icol = st_c1 + ii*sh_sm;
                for (jj=0, jcol=icol; jj<m_res; jj++, jcol += ndata_sm_in)
                {
                  switch (in_sds_info->data_type)
                  {
                    case 20: value = ((int8 *)data_in)[jcol]; break;
                    case 21: value = ((uint8 *)data_in)[jcol]; break;
                    case 22: value = ((int16 *)data_in)[jcol]; break;
                    case 23: value = ((uint16 *)data_in)[jcol]; break;
                    case 24: value = ((int32 *)data_in)[jcol]; break;
                    case 25: value = ((uint32 *)data_in)[jcol]; break;
                  }
                  if ((value != fill_val) && (value < MAX_NUM_CLASS))
		    ++class_cnt[value];
                }  /* for (jj=0; . . . ) */
              }  /* for (ii=0; . . .) */
              for (cl_num=0, max_cl_id=0; cl_num<MAX_NUM_CLASS; cl_num++)
                if (class_cnt[cl_num] > class_cnt[max_cl_id])
                  max_cl_id = cl_num;
	      if (class_cnt[max_cl_id] == 0) max_cl_id = fill_val;
              switch (out_sds_info.data_type)
              {
                case 20: ((int8 *)data_out)[ic] = (int8)max_cl_id; break;
                case 21: ((uint8 *)data_out)[ic] = (uint8)max_cl_id; break;
                case 22: ((int16 *)data_out)[ic] = (int16)max_cl_id; break;
                case 23: ((uint16 *)data_out)[ic] = (uint16)max_cl_id; break;
                case 24: ((int32 *)data_out)[ic] = max_cl_id; break;
                case 25: ((uint32 *)data_out)[ic] = max_cl_id; break;
              }
	    }  /* for (k=0; . . . ) */
	  }  /* for (ic1=0; . .  ) */
	} /* else . . . .  */
        if (SDwritedata(out_sds_info.sds_id, out_start, NULL, out_edge, data_out) == FAIL)
	  fprintf(stderr, "Error writing line of SDS to file in reduce_an_sds_by_class\n");
      }  /* else		*/
    }  /* for (ir=0; . . . 	*/ 
    free(data_in);
    free(data_out);
  }  /* if ((data_in != NULL) . . . 	*/
  free(attr_val);
  SDendaccess(in_sds_info->sds_id);
  SDendaccess(out_sds_info.sds_id);
}  
