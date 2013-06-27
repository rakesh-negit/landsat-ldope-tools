/****************************************************************************
!C

!File: create_sds_ts_stat.c

!Description:
  Compute pixelwise statistics of SDS values in a set of input MODISLand 
  HDF-EOS files. 

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  February 2004. Version 1.0

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

    Yi (Sophyi) Zhang 
    LDOPE                             Science Systems and Applications Inc.
    sophyi@ltpmail.gsfc.nasa.gov      NASA/GSFC Code 922
    phone: 301-614-5497               Greenbelt, MD 20771    

    David Roy 
    LDOPE                             University of Maryland, 
                                      Department of Geography
    droy@kratmos.gsfc.nasa.gov        NASA/GSFC Code 922 (B32)
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
" create_sds_ts_stat - Create a summary statistic file containing one or\n" \
"    more output 2D SDS that describe the mean, standard deviation,\n" \
"    minimum, maximum, sum, and number of observations, computed on\n" \
"    pixelwise basis from a time series of input MODIS Land HDF-EOS data\n"\
"    products.\n" \
" \n" \
"SYNOPSIS \n" \
"    create_sds_ts_stat [-help] [filename] \n" \
" \n" \
"    create_sds_ts_stat\n"\
"             -sds=<sds_name,sds_minval,sds_maxval,f_nop_in,f_nop_out,dt>\n"\
"             -of=<filename>\n"\
"             -param=[avg][,std][,min][,max][,npix][,sum]\n"\
"             f1 [f2 f3. . fn] \n" \
" \n" \
"DESCRIPTION \n" \
"    Compute statistics of SDS values at each pixel from a set of input\n" \
"    MODIS Land HDF-EOS data product. For example, the tool can be used to\n" \
"    compute the time series statistics of SDS values when the input files\n"\
"    are a time series sequence of files for the same geographic area i.e.,\n"\
"    files that have the same tile id. \n" \
" \n" \
"    The output HDF file contains SDS of the following statistics: \n" \
"    mean, standard deviation, minimum, maximum, number of \n" \
"    pixels and sum of SDS values. By default all parameters are output. \n" \
" \n" \
"    The user can specify the range of the SDS values to be considered for \n" \
"    computing statistics and a fill value to be excluded from consideration.\n"\
"    Two or more SDSs are processed by repeating the option -sds option with\n"\
"    the SDS names. An input file is ignored if the user specified SDS is \n" \
"    missing in that file. \n" \
"\n" \
"OPTIONS \n" \
"    -help                Print this help message, If the input filename \n" \
"                         is specified with this option, then the names \n" \
"                         of all the SDS in the file are displayed. \n" \
"    -sds=<SDSname,sds_minval,sds_maxval,f_nop_in,f_nop_out,dt> \n" \
"                         Input SDS name, minimum and maximum of SDS values\n"\
"                         to be considered, no operation fill value in the\n"\
"                         input SDS, no operation fill value in the output\n"\
"                         SDS and output data type.\n" \
"\n" \
"                         Input SDS pixel values equal to 'f_nop_in' are not\n"\
"                         considered. If no values are available for\n"\
"                         computing the statistics at a particular pixel\n"\
"                         then 'f_nop_out' is written to the pixel \n" \
"                         in the output SDS. \n" \
"\n" \
"                         Each of the parameters except the SDS name can be\n" \
"                         replaced by * to assign default values defined as:\n"\
"                             sds_min,sds_max: valid range of input SDS. \n" \
"                             f_nop_in: input SDS fill value \n" \
"                             f_nop_out: f_nop_in \n" \
"                             dt: input SDS data type.\n" \
"                         Valid values for output data type 'dt' are \n" \
"                         INT8, UINT8, INT16, UINT16, INT32, UINT32, \n" \
"                         FLOAT32 \n" \
"\n" \
"                         Output data type of SDS number of pixels  is \n" \
"                         always INT16. \n" \
"\n" \
"                         To process a specific layer of a 3D SDS specify \n" \
"                         the element number of the third dimension as a \n" \
"                         dot extension of the SDS name: sds_name.n (e.g., \n" \
"                         sur_refl_b02.1 = the layer defined by the 1st \n" \
"                         element of the 3rd dimension of the 3D SDS \n" \
"                         sur_refl_b02). \n" \
" \n" \
"                         To process a specific layer of a 4D SDS, specify \n" \
"                         the higher dimension element number(s) as a dot \n" \
"                         extension of the SDS name: sds_name.n.m (e.g., \n" \
"                         Surface_Refl.1.2 = the layer defined by the \n" \
"                         1st element of the 3rd dimension and the 2nd \n" \
"                         element of the 4th dimension of the 4D SDS \n" \
"                         Surface_Refl). \n" \
"\n" \
"                         Note that wildcards and ranges of element values \n" \
"                         may be specified as sds_name.* and as\n" \
"                         sds_name.n1-n2.m respectively. \n" \
"    -param=avg,std,min,max,npix,sum \n" \
"                         Output statistical parameters include average, \n" \
"                         standard deviation, minimum, maximum, number of \n" \
"                         pixels, and sum of SDS values at a pixel. One \n" \
"                         or more of these parameter may be specified in any\n"\
"                         order. If this argument is unspecified then all\n" \
"                         these parameters are output.  \n" \
"    -of=<filename>       Output filename \n" \
" \n" \
"Examples: \n" \
"    create_sds_ts_stat -sds=\"Day_Tile_Snow_Cover,*,*,*,*,*\" \n" \
"           -param=min \n" \
"           MOD10A1.A2002194.h30v11.003.2002199095914.hdf \n" \
"           MYD10A1.A2002194.h30v11.003.2002199015759.hdf \n" \
"           -of=test_ts_stat_MOD10A1.A2002194.h30v11.003.2002199095914.hdf\n"\
" \n" \
"    create_sds_ts_stat -sds=\"Day_Tile_Snow_Cover,*,*,*,*,*\" \n" \
"           -sds=\"Snow_Spatial_QA,*,*,*,*,*\" -param=max \n" \
"           MOD10A1.A2002194.h30v11.003.2002199095914.hdf \n" \
"           MYD10A1.A2002194.h30v11.003.2002199015759.hdf \n" \
"           -of=test_2_ts_stat_MOD10A1.A2002194.h30v11.003.2002199095914.hdf\n"\
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 04/05/04 \n" \

#define USAGE \
"usage: \n" \
"    create_sds_ts_stat [-help] [filename] \n" \
" \n" \
"    create_sds_ts_stat\n"\
"             -sds=<sds_name,sds_min,sds_max,f_nop_in,f_nop_out,dt>\n"\
"             -of=<filename> -param=[avg][,std][,min][,max][,npix][,sum] \n" \
"             f1 f2 f3. . fn \n" \
" \n" \
"OPTIONS \n" \
"    -help                Print this help message, If the input filename \n" \
"                         is specified with this option, then the names \n" \
"                         of all the SDS in the file are displayed. \n" \
"    -sds=<SDSname,sds_minval,sds_maxval,f_nop_in,f_nop_out,dt> \n" \
"                         Input SDS name, minimum and maximum of SDS values\n" \
"                         to be considered, no operation fill value in the\n"\
"                         input SDS, no operation fill value in the output\n"\
"                         SDS and output data type.\n" \
"\n" \
"                         Input SDS pixel values equal to 'f_nop_in' are not\n"\
"                         considered. If no values are available for\n"\
"                         computing the statistics at a particular pixel\n"\
"                         then 'f_nop_out' is written to the pixel in the\n"\
"                         output SDS. \n" \
" \n"\
"                         Each of the parameters except the SDS name can be\n" \
"                         replaced by * to assign default values defined as:\n"\
"                             sds_min,sds_max: valid range of input SDS. \n" \
"                             f_nop_in: input SDS fill value \n" \
"                             f_nop_out: f_nop_in \n" \
"                             dt: input SDS data type.\n" \
"                         Valid values for output data type 'dt' are \n" \
"                         INT8, UINT8, INT16, UINT16, INT32, UINT32, \n" \
"                         FLOAT32 \n" \
" \n" \
"                         Output data type of SDS number of pixels  is \n" \
"                         always INT16. \n" \
" \n" \
"                         To process a specific layer of a 3D SDS specify \n" \
"                         the element number of the third dimension as a \n" \
"                         dot extension of the SDS name: sds_name.n (e.g., \n" \
"                         sur_refl_b02.1 = the layer defined by the 1st \n" \
"                         element of the 3rd dimension of the 3D SDS \n" \
"                         sur_refl_b02). \n" \
" \n" \
"                         To process a specific layer of a 4D SDS, specify \n" \
"                         the higher dimension element number(s) as a dot \n" \
"                         extension of the SDS name: sds_name.n.m (e.g., \n" \
"                         Surface_Refl.1.2 = the layer defined by the \n" \
"                         1st element of the 3rd dimension and the 2nd \n" \
"                         element of the 4th dimension of the 4D SDS \n" \
"                         Surface_Refl). \n" \
"\n" \
"                         Note that wildcards and ranges of element values \n" \
"                         may be specified as sds_name.* and as\n" \
"                         sds_name.n1-n2.m respectively. \n" \
"    -param=avg,std,min,max,npix,sum \n" \
"                         Output statistical parameters include average, \n" \
"                         standard deviation, minimum, maximum, number of \n" \
"                         pixels, and sum of SDS values at a pixel. One \n" \
"                         or more of these parameter may be specified in any\n"\
"                         order. If this argument is unspecified then all\n" \
"                         these parameters are output.  \n" \
"    -of=<filename>       Output filename \n" \
"\n"

#define MAX_NSDS 10

/******************************************************************************
                            Prototypes.
******************************************************************************/

int parse_cmd_create_sds_ts_stat(int argc, char **argv, char **expr, int *nsds, int *fcnt, 
			  char *out_fname, int *param_st);
int read_param(char *expr, char *sds_name, char *range1, char *range2, char *f_nop_in, 
	       char *f_nop_out, char *dt);
void comp_stat(sds_t *in_sds_info, int fcnt, int out_sd_id, char *range1, char *range2, 
	       char *f_nop_in, char *f_nop_out, char *dt, int *param_st);

int main(int argc, char **argv)
/******************************************************************************
!C

!Description:
  Main function for create_sds_ts_stat.

!Input Parameters: (none)
  command line arguments: see help for details.

!Output Parameters: (none)
  return 0 on successful completion of the process

!Revision History:
  See file prologue..

!Team-unique Header:
  See file prologue.

!References and Credits: (see file prologue)

!Design Notes: (none)

!END
********************************************************************************/
{
  char **expr;
  char range1[10], range2[10];
  char sds_name[MAX_SDS_NAME_LEN];
  char out_fname[MAX_PATH_LENGTH];
  char f_nop_in[10], f_nop_out[10], dt[10];
  int param_st[6];
  int i, fcnt, id, fid, isds, nsds, iarg, status;
  int32 out_sd_id;
  sds_t *in_sds_info;

  if (argc == 1)
    {
      status = -1;
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_FAILURE);
    }   
  
  if ((argc==2) && (strcmp(argv[1], "-help")==0))
    {
      fprintf(stderr, "%s\n", HELP);
      exit(EXIT_SUCCESS);
    }
  
  /*  Display SDS names of input HDF file */
  
  if ((argc>=3) && (strcmp(argv[1], "-help")==0))
    {
      /* Print all the SDS names in the input file */
      for (i=2; i<argc; i++)
        {
          if (argv[i][0] != '-')
            {
              display_sds_info_of_file(argv[i]);
            }
        } 
      exit(EXIT_SUCCESS);
    }

  if ((expr = (char **)Calloc2D(MAX_NUM_SDS, MAX_STR_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for expr in main()\n");
  else
  {
    if ((status = parse_cmd_create_sds_ts_stat(argc, argv, expr, &nsds, &fcnt, out_fname, param_st)) == -1)
      fprintf(stderr, "%s\n", USAGE);
    else if (status != 0)
    {
      if ((in_sds_info = (sds_t *)calloc(fcnt, sizeof(sds_t))) == NULL)
        fprintf(stderr, "Cannot allocate memory for in_sds_info in create_sds_ts_stat\n");
      else
      {
        if ((out_sd_id = SDstart(out_fname, DFACC_CREATE)) == FAIL)
	  fprintf(stderr, "Cannot create output HDF file: %s\n", out_fname);
        else
        {
          for (isds=0; isds<nsds; isds++)
          {
            if (read_param(expr[isds], sds_name, range1, range2, f_nop_in, f_nop_out, 
			   dt) == -1)
	      {
		fprintf(stderr, "Invalid argument %s for -sds option. \n", expr[isds]);
		fprintf(stderr, "Argument should be in the form of -sds=<sds_name,sds_min,sds_max,f_nop_in,f_nop_out,dt> \n");
		fprintf(stderr, "Argument is not processed\n");
		exit(EXIT_FAILURE);
	      }
	    else
            {
	      fprintf(stdout, "Processing SDS %s\n", sds_name);
	      for (iarg=1, fid=0; iarg<argc; iarg++)
	      {
	        if (argv[iarg][0] != '-') 
	        {
	          in_sds_info[fid].sd_id = in_sds_info[fid].sds_id = -1;
	          strcpy(in_sds_info[fid].name, sds_name);
	          if (get_sds_info(argv[iarg], &in_sds_info[fid]) == -1)
	            fprintf(stderr, "\tIgnoring input file %s\n", argv[iarg]);
	          else fid++;	
	        }
              } /* for (iarg=1, . . ) */
	      if (fid > 0)
	        comp_stat(in_sds_info, fid, out_sd_id, range1, range2, f_nop_in, 
			  f_nop_out, dt, param_st); 
	      else 
		{
		  fprintf(stderr, "No valid input file. \n");
		  exit(EXIT_FAILURE);
		}
            }
	    for (id=0; id<fid; id++)
	    {
	      SDendaccess(in_sds_info[fid].sds_id);
	      SDend(in_sds_info[fid].sd_id);
	    }
          } /* for (isds=0; . . .) */
          SDend(out_sd_id);
        }
        free(in_sds_info);
      }
    }
    Free2D((void **)expr);
    fprintf(stderr, "Processing done ! \n");
  }
  return 0;
}
  
int parse_cmd_create_sds_ts_stat(int argc, char **argv, char **expr, int *nsds, 
			  int *fcnt, char *out_fname, int *param_st)
/******************************************************************************
!C

!Description:
  Function to parse command line arguments.

!Input Parameters:
  argc: number of input arguments
  argv: string array containing arguments

!Output Parameters:
  expr:      String contains the -SDS option input.
  nsds:      number of SDS
  fcnt:      Number of input file.
  out_fname: output filename
  param_st:  Array contains the flag of -param option input. 
             sum   : param_st[0] 
             avg   : param_st[1]
             std   : param_st[2]
             npix  : param_st[3]
             min   : param_st[4]
             max   : param_st[5]

  return 1 if parsing is succesfull, -1 if not all required parameters input.
 
!Revision History:
    See file prologue.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)

!END
 ********************************************************************************/
{
  int i, p1, p2, len, isds, st;
  char tmp_str[10], param_str[80];

  *fcnt = *nsds = 0;
  out_fname[0] = '\0';
  param_str[0] = '\0';
  for (i=1, st=1, isds=0; i<argc; i++)
  {
    if (is_arg_id(argv[i], "-sds=") == 0)
    {
      get_arg_val(argv[i], expr[isds]);
      isds++;
    }
    else if (is_arg_id(argv[i], "-of=") == 0)
      get_arg_val(argv[i], out_fname);
    else if (is_arg_id(argv[i], "-param=") == 0)
      get_arg_val(argv[i], param_str);
    else if (argv[i][0] != '-') ++*fcnt;
    else fprintf(stderr, "Ignoring unknown option %s\n", argv[i]);
  }
  if (strlen(out_fname) <= 0) {
    st = -1; fprintf(stderr, "Missing output filename\n");
  }
  else 
  {
    if (isds == 0) {
      fprintf(stdout, "No SDS name input. \n");
      st = 0;
    }
    else 
    {
      if (param_str[0] == '\0')
      {
	fprintf(stdout, "No output parameters specified. All parameters output.\n");
	for (i=0; i<6; i++)
	  param_st[i] = 1;
      }
      else 
      {
	for (i=0; i<6; i++)
	  param_st[i] = 0;
	len = (int)strlen(param_str);
	for (i=0, p1=p2=0; i<6; i++)
	{
	  p2 = sd_charpos(param_str, ',', p1);
	  if (p2 == -1) sd_strmid(param_str, p1, len-p1, tmp_str);
	  else sd_strmid(param_str, p1, p2-p1, tmp_str);
	  if (strcmp(tmp_str, "sum") == 0) param_st[0] = 1;
	  else if (strcmp(tmp_str, "avg") == 0) param_st[1] = 1;
	  else if (strcmp(tmp_str, "std") == 0) param_st[2] = 1;
	  else if (strcmp(tmp_str, "npix") == 0) param_st[3] = 1;
	  else if (strcmp(tmp_str, "min") == 0) param_st[4] = 1;
	  else if (strcmp(tmp_str, "max") == 0) param_st[5] = 1;
	  else fprintf(stderr, "Cannot compute parameter %s\n", tmp_str);
	  if (p1 == -1) break;
	  else p1 = p2 + 1;
	}
      }
    }
    *nsds = isds;
  }
  return st;                                          
}
          
int read_param(char *expr, char *sds_name, char *range1, char *range2, 
	       char *f_nop_in, char *f_nop_out, char *dt)
/******************************************************************************
!C

!Description:
  Function to parse the input -sds option.

!Input Parameters:
  expr:      Input sds option string.

!Output Parameters:
  sds_name:  Input SDS name to be processed. SDS name should be input seperately 
             using the -sds option repeatly.
  range1:    Input SDS minimum range.
  range2:    Input SDS maximum range.
  f_nop_in:  User defined SDS attribute fill value for input..
  f_nop_out: User defined SDS attribute fill value for output.
  dt:        Input SDS data type.

  return 1 if parsing is succesfull, -1 if not all required parameters input.
 
!Revision History:
    See file prologue.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)

!END
 ********************************************************************************/
{
  int st = 1;
  int p1, p2;
  int i, nc, len;

  dt[0] = '\0';
  sds_name[0] = '\0';
  f_nop_in[0] = f_nop_out[0] = '\0';
  len = (int)strlen(expr);
  for (i=0, p1=0; i<5; i++)
  {
    p2 = sd_charpos(expr, ',', p1);
    if (p2 != -1)
    {
      nc = p2 - p1;
      if (nc > 0) {
        switch (i) {
          case 0: sd_strmid(expr, p1, nc, sds_name); break;
          case 1: sd_strmid(expr, p1, nc, range1); break;
          case 2: sd_strmid(expr, p1, nc, range2); break;
          case 3: sd_strmid(expr, p1, nc, f_nop_in); break;
          case 4: sd_strmid(expr, p1, nc, f_nop_out); break;
        }
      }
      p1 = p2 + 1;
    }
    else { st = -1; break; }
  }
  if (st == 1) {
    nc = len - p1;
    sd_strmid(expr, p1, nc, dt);
  }
  return st;                                                              
}

void comp_stat(sds_t *in_sds_info, int fcnt, int out_sd_id, char *range1, 
	       char *range2, char *f_nop_in, char *f_nop_out, char *dt, int *param_st)
/******************************************************************************
!C

!Description:
  Function comp_stat to compute pixelwise statistics of input SDS.

!Input Parameters:
  in_sds_info:  Array of input SDS information structure.
  fcnt:         Number of valid input file.
  out_sd_id:    Output file descriptor.
  range1:       Input SDS minimum range.
  range2:       Input SDS maximum range.
  f_nop_in:     User defined SDS attribute fill value for input..
  f_nop_out:    User defined SDS attribute fill value for output.
  dt:           Input SDS data type.
  param_st:     Array contains the flag of -param option input. 
                 sum   : param_st[0] 
                 avg   : param_st[1]
                 std   : param_st[2]
                 npix  : param_st[3]
                 min   : param_st[4]
                 max   : param_st[5]

  return 1 if parsing is succesfull, -1 if not all required parameters input.
 
!Revision History:
    See file prologue.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)

!END
 ********************************************************************************/

{
  int32 out_dt;
  int16 *sds_npix;
  float32 *sds_std;
  int fid, isds, nsds;
  int ic, st_c, offset;
  int in_rank, out_rank;
  sds_t out_sds_info[6];
  int ndata_in, range[2];
  int i, n, m, dim_sz[4];
  int nop_in, nop_out, npix;
  int irow, icol, nrows = 0, ncols = 0;
  int32 in_start[4], in_edge[4];
  int32 out_start[4], out_edge[4];
  double sum, sum2, avg, avg2;
  double min = 0, max = 0, std, sds_val = 0.0;
  void *sds_sum, *sds_mean;
  void **sds_data, *sds_min, *sds_max;
  char out_sds_name[MAX_SDS_NAME_LEN];
  char *sds_str[] = {"Sum", "Mean", "Std", "Npix", "Min", "Max"};

  if (strcmp(dt, "*") == 0)
    out_dt = in_sds_info[0].data_type;
  else {
    if (strcmp(dt, "FLOAT32") == 0) out_dt = 5;
    else if (strcmp(dt, "INT8") == 0) out_dt = 20;
    else if (strcmp(dt, "UINT8") == 0) out_dt = 21;
    else if (strcmp(dt, "INT16") == 0) out_dt = 22;
    else if (strcmp(dt, "UINT16") == 0) out_dt = 23;
    else if (strcmp(dt, "INT32") == 0) out_dt = 24;
    else if (strcmp(dt, "UINT32") == 0) out_dt = 25;
    else {
      fprintf(stderr, "Output data type %s not recognized. Set to default\n", dt);
      out_dt = in_sds_info[0].data_type;
    }                                                                     
  }

  if (strcmp(range1, "*") == 0)
    range[0] = (in_sds_info[0].data_type == 5) ? (int)in_sds_info[0].frange[0] : in_sds_info[0].range[0];
  else range[0] = (int)atoi(range1);
  if (strcmp(range2, "*") == 0)
    range[1] = (in_sds_info[0].data_type == 5) ? (int)in_sds_info[0].frange[1] : in_sds_info[0].range[1];
  else range[1] = (int)atoi(range2);
  if (strcmp(f_nop_in, "*") == 0)
    nop_in = (in_sds_info[0].data_type == 5) ? (int)in_sds_info[0].fill_fval : in_sds_info[0].fill_val;
  else nop_in = (int)atoi(f_nop_in);
  /* if (strcmp(f_nop_out, "*") == 0) */
  nop_out = (strcmp(f_nop_out, "*") == 0) ? nop_in : (int)atoi(f_nop_out);

  get_sds_param(&in_sds_info[0], &n, &m, &out_rank, dim_sz);
 
  nsds = 6;

  for (isds=0; isds<nsds; isds++)
    if (param_st[isds] == 1)
    {
      out_sds_info[isds].sds_id = -1;
      out_sds_info[isds].sd_id = out_sd_id;
      if (isds == 2) out_sds_info[isds].data_type = DFNT_FLOAT32;
      else if (isds == 3) out_sds_info[isds].data_type = DFNT_INT16;
      else out_sds_info[isds].data_type = out_dt;
      out_sds_info[isds].data_size = DFKNTsize(out_dt);;
      out_sds_info[isds].rank=out_rank;
      for (i=0; i<out_rank; i++)
        out_sds_info[isds].dim_size[i] = dim_sz[i];    
      sprintf(out_sds_name, "%s of %s", sds_str[isds], in_sds_info[0].name);
      strcpy(out_sds_info[isds].name, out_sds_name);
      if (open_sds((char *)NULL, &out_sds_info[isds], 'W') != -1)
      {
        if (isds == 2)
          write_attr_fval(out_sds_info[isds].sds_id, out_sds_info[isds].data_type, 
			  0, 0, ATTR_FILL_NAME);
        else
          write_attr_fval(out_sds_info[isds].sds_id, out_sds_info[isds].data_type, 
			  1, nop_out, ATTR_FILL_NAME);
      }
    }

  ndata_in = compute_sds_ndata(&in_sds_info[0]);

  for (isds=0; isds<nsds; isds++)
    if (param_st[isds] == 1) {
      compute_sds_nrows_ncols(&out_sds_info[isds], &nrows, &ncols);
      break;
    }

  if ((sds_data = (void **)Calloc2D(fcnt, ndata_in, in_sds_info[0].data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_data \n");
  
  if ((sds_sum = (void *)calloc(ncols, in_sds_info[0].data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_sum \n");
  
  if ((sds_mean = (void *)calloc(ncols, in_sds_info[0].data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_mean \n");
  
  if ((sds_std = (float32 *)calloc(ncols, sizeof(float32))) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_std \n");
  
  if ((sds_npix = (int16 *)calloc(ncols, sizeof(int16))) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_npix \n");
  
  if ((sds_min = (void *)calloc(ncols, in_sds_info[0].data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_min \n");
  
  if ((sds_max = (void *)calloc(ncols, in_sds_info[0].data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_max \n");

  get_sds_edge(&in_sds_info[0], in_edge); 

  get_sds_edge(&out_sds_info[isds], out_edge); 

  compute_sds_start_offset(&in_sds_info[0], n, m, &st_c, &offset);

  in_rank = in_sds_info[0].rank;

  for (i=0; i<in_rank; i++) 
    in_start[i] = 0;

  for (i=0; i<out_rank; i++) 
    out_start[i] = 0;

  for (irow=0; irow<nrows; irow++)
  {
    if ((in_rank == 2) || (in_sds_info[0].dim_size[0] > in_sds_info[0].dim_size[in_rank-1]))
      in_start[0] = irow;
    else in_start[in_rank-2] = irow;
    if ((out_rank == 2) || (out_sds_info[isds].dim_size[0] > out_sds_info[isds].dim_size[out_rank-1]))
      out_start[0] = irow;                      
    else out_start[out_rank-2] = irow;

    for (fid=0; fid<fcnt; fid++)
      if (SDreaddata(in_sds_info[fid].sds_id, in_start, NULL, in_edge, sds_data[fid]) == FAIL)
	fprintf(stderr, "Error reading data line %d from SDS %s\n", irow, in_sds_info[fid].name);

    ic = st_c;

    for (icol=0; icol<ncols; icol++)
    {
      sum = sum2 = 0.0;
      npix = 0;
      for (fid=0; fid<fcnt; fid++)
      {
        switch(in_sds_info[0].data_type)
        {
	  case 5 : sds_val = ((float32 **)sds_data)[fid][ic]; break;
	  case 20: sds_val = ((int8 **)sds_data)[fid][ic]; break;
	  case 21: sds_val = ((uint8 **)sds_data)[fid][ic]; break;
	  case 22: sds_val = ((int16 **)sds_data)[fid][ic]; break;
	  case 23: sds_val = ((uint16 **)sds_data)[fid][ic]; break;
	  case 24: sds_val = ((int32 **)sds_data)[fid][ic]; break;
	  case 25: sds_val = ((uint32 **)sds_data)[fid][ic]; break;
        }

        if ((sds_val != nop_in) && (sds_val >= range[0]) && (sds_val <= range[1]))
        {
	  sum += sds_val;
          sum2 = sum2 + sds_val*sds_val;
	  if (npix == 0) 
	    min = max = sds_val;
	  else {
	    if (sds_val < min) min = sds_val;
	    if (sds_val > max) max = sds_val;
	  }
	  npix++;
        }
      } /* for (fid=0; . . ) */
      if (npix == 0)
        sum = avg = std = min = max = nop_out;
      else
      {
	avg = sum/(double)npix;
        avg2 = sum2/(double)npix;
        std = sqrt(avg2 - avg*avg); 
      }

      sds_npix[icol] = (int16)npix;
      sds_std[icol] = (float32)std;

      switch(in_sds_info[0].data_type)
      {
	case 5 : ((float32 *)sds_mean)[icol] = (float32)avg;
		 ((float32 *)sds_sum)[icol] = (float32)sum;
		 ((float32 *)sds_min)[icol] = (float32)min;
	         ((float32 *)sds_max)[icol] = (float32)max;
		 break;
	case 20: ((int8 *)sds_mean)[icol] = (int8)avg;
		 ((int8 *)sds_sum)[icol] = (int8)sum;
		 ((int8 *)sds_min)[icol] = (int8)min;
	         ((int8 *)sds_max)[icol] = (int8)max;
		 break;
	case 21: ((uint8 *)sds_mean)[icol] = (uint8)avg;
		 ((uint8 *)sds_sum)[icol] = (uint8)sum;
		 ((uint8 *)sds_min)[icol] = (uint8)min;
	         ((uint8 *)sds_max)[icol] = (uint8)max;
		 break;
	case 22: ((int16 *)sds_mean)[icol] = (int16)avg; 
		 ((int16 *)sds_sum)[icol] = (int16)sum;
	         ((int16 *)sds_min)[icol] = (int16)min;
	         ((int16 *)sds_max)[icol] = (int16)max;
		 break;
	case 23: ((uint16 *)sds_mean)[icol] = (uint16)avg; 
		 ((uint16 *)sds_sum)[icol] = (uint16)sum;
	         ((uint16 *)sds_min)[icol] = (uint16)min;
	         ((uint16 *)sds_max)[icol] = (uint16)max;
		 break;
	case 24: ((int32 *)sds_mean)[icol] = (int32)avg; 
	         ((int32 *)sds_sum)[icol] = (int32)sum;
	         ((int32 *)sds_min)[icol] = (int32)min;
	         ((int32 *)sds_max)[icol] = (int32)max;
		 break;
	case 25: ((uint32 *)sds_mean)[icol] = (uint32)avg; 
	         ((uint32 *)sds_sum)[icol] = (uint32)sum;
	         ((uint32 *)sds_min)[icol] = (uint32)min;
	         ((uint32 *)sds_max)[icol] = (uint32)max;
		 break;
      }    
      ic += offset;
    }

    if (param_st[0] == 1)
    {
      if (SDwritedata(out_sds_info[0].sds_id, out_start, NULL, out_edge, sds_sum) == FAIL)
        fprintf(stderr, "Error writing data line %d to SDS %s\n", irow, out_sds_info[0].name); 
    }

    if (param_st[1] == 1)
    {
      if (SDwritedata(out_sds_info[1].sds_id, out_start, NULL, out_edge, sds_mean) == FAIL)
        fprintf(stderr, "Error writing data line %d to SDS %s\n", irow, out_sds_info[1].name); 
    }

    if (param_st[2] == 1)
    {
      if (SDwritedata(out_sds_info[2].sds_id, out_start, NULL, out_edge, (VOIDP)sds_std) == FAIL)
        fprintf(stderr, "Error writing data line %d to SDS %s\n", irow, out_sds_info[2].name); 
    }

    if (param_st[3] == 1)
    {
      if (SDwritedata(out_sds_info[3].sds_id, out_start, NULL, out_edge, (VOIDP)sds_npix) == FAIL)
        fprintf(stderr, "Error writing data line %d to SDS %s\n", irow, out_sds_info[3].name); 
    }

    if (param_st[4] == 1)
    {
      if (SDwritedata(out_sds_info[4].sds_id, out_start, NULL, out_edge, sds_min) == FAIL)
        fprintf(stderr, "Error writing data line %d to SDS %s\n", irow, out_sds_info[4].name); 
    }

    if (param_st[5] == 1)
    {
      if (SDwritedata(out_sds_info[5].sds_id, out_start, NULL, out_edge, sds_max) == FAIL)
        fprintf(stderr, "Error writing data line %d to SDS %s\n", irow, out_sds_info[5].name); 
    }
  }

  for (isds=0; isds<nsds; isds++)
    if (param_st[isds] == 1)
      SDendaccess(out_sds_info[isds].sds_id);

  Free2D((void **)sds_data);
  free(sds_mean);
  free(sds_std);
  free(sds_npix);
  free(sds_min);
  free(sds_max);
}
