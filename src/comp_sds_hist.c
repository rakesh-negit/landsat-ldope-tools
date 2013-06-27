/****************************************************************************
!C

!File: comp_sds_hist.c

!Description:
   Compute histogram of SDS values in an SDS of an HDF file.  

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original Feburary 2003.

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
    LDOPE                             University of Maryland 
                                      Department of Geography
    droy@kratmos.gsfc.nasa.gov        NASA/GSFC Code 922 (B32)
    phone: 301-614-5571               Greenbelt, MD 20771

!Design Notes: (none)

!END
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mfhdf.h"
#include "sds_rw.h"
#include "qa_tool.h"
#include "main_util.h"
#include "alloc_mem.h"
#include "str_op.h"
#include "l2g.h"


#define MAX_NUM_INT8 256
#define MAX_NUM_INT16 65536 
#define MAX_NUM_INT32 4294967296 

#define HELP \
"NAME \n" \
"    comp_sds_hist - Print histogram of data values in one or more SDS of a\n" \
"                    MODIS Land HDF-EOS data product.\n" \
" \n" \
"SYNOPSIS \n" \
"    comp_sds_hist [-help] [filename]\n" \
"    comp_sds_hist [-sds=<SDS_name1>[,<SDS_name2>. . ]] [-layer]\n" \
"                  [-range=<min,max>] filename\n" \
" \n" \
"DESCRIPTION \n" \
"    Compute histogram of data values in one or more SDS of a MODIS Land\n" \
"    HDF-EOS data product. The histogram may be computed for the user\n" \
"    specified range of SDS values. The output to stdout includes the\n" \
"    contains SDS name, dimension size, fill value and a list of SDS values.\n"\
" \n" \
"    If an SDS is 3D or 4D, then the tool can optionally output the\n" \
"    histogram for each layer/slice of the 3D/4D SDS.\n" \
" \n" \
"    The tool command arguments can be specified in any order.\n" \
" \n" \
"OPTIONS \n" \
"    -help             Display this help message. If the input filename is\n" \
"                      specified with this option, then the names of all\n" \
"                      the SDS in the file are displayed.\n" \
"    <SDS_list>        List of SDSs to read. SDS names are separated by\n" \
"                      commas with no space. By default sds values are\n" \
"                      printed for all SDSs in the input file.\n" \
"    -layer            Compute histogram for every layer/slice separately\n" \
"                      for a 3D/4D SDS.\n" \
"    -range=<min,max>  Histogram range (minimum and maximum values). Default\n"\
"                      is set to valid range of the SDS. Fill value is counted\n" \
"                      separately. If valid range attribute is not available\n"\
"                      the range of the SDS data type is used as the limit.\n" \
"                      The range value for various data type is shown below\n" \
"                      INT8:  (-128, 127)     UINT8: (0, 255)\n" \
"                      INT16: (-32768, 32767) UINT16: (0, 65535)\n" \
"                      INT32: (-2147483648, 2147483647) UINT32: (0, 4294967295)\n" \
"                      FLOAT32: UNDEFINED.\n" \
"                      For float data type the histogram is computed after\n" \
"                      converting the float values to their closest integer.\n" \
"    Filename          input filenames \n" \
" \n" \
"Examples: \n" \
"    comp_sds_hist -sds=sur_refl_b01\n" \
"                  MOD09A1.A2001033.h08v05.001.2001166175830.hdf\n" \
"\n" \
"    comp_sds_hist -layer -sds=Surface_Refl -range=0,10000 -layer\n" \
"                  MODAGAGG.A2000065.h13v02.002.2000075160322.hdf\n" \
"AUTHOR: \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 04/05/2004\n" \

#define USAGE \
"usage:	\n" \
"    comp_sds_hist [-help] [filename]\n" \
"    comp_sds_hist [-sds=<SDS_name1>[,<SDS_name2>. . ]] [-layer]\n" \
"                  [-range=<min,max>] filename\n" \
"\n" \
"OPTIONS \n" \
"    -help             Display this help message. If the input filename is\n" \
"                      specified with this option, then the names of all\n" \
"                      the SDS in the file are displayed.\n" \
"    <SDS_list>        List of SDSs to read. SDS names are separated by\n" \
"                      commas with no space. By default sds values are\n" \
"                      printed for all SDSs in the input file.\n" \
"    -layer            Compute histogram for every layer/slice separately\n" \
"                      for a 3D/4D SDS.\n" \
"    -range=<min,max>  Histogram range (minimum and maximum values). Default\n"\
"                      is set to valid range of the SDS. Fill value is counted\n" \
"                      separately. If valid range attribute is not available\n"\
"                      the range of the SDS data type is used as the limit.\n"\
"                      The range value for various data type is shown below\n"\
"                      INT8:  (-128, 127)     UINT8: (0, 255)\n" \
"                      INT16: (-32768, 32767) UINT16: (0, 65535)\n" \
"                      INT32: (-2147483648, 2147483647) UINT32: (0, 4294967295)\n" \
"                      FLOAT32: UNDEFINED. \n" \
"                      For float data type the histogram is computed after\n" \
"                      converting the float values to their closest integer.\n" \
"    Filename          input filenames \n" \
"\n"

int parse_cmd_comp_sds_hist(int argc, char **argv, int *nd_info, char **sds_names, int *nsds, int *hist_range, int *fcnt);
void compute_comp_sds_hist(char *fname, char **sds_names, int nsds, int nd_info, int *hist_range);
void add_to_hist_from_row(sds_t *sds_info, void *data_in, int ncol, int st_c, int offset, 
			  int *hist_cnt, int *fill_cnt, int *hist_range);
void print_comp_sds_hist(sds_t *sds_info, int **hist_cnt, int *fill_cnt, int n_layer, int *hist_range);

int int8_range[2] = {-128, 127};
int uint8_range[2] = {0, 255};
int int16_range[2] = {-32768, 32767};
int uint16_range[2] = {0, 65535};
int32 int32_range[2] = {(-2147483647 - 1), 2147483647};
uint32 uint32_range[2] = {0, 4294967295U};

int main(int argc, char **argv)
/******************************************************************************
!C

!Description:
  Main function for comp_sds_hist

!Input Parameters: (none)
  command line arguments: see help for details.

!Output Parameters: (none)
  return 0 on successful completion of the process

!Revision History:
  December, 2003 version 1.0 

!Team-unique Header:
  See file prologue.

!References and Credits: (see file prologue)

!Design Notes: (none)

!END
********************************************************************************/

{
  char **sds_names;
  int hist_range[2];
  int nsds, sds_cnt;
  int i, fcnt, nd_info;

  if (argc == 1)
    {
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_FAILURE);
    }   

  if ((argc==2) && ((strcmp(argv[1],"-help")==0) || (strcmp(argv[1], "-h")==0)))
  {
    fprintf(stderr, "%s\n", HELP);
    exit(EXIT_SUCCESS);
  }

  /*  Display SDS names of input HDF file */

if ((argc>=3) && ((strcmp(argv[1], "-help")==0) || (strcmp(argv[1], "-h")==0)))
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


  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for sds_names in sds_range: main()\n");
  else
  {
    if (parse_cmd_comp_sds_hist(argc, argv, &nd_info, sds_names, &nsds, hist_range, &fcnt) == -1)
      {
	fprintf(stderr, "%s\n", USAGE);
	exit(EXIT_FAILURE);
      }	
    else
    {
      for (i=1; i<argc; i++)
      {
        if (argv[i][0] != '-')
        {
	  fprintf(stdout, "Reading from %s\n", argv[i]);
	  if (nsds == 0) 
	    {
	      /* no SDS inputed from command line. process all SDS */
	      sds_cnt = get_sds_names(argv[i], sds_names);
	    }
          else sds_cnt = nsds;
	  compute_comp_sds_hist(argv[i], sds_names, sds_cnt, nd_info, hist_range); 
        }
      }
    }
    Free2D((void **)sds_names);
  }
  fprintf(stderr, "Processing done ! \n");
  return 0;
}

int parse_cmd_comp_sds_hist(int argc, char **argv, int *nd_info, char **sds_names, int *nsds, 
		       int *hist_range, int *fcnt)
/******************************************************************************
!C

!Description:
  Function to parse command line arguments.

!Input Parameters:
  argc: number of input arguments
  argv: string array containing arguments

!Output Parameters:
  nd_info    : flag of if -layer option is specified. Equals to 1 if -layer
		             option inputed.
  sds_names  : input SDS names
  sds_cnt    : number of SDS
  hist_range : User input histogram range.
  fcnt       : number of input files.

  return 1 if parsing is succesfull, -1 if not all required parameters input.
 
!Revision History:
  Original December 2003.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.
    

!Design Notes: (none)

!END
********************************************************************************/

{
  int i, p1, st, len;
  char range_str[50], val_str[25];

  st = 1;
  range_str[0] = '\0';
  *fcnt = *nsds = *nd_info = 0;
  hist_range[0] = hist_range[1] = -111;
  for (i=1; i<argc; i++)
  {
    if (strcmp(argv[i], "-layer") == 0) *nd_info = 1;
    else if (is_arg_id(argv[i], "-sds") == 0)
      get_arg_val_arr(argv[i], sds_names, nsds);
    else if (is_arg_id(argv[i], "-range") == 0)
      get_arg_val(argv[i], range_str);
    else if (argv[i][0] == '-')
      fprintf(stderr, "Ignoring unknown option %s\n", argv[i]);
    else ++*fcnt;
  }
  if (*fcnt == 0) {
    st = -1;
    fprintf(stderr, "Missing input file . . \n");
  }
  if (range_str[0] != '\0')
  {
    p1 = sd_charpos(range_str, ',', 0);
    if (p1 != -1)
    {
      len = (int)strlen(range_str);
      sd_strmid(range_str, 0, p1, val_str);
      hist_range[0] = (int)atoi(val_str);
      sd_strmid(range_str, p1+1, len-p1-1, val_str);
      hist_range[1] = (int)atoi(val_str);
      if (hist_range[0] > hist_range[1]) {
        st = -1;
        fprintf(stderr, "Invalid range option %s\n", range_str);
      }
    }
    else {
      st = -1;
      fprintf(stderr, "Invalid range option %s\n", range_str);
    } 
  }
  else 
    fprintf(stderr, "No range option Input. Using valid range from SDS or default range of the data type. \n");
  if (st == 1)
    if (*nsds == 0) fprintf(stderr, "No SDS name input. Reading all SDSs \n");
  return st;
}

void compute_comp_sds_hist(char *fname, char **sds_names, int nsds, int nd_info, int *hist_range)
/******************************************************************************
!C

!Description:
  Function to compute the SDS histogram.

!Input Parameters:
  fname     : Input filename.
  sds_names : Input SDS names.
  nsds      : Number of input SDS.
  nd_info   : Flag of if -layer option is specified.
  hist_range: Input histogram range.

!Output Parameters:
  None.

!Revision History:
  Original December 2003.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)

!END
********************************************************************************/
{
  void *data_in;
  sds_t sds_info;
  int bsq = 0, rank;
  int n, m, max_m = 0;
  int st_c = 0, offset = 0;
  int irange, n_val;
  int layer_id, n_layer;
  int irow, irank, isds;
  int nrows, ncols, ndata;
  int *fill_cnt, **hist_cnt;
  int32 start[4], edge[4];

  sds_info.sd_id = -1;
  for (isds=0; isds<nsds; isds++)
  {
    sds_info.sds_id = -1;
    strcpy(sds_info.name, sds_names[isds]);
    for (irank=0; irank<4; irank++)
      start[irank] = edge[irank] = 0;
    sds_info.range[0] = sds_info.range[1] = -111;
    if (get_sds_info(fname, &sds_info) != -1)
    {
      if (sds_info.data_type == 5)
      {
	if ((hist_range[0] == -111) && (hist_range[1] == -111))
	  fprintf(stderr, "SDS %s is of float data type. Needs Histogram range input\n", sds_info.name);
      }
      else 
      {
	if ((hist_range[0] == -111) && (hist_range[1] == -111))
        {
	  if ((sds_info.range[0] == -111) && (sds_info.range[1] == -111))
	  {
	    switch(sds_info.data_type) {
	      case 20: hist_range[0] = int8_range[0]; hist_range[1] = int8_range[1]; break; 
	      case 21: hist_range[0] = uint8_range[0]; hist_range[1] = uint8_range[1]; break; 
	      case 22: hist_range[0] = int16_range[0]; hist_range[1] = int16_range[1]; break; 
	      case 23: hist_range[0] = uint16_range[0]; hist_range[1] = uint16_range[1]; break; 
	      case 24: hist_range[0] = int32_range[0]; hist_range[1] = int32_range[1]; break; 
	      case 25: hist_range[0] = uint32_range[0]; hist_range[1] = uint32_range[1]; break; 
	    }
	  }
	  else
	  { 
             hist_range[0] = sds_info.range[0];
             hist_range[1] = sds_info.range[1];
          }
        }
      }
      if ((hist_range[0] != -111) && (hist_range[1] != -111))
      {
        rank = sds_info.rank;
        n_layer = 1;
        if (rank == 1)
        {
          nrows = 1;
	  edge[0] = ncols = ndata = sds_info.dim_size[0];
        }
        else
        {
          bsq = ((rank == 2) || (sds_info.dim_size[0] < sds_info.dim_size[rank - 1])) ? 1 : 0;
          if (bsq == 1)
          {
            nrows = sds_info.dim_size[rank-2];
            ncols = ndata = edge[rank-1] = sds_info.dim_size[rank-1];
            for (irank=0; irank<rank-2; irank++)
            {
	      n_layer *= sds_info.dim_size[irank];
              ndata *= sds_info.dim_size[irank];
              edge[irank] = sds_info.dim_size[irank];
            }
            edge[rank-2] = 1;
	    max_m = sds_info.dim_size[1];
          }
          else
          {
	    nrows = sds_info.dim_size[0];
	    ncols = ndata = sds_info.dim_size[1];
	    for (irank=2; irank<rank; irank++)
	    {
	      n_layer *= sds_info.dim_size[irank];
	      ndata *= sds_info.dim_size[irank];
	      edge[irank] = sds_info.dim_size[irank];
	    }
	    edge[0] = 1;
	    edge[1] = sds_info.dim_size[1];
	    max_m = sds_info.dim_size[3];
          }
        }
  
	if ((nd_info == 0) || (rank <= 2))
	{
	  st_c = 0;
	  offset = 1;
	  n_layer = 1;
	}
  
        if ((data_in = (void *)calloc(ndata, sds_info.data_size)) == NULL)
	  fprintf(stderr, "Cannot allocate memory for data_in in read_sds_range()\n");
        n_val = hist_range[1] - hist_range[0] + 1;
	if ((hist_cnt = (int **)Calloc2D(n_layer, n_val, sizeof(int))) == NULL)
	  fprintf(stderr, "Cannot allocate memory for hist_cnt in read_sds_values\n");	
	if ((fill_cnt = (int *)calloc(n_layer, sizeof(int))) == NULL)
	  fprintf(stderr, "Cannot allocate memory for fill_cnt in read_sds_values\n");	
        if ((data_in != NULL) && (hist_cnt != NULL) && (fill_cnt != NULL)) 
	{
          for (layer_id=0; layer_id<n_layer; layer_id++)
	    fill_cnt[layer_id] = 0;
          for (layer_id=0; layer_id<n_layer; layer_id++)
	    for (irange=0; irange<n_val; irange++)
	      hist_cnt[layer_id][irange] = 0;
          for (irow=0; irow<nrows; irow++)
          {
            if (rank == 1) start[0] = 0;
	    else
	    {
	      if (bsq == 1) start[rank-2] = irow;
	      else start[0] = irow;
	    }
	    if (SDreaddata(sds_info.sds_id, start, NULL, edge, data_in) == FAIL)
	      fprintf(stderr, "Failed to read data row for SDS %s\n", sds_info.name);
	    else
            {
	      for (layer_id=0; layer_id<n_layer; layer_id++)
	      {
	        if (nd_info == 1) 
	        {
		  m = (rank == 3) ? -1 : layer_id%max_m;
		  n = (rank == 3) ? layer_id : layer_id/max_m;
		  compute_sds_start_offset(&sds_info, n, m, &st_c, &offset);
	        }
	        else ncols = ndata;
	        add_to_hist_from_row(&sds_info, data_in, ncols, st_c, offset, hist_cnt[layer_id], &fill_cnt[layer_id], hist_range);
	      }
	    }
          } /* for (irow=0;  . .  ) */
          print_comp_sds_hist(&sds_info, hist_cnt, fill_cnt, n_layer, hist_range);
	}
	Free2D((void **)hist_cnt);
	free(fill_cnt);
        free(data_in);
        if (sds_info.sds_id != -1) SDendaccess(sds_info.sds_id);
      }
    }
  } /* for (isds=0;  . .  ) */
  SDend(sds_info.sd_id);
}
	        
void add_to_hist_from_row(sds_t *sds_info, void *data_in, int ncol, int st_c, 
			  int offset, int *hist_cnt, int *fill_cnt, int *hist_range)
/******************************************************************************
!C

!Description:
  Function to process the SDS histogram row by row.

!Input Parameters:
  sds_info:   The SDS information structure.    
  data_in:    A whole row data
  ncol:       Number of columns in the row.
  st_c:       Starting indices of data to be compute.
  offset:     Offset of the data to be compute.
  hist_cnt:   histogram count.
  fill_cnt:   Fill value count.
  hist_range: Input histogram range.

!Output Parameters:
  None.

!Revision History:
  Original December 2003.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)

!END
********************************************************************************/
{
  float sds_fval;
  int i, id, icol, sds_val;

  switch(sds_info->data_type)
  {
    case 5 : for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_fval = ((float *)data_in)[icol]; 
      	       sds_val = (int)sds_fval; 
      	       if (sds_fval != sds_info->fill_fval) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
	       }	
             }
    	     break;
    case 20: for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_val = ((int8 *)data_in)[icol]; 
      	       if (sds_val == sds_info->fill_val) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
      	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
      	       }
             }
    	     break;
    case 21: for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_val = ((uint8 *)data_in)[icol]; 
      	       if (sds_val == sds_info->fill_val) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
      	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
      	       }
             }
    	     break;
    case 22: for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_val = ((int16 *)data_in)[icol]; 
      	       if (sds_val == sds_info->fill_val) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
      	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
      	       }
             }
    	     break;
    case 23: for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_val = ((uint16 *)data_in)[icol]; 
      	       if (sds_val == sds_info->fill_val) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
      	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
      	       }
             }
    	     break;
    case 24: for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_val = ((int32 *)data_in)[icol]; 
      	       if (sds_val == sds_info->fill_val) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
      	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
      	       }
             }
    	     break;
    case 25: for (i=0, icol=st_c; i<ncol; icol += offset, i++)
    	     {
      	       sds_val = ((uint32 *)data_in)[icol]; 
      	       if (sds_val == sds_info->fill_val) *fill_cnt += 1;
	       else if ((sds_val >= hist_range[0]) && (sds_val <= hist_range[1]))
      	       {
		 id = sds_val - hist_range[0];
		 hist_cnt[id] += 1;
      	       }
             }
    	     break;
  }
}
        
void print_comp_sds_hist(sds_t *sds_info, int **hist_cnt, int *fill_cnt, int n_layer, 
		    int *hist_range)
/******************************************************************************
!C

!Description:
  Function to print out the SDS histogram table.

!Input Parameters:
  sds_info:   The SDS information structure. 
  hist_cnt:   Array containing the histogram.
  fill_cnt:   Number of fill values.
  n_layer:    Number of layer.
  hist_range: Histogram range.

!Output Parameters:
  None.

!Revision History:
  Original December 2003.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)

!END
********************************************************************************/
{
  int i, j = 0, id, sum;
  char fval_str[25], dim_str[80];

  if (sds_info->rank == 2) 
    sprintf(dim_str, "Dimension = (" LONG_INT_FMT " x " LONG_INT_FMT ")", sds_info->dim_size[0], sds_info->dim_size[1]);
  else if (sds_info->rank == 3) 
    sprintf(dim_str, "Dimension = (" LONG_INT_FMT " x " LONG_INT_FMT " x " LONG_INT_FMT ")", sds_info->dim_size[0], sds_info->dim_size[1], sds_info->dim_size[2]);
  else if (sds_info->rank == 4) 
    sprintf(dim_str, "Dimension = (" LONG_INT_FMT " x " LONG_INT_FMT " x " LONG_INT_FMT " x " LONG_INT_FMT ")", sds_info->dim_size[0], sds_info->dim_size[1], sds_info->dim_size[2], sds_info->dim_size[3]);
  if (sds_info->data_type == 5) sprintf(fval_str, "Fill Value = %f", sds_info->fill_fval);
  else sprintf(fval_str, "Fill Value = %ld", sds_info->fill_val);
  fprintf(stdout, "%s:\t%s\t%s\n", sds_info->name, dim_str, fval_str);

  for (i=hist_range[0]; i<=hist_range[1]; i++)
  {
    id = i - hist_range[0];
    for (j=0, sum=0; j<n_layer; j++)
      sum += hist_cnt[j][id];
    if (sum != 0) { 
      fprintf(stdout, "%d", i);
      for (j=0; j<n_layer; j++)
        fprintf(stdout, "\t%d", hist_cnt[j][id]);
      fprintf(stdout, "\n");
    }
  }
  if (fill_cnt[j] != 0) { 
    if (sds_info->data_type == 5)
      fprintf(stdout, "%f", sds_info->fill_fval);
    else 
      fprintf(stdout, "%ld", sds_info->fill_val);
    for (j=0; j<n_layer; j++)
      fprintf(stdout, "\t%d", fill_cnt[j]);
    fprintf(stdout, "\n"); 
  }
}
