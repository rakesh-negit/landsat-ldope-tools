/****************************************************************************
!C

!File: subset_sds.c

!Description:
  Create spatial subset SDSs of one or more SDS from an MODIS Land HDF-EOS 
  data product.

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  December 2003. Version 1.0

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
    LDOPE                             University of Maryland, Department of Geography
    droy@kratmos.gsfc.nasa.gov        NASA/GSFC Code 922 (B32)
    phone: 301-614-5571               Greenbelt, MD 20771   

!Design Notes: (none)

!END
*****************************************************************************/         

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mfhdf.h"
#include "alloc_mem.h"
#include "str_op.h"
#include "sds_rw.h"
#include "main_util.h"
#include "meta.h"

#define HELP \
"NAME \n" \
"    subset_sds - Create spatial subset SDSs of one or more SDS from an\n"\
"                 MODIS Land HDF-EOS data product. \n" \
" \n" \
"SYNOPSIS \n" \
"    subset_sds -help [filename] \n" \
"    subset_sds -of=<output filename> -row=<start,end>  -col=<start,end> \n" \
"               [-sds=<SDS_name1>[,<SDS_name2>,. . ]] filename \n" \
" \n" \
"DESCRIPTION \n" \
"    Create spatial subset SDS(s) of one or more SDS from an input MODIS \n" \
"    Land HDF-EOS data product. Spatial subset is specified by the\n"\
"    row-column range of the data in input product. \n" \
" \n" \
"    This tool supports 2D/3D/4D SDS(s). \n" \
" \n" \
"    The tool command arguments can be specified in any order. \n" \
" \n" \
"OPTIONS \n" \
"    -help                   Display this help message, If the input\n"\
"                            filename is specified with this option, then\n"\
"                            the names of all the SDS in the file are\n"\
"                            displayed. \n" \
"    -sds=<SDS_list>         List of SDS to enlarge. SDS names are separated\n"\
"                            by commas with no space. By default all SDSs\n"\
"                            are processed maintaining the input SDS\n"\
"                            interleaving. \n" \
" \n" \
"                            To process a specific layer of a 3D SDS specify\n"\
"                            the element number of the third dimension as a\n"\
"                            dot extension of the SDS name: sds_name.n\n"\
"                            (e.g., sur_refl_b02.1 = the layer defined by\n"\
"                            the 1st element of the 3rd dimension of the 3D\n"\
"                            SDS sur_refl_b02). \n" \
" \n" \
"                            To process a specific layer of a 4D SDS,\n"\
"                            specify the higher dimension element number(s)\n"\
"                            as a dot extension of the SDS name:\n"\
"                            sds_name.n.m (e.g., Surface_Refl.1.2 = the\n"\
"                            layer defined by the 1st element of the 3rd\n"\
"                            dimension and the 2nd element of the 4th\n"\
"                            dimension of the 4D SDS Surface_Refl). \n" \
" \n" \
"                            Note that wildcards and ranges of element\n"\
"                            values may be specified as sds_name.* and as\n"\
"                            sds_name.n1-n2.m respectively. \n" \
"    -row=<start,end>        Subset row range (start and end are inclusive)\n"\
"    -col=<start,end>        Subset column range (start and end are\n"\
"                            inclusive)\n" \
"    -of=<out filename>      Output filename \n" \
"    filename                Input filename \n" \
" \n" \
"Example: \n" \
"    subset_sds -sds=\"most confident detected fire\" -row=0,10 \n" \
"               -col=1130,1140 MYD14A1.A2003281.h09v05.003.2003296025915.hdf\n"\
"               -of=Subset_MYD14A1.A2003281.h09v05.003.2003296025915.hdf \n" \
" \n" \
"    subset_sds -sds=Surface_Refl -row=10,100 -col=100,200 \n" \
"               MODAGAGG.A2003077.h19v07.004.2003079210733.hdf \n" \
"               -of=Subset_MODAGAGG.A2003077.h19v07.004.2003079210733.hdf \n" \
" \n" \
"    subset_sds -sds=EV_Band26 -row=300,400 -col=1000,1100 \n" \
"               MOD021KM.A2004007.0155.004.2004007095233.hdf \n" \
"               -of=Subset_MOD021KM.A2004007.0155.004.2004007095233.hdf \n" \
"\n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 04/05/2004 \n" \

#define USAGE \
"    subset_sds -help [filename] \n" \
"    subset_sds -of=<output filename> -row=<start,end>  -col=<start,end> \n" \
"               [-sds=<SDS_name1>[,<SDS_name2>,. . ]] filename \n" \
" \n" \
"OPTIONS \n" \
"    -help                   Display this help message, If the input\n"\
"                            filename is specified with this option, then\n"\
"                            the names of all the SDS in the file are\n"\
"                            displayed. \n" \
"    -sds=<SDS_list>         List of SDS to enlarge. SDS names are separated\n"\
"                            by commas with no space. By default all SDSs\n"\
"                            are processed maintaining the input SDS\n"\
"                            interleaving. \n" \
" \n" \
"                            To process a specific layer of a 3D SDS specify\n"\
"                            the element number of the third dimension as a\n"\
"                            dot extension of the SDS name: sds_name.n\n"\
"                            (e.g., sur_refl_b02.1 = the layer defined by\n"\
"                            the 1st element of the 3rd dimension of the 3D\n"\
"                            SDS sur_refl_b02). \n" \
" \n" \
"                            To process a specific layer of a 4D SDS,\n"\
"                            specify the higher dimension element number(s)\n"\
"                            as a dot extension of the SDS name:\n"\
"                            sds_name.n.m (e.g., Surface_Refl.1.2 = the\n"\
"                            layer defined by the 1st element of the 3rd\n"\
"                            dimension and the 2nd element of the 4th\n"\
"                            dimension of the 4D SDS Surface_Refl). \n" \
" \n" \
"                            Note that wildcards and ranges of element\n"\
"                            values may be specified as sds_name.* and as\n"\
"                            sds_name.n1-n2.m respectively. \n" \
"    -row=<start,end>        Subset row range (start and end are inclusive)\n"\
"    -col=<start,end>        Subset column range (start and end are\n"\
"                            inclusive)\n" \
"    -of=<out filename>      Output filename \n" \
"    filename                Input filename \n" \
" \n"


/**************************************************************************************************
                            		Prototypes.
**************************************************************************************************/

int parse_cmd_subset_sds(int argc, char **argv, char **sds_names, int *nsds, int *row_range, 
		int *col_range, char *out_fname, char *in_fname);
/* void subset_an_sds(sds_t *in_sds_info, sds_t *out_sds_info, int *row_range, int *col_range); */
void subset_an_sds(sds_t *in_sds_info, sds_t *out_sds_info, int *row_range, int *col_range);

int main(int argc, char *argv[])
/******************************************************************************
!C

!Description:
  Main function for subset_sds.

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
  sds_t in_sds_info;
  sds_t out_sds_info;
  int i, st, nsds, isds, m, n;
  int row_range[5], col_range[5];
  char **sds_names;
  char in_fname[MAX_PATH_LENGTH];
  char out_fname[MAX_PATH_LENGTH];

  if (argc == 1)
    {
      st = -1;
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_FAILURE);
    }   
  
  if ((argc==2) && (strcmp(argv[1], "-help")==0))
    {
      fprintf(stderr, "%s\n", HELP);
      exit(0);
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

  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for sds_names in subset_sds: main()\n");
      exit(EXIT_FAILURE);
    }
  else
    {
      if ((st = parse_cmd_subset_sds(argc, argv, sds_names, &nsds, row_range, col_range, 
				     out_fname, in_fname)) == -1)
	fprintf(stderr, "%s\n", USAGE);
      else if (st != 0)
	{
	  if (nsds == 0)
	    nsds = get_sds_names(in_fname, sds_names);
	  
	  /* open output HDF file and set global attribute input pointer
	     to contain names of all input files */
	  out_sds_info.sd_id = -1;
	  
	  if ((out_sds_info.sd_id = SDstart(out_fname, DFACC_CREATE)) == FAIL)
	  {
	    fprintf(stderr, "Cannot create the output hdf file %s\n", out_fname);
	    exit(EXIT_FAILURE);
	  }
	  else 
	    {
	      in_sds_info.sd_id = -1;
	      
	      for (isds=0; isds<nsds; ++isds)
		{
		  fprintf(stdout, "Processing SDS %s\n", sds_names[isds]);
		  
		  get_sdsname_dim(sds_names[isds], in_sds_info.name, &n, &m);
		  
		  if ((n != -1) || (m != -1))
		    fprintf(stdout, "A 2D slice of the 3D/4D SDS can't be selected to output\n");
		}
	      
	      in_sds_info.sds_id = -1;
	      
	      if (get_sds_info(in_fname, &in_sds_info) != -1)
		{	    
		  subset_an_sds(&in_sds_info, &out_sds_info, row_range, col_range); 
		}
	      
	      SDend(in_sds_info.sd_id);
	      SDend(out_sds_info.sd_id);
	    }
	}
      Free2D((void **)sds_names);
      fprintf(stderr, "Processing done ! \n");
    }
  return 0;
}

int parse_cmd_subset_sds(int argc, char **argv, char **sds_names, int *nsds, int *row_range,
			 int *col_range, char *out_fname, char *in_fname)
/******************************************************************************
!C

!Description:
  Function to parse command line arguments.

!Input Parameters:
  argc: number of input arguments
  argv: string array containing arguments

!Output Parameters:
  sds_names: input SDS names
  nsds:      number of SDS
  out_fname: output filename
  row_range: Input subseting row range 
  col_range: Input subseting column range 
  out_fname: output filename
  in_fname:  input  filename

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
  int i, cnt;
  int st = 1;
  char **tmp_str;

  if ((tmp_str = (char **)Calloc2D(5, 20, sizeof(char))) == NULL)
  {
    st = 0;
    fprintf(stderr, "Cannot allocate memory for tmp_str in parse_cmd_subset_sds\n");
  }
  else 
  {
    *nsds = 0;
    in_fname[0] = out_fname[0] = '\0';
    row_range[0] = row_range[1] = 0;
    col_range[0] = col_range[1] = 0;
    for (i=1; i<argc; i++)
    {
      cnt = 0;
      if (is_arg_id(argv[i], "-sds=") == 0)
        get_arg_val_arr(argv[i], sds_names, nsds);
      else if (is_arg_id(argv[i], "-of=") == 0)
        get_arg_val(argv[i], out_fname);
      else if (is_arg_id(argv[i], "-row=") == 0) {
        get_arg_val_arr(argv[i], tmp_str, &cnt);
	if (cnt != 2) {
	  st = -1;
	  fprintf(stderr, "Invalid value for option -row: %s\n", argv[i]);
	}
	else {
	  row_range[0] = (int)atoi(tmp_str[0]);
	  row_range[1] = (int)atoi(tmp_str[1]);
        }
      }
      else if (is_arg_id(argv[i], "-col=") == 0) {
        get_arg_val_arr(argv[i], tmp_str, &cnt);
	if (cnt != 2) {
	  st = -1;
	  fprintf(stderr, "Invalid value for option -col: %s\n", argv[i]);
	}
	else {
	  col_range[0] = (int)atoi(tmp_str[0]);
	  col_range[1] = (int)atoi(tmp_str[1]);
        }
      }
      else if (argv[i][0] == '-') 
        fprintf(stderr, "Ignoring unknown option %s\n", argv[i]);
      else strcpy(in_fname, argv[i]);
    }
    if (in_fname[0] == '\0') 
    {
      st = -1;
      fprintf(stderr, "Missing input filename . . \n");
    }
    if (out_fname[0] == '\0') 
    {
      st = -1;
      fprintf(stderr, "Missing output filename . . \n");
    }
    if ((row_range[0] >= row_range[1]) || (row_range[0] < 0) || (row_range[1] < 0)) 
    {
      st = -1;
      fprintf(stderr, "Invalid subset row range . . \n");
    }
    if ((col_range[0] >= col_range[1]) || (col_range[0] < 0) || (col_range[1] < 0)) 
    {
      st = -1;
      fprintf(stderr, "Invalid subset column range . . \n");
    }
  }
  return st;
}

void subset_an_sds(sds_t *in_sds_info, sds_t *out_sds_info, int *row_range, int *col_range)
/*
!C******************************************************************************

!Function: subset_an_sds
       
!Description:

  create a subset of the original SDS to output. Original interleaving is preserved.

!Input Parameters: 
  in_sds_info    Array of input SDS information structure
  out_sds_info   Array of output SDS information structure
  row_range      Row range of input SDS to subset
  col_range      Column range of input SDS to subset

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/
{
  int i, k, irow, ncols; 
  int rank, bsq, st_c, offset;
  int32 in_start[4] = {0, 0, 0, 0};
  int32 in_edge[4] = {0, 0, 0, 0};
  int32 out_start[4] = {0, 0, 0, 0};
  int32 out_edge[4] = {0, 0, 0, 0};
  int32 attr_type, attr_cnt;
  void *data_in, *data_out, *attr_val;
  int n, m, j, nrows, icol, ic;
  int start_row, start_col, end_row, end_col;
  int ndata_sm_in = 0, ndata_sm_out = 0;
  int total_layer = 0;
  rank = in_sds_info->rank;

  n = m = -1;  /* no layer can be selected */
  
  /*create output SDS with proper dimension and other attributes */
  if ((attr_val = (void *)get_sds_attr(in_sds_info->sds_id, "_FillValue", &attr_type, &attr_cnt)) == NULL)
    fprintf(stderr, "Attribute "_FillValue" not defined for output\n");


  bsq = ((rank == 2) || (in_sds_info->dim_size[0] < in_sds_info->dim_size[rank-1])) ? 1 : 0;
  
  /* initialize the dimension size of output sds */
  for (j=0; j<4; j++)
    out_sds_info->dim_size[j] = 0;

  for (j=0; j<rank; j++)
    out_sds_info->dim_size[j] = in_sds_info->dim_size[j];

  start_row = row_range[0] < row_range[1] ? row_range[0] : row_range[1];
  end_row = row_range[1] > row_range[0] ? row_range[1] : row_range[0];
  start_col = col_range[0] < col_range[1] ? col_range[0] : col_range[1];
  end_col = col_range[1] > col_range[0] ? col_range[1] : col_range[0];

  nrows = end_row - start_row + 1;
  ncols = end_col - start_col + 1;

  if (rank == 2)
    {
      in_edge[0] = out_edge[0] = 1;
      in_edge[1] = in_sds_info->dim_size[1];
      out_edge[1] = ncols;
      out_sds_info->dim_size[0] = nrows;
      out_sds_info->dim_size[1] = ncols;
      ndata_sm_out = ncols;
      ndata_sm_in = in_sds_info->dim_size[1];
      total_layer = 1;

      if ((row_range[0] > in_sds_info->dim_size[0]) ||
	  (row_range[1] > in_sds_info->dim_size[0]) ||
	  (col_range[0] > in_sds_info->dim_size[1]) ||
	  (col_range[1] > in_sds_info->dim_size[1]))
	{
	  fprintf(stderr, "Input subsetting range is incorrect\n");
	  exit(EXIT_FAILURE);
	}
    }
  else if (rank > 2)
    {
      for (i=0; i<rank; i++)
	{
	  in_edge[i] = out_edge[i] = in_sds_info->dim_size[i];
	}
      
      if ((n == -1) && (m == -1))
	{
	  if (bsq == 0)
	    {
	      in_edge[0] = out_edge[0] = 1;
	      if ((row_range[0] > in_sds_info->dim_size[0]) ||
		  (row_range[1] > in_sds_info->dim_size[0]) ||
		  (col_range[0] > in_sds_info->dim_size[1]) ||
		  (col_range[1] > in_sds_info->dim_size[1]))
		{
		  fprintf(stderr, "Input subsetting range is incorrect\n");
		  exit(EXIT_FAILURE);
		} 
	      out_edge[1] = ncols;
	      out_sds_info->dim_size[0] = nrows;
	      out_sds_info->dim_size[1] = ncols; 
	      if (rank == 3) 
		{
		  ndata_sm_in = in_sds_info->dim_size[0] * in_sds_info->dim_size[2];
		  ndata_sm_out = ncols * out_sds_info->dim_size[2];
		  total_layer = in_sds_info->dim_size[2];
		}
	      if (rank == 4) 
		{ 
		  ndata_sm_in = in_sds_info->dim_size[0] * in_sds_info->dim_size[2] * in_sds_info->dim_size[3];
		  ndata_sm_out = ncols * out_sds_info->dim_size[2] *out_sds_info->dim_size[3]; 
		  total_layer = in_sds_info->dim_size[2] * in_sds_info->dim_size[3];
		}
	    }
	  else
	    {
	      in_edge[rank-2] = out_edge[rank-2] = 1;
	      if ((row_range[0] > in_sds_info->dim_size[rank-2]) ||
		  (row_range[1] > in_sds_info->dim_size[rank-2]) ||
		  (col_range[0] > in_sds_info->dim_size[rank-1]) ||
		  (col_range[1] > in_sds_info->dim_size[rank-1]))
		{
		  fprintf(stderr, "Input subsetting range is incorrect\n");
		  exit(EXIT_FAILURE);
		}      
	      out_edge[rank-1] = ncols;
	      out_sds_info->dim_size[rank-2] = nrows;
	      out_sds_info->dim_size[rank-1] = ncols;
	      if (rank == 3) 
		{ 
		  ndata_sm_in = in_sds_info->dim_size[0] * in_sds_info->dim_size[2];
		  ndata_sm_out = ncols * out_sds_info->dim_size[0];
		  total_layer = in_sds_info->dim_size[0];	  
		}
	      if (rank == 4) 
		{ 
		  ndata_sm_in = in_sds_info->dim_size[0] * in_sds_info->dim_size[1] * in_sds_info->dim_size[3];
		  ndata_sm_out = ncols * out_sds_info->dim_size[0] * out_sds_info->dim_size[1];
		  total_layer = in_sds_info->dim_size[0] * in_sds_info->dim_size[1];
		}
	    }
	}
    }

  out_sds_info->sds_id = -1;
  out_sds_info->data_type = in_sds_info->data_type;
  out_sds_info->data_size = in_sds_info->data_size;
  out_sds_info->fill_val = in_sds_info->fill_val;
  out_sds_info->rank = in_sds_info->rank;
  sprintf(out_sds_info->name, "Subset_%s", in_sds_info->name);
      
  /* compute start and offset for the requested SDS layer */
  
  compute_sds_start_offset(in_sds_info, n, m, &st_c, &offset);
  
  if ((data_in = (void *)calloc(ndata_sm_in, out_sds_info->data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_in in subset_an_sds \n");
  if ((data_out = (void *)calloc(ndata_sm_out, out_sds_info->data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_in in subset_an_sds \n");

  /*   if (open_sds(out_fname, &out_sds_info, 'W') != -1)   */
  if (open_sds((char *)NULL, out_sds_info, 'W') != -1)
    {
      if (attr_val != NULL)
	{
	  if (SDsetattr(out_sds_info->sds_id, "_FillValue", attr_type, 1, (VOIDP)attr_val) == FAIL)
	    fprintf(stderr, " Cannot write sds attribute _FillValue in subset_an_sds \n");
	  free(attr_val);
	}

      if ((data_in != NULL) && (data_out != NULL))
	{
	  for (irow = start_row; irow <end_row +1; irow++)
	    {
	      if ((rank == 2) || (bsq == 0))
		{
		  in_start[0] = irow;
		}
	      else
		{
		  in_start[rank - 2] = irow;
		}
	      
	      if (SDreaddata( in_sds_info->sds_id, in_start, NULL, in_edge, data_in) == FAIL)
		{
		  fprintf(stderr, "Error reading SDS data line in subset_sds: %s", 
			  in_sds_info->name);
		  break;
		}
	      
	      if ((rank == 2) || (bsq == 0))
		{
		  for (icol=start_col, ic=0, i=0; icol < end_col + 1; icol++, i++)
		    {
		      j = icol*total_layer + st_c;
		      for (k=0; k<total_layer; k++, j++, ic++)
			{
			  switch(in_sds_info->data_type)
			    {
			    case 5:
			      ((float32 *)data_out)[ic] = ((float32 *)data_in)[j];
			      break;
			    case 20:
			      ((int8 *)data_out)[ic]= ((int8 *)data_in)[j];
			      break;
			    case 21:
			      ((uint8 *)data_out)[ic] = ((uint8 *)data_in)[j];
			      break;
			    case 22:
			      ((int16 *)data_out)[ic] = ((int16 *)data_in)[j];
			      break;
			    case 23:
			      ((uint16 *)data_out)[ic] = ((uint16 *)data_in)[j];
			      break;
			    case 24:
			      ((int32 *)data_out)[ic] = ((int32 *)data_in)[j];
			      break;
			    case 25:
			      ((uint32 *)data_out)[ic] = ((uint32 *)data_in)[j];
			      break;
			    default: fprintf(stdout,
                                        "HDF datatype " LONG_INT_FMT " not supported", 
                                        in_sds_info->data_type);
			    }
			}
		    }
		}
	      else  /* (bsq == 1) */
		{
		  for (k=0, ic=0; k<total_layer; k++)
		    {
		      for (icol=start_col; icol < end_col + 1; icol++, ic++)
			{
			  j = k*1200 + icol + st_c;
			  
			  switch(in_sds_info->data_type)
			    {
			    case 5:
			      ((float32 *)data_out)[ic] = ((float32 *)data_in)[j];
			      break;
			    case 20:
			      ((int8 *)data_out)[ic]= ((int8 *)data_in)[j];
			      break;
			    case 21:
			      ((uint8 *)data_out)[ic] = ((uint8 *)data_in)[j];
			      break;
			    case 22:
			      ((int16 *)data_out)[ic] = ((int16 *)data_in)[j];
			      break;
			    case 23:
			      ((uint16 *)data_out)[ic] = ((uint16 *)data_in)[j];
			      break;
			    case 24:
			      ((int32 *)data_out)[ic] = ((int32 *)data_in)[j];
			      break;
			    case 25:
			      ((uint32 *)data_out)[ic] = ((uint32 *)data_in)[j];
			      break;
			    default: fprintf(stdout,
                                        "HDF datatype " LONG_INT_FMT " not supported", 
                                        in_sds_info->data_type);
			    }
			}
		    }
		}

	      	      
	      if (SDwritedata(out_sds_info->sds_id, out_start, NULL, out_edge, data_out) == FAIL)
		{
		  fprintf(stderr,
                     "Error writing data line " LONG_INT_FMT " to SDS %s in subset_sds",
                     out_start[0],out_sds_info->name);
		  break;
		}
	      if ((rank > 2) && (bsq == 1))
		++out_start[rank-2];
	      else ++out_start[0];
	    } /* for (irow */
	}
      if (data_in != NULL) free(data_in);
      if (data_out != NULL) free(data_out);
    }
  else 
    {
      if (attr_val != NULL) free(attr_val);
    }
  if (out_sds_info->sds_id != -1) SDendaccess(out_sds_info->sds_id);
}
