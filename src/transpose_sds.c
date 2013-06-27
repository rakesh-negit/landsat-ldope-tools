/****************************************************************************
!C

!File: transpose_sds.c

!Description:
  Transpose one or more SDS in a MODIS Land HDF-EOS data product by rotating 
  the SDS 180 degrees in a clockwise direction. 

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

#include "meta.h"
#include "mfhdf.h"
#include "sds_rw.h"
#include "str_op.h"
#include "qa_tool.h"
#include "main_util.h"
#include "alloc_mem.h"

#define HELP \
"NAME\n" \
"    transpose_sds - Transpose one or more SDS in a MODIS Land HDF-EOS data\n"\
"                    product by rotating the SDS 180 degrees in a clockwise\n"\
"                    direction. \n" \
"\n" \
"SYNOPSIS\n" \
"    transpose_sds help [filename] \n" \
"    transpose [-sds=<SDS_name>1[,<SDS_name2>. . ]] -of=<output_file> \n" \
"              [-meta] filename \n" \
"\n" \
"DESCRIPTION \n" \
"    Transpose one or more SDS in a MODIS Land HDF-EOS data product by \n" \
"    rotating the SDS 180 degrees in a clockwise direction. This tool when \n" \
"    applied on an Aqua L2 granule, data lines in the output granule will \n" \
"    be aligned as seen in a Terra L2 granule. This tool enables qualitative\n"\
"    comparison of MODIS Aqua and Terra Level 2 or Level 1 granules. \n" \
" \n" \
"    This tool supports 2D/3D/4D SDSs. \n" \
" \n" \
"    The tool command arguments can be specified in any order. \n" \
"\n" \
"OPTIONS\n" \
"    -help            Display this help message, If the input filename is \n" \
"                     specified with this option, then the names of all the\n" \
"                     SDS in the file are displayed. \n" \
"    -sds=<SDS_list>  List of SDS to enlarge. SDS names are separated by \n" \
"                     commas with no space. By default all SDSs are\n"\
"                     processed maintaining the input SDS interleaving. \n" \
"                     To process a specific layer of a 3D SDS specify the \n" \
"                     element number of the third dimension as a dot\n"\
"                     extension of the SDS name: sds_name.n (e.g.,\n"\
"                     sur_refl_b02.1 = the layer defined by the 1st element\n"\
"                     of the 3rd dimension of the 3D SDS sur_refl_b02). \n" \
"                     To process a specific layer of a 4D SDS, specify the\n"\
"                     higher dimension element number(s) as a dot extension\n"\
"                     of the SDS name: sds_name.n.m (e.g.,\n"\
"                     Surface_Refl.1.2 = the layer defined by the 1st\n"\
"                     element of the 3rd dimension and the 2nd element of \n" \
"                     the 4th dimension of the 4D SDS Surface_Refl). \n" \
"                     Note that wildcards and ranges of element values may\n"\
"                     be specified as sds_name.* and as sds_name.n1-n2.m\n"\
"                     respectively. \n" \
"    -of=filename     Output filename \n" \
"    -meta            Copy metadata from the input file to the output. \n" \
"    filename         Input filename \n" \
" \n" \
"Examples: \n" \
"    transpose_sds -sds='500m Surface Reflectance Band 1,500m Surface\n"\
"                   Reflectance Band 3,500m Surface Reflectance Band 4'\n"\
"                   MOD09.A2002123.0040.003.2002125174437.hdf \n" \
"                   -of=test.hdf \n" \
"    transpose_sds -sds=Cloud_Mask,Quality_Assurance \n" \
"                  -of=test.hdf MOD35_L2.A2002123.0040.003.2002124023706.hdf\n"\
"                  -meta \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 04/05/2004 \n" \

#define USAGE \
"Usage: \n" \
"    transpose_sds help [filename] \n" \
"    transpose [-sds=<SDS_name>1[,<SDS_name2>. . ]] -of=<output_file> \n" \
"              [-meta] filename \n" \
" \n" \
"OPTIONS \n" \
"    -help            Display this help message, If the input filename is \n" \
"                     specified with this option, then the names of all the\n" \
"                     SDS in the file are displayed. \n" \
"    -sds=<SDS_list>  List of SDS to enlarge. SDS names are separated by \n" \
"                     commas with no space. By default all SDSs are\n"\
"                     processed maintaining the input SDS interleaving. \n" \
"                     To process a specific layer of a 3D SDS specify the \n" \
"                     element number of the third dimension as a dot\n"\
"                     extension of the SDS name: sds_name.n (e.g.,\n"\
"                     sur_refl_b02.1 = the layer defined by the 1st element\n"\
"                     of the 3rd dimension of the 3D SDS sur_refl_b02). \n" \
"                     To process a specific layer of a 4D SDS, specify the\n"\
"                     higher dimension element number(s) as a dot extension\n"\
"                     of the SDS name: sds_name.n.m (e.g.,\n"\
"                     Surface_Refl.1.2 = the layer defined by the 1st\n"\
"                     element of the 3rd dimension and the 2nd element of \n" \
"                     the 4th dimension of the 4D SDS Surface_Refl). \n" \
"                     Note that wildcards and ranges of element values may\n"\
"                     be specified as sds_name.* and as sds_name.n1-n2.m\n"\
"                     respectively. \n" \
"    -of=filename     Output filename \n" \
"    -meta            Copy metadata from the input file to the output. \n" \
"    filename         Input filename \n" \
"\n"


int parse_cmd_transpose_sds(int argc, char **argv, char **sds_names, int *nsds, char *in_fname, char *out_fname, int *m_opt);
void transpose_an_sds(sds_t *in_sds_info, sds_t *out_sds_info);

int main(int argc, char *argv[])
/******************************************************************************
!C

!Description:
  Main function for transpose_sds.

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
  int isds, nsds, m_opt, i;
  sds_t in_sds_info, out_sds_info;
  char **sds_names;
  char in_fname[MAX_PATH_LENGTH];
  char out_fname[MAX_PATH_LENGTH];

  if (argc == 1)
    {
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

  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Error allocating memory for sds_names in transpose_sds: main()\n");
  else 
    {
      if (parse_cmd_transpose_sds(argc, argv, sds_names, &nsds, in_fname, out_fname, &m_opt) == -1)
        fprintf(stderr, "%s\n", USAGE);
      else if (nsds == 0)
	fprintf(stderr, "No SDS to process . . . \n");
      else
	{
	  if ((in_sds_info.sd_id = SDstart(in_fname, DFACC_READ)) == FAIL)
	    {
	      fprintf(stderr, "Cannot open input HDF file %s\n", in_fname);
	      exit(EXIT_FAILURE);
	    }
	  else
	    {
	      if ((out_sds_info.sd_id = SDstart(out_fname, DFACC_CREATE)) == FAIL)
		{
		  fprintf(stderr, "Cannot open output HDF file %s\n", out_fname);
		  exit(EXIT_FAILURE);
		}	      
	      else
		{
		  fprintf(stdout, "Processing file %s\n", in_fname);
		  for (isds=0; isds<nsds; isds++) 
		    {
		      strcpy(in_sds_info.name, sds_names[isds]);
		      in_sds_info.sds_id = -1;
		      if (get_sds_info(in_fname, &in_sds_info) != -1)
			{
			  strcpy(out_sds_info.name, in_sds_info.name);
			  out_sds_info.rank = in_sds_info.rank;
			  for (i=0; i<in_sds_info.rank; i++)
			    out_sds_info.dim_size[i] = in_sds_info.dim_size[i];
			  out_sds_info.data_type = in_sds_info.data_type;
			  out_sds_info.data_size = in_sds_info.data_size;
			  out_sds_info.sds_id = -1;
			  if (open_sds((char *)NULL, &out_sds_info, 'W') != -1)
			    {
			      transpose_an_sds(&in_sds_info, &out_sds_info);
			      write_all_sds_attrs(in_sds_info.sds_id, out_sds_info.sds_id, in_sds_info.nattr);
			      SDendaccess(out_sds_info.sds_id);
			    }
			  SDendaccess(in_sds_info.sds_id);
			}
		    } /* for (isds=0; . . .) */
		  if (m_opt == 1) 
		    write_metadata(in_sds_info.sd_id, out_sds_info.sd_id);
		  SDend(out_sds_info.sd_id);
		}
	      SDend(in_sds_info.sd_id);
	    }
	}
    }
  fprintf(stderr, "Processing done ! \n");
  return 0;
}

int parse_cmd_transpose_sds(int argc, char **argv, char **sds_names, int *nsds, 
			    char *in_fname, char *out_fname, int *m_opt)
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
  in_fname:  input  filename
  out_fname: output filename
  m_opt:     flag of whether the -meta option is selected.

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
  int i, st;

  st = 1;
  *m_opt = *nsds = 0;
  in_fname[0] = '\0';
  out_fname[0] = '\0';
  for (i=0, st=1; i<argc; i++)
  {
    if (is_arg_id(argv[i], "-of") == 0) 
      get_arg_val(argv[i], out_fname);
    else if (is_arg_id(argv[i], "-sds") == 0) 
      get_arg_val_arr(argv[i], sds_names, nsds);
    else if (strcmp(argv[i], "-meta") == 0)
      *m_opt = 1;
    else if (argv[i][0] == '-')
      fprintf(stderr, "Ignoring invalid input argument %s\n", argv[i]);
    else strcpy(in_fname, argv[i]);
  }
  if (in_fname[0] == '\0') { fprintf(stderr, "Missing input file \n"); st = -1; }
  if (out_fname[0] == '\0') { fprintf(stderr, "Missing output file \n"); st = -1; }
  if (st == 1)
  {
    if (*nsds == 0) {
      fprintf(stdout, "No SDS names specified. Processing all SDS\n");
      *nsds = get_sds_names(in_fname, sds_names);
    }
  }
  return st;
}

void transpose_an_sds(sds_t *in_sds_info, sds_t *out_sds_info)
/******************************************************************************
!C

!Description:
  Function transpose a single input SDS.

!Input Parameters:
  in_sds_info: Array of input SDS information structure

!Input/Output Parameters:
  out_sds_info: Array of output SDS information structure

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
  int ic_in, ic_out;
  int irow, icol, iobs;
  int nrows, ncols, nobs;
  int i, bsq, rank, ndata;
  int32 in_edge[4], out_edge[4]; 
  int32 in_start[4], out_start[4]; 
  void *data_in, *data_out;

  fprintf(stdout, "\tProcessing SDS %s\n", in_sds_info->name);
  rank = in_sds_info->rank;
  bsq = ((rank == 2) || (in_sds_info->dim_size[0] < in_sds_info->dim_size[rank-1])) ? 1 : 0;
  if (bsq == 1)
  {
    nrows = in_sds_info->dim_size[rank-2];
    ncols = ndata = in_sds_info->dim_size[rank-1];
    for (i=0, nobs=1; i<rank-2; i++)
    {
      ndata *= in_sds_info->dim_size[i];
      nobs *= in_sds_info->dim_size[i];
    }
  }
  else
  {
    nrows = in_sds_info->dim_size[0];
    ncols = ndata = in_sds_info->dim_size[1];
    for (i=2, nobs=1; i<rank; i++)
    {
      ndata *= in_sds_info->dim_size[i];
      nobs *= in_sds_info->dim_size[i];
    }
  }

  if ((data_in = (void *)calloc(ndata, in_sds_info->data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_in in transpose_an_sds\n");
  if ((data_out = (void *)calloc(ndata, out_sds_info->data_size)) == NULL)
    fprintf(stderr, "Cannot allocate memory for data_out in transpose_an_sds\n");
  if ((data_in != NULL) && (data_out != NULL))
  {
    for (i=0; i<rank; i++)
    {
      in_edge[i] = out_edge[i] = in_sds_info->dim_size[i];
      in_start[i] = out_start[i] = 0;
    }
    if (bsq == 0) in_edge[0] = out_edge[0] = 1;
    else in_edge[rank-2] = out_edge[rank-2] = 1;

    for (irow=0; irow<nrows; irow++)
    {
      if (bsq == 1) { 
	in_start[rank-2] = irow;
	out_start[rank-2] = nrows - 1 - irow;
      }
      else {
	in_start[0] = irow;
	out_start[0] = nrows - 1- irow;
      }
      if (SDreaddata(in_sds_info->sds_id, in_start, NULL, in_edge, data_in) == FAIL)
      {
        fprintf(stderr, "Cannot read data line %d from SDS %s in transpose_an_sds\n", irow, in_sds_info->name);
        break;
      }
      else
      {
	if (bsq == 1)
	{
	  for (iobs=0, ic_out=0; iobs<nobs; iobs++)
	  {
	    ic_in = iobs*ncols + ncols - 1; 
	    switch(in_sds_info->data_type)
	    {
	      case 20: for (icol=0; icol<ncols; icol++, ic_in--, ic_out++)
			 ((int8 *)data_out)[ic_out] = ((int8 *)data_in)[ic_in]; 
		       break;
	      case 21: for (icol=0; icol<ncols; icol++, ic_in--, ic_out++)
			 ((uint8 *)data_out)[ic_out] = ((uint8 *)data_in)[ic_in];
		       break;
	      case 22: for (icol=0; icol<ncols; icol++, ic_in--, ic_out++)
			 ((int16 *)data_out)[ic_out] = ((int16 *)data_in)[ic_in];
		       break;
	      case 23: for (icol=0; icol<ncols; icol++, ic_in--, ic_out++)
			 ((uint16 *)data_out)[ic_out] = ((uint16 *)data_in)[ic_in];
		       break;
	      case 24: for (icol=0; icol<ncols; icol++, ic_in--, ic_out++)
			 ((int32 *)data_out)[ic_out] = ((int32 *)data_in)[ic_in];
		       break;
	      case 25: for (icol=0; icol<ncols; icol++, ic_in--, ic_out++)
			 ((uint32 *)data_out)[ic_out] = ((uint32 *)data_in)[ic_in];
		       break;
	    }
	  } /* for (iobs=0; . . ) */
	}
	else
	{
	  for (icol=0, ic_out=0; icol<ncols; icol++)
	  {
	    ic_in = (ncols - 1 - icol)*nobs;; 
	    switch(in_sds_info->data_type)
	    {
	      case 20: for (iobs=0; iobs<nobs; iobs++, ic_out++, ic_in++)
			 ((int8 *)data_out)[ic_out] = ((int8 *)data_in)[ic_in]; 
		       break;
	      case 21: for (iobs=0; iobs<nobs; iobs++, ic_out++, ic_in++)
			 ((uint8 *)data_out)[ic_out] = ((uint8 *)data_in)[ic_in]; 
		       break;
	      case 22: for (iobs=0; iobs<nobs; iobs++, ic_out++, ic_in++)
			 ((int16 *)data_out)[ic_out] = ((int16 *)data_in)[ic_in]; 
		       break;
	      case 23: for (iobs=0; iobs<nobs; iobs++, ic_out++, ic_in++)
			 ((uint16 *)data_out)[ic_out] = ((uint16 *)data_in)[ic_in]; 
		       break;
	      case 24: for (iobs=0; iobs<nobs; iobs++, ic_out++, ic_in++)
			 ((int32 *)data_out)[ic_out] = ((int32 *)data_in)[ic_in]; 
		       break;
	      case 25: for (iobs=0; iobs<nobs; iobs++, ic_out++, ic_in++)
			 ((uint32 *)data_out)[ic_out] = ((uint32 *)data_in)[ic_in]; 
		       break;
	    }
	  } /* for (icol=0; . .  ) */
	}
        if (SDwritedata(out_sds_info->sds_id, out_start, NULL, out_edge, data_out) == FAIL)
        {
          fprintf(stderr, "Cannot write data line %d to SDS %s in transpose_an_sds\n", irow, out_sds_info->name);
          break;
        }
      }
    } /* for (irow=0; . . .) */
  }
  if (data_in != NULL) free(data_in);
  if (data_out != NULL) free(data_out);
}
