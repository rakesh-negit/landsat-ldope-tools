/****************************************************************************
!C

!File: sds2bin.c

!Description:

  This file contains the routines for dumping an SDS to a binary file. 

!Input Parameters: (none)

!Output Parameters: (none)

!Revision History:

    Version 1.0    October, 2002

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
#include <string.h>
#include <stdlib.h>
#include "mfhdf.h"

#include "str_op.h"
#include "alloc_mem.h"
#include "qa_tool.h"
#include "meta.h"
#include "main_util.h"
#include "sds_rw.h"

typedef unsigned char uchar;

#define HELP \
"NAME \n" \
"    sds2bin - Convert an SDS of a MODIS Land HDF-EOS data product to binary\n"\
"              format.\n" \
" \n" \
"SYNOPSIS \n" \
"    sds2bin -help [filename] \n" \
"    sds2bin -of=<output filename> -sds=<SDSname> filename \n" \
" \n" \
"DESCRIPTION \n" \
"    Convert a user specified SDS from a MODIS Land HDF-EOS data product to\n"\
"    an output binary format file. \n" \
" \n" \
"    The tool supports 2D/3D/4D SDS. \n" \
" \n" \
"    The tool command arguments can be specified in any order. \n" \
" \n" \
"ARGUMENTS \n" \
"    -help [filename]         Display this help message. If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDS in the file are\n"\
"                             displayed. \n" \
"    -of=<filename>           Output filename. \n" \
"    -sds=<SDS name>          SDS to converted. \n" \
"                             To process a specific layer of a 3D SDS\n"\
"                             specify the element number of the third\n"\
"                             dimension as a dot extension of the SDS name:\n"\
"                             sds_name.n (e.g., sur_refl_b02.1 = the layer\n"\
"                             defined by the 1st element of the 3rd\n"\
"                             dimension of the 3D SDS sur_refl_b02). \n" \
"                             To process a specific layer of a 4D SDS,\n"\
"                             specify the higher dimension element number(s)\n"\
"                             as a dot extension of the SDS name: \n" \
"                             sds_name.n.m (e.g., Surface_Refl.1.2 = the\n"\
"                             layer defined by the 1st element of the 3rd\n"\
"                             dimension and the 2nd element of the 4th\n"\
"                             dimension of the 4D SDS Surface_Refl). \n" \
"                             Note that wildcards and ranges of element\n"\
"                             values may be specified as sds_name.* and as\n"\
"                             sds_name.n1-n2.m respectively. \n" \
"    filename                 Input filename. \n" \
" \n" \
"EXAMPLES \n" \
"    sds2bin -sds=sur_refl_b01 -of=sur_refl_b01.img \n" \
"            MOD09A1.A2001145.h20v10.003.2001214125825.hdf \n" \
" \n" \
"    sds2bin -sds=EV_1KM_Emissive -of=ev_1km_emissive.img \n" \
"            MYD021KM.A2002189.2040.003.2002191123800.hdf \n" \
" \n" \
"    sds2bin -sds=\"BRDF_Albedo_Parameters.1.2\" -of=brdf_albedo.img  \n" \
"            MYD43B1.A2002177.h11v11.003.2002210233848.hdf \n" \
"        {Note: This examples converts the layer defined by the 1st element\n"\
"               of the 3rd dimension and the 2nd element of the 4th\n"\
"               dimension of the SDS BRDF_Albedo_parameters to the binary\n"\
"               format file brdf_albedo.img. The output is a 2D binary image}\n"\
" \n" \
"    sds2bin -sds=\"BRDF_Albedo_Parameters.1-7.1\" -of=brdf_albedo.img \n" \
"            MYD43B1.A2002177.h11v11.003.2002210233848.hdf  \n" \
"        {Note: This examples converts the 7 layers defined by the 1st seven\n"\
"               elements of the 3rd dimension and the 1st element of the 4th\n"\
"               dimension of the SDS BRDF_Albedo_parameters to the binary\n"\
"               format file brdf_albedo.img. The 7 layers are output as a\n"\
"               single 3D binary image where the number of elements in the\n"\
"               3rd dimension is 7} \n" \
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 08/08/2002 \n" \

#define USAGE \
"usage: \n" \
"    sds2bin -help [filename] \n" \
"    sds2bin -of=<output filename> -sds=<SDSname> filename \n" \
" \n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDS in the file are\n"\
"                             displayed. \n" \
"    -of=<filename>           Output filename. \n" \
"    -sds=<SDS name>          SDS to converted. \n" \
"                             To process a specific layer of a 3D SDS\n"\
"                             specify the element number of the third\n"\
"                             dimension as a dot extension of the SDS name:\n"\
"                             sds_name.n (e.g., sur_refl_b02.1 = the layer\n"\
"                             defined by the 1st element of the 3rd\n"\
"                             dimension of the 3D SDS sur_refl_b02). \n" \
"                             To process a specific layer of a 4D SDS,\n"\
"                             specify the higher dimension element number(s)\n"\
"                             as a dot extension of the SDS name: \n" \
"                             sds_name.n.m (e.g., Surface_Refl.1.2 = the\n"\
"                             layer defined by the 1st element of the 3rd\n"\
"                             dimension and the 2nd element of the 4th\n"\
"                             dimension of the 4D SDS Surface_Refl). \n" \
"                             Note that wildcards and ranges of element\n"\
"                             values may be specified as sds_name.* and as\n"\
"                             sds_name.n1-n2.m respectively. \n" \
"    filename                 Input filename. \n" \
" \n"



/*****************************************************************************
                            Prototypes. 
*****************************************************************************/

int main(int argc, char **argv);
int parse_cmd_sds2bin(int argc, char **argv, char *in_fname, char *out_fname, char **sds_names);
int sds2bin(char *in_fname, char **sds_names, int isds, char *out_fname);

/*************************************************************************************/

int main(int argc, char **argv)
/*
!C************************************************************************************
       
!Description:
 The sds2bin routine takes 3 argument: SDS name, output file name and input
 HDF file name. It dump the user specified SDS as binary to an output file.

!Input Parameters: (none)

!Output Parameters: 
  (returns)  Completion status: 
             0 - successful completion
	     1 - error exit
    
!Revision History: (none)

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
**************************************************************************************/
{
  int status = 0, i, k; 
  int isds;
  int32 msds, nattr, rank, dim_size[4];
  int32 sd_id, sds_id, dt;
  char in_fname[MAX_PATH_LENGTH];
  char out_fname[MAX_PATH_LENGTH];
  char dim_str[MAX_STR_LEN];  
  char **tmp_sds_names, **sds_names;
  char name[MAX_SDS_NAME_LEN], sds_name[MAX_SDS_NAME_LEN];

  if (argc == 1)
    {
      status = -1;
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      return EXIT_FAILURE;
    }   
  
  /* display help */
  if ((argc == 2) && (strcmp(argv[1], "-help") == 0))
    {
      fprintf(stderr, "%s\n", HELP);
      return 0;
    }
  
  /*  Display SDS names of input HDF file */
  
  if ((argc>=3) && (strcmp(argv[1], "-help")==0))
    {
      /* Print all the SDS names in the input file */
      msds = 0;
      for (i=2; i<argc; i++)
	{
	  if (argv[i][0] != '-')
	    {
	      if ((sd_id = SDstart(argv[i], DFACC_READ)) == FAIL)
		fprintf(stderr, "Cannot open the HDF file %s\n", argv[i]);
	      else
		{
		  if (SDfileinfo(sd_id, &msds, &nattr) == FAIL)
		    fprintf(stderr, "Error reading information for HDF file %s\n",argv[i]);
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
      return EXIT_SUCCESS;
    }

  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for sds_names in sds2bin (main)\n");
      return EXIT_FAILURE;
    }
  if ((tmp_sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for tmp_sds_names in sds2bin (main)\n");
      return EXIT_FAILURE;
    }
  
  if ((tmp_sds_names != NULL) &&(sds_names !=NULL))
    {
      status = parse_cmd_sds2bin(argc, argv, in_fname, out_fname, sds_names);
      
      if (status == -1)
	{
	  fprintf(stderr, "%s\n", USAGE);
	  return EXIT_FAILURE;
	}
      else if (status != 0)
	{	 
	  isds = 1;
	  strcpy(tmp_sds_names[0], sds_names[0]);
	  update_nd_sdsnames(tmp_sds_names, &isds, in_fname);
	  status = sds2bin(in_fname, tmp_sds_names, isds, out_fname);
	}
    }
  Free2D((void **)tmp_sds_names); 
  Free2D((void **)sds_names);
  fprintf(stderr, "Processing done ! \n");
  return status;
}	  

int sds2bin(char *in_fname, char **sds_names, int sds_cnt, char *out_fname )
/*
!C******************************************************************************

!Function: sds2bin
       
!Description:

  Convert user specified SDS to binary output. Original interleaving is preserved.

!Input Parameters: 
  in_fname     Input HDF filename.
  sds_names    User specified sds names. Please note that for 3D or 4D sds, if the
               user specified layer range, the sds_names is in the format as 
	       sds_name.n.m
  sds_cnt      Number of input sds layers.
  out_fname    Output filename. Please note output file is in binaray format.

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/

{
  int bsq;
  int n, m, isds;
  size_t i, ic;
  int rank,st_c, offset, k;
  int nrow = 0, irow;
  size_t ndata_in = 0;
  size_t ndata_out = 0;
  FILE *fp;  
  int32 in_edge[4] = {0,0,0,0};  
  int32 in_start[4] = {0,0,0,0};   
  void *data_in, *data_out; /*, *attr_val; */
  sds_t in_sds_info;
  
  in_sds_info.sd_id = in_sds_info.sds_id = -1;
  
  if ((fp = fopen(out_fname, "w")) == NULL)
    {
      fprintf(stderr, "Cannot create output file: %s\n", out_fname);  
      exit(EXIT_FAILURE);
    }
  
  if ((in_sds_info.sd_id = SDstart(in_fname, DFACC_READ)) == FAIL)
    {
      fprintf(stderr, "Cannot open input HDF file %s \n", in_fname);
      exit(EXIT_FAILURE);
    }
  
  for (isds=0; isds<sds_cnt; isds++)
  {
    get_sdsname_dim(sds_names[isds], in_sds_info.name, &n, &m);
    fprintf(stderr, "	Processing SDS %s\n", sds_names[isds]);  
    in_sds_info.sds_id = -1;

    if (get_sds_info((char *)NULL, &in_sds_info) == -1)
    {
      if (in_sds_info.sds_id != -1) 
        SDendaccess(in_sds_info.sds_id);
     
      continue;
    } 

    rank = in_sds_info.rank;
    bsq = ((rank == 2) || (in_sds_info.dim_size[0] < in_sds_info.dim_size[rank-1])) ? 1 : 0;

    if (rank == 2)
      ndata_in = ndata_out = in_sds_info.dim_size[1];
    else if (rank > 2)
    {
      if (in_sds_info.dim_size[0] > in_sds_info.dim_size[rank-1])
      {
        ndata_in = in_sds_info.dim_size[1];
        for (k=2; k<rank; k++)
          ndata_in *= in_sds_info.dim_size[k];
        ndata_out = ((n == -1) && (m == -1)) ? ndata_in : (size_t)in_sds_info.dim_size[1];
      }
      else
      {
	ndata_in = in_sds_info.dim_size[rank-1];
        for (k=0; k<rank-2; k++)
          ndata_in *= in_sds_info.dim_size[k];
        ndata_out = ((n == -1) && (m == -1)) ? ndata_in : (size_t)in_sds_info.dim_size[rank-1];
      }
    }

    if ((data_in = (void *)calloc(ndata_in, in_sds_info.data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_in in sds2bin\n");
    
    if ((data_out = (void *)calloc(ndata_out, in_sds_info.data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_out in sds2bin\n");
    
    if ((data_in != NULL) && (data_out != NULL))
      {
	if (rank == 2)
	  {
	    in_edge[0] = 1;
	    in_edge[1] = in_sds_info.dim_size[1];
	    nrow = in_sds_info.dim_size[0];
	  }
	else if (rank > 2)
	  {
	    nrow = (bsq == 1) ? in_sds_info.dim_size[rank-2] : in_sds_info.dim_size[0];
	    for (k=0; k<rank; k++)
	      in_edge[k] = in_sds_info.dim_size[k];
	    if ((n == -1) && (m == -1))
	      {
		if (bsq == 0) in_edge[0] = 1;
		else in_edge[rank-2] = 1;
	      }
	    else
	      {
		if (bsq == 0) in_edge[0] = 1;
		else in_edge[rank-2] = 1;
	      }
	  }
	
	compute_sds_start_offset(&in_sds_info, n, m, &st_c, &offset);
	
	
	for (irow=0; irow<nrow; irow++)
	  {
	    if ((rank == 2) || (in_sds_info.dim_size[0] > in_sds_info.dim_size[rank-1]))
	      in_start[0] = irow;
	    else
	      in_start[rank-2] = irow;
	    
	    if (SDreaddata(in_sds_info.sds_id, in_start, NULL, in_edge, data_in) == FAIL)
	      {
		fprintf(stderr, "Cannot read data line from SDS %s in sds2bin\n", in_sds_info.name);
		break;
	      }
	    
	    for (i=0, ic = st_c; i<ndata_out; ic += offset, i++)
	      {
		switch(in_sds_info.data_type)
		  {
		  case 5:
		    ((float32 *)data_out)[i] = ((float32 *)data_in)[ic];
		    break;
		  case 6:
		    ((float64 *)data_out)[i] = ((float64 *)data_in)[ic];
		    break;
		  case 20:
		    ((int8 *)data_out)[i] = ((int8 *)data_in)[ic];
		    break;
		  case 21:
		    ((uint8 *)data_out)[i] = ((uint8 *)data_in)[ic];
		    break;
		  case 22:
		    ((int16 *)data_out)[i] = ((int16 *)data_in)[ic];
		    break;
		  case 23:
		    ((uint16 *)data_out)[i] = ((uint16 *)data_in)[ic];
		    break;
		  case 24:
		    ((int32 *)data_out)[i] = ((int32 *)data_in)[ic];
		    break;
		  case 25:
		    ((uint32 *)data_out)[i] = ((uint32 *)data_in)[ic];
		    break;
		  default: fprintf(stdout,
                             "HDF datatype " LONG_INT_FMT " not supported",
                             in_sds_info.data_type);
		  }
	      }
	    
	    if (fwrite(data_out, in_sds_info.data_size, ndata_out, fp ) != ndata_out)
	      {
		fprintf(stderr, "Error writing data to file %s\n", out_fname);
	      }
	  }
      }
    if (data_in != NULL) free(data_in);
    if (data_out != NULL) free(data_out);
    SDendaccess(in_sds_info.sds_id);
  }	/* for (isds = 0; . . . */
  
  SDend(in_sds_info.sd_id);
  
  fclose(fp);
  return 0; 
}

int parse_cmd_sds2bin(int argc, char **argv, char *in_fname, char *out_fname, 
                      char **sds_names)
/*
!C*********************************************************************************

!Function: parse_cmd_sds2bin
       
!Description:
   function parse_cmd_sds2bin parse the argument string to the corresponding
   input filename (in_fname), output filename (out_fname) and SDS name
   (sds_name).
   It print out a error message if there the input file, output file or 
   SDS name is missing. 

!Input Parameters: 
  argc      input argument count
  argv      input argument vector

!Input/output Parameters:

  in_fname  Input HDF file name.
  out_fname Output file name.
  sds_names User input sds_name. Should only input one sds name since sds2bin is 
            intend to write one SDS(all or some layers) out to output file.

!Output Parameters:

  ret   :    1  - successful completion
	     -1 - error exit   

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/
{
  int i;
  int ret;            /* return value */
  int sds_cnt;
  
  in_fname[0] = '\0';
  out_fname[0] = '\0';
  
  for (i=1, ret=1; i<argc; i++)
    {
      if (is_arg_id(argv[i], "-sds") == 0) 
	{
	  sds_cnt = 0;
	  get_arg_val_arr(argv[i], sds_names, &sds_cnt);
	  if (sds_cnt != 1)
	    {
	      fprintf(stderr, "Error in input SDS name. Should input one and only one SDS names. \n");
	      exit(EXIT_FAILURE);
	    }
	}
      else if (is_arg_id(argv[i], "-of") == 0) get_arg_val(argv[i], out_fname);
      else if (argv[i][0] == '-') fprintf(stderr, "Ignoring invalid option %s\n", argv[i]);
      else 
	{
	  strcpy(in_fname, argv[i]);
	}
    }
  if (strlen(in_fname) <= 0) 
    { 
      fprintf(stderr, "Missing input file \n"); 
      ret = -1; 
    }
  if (strlen(out_fname) <= 0) 
    { 
      fprintf(stderr, "Missing output file \n"); 
      ret = -1; 
    }
  if (strlen(sds_names[0]) <= 0) 
    { 
      fprintf(stderr, "Missing input SDS name \n"); 
      ret = -1; 
    }
  return ret;
}

