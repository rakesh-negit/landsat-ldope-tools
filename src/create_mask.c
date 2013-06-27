/****************************************************************************
!C

!File: create_mask.c
!Description:

  This file contains the routines for creating a 2D mask SDS 

!Input Parameters: (none)

!Output Parameters: (none)

!Revision History:

    Version 1.0    April, 2004

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
#include <stdlib.h>
#include <string.h>

#include "mfhdf.h"
#include "sds_rw.h"
#include "qa_tool.h"
#include "main_util.h"
#include "alloc_mem.h"
#include "str_op.h"
#include "mask_sds_lib.h"

#define OFF_VAL 0
#define ON_VAL 255

#define HELP \
"NAME \n" \
"    create_mask --- Create an output masking SDS containing two values, an\n"\
"                    user defined 'ON' value at pixels where the masking\n"\
"                    criteria are satisfied and an 'OFF' value elsewhere.\n"\
"                    The mask criteria are specified using relational and\n"\
"                    logical operators applied to the SDS of the same or\n" \
"                    different L2/L3/L4 MODIS Land HDF-EOS data products. \n" \
" \n" \
"SYNOPSIS \n" \
"    create_mask -help [filename] \n" \
"    create_mask -of=<output filename>  \n" \
"                -mask=<mask1>[,AND|OR,<mask2>][,?] [-on=<output ON value>]\n" \
"                [-off=<output OFF value>] \n" \
"      where maskn = <filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
" \n" \
"DESCRIPTION \n" \
"    Create an output masking SDS containing two values: an user defined\n"\
"    'ON' value at pixels where the masking criteria are satisfied and an\n"\
"    'OFF' value elsewhere.\n" \
" \n" \
"    The mask criteria are specified using relational and logical operators\n" \
"    applied to the SDS of the same or different L2/L3/L4 MODIS Land HDF-EOS\n"\
"    data products. SDS(s) used to define the masking criteria must have the\n"\
"    same or lower resolution as the input file SDS(s) to be masked.\n" \
" \n" \
"    The mask criteria are defined by a combination of one or more\n"\
"    individual masks. Each mask is defined by testing SDS bits against bit\n"\
"    values using a relational operator. Testing using a decimal value is\n"\
"    also supported.  Different masks are combined using the logical \"AND\"\n"\
"    or \"OR\" operators.\n" \
" \n" \
"    Masking criteria cannot be applied at pixels where one or more of the\n"\
"    mask SDS(s) have fill values. A mask fill value will be output at these\n"\
"    pixels. The mask fill value may be optionally specified or will be set\n"\
"    to 255 by default. \n" \
" \n" \
"    This tool supports 2D/3D/4D SDSs. Note, only a two dimensional (2D) SDS\n"\
"    or a 2D layer of a 3D/4D SDS can be used to make a mask. \n" \
" \n" \
"    The tool command arguments can be specified in any order.\n" \
" \n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If a filename is\n" \
"                             specified with this option, then the names of\n"\
"                             all the SDSs in the file are displayed.\n" \
"    -of=<filename>           Output filename.\n" \
"    -meta                    Copy metadata from the input file to the\n"\
"                             output file.\n" \
"    -mask=<mask1>[,AND|OR,<mask2>[,..]] \n" \
"       where maskn=< filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
"                             Define a mask from one or more individual\n"\
"                             masks combined using the logical operators\n"\
"                             \"AND\" or \"OR\". \n" \
"                             Each individual mask consists of: \n" \
"    -filename=               MODIS Land product file \n" \
"    -SDSname=                name of an SDS in the file \n" \
"    -bit_numbers=            A list or range of SDS bits \n" \
"    -operator=               relational operator (>, <, <=, >=, ==, !=) \n" \
"    -bit_values=             bit values that are tested against \n" \
"                             The bits in bit_numbers are specified by the\n"\
"                             lower bit followed by the higher bit and the\n"\
"                             bit_values are specified in the reverse\n"\
"                             order.  For example, \n" \
"                             0-2,4==0101 signifies bits 4,2,1,0==0101. \n" \
"                             If the bit_numbers are omitted, then the\n"\
"                             bit_values are parsed as a decimal value.\n"\
"                             This provides a convenient way to refer to a\n"\
"                             specific value, instead of a list of bits.\n"\
"                             For example, -mask=file,SDS,>=200 makes a\n"\
"                             mask where only the SDS values in the file\n"\
"                             greater than or equal to 200 are considered.\n"\
" \n" \
"                             If several masks are combined together then\n"\
"                             '*' may be used in place of the filename\n"\
"                             and/or SDS name to specify the same filename\n"\
"                             and/or SDS name used in the previous mask.\n"\
"                             For example,\n" \
"                                -mask=file1,SDS1,0-2,4==0101,AND,*,*,4-5==10\n"\
" \n" \
"                             To specify a 3D SDS layer write the element\n"\
"                             number of the third dimension as a dot\n"\
"                             extension of the SDS name: sds_name.n (e.g.,\n"\
"                             sur_refl_b02.1 = the layer defined by the 1st\n"\
"                             element of the 3rd dimension of the 3D SDS\n"\
"                             sur_refl_b02). \n" \
" \n" \
"                             To specify a 4D SDS layer write the higher\n"\
"                             dimension element number(s) as a dot extension\n"\
"                             of the SDS name: sds_name.n.m (e.g.,\n"\
"                             Surface_Refl.1.2 = the layer defined by the\n"\
"                             1st element of the 3rd dimension and 2nd\n"\
"                             element of the 4th dimesnsion of the 4D SDS\n"\
"                             Surface_Refl).  \n" \
" \n" \
"    -on=<ON value>           User defined output ON value. \n" \
" \n" \
"    -off=<OFF value>         User defined output OFF value. \n" \
" \n" \
"Examples: \n" \
"    create_mask -of=land_mask.hdf -on=255 -off=0 \n" \
"         -mask=\"MOD09A1.A1996214.h12v04.002.hdf,sur_refl_state_500m,\n"\
"                3-5==001\"\n"\
" \n" \
"    create_mask -of=clear_land.hdf -on=100 -off=0 \n" \
"         -mask=\"MOD09A1.A1996214.h12v04.002.hdf,sur_refl_state_500m,\n"\
"                3-5==001,AND,*,*,01==00\" \n" \
" \n" \
"    create_mask -of=mod35_cloudy_land.hdf -on=1 -off=0 \n" \
"         -mask=\"MOD35_L2.A1996213.1024.002.hdf,\n"\
"               Cloud_Mask.1,1-2==00,AND,*,*,6-7==11\" \n" \
" \n" \
"    create_mask -of=agg_b01_obs1_land_qc.hdf -on=1 -off=0 \n" \
"         -mask=\"MODAGAGG.A1996214.h12v04.001.hdf,Band_QC.1.1,2-5==1100, \n" \
"               AND,*,Aggregate_QC.1.1,3-5==001\" \n" \
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"    Version 1.0, 04/05/2004 \n" \

#define USAGE \
"usage:  \n" \
"    create_mask -help [filename] \n" \
"    create_mask -of=<output filename> \n" \
"                -mask=<mask1>[,AND|OR,<mask2>][,?] [-on=<output ON value>]\n" \
"                [-off=<output OFF value>] \n" \
"       where maskn = <filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
" \n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If a filename is \n" \
"                             specified with this option, then the names of \n"\
"                             all the SDSs in the file are displayed.\n" \
"    -of=<filename>           Output filename.\n" \
"    -meta                    Copy metadata from the input file to the\n"\
"                             output file.\n" \
"    -mask=<mask1>[,AND|OR,<mask2>[,..]]\n" \
"       where maskn=< filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
"                             Define a mask from one or more individual\n"\
"                             masks combined using the logical operators\n"\
"                             \"AND\" or \"OR\".\n" \
"                             Each individual mask consists of: \n" \
"    -filename                MODIS Land product file \n" \
"    -SDSname                 Name of an SDS in the file \n" \
"    -bit_numbers             A list or range of SDS bits \n" \
"    -operator                Relational operator (>, <, <=, >=, ==, !=) \n" \
"    -bit_values              Bit values that are tested against \n" \
"                             The bits in bit_numbers are specified by the\n"\
"                             lower bit followed by the higher bit and the\n"\
"                             bit_values are specified in the reverse order.\n"\
"                             For example, 0-2,4==0101 signifies bits\n"\
"                             4,2,1,0==0101. If the bit_numbers are omitted,\n"\
"                             then the bit_values are parsed as a decimal\n"\
"                             value. This provides a convenient way to refer\n"\
"                             to a specific value, instead of a list of\n"\
"                             bits. For example, -mask=file,SDS,>=200 makes\n"\
"                             a mask where only the SDS values in the file\n"\
"                             greater than or equal to 200 are considered. \n"\
" \n" \
"                             If several masks are combined together then\n"\
"                             '*' may be used in place of the filename\n"\
"                             and/or SDS name to specify the same filename\n"\
"                             and/or SDS name used in the previous mask. For\n"\
"                             example,  \n" \
"                               -mask=file1,SDS1,0-2,4==0101,AND,*,*,4-5==10\n"\
" \n" \
"                             To specify a 3D SDS layer write the element\n"\
"                             number of the third dimension as a dot\n"\
"                             extension of the SDS name: sds_name.n (e.g.,\n"\
"                             sur_refl_b02.1 = the layer defined by the 1st\n"\
"                             element of the 3rd dimension of the 3D SDS\n"\
"                             sur_refl_b02). \n" \
"\n" \
"                             To specify a 4D SDS layer write the higher\n"\
"                             dimension element number(s) as a dot extension\n"\
"                             of the SDS name: sds_name.n.m (e.g.,\n"\
"                             Surface_Refl.1.2 = the layer defined by the\n"\
"                             1st element of the 3rd dimension and 2nd\n"\
"                             element of the 4th dimesnsion of the 4D SDS\n"\
"                             Surface_Refl).  \n" \
" \n" \
"    -on=<ON value>           User defined output ON value. \n" \
" \n" \
"    -off=<OFF value>         User defined output OFF value. \n" \
" \n"


/*****************************************************************************
                            Prototypes. 
*****************************************************************************/

int parse_cmd_create_mask(int argc, char **argv, char *mask_str, char *out_fname, 
		int *on_val, int *off_val);
void generate_mask(char *m_str, char *out_fname, int on_val, int off_val);

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
  int st,  i, k, isds;
  int on_val, off_val;
  char mask_str[MAX_STR_LEN];
  char out_fname[MAX_PATH_LENGTH];
  int32 msds, nattr, rank, dim_size[4];
  int32 sd_id, sds_id, dt;
  char dim_str[MAX_STR_LEN];
  char name[MAX_SDS_NAME_LEN], sds_name[MAX_SDS_NAME_LEN];
		
  if (argc == 1)
    {
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_FAILURE);
    }   
  
  /* display help */
  if ((argc==2) && (strcmp(argv[1], "-help")==0))
  {
    fprintf(stderr, "%s\n", HELP);
    exit(0);
  }

  /* display the SDS information for user */
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
      exit(EXIT_SUCCESS);
    }

  st = parse_cmd_create_mask(argc, argv, mask_str, out_fname, &on_val, &off_val);

  if (st == -1)
    fprintf(stderr, "%s\n", USAGE);
  else if (st != 0)
    generate_mask(mask_str, out_fname, on_val, off_val);
		
	 fprintf(stderr, "Processing done ! \n");
  return 0;
}

int parse_cmd_create_mask(int argc, char **argv, char *mask_str, char *out_fname, int *on_val, 
	int *off_val)

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
  int i, st;
  char on_val_str[20];
  char off_val_str[20];
	char fill_val_str[20];
		
  *on_val = ON_VAL;
  *off_val = OFF_VAL;
  mask_str[0] = '\0';
  out_fname[0] = '\0';
  on_val_str[0] = '\0';
  off_val_str[0] = '\0';
		fill_val_str[0] = '\0';

  st = 1;
  
  for (i=1; i<argc; i++) {
    if (strstr(argv[i], "-mask=") != NULL)
	get_arg_val(argv[i], mask_str);
    else if (strstr(argv[i], "-on=") != NULL) get_arg_val(argv[i], on_val_str);
    else if (strstr(argv[i], "-off=") != NULL) get_arg_val(argv[i], off_val_str);
    else if (strstr(argv[i], "-of=") != NULL) get_arg_val(argv[i], out_fname);
				else if (strstr(argv[i], "-fill=") != NULL) get_arg_val(argv[i], fill_val_str);
    else fprintf(stderr, "Igonoring invalid argument %s\n", argv[i]);
  }
  if ((mask_str[0] == '\0') || (out_fname[0] == '\0')) st = -1;
  if (mask_str[0] == '\0') fprintf(stderr, "Missing input mask_str . . \n");
  if (out_fname[0] == '\0') fprintf(stderr, "Missing output filename . . \n");

  if (st == 1) 
  {
    *off_val = (int)atoi(off_val_str);
    *on_val = (int)atoi(on_val_str);
    if (*off_val == *on_val)
    {
      *off_val = OFF_VAL; *on_val = ON_VAL;
      if ((strlen(on_val_str) > 0) || (strlen(off_val_str) > 0))
        fprintf(stderr, "Output mask SDS values (ON and OFF) invalid: %s %s\n", on_val_str, off_val_str);
      fprintf(stderr, "Ouput SDS values (ON and OFF) set to default: %d %d\n", *on_val, *off_val);
    }
    else if ((*off_val < OFF_VAL) || (*off_val > ON_VAL) || (*on_val < OFF_VAL) || (*on_val > ON_VAL))
    {
      *off_val = OFF_VAL; *on_val = ON_VAL;
      fprintf(stderr, "Output mask SDS values (ON and OFF) invalid: %s %s\n", on_val_str, off_val_str);
      fprintf(stderr, "Ouput SDS values (ON and OFF) set to default: %d %d\n", *on_val, *off_val);
    }
  }
  return st;
}

void generate_mask(char *m_str, char *out_fname, int on_val, int off_val)
/*
!C******************************************************************************

!Function: generate_mask
       
!Description:

  Create user specified masking SDS.

!Input Parameters: 
  m_str        String contains the user specified masking logic.
  out_fname    Output file containing the newly created masking SDS.
  on_val       Output SDS value where mask is true (default 255)
  off_val      Output SDS value whee mask is false (default 0)

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*******************************************************************************/

{
  int st = 1;
  uint8 *mask_row;
  int i, j, p1;
  int rank, len, irow;
  int n_op = 0, i_op, j_op;
  int obs_num[MAX_NUM_OP];
  void *data_qa[MAX_NUM_OP];
  int32 edge[4] = {0, 0, 0, 0};
  int32 start[4] = {0, 0, 0, 0}; 
  
  /* int32 edge[4] = {0, 0};
     int32 start[4] = {0, 0}; */
  
  int32 *data_qa_nadd[MAX_NUM_OP];
  char sdsi_name[MAX_SDS_NAME_LEN];
  char sdsj_name[MAX_SDS_NAME_LEN];
  char num_str[5], **mask_str, **qa_fnames;
  sds_t out_sds_info, qa_sds_nobs_info[MAX_NUM_OP];
  sds_t qa_sds_info[MAX_NUM_OP], qa_sdsc_info[MAX_NUM_OP];
  int val_opt[MAX_NUM_OP], sel_qa_op[MAX_NUM_OP], fqa_l2g[MAX_NUM_OP];
  int res_s[MAX_NUM_OP], res_l[MAX_NUM_OP], rel_op[MAX_NUM_OP];
  unsigned long bit_mask_arr[MAX_NUM_OP], mask_val_arr[MAX_NUM_OP];

  if ((mask_str = (char **)Calloc2D(4*MAX_NUM_OP, 2*MAX_PATH_LENGTH, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for mask_str generate_mask\n");
  else n_op = get_mask_string(m_str, mask_str, val_opt, fqa_l2g);

  if ((qa_fnames = (char **)Calloc2D(MAX_NUM_OP, MAX_PATH_LENGTH, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for qa_fnames generate_mask\n");
  else 
    st = get_parameters(mask_str, n_op, sel_qa_op, qa_fnames, qa_sds_info, bit_mask_arr, 
			mask_val_arr, val_opt, rel_op);
  if ((n_op != -1) && (st != -1))
    st = get_qa_sds_info(qa_fnames, qa_sds_info, qa_sdsc_info, fqa_l2g, n_op);
  if (st != -1) 
  {
    if ((st = get_res_factors(&qa_sds_info[0], qa_sds_info, n_op, res_l, res_s)) != -1)
    {
      strcpy(out_sds_info.name, "Mask_sds");
      out_sds_info.data_type = DFNT_UINT8;
      out_sds_info.sd_id = out_sds_info.sds_id = -1;
      rank = qa_sds_info[0].rank;
      out_sds_info.rank = 2; 
      if ((rank == 2) || (qa_sds_info[0].dim_size[0] > qa_sds_info[0].dim_size[rank-2]))
      {
        out_sds_info.dim_size[0] = qa_sds_info[0].dim_size[0];
        out_sds_info.dim_size[1] = qa_sds_info[0].dim_size[1];
      }
      else
      {
	out_sds_info.dim_size[0] = qa_sds_info[0].dim_size[rank-2];
        out_sds_info.dim_size[1] = qa_sds_info[0].dim_size[rank-1];
      }
      if (open_sds(out_fname, &out_sds_info, 'W') != -1)
      {
        if ((mask_row = (uint8 *)calloc(out_sds_info.dim_size[1], sizeof(uint8))) == NULL)
          fprintf(stderr, "Cannot allocate memory for mask_row in generate_mask\n");
        st = open_qa_sds_nsds((char *)NULL, (sds_t *)NULL, (sds_t *)NULL, (sds_t *)NULL, 1, 
			qa_fnames, qa_sds_info, qa_sdsc_info, qa_sds_nobs_info, fqa_l2g, n_op);
        if (st != -1)
          st = malloc_qa_sds(qa_sds_info, n_op, fqa_l2g, data_qa, data_qa_nadd);

        if ((st != -1) && (mask_row != NULL))
        {
          for (i_op=0; i_op<=n_op; i_op++)
            if (fqa_l2g[i_op] == 1)
            {
              len = (int)strlen(qa_sdsc_info[i_op].name);
              p1 = sd_charpos(qa_sdsc_info[i_op].name, '.', 0);
              sd_strmid(qa_sdsc_info[i_op].name, p1+1, len-p1-1, num_str);
              obs_num[i_op] = (int)atoi(num_str);
            }
            else obs_num[i_op] = 1;

	  edge[0] = 1; edge[1] = out_sds_info.dim_size[1];
          for (irow=0; irow<out_sds_info.dim_size[0]; irow++)
          {
	    start[0] = irow;
            read_qa_sds(qa_sds_info, qa_sdsc_info, qa_sds_nobs_info, n_op, data_qa, data_qa_nadd,
                            irow, res_l, fqa_l2g, obs_num);
            process_mask_data(data_qa, out_sds_info.dim_size[1], qa_sds_info, n_op, sel_qa_op,
		bit_mask_arr, mask_val_arr, rel_op, res_s, mask_row, on_val, off_val, MASK_FILL);
	    if (SDwritedata(out_sds_info.sds_id, start, NULL, edge, (VOIDP)mask_row) == FAIL)
	      fprintf(stderr, "Error writing a line of data to output SDS in generate_mask\n");
          } /* for (irow=0; . . . */
        }

        /* close all SDS and HDF files */
        for (i_op=0; i_op<=n_op; i_op++)
          if (fqa_l2g[i_op] == 1)
          {
            for (j_op=0; j_op<i_op; j_op++)
              /* The following if statement does nothing except create a
               *  compile warning...
              if ((qa_sds_nobs_info[i_op].sd_id == qa_sds_nobs_info[j_op].sd_id) &&
                              (qa_sds_nobs_info[i_op].sds_id == qa_sds_nobs_info[j_op].sds_id)) ;
              */
              if (j_op >= i_op) SDendaccess(qa_sds_nobs_info[i_op].sds_id);
          }

        for (i_op=0; i_op<=n_op; i_op++)
          if (fqa_l2g[i_op] == 1)
          {
            for (j_op=0; j_op<i_op; j_op++)
              if (qa_sds_info[i_op].sd_id == qa_sds_info[j_op].sd_id) break;
            if (j_op >= i_op) free(data_qa_nadd[i_op]);
          }

        close_qa_hdf((char *)NULL, (sds_t *)NULL, qa_fnames, qa_sds_info, n_op);
	if (out_sds_info.sds_id != -1) SDendaccess(out_sds_info.sds_id);
	if (out_sds_info.sd_id != -1) SDend(out_sds_info.sd_id);
        if (mask_row != NULL) free(mask_row);
        if (data_qa[0] != NULL) free(data_qa[0]);
        for (i=1; i<=n_op; i++)
        {
          for (j=0; j<i; j++)
	  {
            if ((p1 = sd_charpos(qa_sds_info[i].name, '.', 0)) == -1)
              strcpy(sdsi_name, qa_sds_info[i].name);
            else
              sd_strmid(qa_sds_info[i].name, 0, p1, sdsi_name);
            if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
              strcpy(sdsj_name, qa_sds_info[j].name);
            else
              sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
            if ((strcmp(qa_fnames[i], qa_fnames[j]) == 0) && (strcmp(sdsi_name, sdsj_name) == 0)) 
	      break;
	  }
          if (j >= i)
          {
            if (data_qa[i] != NULL) free(data_qa[i]);
          }
        } /* for (i=1; . . . */
      }
    }
  }
  Free2D((void **)mask_str);
  Free2D((void **)qa_fnames);
}
