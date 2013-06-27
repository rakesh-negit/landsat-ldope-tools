/****************************************************************************
!C

!File: unpack_sds_bits.c

!Description:
  Contains routines for unpacking selected SDS bit fields

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original May 1998. Version 1.1
  Modified June 1998 Version 1.1 (process multiple input files)
  Modified Dec 1998 Version 1.1 (update for 3D/4D sds)
  Modifeid Jan 1999 Version 1.1 (update for range bits)
  Modified June 1999 Version 1.2 (support specific value of higher dimension)
  Modified September 1999 (included -meta option)
  Modified December 1999 (fixed bug: probem with 3D SDS in bit sequential)
  Modified March 2012 Version 2.0 by Xiaoping Zhang, Sigma Space Corp.
  	(to add fill value option on command line)
  Modified Jun 2012 Version 2.1  by Xiaoping Zhang, Sigma Space Corp.
  	(fixed bug: probem with like MOD13C1/MOD13C2 SDS names containing "." character)
  Modified Jun 2012 Version 2.2  by Xiaoping Zhang, Sigma Space Corp.
  	(fixed bug: problem with the dimension layer number(n and m) not set to -1 initially)

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

!Design Notes: (none)

!END
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "isoc.h"
#include "mfhdf.h"
#include "qa_tool.h"
#include "sds_rw.h"
#include "str_op.h"
#include "meta.h"
#include "alloc_mem.h"
#include "main_util.h"

#define HELP \
"NAME \n" \
"    unpack_sds_bits - unpack bit fields of sds data sets \n" \
" \n" \
"SYNOPSIS \n" \
"    unpack_sds_bits -help\n" \
"    unpack_sds_bits [-sds=SDSname1[,SDSname2...]] -of=output_file\n"\
"                    -bn=Bitnumbers -meta filename\n"\
"                    -fill_value=<fill value> \n" \
" \n" \
"    echo [-sds=SDSname1[,SDSname2...]] filename -of=output_file \n" \
"         -bn=Bitnumbers -meta -fill_value=<fill value> | unpack_sds_bits \n" \
" or unpack_sds_bits < argument_file \n" \
" \n" \
"DESCRIPTION \n" \
"    The MODIS Land per-pixel QA information and other information, such as\n"\
"    for example, the land-sea mask, logical criteria used by the algorithm,\n"\
"    and cloud state, are stored in an efficient bit encoded manner. This\n"\
"    tool decodes requested bit fields and writes them to an output HDF\n"\
"    file. The output SDS data type is uint8, uint16 or uint32 depending on\n"\
"    the number of unpacked bits. Note that the unpacked bits are stored in\n"\
"    the least significant bits of the output SDS. Refer to the MODIS\n"\
"    product file specifications for information on which bits to select for\n"\
"    unpacking. \n" \
" \n" \
"    This tool supports 2D/3D/4D SDS.\n" \
" \n" \
"    The tool command arguments can be specified in any order. \n" \
" \n" \
"OPTIONS \n" \
"    -help            Display this help message\n" \
"    -of=filename     Output file \n" \
"    -sds=<SDS list>  List of SDS to be unpacked (separated by commas). If\n"\
"                     the SDS is 3D enter each SDS name in the list as\n"\
"                     sdsname.n and if 4D enter sds_name.n.m where n and m\n"\
"                     are the specific index (1-based) of the higher\n"\
"                     dimension to unpack. sds_name.*.m will unpack all the\n"\
"                     layers in the 3rd dimesion for layer m in the fourth\n"\
"                     dimension. range of dimession value can be specified\n"\
"                     as sds_name.n1-n2.m1-m2. Note that dimension values\n"\
"                     cannot be separated by comma. Default is all SDS and\n"\
"                     all dimensions; and retain original interleaving. \n" \
"    -bn=<Bitnumbers> List of bit numbers separated by commas. Range of \n" \
"                     continuous bit numbers are specified by '-'. \n" \
"                     (e.g. -bn=4,10 -bn=4-8  -bn='3-5, 9-11, 15, 18') \n" \
"    -meta            Copy metadata from input file to output \n" \
"    -fill_value      Specify the fill value. When used with user specified\n"\
"                     value for the fill value this option will override the\n"\
"                     fill value in the input HDF file.\n" \
"    filename         Input filename. \n" \
" \n" \
"Examples: \n" \
"    unpack_sds_bits -sds=Cloud_Mask.1 -bn=1-2 -of=cloud_bits.hdf \n" \
"                    MYD35_L2.A2002189.2040.003.2002191125354.hdf \n" \
"        {Note: This example unpacks bits 1 and 2 of the layer defined by\n"\
"               the 1st element of the 3rd dimension of the SDS Cloud_Mask.} \n"\
" \n" \
"    unpack_sds_bits -sds=sur_refl_qc_500m -bn=10-13 -of=srefl_qc_bits.hdf\n" \
"                     MOD09A1.A1996222.h12v04.8days.002.hdf\n" \
" \n" \
"    unpack_sds_bits -sds=sur_refl_qc_500m -bn=10-13,14-17,18-21\n"\
"                    -of=srefl_qc_bits.hdf\n" \
"                    MOD09A1.A1996222.h12v04.8days.002.hdf\n" \
" \n" \
"    unpack_sds_bits -sds=\"most confident detected fire\" -bn=0-3\n"\
"                    -of=fire_bits.hdf  \n" \
"                    MOD14A1.A2002185.h30v11.003.2002204204451.hdf  \n" \
" \n" \
"    unpack_sds_bits -sds=Band_QC.1.1-2 -bn=2-5,6-7 -of=qc_obs1_b1b2.hdf\n" \
"                    MODAGAGG.A1996214.h12v04.001.hdf\n" \
" \n" \
"    unpack_sds_bits -sds=\"Band_QC.1.*\" -sds=Band_QC.1.2 -bn=2-5,6-7\n"\
"                    -of=qc_obs1_b1b2.hdf\n" \
"                     MODAGAGG.A1996214.h12v04.001.hdf\n" \
" \n" \
"    unpack_sds_bits -sds=Cloud_Mask.1 -bn=1-2 -of=clouds.hdf\n" \
"                    MOD35_L2.A1997223.1034.002.1999141233243.hdf\n" \
" \n" \
"    unpack_sds_bits -sds=QC_Day\n"\
"                    -of=MOD11A2.A2012065.h17v07.005.DayMQA_fillOption.hdf\n"\
"                    -bn=0-1 -fill_value=255\n" \
"                    MOD11A2.A2012065.h17v07.005.2012075043424.hdf\n" \
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Y. Zhang \n" \
"    Documentation: S.Devadiga and D. Roy \n" \
" \n" \
"Version 2.2, 06/15/2012\n" \
"Please report problems to Sadashiva Devadiga (devadiga@ltpmail.gsfc.nasa.gov)\n" 

#define USAGE \
"usage:	unpack_sds_bits [-sds=SDSname1[,SDSname2...]] -of=output_file\n" \
"                       [-bn=Bitnumbers] -meta -fill_value=<fill value>\n"\
"                       filename\n" \
" \n" \
"OPTIONS \n" \
"    -help            Display this help message\n" \
"    -of=filename     Output file \n" \
"    -sds=<SDS list>  List of SDS to be unpacked (separated by commas). If\n"\
"                     the SDS is 3D enter each SDS name in the list as\n"\
"                     sdsname.n and if 4D enter sds_name.n.m where n and m\n"\
"                     are the specific index (1-based) of the higher\n"\
"                     dimension to unpack. sds_name.*.m will unpack all the\n"\
"                     layers in the 3rd dimesion for layer m in the fourth\n"\
"                     dimension. range of dimession value can be specified\n"\
"                     as sds_name.n1-n2.m1-m2. Note that dimension values\n"\
"                     cannot be separated by comma. Default is all SDS and\n"\
"                     all dimensions; and retain original interleaving. \n" \
"    -bn=<Bitnumbers> List of bit numbers separated by commas. Range of \n" \
"                     continuous bit numbers are specified by '-'. \n" \
"                     (e.g. -bn=4,10 -bn=4-8  -bn='3-5, 9-11, 15, 18') \n" \
"    -meta            Copy metadata from input file to output \n" \
"    -fill_value      Specify the fill value. When used with user specified\n"\
"                     value for the fill value this option will override the\n"\
"                     fill value in the input HDF file.\n" \
"    filename         Input filename. \n" \
" \n"


/*******************************************************************************
                            Static variable.
********************************************************************************/

static unsigned int BIT[] = {
  0x1,        0x3,        0x7,        0xf,
  0x1f,       0x3f,       0x7f,       0xff,
  0x1ff,      0x3ff,      0x7ff,      0xfff,
  0x1fff,     0x3fff,     0x7fff,     0xffff,
  0x1ffff,    0x3ffff,    0x4ffff,    0xfffff,
  0x1fffff,   0x3fffff,   0x7fffff,   0xffffff,
  0x1ffffff,  0x3ffffff,  0x7ffffff,  0xfffffff,
  0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff
};

/*******************************************************************************
                            Prototypes.
********************************************************************************/

void unpack_sds(char *in_fname, char **sds_names, char **bn_str, int bn_cnt,
	int32 sd_id, int sds_cnt, int if_cnt, int m_opt, char *fillVal);
void unpack_bits(void *data_in, void *data_out, sds_t *sds_info, int bn_cnt, int *bn_arr,
	int ndata, int st_c, int offset);
int parse_cmd_unpack_sds_bits(int argc, char **argv, int *if_cnt, char **sds_names, 
	int *sds_cnt, char **bn_str, int *bn_cnt, char *out_fname, int *m_opt,  char *fillVal);
int getSdsnameDim(char *sdsname_str, int sd_id, int *n, int *m); 

int main(int argc, char **argv)
{
  int m_opt;
  int i, status;
  int stdin_argc = 0;
  int nsds, sds_cnt;
  int fid, if_cnt, bn_cnt;
  int32 out_sd_id;
  char **bn_str, **sds_names;
  char **stdin_argv = NULL, **tmp_sds_names;
  char out_fname[MAX_PATH_LENGTH], in_fname[MAX_PATH_LENGTH];
  char fillVal[10];

  int k, isds;
  int32 sd_id, sds_id, dt;
  int32 msds, nattr, rank, dim_size[4];
  char dim_str[MAX_STR_LEN];
  char name[MAX_SDS_NAME_LEN], sds_name[MAX_SDS_NAME_LEN];

  if ((argc==2) && ((strcmp(argv[1], "-h")==0) || (strcmp(argv[1], "-help")==0)))
  {
    fprintf(stderr, "%s\n", HELP);
    exit(0);
  }

  if ((argc>=3) && ((strcmp(argv[1], "-h")==0) || (strcmp(argv[1], "-help")==0)))
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

  if (argc == 1)
  {
    if ((stdin_argv = (char **)Calloc2D(MAX_NUM_PARAM, MAX_PARAM_LENGTH, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for stdin_argv in data_fedex(main)\n");
      exit(1);
    }
    if (parse_stdin(&stdin_argc, stdin_argv) == -1)
    {
      Free2D((void *)stdin_argv);
      exit(0);
    }
  }

  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
      fprintf(stderr, "Cannot allocate memory for sds_names in unpack_sds_bits(main)\n");
  if ((tmp_sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
      fprintf(stderr, "Cannot allocate memory for tmp_sds_names in unpack_sds_bits(main)\n");
  if ((bn_str = (char **)Calloc2D(MAX_NUM_BITS, MAX_NUM_BITS, sizeof(char))) == NULL)
      fprintf(stderr, "Cannot allocate memory for bn_str in unpack_sds_bits(main)\n");
  if ((sds_names != NULL) && (tmp_sds_names != NULL) && (bn_str != NULL))
    {
    if (argc == 1)
      status = parse_cmd_unpack_sds_bits(stdin_argc, stdin_argv, &if_cnt, sds_names, &sds_cnt, 
		bn_str, &bn_cnt, out_fname, &m_opt, fillVal);
    else 
      status = parse_cmd_unpack_sds_bits(argc, argv, &if_cnt, sds_names, &sds_cnt, 
        bn_str, &bn_cnt, out_fname, &m_opt, fillVal);
  if (status == -1)
    fprintf(stderr, "%s\n", USAGE);
  else if (status != 0)
    {
      if ((out_sd_id = SDstart(out_fname, DFACC_CREATE)) == FAIL)
        fprintf(stderr, "Cannot open output HDF file %s\n", out_fname);
      else 
	{
      if (argc != 1) stdin_argc = argc;
      for (fid=1; fid<stdin_argc; fid++)
    {
	  if (argc == 1) 
	    strcpy(in_fname, stdin_argv[fid]);
	  else 
        strcpy(in_fname, argv[fid]);
	  if (in_fname[0] == '-') continue;
	  for (i=0; i<sds_cnt; i++)
         strcpy(tmp_sds_names[i], sds_names[i]);
         fprintf(stdout, "%d %s\n", sds_cnt, tmp_sds_names[0]);
	     if ((sds_cnt == 1) && (strcmp(tmp_sds_names[0], "all") == 0))
		   nsds = get_sds_names(in_fname, tmp_sds_names);
	     else nsds = sds_cnt;
	     fprintf(stderr, "\nProcessing input hdf file: %s\n", in_fname);
	     fprintf(stderr, "---------------------------------------------------------------\n");
	     update_nd_sdsnames(tmp_sds_names, &nsds, in_fname);
	     unpack_sds(in_fname, tmp_sds_names, bn_str, bn_cnt, out_sd_id, nsds, if_cnt, m_opt, fillVal); 
	    }
	  SDend(out_sd_id);
	  }
    }
  }
  Free2D((void **)bn_str);
  Free2D((void **)sds_names);
  Free2D((void **)tmp_sds_names);
  if (argc == 1)
    Free2D((void **)stdin_argv);
  fprintf(stderr, "Done!\n");
  return 0;
}

int parse_cmd_unpack_sds_bits(int argc, char **argv, int *if_cnt, char **sds_names, 
		int *sds_cnt, char **bn_str, int *bn_cnt, char *out_fname, int *m_opt, char *fillVal)
{
  int i, status;

  status = 1;
  out_fname[0] = fillVal[0] = '\0';
  *bn_cnt = *if_cnt = 0;
  *m_opt = *sds_cnt = 0;

  for (i=1; i<argc; i++)
  {
    if (is_arg_id(argv[i], "-sds=") == 0) 
      get_arg_val_arr(argv[i], sds_names, sds_cnt);

    else if (is_arg_id(argv[i], "-sds") == 0) ;

    else if ((is_arg_id(argv[i], "-of=") == 0) || (is_arg_id(argv[i], "-o=") == 0)) 
      get_arg_val(argv[i], out_fname);

    else if ((is_arg_id(argv[i], "-bit=") == 0) || (is_arg_id(argv[i], "-bn=") == 0))
      get_arg_val_arr(argv[i], bn_str, bn_cnt);

    else if (strcmp(argv[i], "-meta") == 0) *m_opt = 1;

    else if (is_arg_id(argv[i], "-fill_value=") == 0) get_arg_val(argv[i], fillVal);

    else if (argv[i][0] == '-')
      fprintf(stderr, "Unknown option %s\n", argv[i]);

    else ++*if_cnt;
  }

  if ((out_fname[0] == '\0') || (*bn_cnt == 0) || (*if_cnt == 0)) status = -1;
  if (*bn_cnt == 0) fprintf(stderr, "No bit numbers specified \n");
  if (*if_cnt == 0) fprintf(stderr, "No input filename specified\n");
  if (out_fname[0] == '\0') fprintf(stderr, "No output filename specified \n");

  if (status != -1)
  {
    if (*sds_cnt == 0)
    {
      *sds_cnt = 1;
      strcpy(sds_names[0], "all");
      fprintf(stderr, "No SDS name input. Reading all SDS . . \n");
    }
  }
  return status;
}

void unpack_sds(char *in_fname, char **sds_names, char **bn_str, int bn_cnt, int32 out_sd_id, 
	int sds_cnt, int if_cnt, int m_opt, char *fillVal)
{
  int bsq;
  int i, j, n = -1, m = -1;
  int st_c, offset;
  int rank, out_rank;
  int out_dim_size[4];
  int nrow = 0, irow, isds;
  int ndata_in = 0, ndata_out = 0;
  int max_nbit, *nbit, **bn_arr;
  int32 attr_type, attr_cnt;
  int32 in_edge[4] = {0, 0, 0, 0}; 
  int32 in_start[4] = {0, 0, 0, 0};
  int32 out_edge[4] = {0, 0, 0, 0}; 
  int32 out_start[4] = {0, 0, 0, 0};
  char tmp_str[MAX_PATH_LENGTH];
  void *data_in, *data_out, *attr_val = NULL;
  sds_t in_sds_info, out_sds_info[MAX_NUM_SDS];
  char attr_name[20];

/*  int p1, p2, p11, len;
  char nd_ext[20], md_ext[20];
  char sds_name[MAX_SDS_NAME_LEN]; */

  if ((bn_arr = (int **)Calloc2D(MAX_NUM_BITS, MAX_NUM_BITS, sizeof(int))) == NULL)
    fprintf(stderr, "Cannot allocate memory for bn_arr in unpack_sds\n");
  if ((nbit = (int *)calloc(bn_cnt, sizeof(int))) == NULL)
    fprintf(stderr, "Cannot allocate memory for nbit in unpack_sds\n");
  if ((bn_arr != NULL) && (nbit != NULL))
     for (i=0; i<bn_cnt; i++)
	    get_bit_num_arr(bn_str[i], &nbit[i], bn_arr[i]);
  else
  {
    Free2D((void **)bn_arr);
    if (nbit != NULL) free(nbit);
    return;
  }
  for (i=1, max_nbit=nbit[0]; i<bn_cnt; i++)
    if (max_nbit < nbit[i])
      max_nbit = nbit[i];
 
  if ((sds_cnt == 1) && (strcmp(sds_names[0], "all") == 0))
    sds_cnt = get_sds_names(in_fname, sds_names);
  if (sds_cnt <= 0) return;
  if ((in_sds_info.sd_id = SDstart(in_fname, DFACC_READ)) == FAIL)
  {
    fprintf(stderr, "Cannot open input HDF file %s\n", in_fname);
    return;
  }
  if (if_cnt > 1)
  {
    strcpy(tmp_str, in_fname);
    rm_path(tmp_str);
    tmp_str[24] = '\0';
  }

  for (isds=0; isds<sds_cnt; isds++)
  {
//  get_sdsname_dim(sds_names[isds], in_sds_info.name, &n, &m);
    fprintf(stderr, "	Processing SDS %s\n", sds_names[isds]);  
    in_sds_info.sds_id = -1;
    strcpy(in_sds_info.name, sds_names[isds]);
    if (get_sds_info((char *)NULL, &in_sds_info) == -1)
    {
      if (in_sds_info.sds_id != -1) 
        SDendaccess(in_sds_info.sds_id);
      continue;
    }

    getSdsnameDim(sds_names[isds], in_sds_info.sd_id, &n, &m);

    if (fillVal[0] != '\0')       
      switch(in_sds_info.data_type)
      {
	case  5: in_sds_info.fill_val = (long)(float32)atof(fillVal); break;
	case  6: in_sds_info.fill_val = (long)(float64)atof(fillVal); break;
	case 20: in_sds_info.fill_val = (int8)atoi(fillVal); break;
	case 21: in_sds_info.fill_val = (uint8)atoi(fillVal); break;
	case 22: in_sds_info.fill_val = (int16)atoi(fillVal); break;
	case 23: in_sds_info.fill_val = (uint16)atoi(fillVal); break;
	case 24: in_sds_info.fill_val = (int32)atoi(fillVal); break;
	case 25: in_sds_info.fill_val = (uint32)atoi(fillVal); break;
      } 
    else
    {  
      if ((attr_val = (void *)get_sds_attr(in_sds_info.sds_id, "_FillValue", &attr_type,
		&attr_cnt)) == NULL)
      fprintf(stderr, "Attribute "_FillValue" not defined for output\n");
    else
      switch(attr_type)
      {
          case 20: in_sds_info.fill_val = ((int8 *)attr_val)[0]; break;
          case 21: in_sds_info.fill_val = ((uint8 *)attr_val)[0]; break;
          case 22: in_sds_info.fill_val = ((int16 *)attr_val)[0]; break;
          case 23: in_sds_info.fill_val = ((uint16 *)attr_val)[0]; break;
          case 24: in_sds_info.fill_val = ((int32 *)attr_val)[0]; break;
          case 25: in_sds_info.fill_val = ((uint32 *)attr_val)[0]; break;
      }
    }
    rank = in_sds_info.rank;
    bsq = ((rank == 2) || (in_sds_info.dim_size[0] < in_sds_info.dim_size[rank-1])) ? 1 : 0;
    out_rank = ((n == -1) && (m == -1)) ? in_sds_info.rank : 2;
    if (in_sds_info.rank == out_rank)
      for (j=0; j<rank; j++)
        out_dim_size[j] = in_sds_info.dim_size[j];
    else
    {
      if (bsq == 1)
      {
	out_dim_size[0] = in_sds_info.dim_size[rank-2];
	out_dim_size[1] = in_sds_info.dim_size[rank-1];
      }
      else 
      {
	out_dim_size[0] = in_sds_info.dim_size[0];
	out_dim_size[1] = in_sds_info.dim_size[1];
      }
    }

    for (j=0; j<bn_cnt; j++)
    { 
      out_sds_info[j].sds_id = -1;
      out_sds_info[j].sd_id = out_sd_id;
      out_sds_info[j].data_size = in_sds_info.data_size;
      if (max_nbit < 8) out_sds_info[j].data_type = DFNT_UINT8;
      else if (max_nbit < 16) out_sds_info[j].data_type = DFNT_UINT16;
      else out_sds_info[j].data_type = DFNT_UINT32;
      out_sds_info[j].data_size = DFKNTsize(out_sds_info[j].data_type);
      out_sds_info[j].rank = out_rank;
      for (i=0; i<out_rank; i++)
        out_sds_info[j].dim_size[i] = out_dim_size[i];
if (if_cnt == 1)
        sprintf(out_sds_info[j].name, "%s_bits_%s", sds_names[isds], bn_str[j]);
      else
        sprintf(out_sds_info[j].name, "%s_bits_%s_%s", sds_names[isds], bn_str[j], tmp_str);
      if (open_sds((char *)NULL, &out_sds_info[j], 'W') == -1) continue;

      strcpy(attr_name, "_FillValue");
      if (max_nbit < 8)
        write_attr_fval(out_sds_info[j].sds_id, DFNT_UINT8, 1, in_sds_info.fill_val, attr_name);
      else if(max_nbit < 16)
        write_attr_fval(out_sds_info[j].sds_id, DFNT_UINT16, 1, in_sds_info.fill_val, attr_name);
      else
        write_attr_fval(out_sds_info[j].sds_id, DFNT_UINT32, 1, in_sds_info.fill_val, attr_name);
    }

    if (rank == 2)
      ndata_in = ndata_out = in_sds_info.dim_size[1];
    else if (rank > 2)
    {
      if (in_sds_info.dim_size[0] > in_sds_info.dim_size[rank-1])
      {
        ndata_in = in_sds_info.dim_size[1];
        for (i=2; i<rank; i++)
          ndata_in *= in_sds_info.dim_size[i];
        ndata_out = ((n == -1) && (m == -1)) ? ndata_in : in_sds_info.dim_size[1];
      }
      else
      {
	ndata_in = in_sds_info.dim_size[rank-1];
        for (i=0; i<rank-2; i++)
          ndata_in *= in_sds_info.dim_size[i];
        ndata_out = ((n == -1) && (m == -1)) ? ndata_in : in_sds_info.dim_size[rank-1];
      }
    }

    if ((data_in = (void *)calloc(ndata_in, in_sds_info.data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_in in unpack_sds\n");
    if ((data_out = (void *)calloc(ndata_out, out_sds_info[0].data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_out in unpack_sds\n");

    if ((data_in != NULL) && (data_out != NULL))
    {
      if (rank == 2)
      {
        in_edge[0] = out_edge[0] = 1;
	in_edge[1] = out_edge[1] = in_sds_info.dim_size[1];
	nrow = in_sds_info.dim_size[0];
      }
      else if (rank > 2)
      {
	nrow = (bsq == 1) ? in_sds_info.dim_size[rank-2] : in_sds_info.dim_size[0];
	for (i=0; i<rank; i++)
	  in_edge[i] = out_edge[i] = in_sds_info.dim_size[i];
	if ((n == -1) && (m == -1))
	{
	  if (bsq == 0) in_edge[0] = out_edge[0] = 1;
	  else in_edge[rank-2] = out_edge[rank-2] = 1;
	}
	else
	{
	  out_edge[0] = 1;
	  out_edge[1] = out_dim_size[1];
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
	if ((rank>2) && (n == -1) && (m == -1) && (in_sds_info.dim_size[rank-1]>in_sds_info.dim_size[0]))
	  out_start[rank-2] = irow;
        else 
	  out_start[0] = irow;
        if (SDreaddata(in_sds_info.sds_id, in_start, NULL, in_edge, data_in) == FAIL)
	{
	  fprintf(stderr, "Cannot read data line from SDS %s in unpack_sds\n", in_sds_info.name);
	  break;
	}
	for (i=0; i<bn_cnt; i++)
        {
          unpack_bits(data_in, data_out, &in_sds_info, nbit[i], bn_arr[i], ndata_out, st_c, offset);
          if (SDwritedata(out_sds_info[i].sds_id, out_start, NULL, out_edge, data_out) == FAIL)
	    fprintf(stderr, "Cannot write data line to SDS %s in unpack_sds\n", out_sds_info[i].name);
	}
      }
    }
    if (data_in != NULL) free(data_in);
    if (data_out != NULL) free(data_out);
    if (attr_val != NULL) free(attr_val);
    SDendaccess(in_sds_info.sds_id);
    for (i=0; i<bn_cnt; i++)
      SDendaccess(out_sds_info[i].sds_id);
  }	/* for (isds = 0; . . . */
  free(nbit);
  Free2D((void **)bn_arr);
  if (m_opt == 1)
    copy_metadata(in_sds_info.sd_id, out_sd_id);
  SDend(in_sds_info.sd_id);
}

void unpack_bits(void *data_in, void *data_out, sds_t *sds_info, int bn_cnt, int *bn_arr,
	int ndata, int st_c, int offset)
{
  int val_at_ij = 0;
  int i, ic;

  for (i=0, ic=st_c; i<ndata; ic += offset, i++)
  {
    switch(sds_info->data_type)
    {
      case 20: val_at_ij = ((int8 *)data_in)[ic]; break;
      case 21: val_at_ij = ((uint8 *)data_in)[ic]; break;
      case 22: val_at_ij = ((int16 *)data_in)[ic]; break;
      case 23: val_at_ij = ((uint16 *)data_in)[ic]; break;
      case 24: val_at_ij = ((int32 *)data_in)[ic]; break;
      case 25: val_at_ij = ((uint32 *)data_in)[ic]; break;
    }
    if (val_at_ij == sds_info->fill_val)
    {
      if (bn_cnt < 8) ((uint8 *)data_out)[i] = (uint8)FILL_VALUE_UINT8;
      else if (bn_cnt < 16) ((uint16 *)data_out)[i] = (uint16)FILL_VALUE_UINT16;
      else ((uint32 *)data_out)[i] = (uint32)FILL_VALUE_UINT32;
    }
    else
    {
      val_at_ij = val_at_ij >> bn_arr[0];
      val_at_ij = val_at_ij & BIT[bn_cnt-1];
      if (bn_cnt < 8) ((uint8 *)data_out)[i] = (uint8)val_at_ij;
      else if (bn_cnt < 16) ((uint16 *)data_out)[i] = (uint16)val_at_ij;
      else ((uint32 *)data_out)[i] = (uint32)val_at_ij;
    }
  } 
}

int getSdsnameDim(char *sdsname_str, int sd_id, int *n, int *m) 
{  
  int st = 1;
  int p1, p2, p11, len;
  char nd_ext[20], md_ext[20];
  char sds_name[MAX_SDS_NAME_LEN];
  sds_t sds_info;

  memset( &sds_info, 0, sizeof(sds_t));
  *n = *m = -1;
  nd_ext[0] = md_ext[0] = '\0';
  if (sd_charpos(sdsname_str, '(', 0) == -1)
  {      
    if ((p1 = sd_charpos(sdsname_str, '.', 0)) != -1)
    {		  
      sds_info.sd_id = sd_id;
      sds_info.sds_id = -1;		  
      sd_strmid(sdsname_str, 0, p1, sds_name);
      strcpy(sds_info.name, sds_name);
      if (get_sds_info((char *)NULL, &sds_info) != -1)
      {
	p1++;
	len = (int)strlen(sdsname_str);
	if ((p2 = sd_charpos(sdsname_str, '.', p1)) != -1) 
	{
	  sd_strmid(sdsname_str, p1, p2-p1, nd_ext); 
	  *n = (int)atoi(nd_ext);
	  p2++;
	  sd_strmid(sdsname_str, p2, len-p2, md_ext); 
	  *m = (int)atoi(md_ext);	
	}
	else
	{
	  sd_strmid(sdsname_str, p1, len-p1, nd_ext);
	  *n = (int)atoi(nd_ext);	
	}
      } 
      else
      {
	p1++;
	if ((p11 = sd_charpos(sdsname_str, '.', p1)) != -1) 
	{		      		      		      
	  sd_strmid(sdsname_str, 0, p11, sds_name);
	  strcpy(sds_info.name, sds_name);	    
	  len = (int)strlen(sdsname_str);
	  if (get_sds_info((char *)NULL, &sds_info) != -1)
	  {
	    p11++;
	    if ((p2 = sd_charpos(sdsname_str, '.', p11)) != -1) 
	    {
	      sd_strmid(sdsname_str, p11, p2-p11, nd_ext); 
	      *n = (int)atoi(nd_ext);
	      p2++;
	      sd_strmid(sdsname_str, p2, len-p2, md_ext); 
	      *m = (int)atoi(md_ext);	
	    }
	    else
	    {
	      sd_strmid(sdsname_str, p11, len-p11, nd_ext);
	      *n = (int)atoi(nd_ext);	
	    }
	  }		      
	}	
	else
	{
	  strcpy(sds_info.name, sdsname_str);
      	  if (get_sds_info((char *)NULL, &sds_info) != -1)
	  {
	    *n = -1;
	    *m = -1;
	  }
	  else
	  {
	    fprintf(stderr, "Cannot find the SDS %s \n", sds_info.name);
	    st = -1;
	  }  
	}	    		    		    
      }		   
    }
    else 
    {
      *n = -1;
      *m = -1;
    }
  }
  else 
    st = -1;
  if ((*n == 0) || (*m == 0))
  {
    fprintf(stderr, "Invalid layer number in SDS name: %s\n", sdsname_str);
    st = -1;
  }
  else
  /* number of layers retrieved for 3D and 4D sds */		    
  {
    if (*n != -1) --*n;
    if (*m != -1) --*m;
  }    
  if (sds_info.sds_id != -1) SDendaccess(sds_info.sds_id);
  if (sds_info.sd_id != -1) SDendaccess(sds_info.sd_id);
  return st;
}

