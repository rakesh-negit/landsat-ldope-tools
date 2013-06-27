/****************************************************************************
!C

!File: math_sds.c

!Description:
  Perform a simple arithmetic on two input SDSs of the same or different 
  MODIS Land HDF-EOS data products. 

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
"    math_sds - Perform simple arithmatic on two input SDSs of the same or\n" \
"              different MODIS Land HDF-EOS data products. \n" \
" \n" \
"SYNOPSIS \n" \
"    math_sds -help [filename] \n" \
"    math_sds -of=filename \n" \
"             -math=<arithmetic expression>,dt,f_nop1,f_nop2,f_nop3,f_ovf \n" \
"       where arithmetic_expression = <SDS_name1,f1>,<op>,<SDS_name2,f2> \n" \
" \n" \
"DESCRIPTION \n" \
"    Perform simple arithmetic on two input SDSs of the same or different \n" \
"    MODIS Land HDF-EOS data products. Output is an HDF file containing  \n" \
"    the result of the pixelwise arithmetic operation. \n" \
" \n" \
"    The math option may be repeated with different arithmetic expressions \n" \
"    to perform additional arithmetic on the same or different SDSs. \n" \
"    Resulting SDSs of all arithmetic expressions are output to the same\n"\
"    output HDF file.\n" \
" \n" \
"    If the input SDSs in an arithmetic expression are of different\n"\
"    resolution then the resolution of one SDS must be an integral multiple\n"\
"    of the other and the output SDS will be of the higher of the two input\n"\
"    resolutions. \n" \
" \n" \
"   This tool supports 2D/3D/4D SDSs.\n" \
" \n" \
"   The tool command arguments can be specified in any order. \n" \
"\n" \
"OPTIONS \n" \
"    -help              Print this help message, If the input filename is\n"\
"                       specified with this option, then the names of all\n"\
"                       the SDS in the file are displayed. \n" \
"    -math=<arithmetic expression>,dt,f_nop1,f_nop2,f_nop3,f_ovf \n" \
" \n" \
"                       The output data type 'dt' can be one of INT8, UINT8,\n"\
"                       INT16, UINT16, INT32, UINT32, or FLOAT32.\n" \
" \n" \
"                       f_nop1 and f_nop2 are no operation fill values for\n" \
"                       the two input SDS. The arithmetic operation is not \n" \
"                       performed at a pixel if the SDS value at the pixel\n" \
"                       in the input SDS1 is f_nop1 or SDS2 is f_nop2. If\n"\
"                       these arguments are unspecified then the input SDS\n"\
"                       fill values are used.\n" \
" \n" \
"                       If the math operation cannot be performed at a pixel\n"\
"                       then the fill value f_nop3 is written to the output \n"\
"                       SDS pixel. If this argument is unspecified then the\n"\
"                       SDS 1 fill value is used.\n" \
" \n" \
"                       The overflow fill value f_ovf is written to the\n"\
"                       pixel in the output SDS if the arithmetic operation\n"\
"                       at the pixel results in an overflow. If this\n"\
"                       argument is unspecified then the SDS1 fill value is\n"\
"                       used. \n" \
" \n" \
"                       To set any fill values to default the * symbol may\n" \
"                       be specified in place of the actual value. \n" \
" \n" \
"                       where arithmetic expression=\n"\
"                           <SDS_name1,f1>,<op>,<SDS_name2,op2> \n" \
" \n" \
"                       The arithmetic expression contains two operands\n"\
"                       separated by an operator. Each operand contains an\n"\
"                       SDS and file containing the SDS separated by a\n"\
"                       comma. The mathematical operator can be one of\n"\
"                       (+, -, *, /). Expressions are evaluated from left\n"\
"                       to right. \n" \
" \n" \
"                       To process a specific layer of a 3D SDS specify the\n"\
"                       element number of the third dimension as a dot\n"\
"                       extension of the SDS name: sds_name.n (e.g.,\n"\
"                       sur_refl_b02.1 = the layer defined by the 1st\n"\
"                       element of the 3rd dimension of the 3D SDS\n"\
"                       sur_refl_b02). \n" \
" \n" \
"                       To process a specific layer of a 4D SDS, specify the\n"\
"                       higher dimension element number(s) as a dot\n"\
"                       extension of the SDS name: sds_name.n.m (e.g.,\n"\
"                       Surface_Refl.1.2 = the layer defined by the 1st\n"\
"                       element of the 3rd dimension and the 2nd element of\n" \
"                       the 4th dimension of the 4D SDS Surface_Refl). \n" \
" \n" \
"                       Note that wildcards and ranges of element values may\n"\
"                       be specified as sds_name.* and as sds_name.n1-n2.m\n"\
"                       respectively. \n" \
" \n" \
"                       Output data type 'dt' can be one of INT8, UINT8,\n"\
"                       INT16, UINT16, INT32, UINT32, FLOAT32. \n" \
" \n" \
"                       f_nop1, f_nop2 and f_nop3 are no operation fill\n"\
"                       values in the two input SDS. The arithmetic\n"\
"                       operation is not performed at a pixel if the SDS\n"\
"                       value at the pixel in the input SDS1 is f_nop1 or\n"\
"                       SDS2 is f_nop2. If these arguments are not input by\n"\
"                       the user then the SDS fill values are used as no\n"\
"                       operation fill values. \n" \
" \n" \
"                       If the math operation cannot be performed at a pixel\n"\
"                       then the fill value f_nop3 is written to the pixel\n"\
"                       in the SDS. If this argument is not input by the\n"\
"                       user then the fill value of SDS1 is used.\n" \
" \n" \
"                       Overflow fill value f_ovf is written to the pixel in\n"\
"                       the output SDS if the arithmetic operation at the\n"\
"                       pixel results in an overflow. If this argument is\n"\
"                       not input by the user then the fill value of SDS1 is\n"\
"                       used. \n" \
" \n" \
"                       To set any of these fill values to default * can\n"\
"                       used in place of the actual value. \n" \
" \n" \
"    -of=<filename>     Output filename \n" \
" \n" \
"Examples: \n" \
"    math_sds -of=diff_temp.hdf \n" \
"        -math=LST_Day_1km,MOD11A1.A2001044.h09v05.002.2001154214604.hdf,-,\n"\
"        LST_Night_1km,MOD11A1.A2001044.h09v05.002.2001154214604.hdf,INT16,\n"\
"        0,0,-1,-2 \n" \
" \n" \
"    math_sds -of=MOD09A1.h12v03.b03.49-53.hdf\n"\
"        -math=sur_refl_b02,MOD09A1.A1999049.h12v03.001.1999265195217.hdf,/,\n"\
"        sur_refl_b03,MOD09A1.A1999053.h12v03.001.1999250150605.hdf,INT16,\n"\
"        *,*,*,*, \n" \
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 04/05/2004 \n" \

#define USAGE \
"usage: \n" \
"    math_sds -help [filename] \n" \
"    math_sds -of=filename \n" \
"             -math=<arithmetic expression>,dt,f_nop1,f_nop2,f_nop3,f_ovf \n" \
"        where arithmetic_expression = <SDS_name1,f1>,<op>,<SDS_name2,f2> \n" \
" \n" \
"OPTIONS \n" \
"    -help              Print this help message, If the input filename is\n"\
"                       specified with this option, then the names of all\n"\
"                       the SDS in the file are displayed. \n" \
"    -math=<arithmetic expression>,dt,f_nop1,f_nop2,f_nop3,f_ovf \n" \
" \n" \
"                       The output data type 'dt' can be one of INT8, UINT8,\n"\
"                       INT16, UINT16, INT32, UINT32, or FLOAT32.\n" \
" \n" \
"                       f_nop1 and f_nop2 are no operation fill values for\n" \
"                       the two input SDS. The arithmetic operation is not \n" \
"                       performed at a pixel if the SDS value at the pixel\n" \
"                       in the input SDS1 is f_nop1 or SDS2 is f_nop2. If\n"\
"                       these arguments are unspecified then the input SDS\n"\
"                       fill values are used.\n" \
" \n" \
"                       If the math operation cannot be performed at a pixel\n"\
"                       then the fill value f_nop3 is written to the output \n"\
"                       SDS pixel. If this argument is unspecified then the\n"\
"                       SDS 1 fill value is used.\n" \
" \n" \
"                       The overflow fill value f_ovf is written to the\n"\
"                       pixel in the output SDS if the arithmetic operation\n"\
"                       at the pixel results in an overflow. If this\n"\
"                       argument is unspecified then the SDS1 fill value is\n"\
"                       used. \n" \
" \n" \
"                       To set any fill values to default the * symbol may\n" \
"                       be specified in place of the actual value. \n" \
" \n" \
"                       where arithmetic expression=\n"\
"                           <SDS_name1,f1>,<op>,<SDS_name2,op2> \n" \
" \n" \
"                       The arithmetic expression contains two operands\n"\
"                       separated by an operator. Each operand contains an\n"\
"                       SDS and file containing the SDS separated by a\n"\
"                       comma. The mathematical operator can be one of\n"\
"                       (+, -, *, /). Expressions are evaluated from left\n"\
"                       to right. \n" \
" \n" \
"                       To process a specific layer of a 3D SDS specify the\n"\
"                       element number of the third dimension as a dot\n"\
"                       extension of the SDS name: sds_name.n (e.g.,\n"\
"                       sur_refl_b02.1 = the layer defined by the 1st\n"\
"                       element of the 3rd dimension of the 3D SDS\n"\
"                       sur_refl_b02). \n" \
" \n" \
"                       To process a specific layer of a 4D SDS, specify the\n"\
"                       higher dimension element number(s) as a dot\n"\
"                       extension of the SDS name: sds_name.n.m (e.g.,\n"\
"                       Surface_Refl.1.2 = the layer defined by the 1st\n"\
"                       element of the 3rd dimension and the 2nd element of\n" \
"                       the 4th dimension of the 4D SDS Surface_Refl). \n" \
" \n" \
"                       Note that wildcards and ranges of element values may\n"\
"                       be specified as sds_name.* and as sds_name.n1-n2.m\n"\
"                       respectively. \n" \
" \n" \
"                       Output data type 'dt' can be one of INT8, UINT8,\n"\
"                       INT16, UINT16, INT32, UINT32, FLOAT32. \n" \
" \n" \
"                       f_nop1, f_nop2 and f_nop3 are no operation fill\n"\
"                       values in the two input SDS. The arithmetic\n"\
"                       operation is not performed at a pixel if the SDS\n"\
"                       value at the pixel in the input SDS1 is f_nop1 or\n"\
"                       SDS2 is f_nop2. If these arguments are not input by\n"\
"                       the user then the SDS fill values are used as no\n"\
"                       operation fill values. \n" \
" \n" \
"                       If the math operation cannot be performed at a pixel\n"\
"                       then the fill value f_nop3 is written to the pixel\n"\
"                       in the SDS. If this argument is not input by the\n"\
"                       user then the fill value of SDS1 is used.\n" \
" \n" \
"                       Overflow fill value f_ovf is written to the pixel in\n"\
"                       the output SDS if the arithmetic operation at the\n"\
"                       pixel results in an overflow. If this argument is\n"\
"                       not input by the user then the fill value of SDS1 is\n"\
"                       used. \n" \
" \n" \
"                       To set any of these fill values to default * can\n"\
"                       used in place of the actual value. \n" \
" \n" \
"    -of=<filename>     Output filename \n" \
" \n"

#define MAX_NSDS 10

/******************************************************************************
                            Prototypes.
******************************************************************************/

int parse_cmd_math_sds(int argc, char **argv, char **expr, int *n_op, char *f3);
int read_param(char *expr, char *sds1, char *sds2, char *f1, char *f2, char *op_t, 
	       char *dt, char *f_nop1, char *f_nop2, char *f_nop3, char *f_ovf);
void compute_math_sds(sds_t *sds1_info, sds_t *sds2_info, sds_t *sds3_info, char op_t, 
		      char *dt, char *f_nop, char *f_ovf);
void get_sds_param(sds_t *sds_info, int *n, int *m, int *rank, int *dim_size);
int check_sds_param(int rank1, int rank2, int *dim_size1, int *dim_size2, int *sc_dim, 
		    int *bd);
void check_fsds_id(char *sds1, char *sds2, char *f1, char *f2, int *st_sds, int *st_f);
void check_sds_name(char *sds_name);
/******************************************************************************
                            Global Varible Initialization.
******************************************************************************/

int int8_min = -128;
int int8_max = 127;
int uint8_min = 0;
int uint8_max = 255;
int int16_min = -32768;
int int16_max = 32767;
int uint16_min = 0;
int uint16_max = 65535;
int int32_min = (-2147483647 - 1);
int int32_max = 2147483647;
int uint32_min = 0;
#ifdef WIN32
uint32 uint32_max = 4294967295;
#else
uint uint32_max = 4294967295UL;
#endif

int main(int argc, char **argv)
/******************************************************************************
!C

!Description:
  Main function for math_sds.

!Input Parameters: (none)
  command line arguments: see help for details.

!Output Parameters: (none)
  return 0 on successful completion of the process.

!Revision History:
  See file prologue..

!Team-unique Header:
  See file prologue.

!References and Credits: (see file prologue)

!Design Notes: (none)

!END
********************************************************************************/
{
  int i, i_op, n_op;
  int st1, st2 = 0;
  int st_f = 0, st_sds;
  int status;
  char f_nop1[10], f_nop2[10];
  char dt[10], f_nop3[10], f_ovf[10];
  char op_t, **expr;
  char sds1[MAX_SDS_NAME_LEN], sds2[MAX_SDS_NAME_LEN]; 
  char f1[MAX_PATH_LENGTH], f2[MAX_PATH_LENGTH], f3[MAX_PATH_LENGTH];
  sds_t sds1_info, sds2_info, sds3_info;

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

  memset( &sds1_info, 0, sizeof(sds_t) );
  memset( &sds2_info, 0, sizeof(sds_t) );

  if ((expr = (char **)Calloc2D(MAX_NUM_OP, MAX_STR_LEN, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for expr in main()\n"); 
      exit(EXIT_FAILURE);
    }
  else 
  {
    status = parse_cmd_math_sds(argc, argv, expr, &n_op, f3);
    if (status == -1)
      {
	fprintf(stderr, "%s\n", USAGE);
	exit(EXIT_FAILURE);
      }
    else if (status != 0)
      {
	if ((sds3_info.sd_id = SDstart(f3, DFACC_CREATE)) == -1)
	  fprintf(stderr, "Cannot open output HDF file %s\n", f3);
	else
	  { 
	    for (i_op=0; i_op<n_op; i_op++) 
	      {
		if (read_param(expr[i_op], sds1, sds2, f1, f2, &op_t, dt, f_nop1, 
			       f_nop2, f_nop3, f_ovf) == -1)
		  fprintf(stderr, "Cannot process operation %s\n", expr[i_op]);
		else 
		  {
		    sds1_info.sd_id = sds2_info.sd_id = -1;
		    strcpy(sds1_info.name, sds1);
		    strcpy(sds2_info.name, sds2);
		    check_fsds_id(sds1, sds2, f1, f2, &st_sds, &st_f);
		    sds1_info.sds_id = sds2_info.sds_id = -1;
		    st1 = get_sds_info(f1, &sds1_info);
		    if (st_f == 1)
		      {
			sds2_info.sd_id = sds1_info.sd_id;
			if (st_sds == 1)
			  {
			    sds2_info.sds_id = sds1_info.sds_id;
			    sds2_info.rank = sds1_info.rank;
			    for (i=0; i<sds1_info.rank; i++)
			      sds2_info.dim_size[i] = sds1_info.dim_size[i];
			    sds2_info.data_type = sds1_info.data_type;
			    sds2_info.data_size = sds1_info.data_size;
			  }
			else
			  st2 = get_sds_info(f2, &sds2_info);
		      }
		    else
		      st2 = get_sds_info(f2, &sds2_info);
		    if ((st1 != -1) && (st2 != -1))
		      {
			if (f_nop1[0] != '\0')
			  switch(sds1_info.data_type)
			    {
			    case 5: sds1_info.fill_fval = (float32)atof(f_nop1);
			    case 20: sds1_info.fill_val = (int8)atoi(f_nop1);
			    case 21: sds1_info.fill_val = (uint8)atoi(f_nop1);
			    case 22: sds1_info.fill_val = (int16)atoi(f_nop1);
			    case 23: sds1_info.fill_val = (uint16)atoi(f_nop1);
			    case 24: sds1_info.fill_val = (int32)atoi(f_nop1);
			    case 25: sds1_info.fill_val = (uint32)atoi(f_nop1);
			    }
			if (f_nop2[0] != '\0')
			  switch(sds2_info.data_type)
			    {
			    case 5: sds2_info.fill_fval = (float32)atof(f_nop2);
			    case 20: sds2_info.fill_val = (int8)atoi(f_nop2);
			    case 21: sds2_info.fill_val = (uint8)atoi(f_nop2);
			    case 22: sds2_info.fill_val = (int16)atoi(f_nop2);
			    case 23: sds2_info.fill_val = (uint16)atoi(f_nop2);
			    case 24: sds2_info.fill_val = (int32)atoi(f_nop2);
			    case 25: sds2_info.fill_val = (uint32)atoi(f_nop2);
			    }
			compute_math_sds(&sds1_info, &sds2_info, &sds3_info, op_t, dt, f_nop3, f_ovf);
		      }
		    if (sds1_info.sds_id != -1)
		      SDendaccess(sds1_info.sds_id);
		    if (sds2_info.sds_id != -1)
		      {
			if  (st_f != 1)
			  SDendaccess(sds2_info.sds_id);
			else if (st_sds != 1)
			  SDendaccess(sds2_info.sds_id);
		      }
		  }
		if (sds1_info.sd_id != -1)
		  SDend(sds1_info.sd_id);
		if ((sds2_info.sd_id != -1) && (st_f != 1))
		  SDend(sds2_info.sd_id);
	      }
	    SDend(sds3_info.sd_id);
	  }
      }
  }
  fprintf(stderr, "Processing done ! \n");
  return 0; 
}

int parse_cmd_math_sds(int argc, char **argv, char **expr, int *n_op, char *f3) 
/******************************************************************************
!C

!Description:
  Function to parse command line arguments.

!Input Parameters:
  argc: number of input arguments
  argv: string array containing arguments

!Output Parameters:
  expr:      String contains the -SDS option input.
  n_op:      Operation flag.
  f3:        Output file name.

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
  int i, i_op, st;

  st = 1;
  *n_op = 0;
  f3[0] = '\0';
  for (i=1, i_op=0; i<argc; i++)
  {
    if (is_arg_id(argv[i], "-math") == 0)
    {
      get_arg_val(argv[i], expr[i_op]);
      i_op++;
    }
    else if (is_arg_id(argv[i], "-of") == 0)
      get_arg_val(argv[i], f3);
    else fprintf(stderr, "Ignoring unknown option %s\n", argv[i]);
  }
  if (strlen(f3) <= 0) {
    st = -1; fprintf(stderr, "Missing output filename\n");
  }
  if (i_op == 0) {
    st = -1; fprintf(stderr, "No math operation specified\n");
  }
  if (st == 1) *n_op = i_op;
  return st;
}
	              
int read_param(char *expr, char *sds1, char *sds2, char *f1, char *f2, char *op_t, 
	       char *dt, char *f_nop1, char *f_nop2, char *f_nop3, char *f_ovf)

/******************************************************************************
!C

!Description:
  Function to parse the -math option.

!Input Parameters:
  expr:      String contains the -math option input.

!Output Parameters:
  sds1:      The left operand SDS.
  sds2:      The right operand SDS.
  f1:        The input HDF file containing the left operand SDS.
  f2:        The input HDF file containing the right operand SDS.
  op_t:      The string containing the arithmetic operator.
  dt:        Output data type. Default is the data type of the left oeprand SDS.
  f_nop1:    Input SDS fill value for left operand SDS. The arithmetic is not 
             performed if SDS values in left operand SDS equals f_nop1.
  f_nop2:    Input SDS fill value for right oeprand SDS. The arithmetic is not 
             performed if SDS values in right operand SDS equals f_nop2.
  f_nop3:    Input SDS fill value for output SDS.
  f_ovf:     SDS fill value for output SDS if the arithmetic operation results in 
             overflow. The overflow fill option is not implemented for float data
             type.

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
  char op[80];

  op_t[0] = '\0';
  dt[0] = '\0';
  f1[0] = f2[0] = '\0';
  sds1[0] = sds2[0] = '\0';
  f_nop3[0] = f_ovf[0] = '\0';
  f_nop1[0] = f_nop2[0] = '\0';

  len = (int)strlen(expr);
  for (i=0, p1=0; i<9; i++)
  {
    p2 = sd_charpos(expr, ',', p1);
    if (p2 != -1)
    {
      nc = p2 - p1;
      if (nc > 0) {
        switch (i) {
          case 0: sd_strmid(expr, p1, nc, sds1); break;
          case 1: sd_strmid(expr, p1, nc, f1); break;
          case 2: sd_strmid(expr, p1, nc, op);
		  if (strlen(op) == 1) *op_t = op[0];
		  else st = -1;
		  break;
          case 3: sd_strmid(expr, p1, nc, sds2); break;
          case 4: sd_strmid(expr, p1, nc, f2); break;
          case 5: sd_strmid(expr, p1, nc, dt); 
          case 6: sd_strmid(expr, p1, nc, f_nop1); break;
          case 7: sd_strmid(expr, p1, nc, f_nop2); break;
          case 8: sd_strmid(expr, p1, nc, f_nop3); break;
	}
      }
      p1 = p2 + 1;
    }
    else { st = -1; break; }
  }
  if (st == 1) {
    nc = len - p1;
    sd_strmid(expr, p1, nc, f_ovf);
    if (strcmp(f_ovf, "*") == 0) f_ovf[0] = '\0';
    if (strcmp(f_nop1, "*") == 0) f_nop1[0] = '\0';
    if (strcmp(f_nop2, "*") == 0) f_nop2[0] = '\0';
    if (strcmp(f_nop3, "*") == 0) f_nop3[0] = '\0';
    if (strcmp(dt, "*") == 0) dt[0] = '\0';
  }
  return st;
}

void compute_math_sds(sds_t *sds1_info, sds_t *sds2_info, sds_t *sds3_info, 
		      char op_t, char *dt, char *f_nop, char *f_ovf)
/******************************************************************************
!C

!Description:
  Function compute_math_sds to perform the arithmetic operation on two input SDSs.

!Input Parameters:
  sds1_info: The input SDS information structure of the left operand SDS.
  sds2_info: The input SDS information structure of the right operand SDS.
  sds3_info: The SDS information structure of the output SDS.
  op_t:      The string containing the arithmetic operator.
  dt:        Output data type. Default is the data type of the left oeprand SDS.
  f_nop:     Input SDS fill value for output SDS.
  f_ovf:     SDS fill value for output SDS if the arithmetic operation results in 
             overflow. The overflow fill option is not implemented for float data
             type.

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
  int k1, k2;
  int val1 = 0, val2 = 0;
  int nrows, ncols;
  int st_c1, st_c2;
  int ic1, ic2, ic3;
  int i, ir, sc_dim;
  int offset1, offset2;
  int n1, n2, m1, m2, bd;
  int rank1, rank2, rank3;
  int dim_sz1[4], dim_sz2[4];
  int ndata1, ndata2, ndata3;
  void *data1, *data2, *data3;
  int32 edge1[4] = {0, 0, 0, 0};
  int32 edge2[4] = {0, 0, 0, 0};
  int32 edge3[4] = {0, 0, 0, 0};
  int32 start1[4] = {0, 0, 0, 0};
  int32 start2[4] = {0, 0, 0, 0};
  int32 start3[4] = {0, 0, 0, 0};
  int fval_st;
  float fval1 = 0.0, fval2 = 0.0;
  double of_ovf, of_nop;
  double dval1, dval2, dval3 = 0.0;

  fprintf(stdout, "Processing SDS: %s %c %s\n", sds1_info->name, op_t, sds2_info->name);

  check_sds_name(sds1_info->name);
  check_sds_name(sds2_info->name);

  get_sds_param(sds1_info, &n1, &m1, &rank1, dim_sz1);
  get_sds_param(sds2_info, &n2, &m2, &rank2, dim_sz2);
  
  if (check_sds_param(rank1, rank2, dim_sz1, dim_sz2, &sc_dim, &bd) == -1)
    return;
  
  if (sds1_info->data_type != sds2_info->data_type)
    fprintf(stderr, "Input SDSs are of different data type: Continues anyway . . \n");

  sds3_info->sds_id = -1;
  sds3_info->rank = rank3 = rank1;
  for (i=0; i<rank3; i++)
    sds3_info->dim_size[i] = (bd == 1) ? dim_sz1[i] : dim_sz2[i];
  if (dt[0] == '\0') sds3_info->data_type = sds1_info->data_type;
  else {
    if (strcmp(dt, "FLOAT32") == 0) sds3_info->data_type = 5;
    else if (strcmp(dt, "INT8") == 0) sds3_info->data_type = 20;
    else if (strcmp(dt, "UINT8") == 0) sds3_info->data_type = 21;
    else if (strcmp(dt, "INT16") == 0) sds3_info->data_type = 22;
    else if (strcmp(dt, "UINT16") == 0) sds3_info->data_type = 23;
    else if (strcmp(dt, "INT32") == 0) sds3_info->data_type = 24;
    else if (strcmp(dt, "UINT32") == 0) sds3_info->data_type = 25;
    else {
      fprintf(stderr, "Output data type %s not recognized. Set to default\n", dt);
      sds3_info->data_type = sds1_info->data_type;
    }
  }
  sds3_info->data_size = DFKNTsize(sds3_info->data_type);
  sprintf(sds3_info->name, "%s%c%s", sds1_info->name, op_t, sds2_info->name);
  if (f_nop[0] == '\0') {
    if (sds1_info->data_type == 5) of_nop = (double)sds1_info->fill_fval; 
    else of_nop = (double)sds1_info->fill_val;
  }
  else of_nop = (double)atof(f_nop);
  if (f_ovf[0] == '\0') of_ovf = of_nop; 
  else of_ovf = (double)atof(f_ovf);

  if (open_sds((char *)NULL, sds3_info, 'W') != -1)
  {
    ndata1 = compute_sds_ndata(sds1_info);
    ndata2 = compute_sds_ndata(sds2_info);
    ndata3 = compute_sds_ndata(sds3_info);
    if (bd == 1) ndata1 *= sc_dim;
    else ndata2 *= sc_dim;
    ndata3 *= sc_dim;
    if ((data1 = (void *)calloc(ndata1, sds1_info->data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data1 in math_sds: compute_math_sds()\n");
    if ((data2 = (void *)calloc(ndata2, sds2_info->data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data2 in math_sds: compute_math_sds()\n");
    if ((data3 = (void *)calloc(ndata3, sds3_info->data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data3 in math_sds: compute_math_sds()\n");
    if ((data1 != NULL) && (data2 != NULL) && (data3 != NULL))
    {
      rank1 = sds1_info->rank;
      rank2 = sds2_info->rank;
      get_sds_edge(sds1_info, edge1);
      get_sds_edge(sds2_info, edge2);
      get_sds_edge(sds3_info, edge3);
      compute_sds_start_offset(sds1_info, n1, m1, &st_c1, &offset1);
      compute_sds_start_offset(sds2_info, n2, m2, &st_c2, &offset2);
      compute_sds_nrows_ncols(sds3_info, &nrows, &ncols);
      for (ir=0; ir<nrows; ir++)
      {
        if ((rank1 == 2) || (sds1_info->dim_size[0] > sds1_info->dim_size[rank1-1]))
          start1[0] = (bd == 1) ? ir : ir/sc_dim;
        else 
	  start1[rank1-2] = (bd == 1) ? ir : ir/sc_dim;
        if ((rank2 == 2) || (sds2_info->dim_size[0] > sds2_info->dim_size[rank2-1]))
          start2[0] = (bd == 2) ? ir : ir/sc_dim;
        else 
	  start2[rank2-2] = (bd == 2) ? ir : ir/sc_dim;
        if ((rank3 == 2) || (sds3_info->dim_size[0] > sds3_info->dim_size[rank3-1]))
          start3[0] = ir;
        else 
	  start3[rank3-2] = ir;
	if (SDreaddata(sds1_info->sds_id, start1, NULL, edge1, data1) == FAIL)
	{
	  fprintf(stderr, "Cannot read dataline from SDS %s in compute_math_sds()\n", sds1_info->name);
	  break;
	}
	if (SDreaddata(sds2_info->sds_id, start2, NULL, edge2, data2) == FAIL)
	{
	  fprintf(stderr, "Cannot read dataline from SDS %s in compute_math_sds()\n", sds2_info->name);
	  break;
	}
        ic1 = st_c1; ic2 = st_c2;
	k1 = k2 = 0;
        for (ic3=0; ic3<ncols; ic3++)
        {
          switch(sds1_info->data_type) 
	  {
            case 5 : fval1 = ((float32 *)data1)[ic1]; break;
            case 20: val1 = ((int8 *)data1)[ic1]; break;
            case 21: val1 = ((uint8 *)data1)[ic1]; break;
            case 22: val1 = ((int16 *)data1)[ic1]; break;
            case 23: val1 = ((uint16 *)data1)[ic1]; break;
            case 24: val1 = ((int32 *)data1)[ic1]; break;
            case 25: val1 = ((uint32 *)data1)[ic1]; break;
	  }
          switch(sds2_info->data_type) 
	  {
            case 5 : fval2 = ((float32 *)data2)[ic2]; break;
            case 20: val2 = ((int8 *)data2)[ic2]; break;
            case 21: val2 = ((uint8 *)data2)[ic2]; break;
            case 22: val2 = ((int16 *)data2)[ic2]; break;
            case 23: val2 = ((uint16 *)data2)[ic2]; break;
            case 24: val2 = ((int32 *)data2)[ic2]; break;
            case 25: val2 = ((uint32 *)data2)[ic2]; break;
	  }
	  fval_st = 0;
	  if (sds1_info->data_type == 5) {
	    if (fval1 == sds1_info->fill_fval) fval_st = 1;
	  }
	  else if (val1 == sds1_info->fill_val) fval_st = 1;
	  if (sds2_info->data_type == 5) {
	    if (fval2 == sds2_info->fill_fval) fval_st = 1;
	  }
	  else if (val2 == sds2_info->fill_val) fval_st = 1;

	  if (fval_st == 0)
	  {
	    if (sds1_info->data_type == 5) dval1 = (double)fval1;
	    else dval1 = (double)val1;
	    if (sds2_info->data_type == 5) dval2 = (double)fval2;
	    else dval2 = (double)val2;

	    if (op_t == '-') dval3 = dval1 - dval2;
	    else if (op_t == '|') dval3 = abs((int)(dval1 - dval2));
	    else if (op_t == '+') dval3 = dval1 + dval2;
	    else if (op_t == '/') dval3 = (float)dval1 / (float)dval2 * 10000;
	    else if (op_t == '*') dval3 = dval1 * dval2;

            switch(sds3_info->data_type) 
	    {
              case 5 : ((float32 *)data3)[ic3] = (float32)dval3;
		       break;
              case 20: if ((dval3 < int8_min) || (dval3 > int8_max))
			 ((int8 *)data3)[ic3] = (int8)of_ovf;
		       else ((int8 *)data3)[ic3] = (int8)dval3;
		       break;
              case 21: if ((dval3 < uint8_min) || (dval3 > uint8_max))
			 ((uint8 *)data3)[ic3] = (uint8)of_ovf;
		       else ((uint8 *)data3)[ic3] = (uint8)dval3;
		       break;
              case 22: if ((dval3 < int16_min) || (dval3 > int16_max))
			 ((int16 *)data3)[ic3] = (int16)of_ovf;
		       else ((int16 *)data3)[ic3] = (int16)dval3;
		       break;
              case 23: if ((dval3 < uint16_min) || (dval3 > uint16_max))
			 ((uint16 *)data3)[ic3] = (uint16)of_ovf;
		       else ((uint16 *)data3)[ic3] = (uint16)dval3;
		       break;
              case 24: if ((dval3 < int32_min) || (dval3 > int32_max))
			 ((int32 *)data3)[ic3] = (int32)of_ovf;
		       else ((int32 *)data3)[ic3] = (int32)dval3;
		       break;
              case 25: if ((dval3 < uint32_min) || (dval3 > uint32_max))
			 ((uint32 *)data3)[ic3] = (uint32)of_ovf;
		       else ((uint32 *)data3)[ic3] = (uint32)dval3;
		       break;
	    }
	  }
	  else
	  {
            switch(sds3_info->data_type) 
	    {
              case 5 : ((float32 *)data3)[ic3] = (float32)of_nop; break;
              case 20: ((int8 *)data3)[ic3] = (int8)of_nop; break;
              case 21: ((uint8 *)data3)[ic3] = (uint8)of_nop; break;
              case 22: ((int16 *)data3)[ic3] = (int16)of_nop; break;
              case 23: ((uint16 *)data3)[ic3] = (uint16)of_nop; break;
              case 24: ((int32 *)data3)[ic3] = (int32)of_nop; break;
              case 25: ((uint32 *)data3)[ic3] = (uint32)of_nop; break;
	    }
	  }
          if (bd == 0) {
            ic1 += offset1; ic2 += offset2;
          }
          else if (bd == 1) 
	  {
	    ic1 += offset1;
	    k2++;
	    if (k2 == sc_dim) {
	      k2 = 0; ic2 += offset2;
	    }
	  }
          else 
	  {
	    ic2 += offset2;
	    k1++;
	    if (k1 == sc_dim) {
	      k1 = 0; ic1 += offset1;
	    }
	  }
        }
	if (SDwritedata(sds3_info->sds_id, start3, NULL, edge3, data3) == FAIL)
	{
	  fprintf(stderr, "Cannot write dataline for SDS %s in compute_math_sds()\n", sds3_info->name);
	  break;
	}
      } /* for (ir=0; . . .) */
      SDendaccess(sds3_info->sds_id);
    }
    if (data1 != NULL) free(data1);
    if (data2 != NULL) free(data2);
    if (data3 != NULL) free(data3);
  }
}

int check_sds_param(int rank1, int rank2, int *dim_size1, int *dim_size2, 
		    int *sc_dim, int *bd)
/******************************************************************************
!C

!Description:
  Function compute_math_sds to perform the arithmetic operation on two input SDSs.

!Input Parameters:
  rank1:     The rank of the left operand SDS.
  rank2:     The rank of the right operand SDS.
  dim_size1: The dimension sizes of the left operand SDS.
  dim_size2: The dimension sizes of the right operand SDS.   

!Output Parameters:
  sc_dim:    Output dimension size.
  bd:        Flag of if the Input SDSs dimensions are intergral multiples of the other.
             bd = 0   equal.
             bd = 1   Dim SDS1 is greater than Dim SDS2.
             bd = 2   Dim SDS1 is smaller than Dim SDS2.

!Return Value:
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
  int sc_sz[2];
  int sz1[2], sz2[2];

  st = 1;
  if (rank1 != rank2)
  {
    fprintf(stderr, "Input SDS are of different rank: %d \t %d\n", rank1, rank2);
    st = -1;
  }
  else
  {
    if (rank1 == 2)
    {
      sz1[0] = dim_size1[0]; sz1[1] = dim_size1[1]; 
      sz2[0] = dim_size2[0]; sz2[1] = dim_size2[1]; 
    }
    else
    {
      if (dim_size1[0] < dim_size1[rank1-1])
      {
        sz1[0] = dim_size1[rank1-2]; sz1[1] = dim_size1[rank1-1]; 
        sz2[0] = dim_size2[rank1-2]; sz2[1] = dim_size2[rank1-1]; 
      } 
      else
      {
        sz1[0] = dim_size1[0]; sz1[1] = dim_size1[1]; 
        sz2[0] = dim_size2[0]; sz2[1] = dim_size2[1]; 
      } 
    }
    if (sz1[0] == sz2[0]) *bd = 0;
    else if (sz1[0] > sz2[0]) *bd = 1;
    else *bd = 2;
    for (i=0; i<2; i++)
    {
      sc_sz[i] = (*bd == 1) ? sz1[i]%sz2[i] : sz2[i]%sz1[i];
      if (sc_sz[i] != 0) 
      {
	fprintf(stderr, "Input SDSs dimensions are not integral multiples\n");
	st = -1;
	break;
      }
    }
    for (i=0; i<2; i++)
      sc_sz[i] = (*bd == 1) ? sz1[i]/sz2[i] : sz2[i]/sz1[i];
    if (sc_sz[0] != sc_sz[1])
    {
      fprintf(stderr, "All dimensions of input SDSs are not of same mutliples\n");
      st  = -1;
    }
  }
  *sc_dim = sc_sz[0];
  return st;
}
    
void check_fsds_id(char *sds1, char *sds2, char *f1, char *f2, int *st_sds, int *st_f)
/******************************************************************************
!C

!Description:
  Function check_fsds_id to check whether the input files are the same or the 
  input SDSs names are the same.

!Input Parameters:
  sds1:     First input SDS name.         
  sds2:     Second input SDS name.
  f1:       Input HDF file name of the first SDS.
  f2:       Input HDF file name of the second SDS.

!Output Parameters:
  st_sds:   Flag of whether the input SDS names are same.
  st_f:     Flag of whether the input HDF files are same.

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
  int p;
  char org_sds1[MAX_SDS_NAME_LEN];
  char org_sds2[MAX_SDS_NAME_LEN];

  *st_f = (strcmp(f1, f2) == 0) ? 1 : 0; 
  p = sd_charpos(sds1, '.', 0);
  if (p == -1) strcpy(org_sds1, sds1);
  else sd_strmid(sds1, 0, p, org_sds1); 
  p = sd_charpos(sds2, '.', 0);
  if (p == -1) strcpy(org_sds2, sds2);
  else sd_strmid(sds2, 0, p, org_sds2); 
  *st_sds = (strcmp(org_sds1, org_sds2) == 0) ? 1 : 0;
}
  
void check_sds_name(char *sds_name)
/******************************************************************************
!C

!Description:
  Function check_sds_name check if the input SDS name is valid for math_sds.
  If user use * or - to specify layers in 3D or 4D SDS, it's invalid for math_sds.

!Input Parameters:
  sds_name: Input SDS name.

!Revision History:
    See file prologue.

!Team-unique Header:
    See file prologue

!References and Credits:
    See file prologue.

!Design Notes: (none)
    We don't support the use of * or - to specify layer ranges in 3D/4D SDS since
    it will generate every combination of the the layers using the arithmatic 
    opertion. It may not be what the user intend to do and the result may seem 
    to be confusing to the user.

!END
********************************************************************************/
{
  if (sd_charpos(sds_name, '*', 0) != -1)
    {
      fprintf(stderr, "Sorry, Use of '*' option for layer is not valid for math_sds. \n");
      exit(EXIT_FAILURE);
    }
  if (sd_charpos(sds_name, '-', 0) != -1)
    {
      fprintf(stderr, "Sorry, Use of '-' option for layer is not valid for math_sds. \n");
      exit(EXIT_FAILURE);
    }
}
