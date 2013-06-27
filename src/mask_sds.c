
/****************************************************************************
!C

!File: mask_sds.c

!Description:

  This file contains routines for masking sds.

!Input Parameters: (none)

!Output Parameters: (none)

!Revision History:

    Version 1.0    May, 2003

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

    Uses library routines defined in mask_sds_lib.h

!END
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mfhdf.h"
#include "qa_tool.h"
#include "sds_rw.h"
#include "alloc_mem.h"
#include "main_util.h"
#include "meta.h"
#include "str_op.h"
#include "mask_sds_lib.h"

#define MAX_NSDS 20

#define HELP \
"NAME \n" \
" mask_sds - Mask one of more SDS of a MODIS Land HDF-EOS data product file\n" \
"    and output the SDS values at pixels where the mask criteria are met and\n"\
"    output fill values elsewhere. The mask criteria are specified using\n"\
"    relational and logical operators applied to the SDS of the same or\n"\
"    different L2/L3/L4 MODIS Land HDF-EOS data products. \n" \
"\n" \
"SYNOPSIS \n" \
"    mask_sds -help [filename] \n" \
"\n" \
"    mask_sds -of=<output filename> -sds=<SDSname1>[,<SDSname2>[,...]]> \n" \
"             [-fill=<mask fill value>] -mask=<mask1>[,AND|OR,<mask2>[,...]]\n"\
"             [-meta] filename \n" \
"       where maskn=< filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
"\n" \
"DESCRIPTION \n" \
"    Mask one of more SDS of a MODIS Land HDF-EOS data product file and\n"\
"    output the SDS values at pixels where the mask criteria are met and\n"\
"    output mask fill values elsewhere. \n" \
"\n" \
"    The mask criteria are specified using relational and logical operators\n"\
"    applied to the SDS of the same or different L2/L3/L4 MODIS Land HDF-EOS\n"\
"    data products.  The SDS(s) used to define the masking criteria must\n"\
"    have the same or lower spatial resolution as the input file SDS(s) to\n"\
"    be masked. \n" \
"\n" \
"    The mask criteria are defined by a combination of one or more\n"\
"    individual masks. Each mask is defined by testing SDS bits against bit\n"\
"    values using a relational operator. Testing using a decimal value is\n"\
"    also supported. Different masks are combined using the logical AND or\n"\
"    OR operators.\n" \
"\n" \
"    Pixel values that do not meet the masking criteria are assigned a mask\n"\
"    fill value. The mask fill value may be optionally specified or will be\n"\
"    set automatically. If the specified mask fill value is set equal to the\n"\
"    input SDS fill value or to a valid input file SDS value then the tool\n"\
"    will issue a warning message and request another mask fill value. \n" \
"\n" \
"    If pixels in the input file SDS(s) have fill values they cannot be\n"\
"    masked. The input file SDS fill value will be output at these pixels. \n" \
"\n" \
"    If pixels in the SDS(s) used to define the masking criteria have fill\n"\
"    values then the masking cannot be performed. The mask fill value will\n"\
"    be output at these pixels. \n" \
"\n" \
"    This tool supports 2D/3D/4D SDSs. Note, only a two dimensional (2D) SDS\n"\
"    or a 2D layer of a 3D/4D SDS can be used to make a mask. \n" \
"\n" \
"    The tool command arguments can be specified in any order. \n" \
"\n" \
"OPTIONS \n" \
"    -help [filename]        Display this help message. If the input\n"\
"                            filename is specified with this option, then\n"\
"                            the names of all the SDSs in the file are\n" \
"                            displayed. \n" \
"    -of=<filename>          Output filename.\n" \
"    -sds=<SDS list>         List of SDSs present in the input file to be\n"\
"                            masked and written to the output file.  SDS\n"\
"                            names must be specified separated by commas\n" \
"                            with no space. \n" \
"                            To mask a specific layer of a 3D SDS specify\n"\
"                            the element number of the third dimension as a\n"\
"                            dot extension of the SDS name: sds_name.n\n" \
"                            (e.g., sur_refl_b02.1 = the layer defined by\n"\
"                            the 1st element of the 3rd dimension of the 3D\n"\
"                            SDS sur_refl_b02). A wildcard may be used as\n" \
"                            sds_name.* \n" \
"                            To mask a specific layer of a 4D SDS, specify\n"\
"                            the higher dimension element number(s) as a dot\n"\
"                            extension of the SDS name: sds_name.n.m \n" \
"                            (e.g., Surface_Refl.1.2 = the layer defined by\n"\
"                            the 1st element of the 3rd dimension and 2nd\n"\
"                            element of the 4th dimesnsion of the 4D SDS\n" \
"                            Surface_Refl). A range of element values may be\n"\
"                            specified as sds_name.n1-n2.m \n" \
"    -fill=<mask fill value> User specified mask fill value. The fill value\n"\
"                            is stored as the SDS attribute Mask_FillValue\n"\
"                            in the output file and is printed to stdout\n"\
"                            when the tool runs. If a mask fill value is not\n"\
"                            specified then an arbitrary value not equal to\n"\
"                            the input file SDS fill value and not equal to\n"\
"                            a valid input file SDS value is assigned. \n" \
"    -meta                   Copy metadata from the input file to the output\n"\
"                            file. \n" \
"    -mask=<mask1>[,AND|OR,<mask2>[,..]] \n" \
"        where maskn=< filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
"                            Define a mask from one or more individual masks\n"\
"                            combined using the logical operators AND or OR.\n"\
"                            Each individual mask consists of: \n" \
"                                -filename= MODLAND product file \n" \
"                                -SDSname= name of an SDS in the file \n" \
"                                -bit_numbers= a list or range of SDS bits \n" \
"                                -operator= relational operator (>, <, <=,\n"\
"                                    >=, ==, !=) \n" \
"                                -bit_values= bit values that are tested\n"\
"                                    against \n" \
"                            The bits in bit_numbers are specified by the\n"\
"                            lower bit followed by the higher bit and the\n"\
"                            bit_values are specified in the reverse order.\n" \
"                            For example, 0-2,4==0101 signifies bits\n"\
"                            4,2,1,0==0101. If the bit_numbers are omitted,\n"\
"                            then the bit_values are parsed as a decimal\n"\
"                            value. This provides a convenient way to refer\n"\
"                            to a specific value, instead of a list of bits.\n"\
"                            For example, -mask=file,SDS,>=200 makes a mask\n"\
"                            where only the SDS values in the file greater\n"\
"                            than or equal to 200 are considered. If several\n"\
"                            masks are combined together then * may be used\n"\
"                            in place of the filename and/or SDS name to\n"\
"                            specify the same filename and/or SDS name used\n"\
"                            in the previous mask. For example, \n" \
"                                -mask=file1,SDS1,0-2,4==0101,AND,*,*,4-5==10\n"\
"                            To specify a 3D SDS layer write the element\n"\
"                            number of the third dimension as a dot\n"\
"                            extension of the SDS name: sds_name.n (e.g.,\n"\
"                            sur_refl_b02.1 = the layer defined by the 1st\n"\
"                            element of the 3rd dimension of the 3D SDS \n" \
"                            sur_refl_b02). \n" \
"                            To specify a 4D SDS layer write the higher\n"\
"                            dimension element number(s) as a dot extension\n"\
"                            of the SDS name: sds_name.n.m (e.g.,\n"\
"                            Surface_Refl.1.2 = the layer defined by the 1st\n"\
"                            element of the 3rd dimension and 2nd element \n" \
"                            of the 4th dimesnsion of the 4D SDS\n"\
"                            Surface_Refl). \n" \
"    filename                Input filename \n" \
"\n" \
"EXAMPLE \n" \
"    mask_sds -sds=Snow_Cover -of=mod10_mask.hdf \n" \
"    -mask=MOD10A1.A2002122.h30v11.004.2002101010341.hdf,Snow_Spatial_QA,1-3==000\n" \
"    MOD10A1.A2002122.h30v11.004.2002101010341.hdf\n" \
"\n"\
"    The above example outputs to mod10_mask.hdf the Snow_ Cover values\n" \
"    stored in MOD10A1.A2003122.h30v11.004.2003134090541.hdf where\n"\
"    Snow_Spatial_QA bits 1-3 equal 000. \n" \
"    Pixels where Snow_Spatial_QA bits 1-3 do not equal 000 are assigned an\n"\
"    automatically generated mask fill value. \n" \
"\n" \
"    mask_sds -sds=Snow_Cover -of=mod10_mask.hdf \n" \
"    -mask=MOD10A1.A2002122.h30v11.004.2002101010341.hdf,Snow_Spatial_QA,1-3==000\n"\
"    -fill=251 MOD10A1.A2002122.h30v11.004.2002101010341.hdf \n" \
"\n"\
"    The above example is the same as the first but pixels where\n"\
"    Snow_Spatial_QA bits 1-3 do not equal 000 are assigned a mask fill\n"\
"    value of 251. \n" \
"\n" \
"    mask_sds -sds=sur_refl_b01,sur_refl_b04,sur_refl_b03\n"\
"    -of=mod09_RGB_mask.hdf \n" \
"    -mask='MOD09A1.A2003113.h24v04.003.2003134160954.hdf,sur_refl_b01,>1500'\n" \
"    MOD09A1.A2003113.h24v04.003.2003134160954.hdf \n" \
"\n"\
"    The above example outputs to mod09_RGB_mask.hdf the SDSs sur_refl_b01,\n"\
"    sur_refl_b04, sur_refl_b03 in\n"\
"    MOD09A1.A2003113.h24v04.003.2003134160954.hdf where sur_refl_b01 is\n"\
"    greater than 1500. Pixels where sur_refl_b01 <= 1500 are assigned an\n"\
"    automatically generated mask fill value. \n" \
" \n" \
"    mask_sds -sds=Snow_Cover -of=howzat.hdf \n" \
"    -mask='MOD09A1.A2003113.h24v04.003.2003134160954.hdf,sur_refl_b01,>1500'\n"\
"    MOD10A1.A2003113.h24v04.004.2003101010541.hdf \n" \
"\n"\
"    The above example outputs to howzat.hdf the Snow_Cover values stored \n" \
"    in MOD10A1.A2003113.h24v04.004.2003101010541.hdf at those pixels in \n" \
"    MOD09A1.A2003113.h24v04.003.2003134160954.hdf where sur_refl_b01 is\n"\
"    greater than 1500. \n" \
"    Pixels where sur_refl_b01 <= 1500 are assigned an automatically\n"\
"    generated mask fill value. \n" \
"\n" \
"    mask_sds -sds=sur_refl_b03 -of=mod09_qa_mask.hdf \n" \
"    -mask='MOD09A1.A2003113.h24v04.003.2003134160954.hdf,sur_refl_qc_500m,\n"\
"    10-13==1100,AND,*,sur_refl_state_500m,3-5==001'\n"\
"    MOD09A1.A2003113.h24v04.003.2003134160954.hdf \n" \
"\n"\
"    The above example outputs to mod09_qa_mask.hdf the SDS sur_rel_b03 at\n"\
"    those pixels in MOD09A1.A2003113.h24v04.003.2003134160954.hdf where the\n"\
"    sur_refl_qc_500m bits 10-13 equal 1100 AND sur_refl_state_500m bits 3-5\n"\
"    equal 001. Pixels where sur_refl_qc_500m bits 10-13 do not equal 1100\n"\
"    or sur_refl_state_500m bits 3-5 do not equal 001 are assigned an \n" \
"    automatically generated  mask fill value. \n" \
"\n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
"\n" \
"Version 1.0, 05/15/2003 \n" 

#define USAGE \
"usage: \n" \
"    mask_sds -help [filename] \n" \
"\n" \
"    mask_sds -of=<output filename> -sds=<SDSname1>[,<SDSname2>[,...]]> \n" \
"    [-fill=<mask fill value>] -mask=<mask1>[,AND|OR,<mask2>[,...]] [-meta]\n"\
"    filename \n" \
"       where maskn=< filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
"\n" \
"OPTIONS \n" \
"    -help [filename]        Display this help message. If the input\n"\
"                            filename is specified with this option, then\n"\
"                            the names of all the SDSs in the file are\n" \
"                            displayed. \n" \
"    -of=<filename>          Output filename.\n" \
"    -sds=<SDS list>         List of SDSs present in the input file to be\n"\
"                            masked and written to the output file.  SDS\n"\
"                            names must be specified separated by commas\n"\
"                            with no space. \n" \
"                            To mask a specific layer of a 3D SDS specify\n"\
"                            the element number of the third dimension as a\n"\
"                            dot extension of the SDS name: sds_name.n\n" \
"                            (e.g., sur_refl_b02.1 = the layer defined by\n"\
"                            the 1st element of the 3rd dimension of the 3D\n"\
"                            SDS sur_refl_b02). A wildcard may be used as\n" \
"                            sds_name.* \n" \
"                            To mask a specific layer of a 4D SDS, specify\n"\
"                            the higher dimension element number(s) as a dot\n"\
"                            extension of the SDS name: sds_name.n.m \n" \
"                            (e.g., Surface_Refl.1.2 = the layer defined by\n"\
"                            the 1st element of the 3rd dimension and 2nd\n"\
"                            element of the 4th dimesnsion of the 4D SDS\n" \
"                            Surface_Refl). A range of element values may be\n"\
"                            specified as sds_name.n1-n2.m \n" \
"    -fill=<mask fill value> User specified mask fill value. The fill value\n"\
"                            is stored as the SDS attribute in the output\n"\
"                            file and is printed to stdout when the tool\n"\
"                            runs. If a mask fill value is not specified\n" \
"                            then an arbitrary value not equal to the input\n"\
"                            file SDS fill value and not equal to a valid\n"\
"                            input file SDS value is assigned. \n" \
"    -meta                   Copy metadata from the input file to the output\n"\
"                            file. \n" \
"    -mask=<mask1>[,AND|OR,<mask2>[,..]] \n" \
"        where maskn=< filename>,<SDSname>,<bit_numbers operator bit_values>\n"\
"                            Define a mask from one or more individual masks\n"\
"                            combined using the logical operators AND or OR.\n"\
"                            Each individual mask consists of: \n" \
"                                -filename= MODLAND product file \n" \
"                                -SDSname= name of an SDS in the file \n" \
"                                -bit_numbers= a list or range of SDS bits \n" \
"                                -operator= relational operator (>, <, <=,\n"\
"                                    >=, ==, !=) \n" \
"                                -bit_values= bit values that are tested\n"\
"                                    against \n" \
"                            The bits in bit_numbers are specified by the\n"\
"                            lower bit followed by the higher bit and the\n"\
"                            bit_values are specified in the reverse order.\n" \
"                            For example, 0-2,4==0101 signifies bits\n"\
"                            4,2,1,0==0101. If the bit_numbers are omitted,\n"\
"                            then the bit_values are parsed as a decimal\n"\
"                            value. This provides a convenient way to refer\n"\
"                            to a specific value, instead of a list of bits.\n"\
"                            For example, -mask=file,SDS,>=200 makes a mask\n"\
"                            where only the SDS values in the file greater\n"\
"                            than or equal to 200 are considered. If several\n"\
"                            masks are combined together then * may be used\n"\
"                            in place of the filename and/or SDS name to\n"\
"                            specify the same filename and/or SDS name used\n" \
"                            in the previous mask. For example, \n" \
"                                -mask=file1,SDS1,0-2,4==0101,AND,*,*,4-5==10\n"\
"                            To specify a 3D SDS layer write the element\n"\
"                            number of the third dimension as a dot\n"\
"                            extension of the SDS name: sds_name.n (e.g.,\n"\
"                            sur_refl_b02.1 = the layer defined by the 1st\n"\
"                            element of the 3rd dimension of the 3D SDS \n" \
"                            sur_refl_b02). \n" \
"                            To specify a 4D SDS layer write the higher\n"\
"                            dimension element number(s) as a dot extension\n"\
"                            of the SDS name: sds_name.n.m (e.g.,\n"\
"                            Surface_Refl.1.2 = the layer defined by the 1st\n"\
"                            element of the 3rd dimension and 2nd element of\n"\
"                            the 4th dimesnsion of the 4D SDS Surface_Refl).\n"\
"    filename                Input filename \n" \
"\n" 

/****************************************************************************************
                              Prototypes
****************************************************************************************/

int parse_cmd_mask_sds(int argc, char **argv, int *arg_cnt, char **arg_list, int *opt, 
		int *fqa_l2g, char *m_str, char **sds_names, int *sds_cnt, char *in_fname, 
		int *m_opt, int *f_opt, int *fill_val);
int mask_nsds(int fin_l2g, int *fqa_l2g, char *m_str, char **arg_list, int arg_cnt, int *opt_arr, 
		char **sds_names, int nsds, int32 out_sd_id, int m_opt, int f_opt, int fill_val);
void mask_nsds_data_row(void **data_in, void **data_out, uint8 *data_mask, int ndata_out, 
		int ndata_mask, int nsds, int bsq, int *st_c, int *offset, sds_t *sds_info, 
		int *mask_fill);

/*************************************************************************************/

int main(int argc, char *argv[])

/*
!C************************************************************************************
       
!Description:
 The mask_sds routine takes a list of arguments: input MODLand L3 product file name, 
 mask, a combination of SDS from the same or many different HDF files, and the output file 
 name. Then it output to the output file the SDSs containing data of pixels from the 
 specified L3 SDSs of the input file that verify the given mask.

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
  int32 out_sd_id;
  int status, st;
  int mm, dd, yy;
  int i, k, p1, p2, nday, day_id;
  int m_opt, nobs = 0, f_opt, fill_val;
  int sds_cnt, nsds, meta_cnt, arg_cnt;
  int len, isds, fin_l2g;
  int32 msds;
  int32 nattr, rank, dim_size[4], sd_id, sds_id, dt; 
  int opt[MAX_NUM_OP], fqa_l2g[MAX_NUM_OP];
  char tmp_str[10];
  char meta_name[MAX_META_NAME_LEN];
  char mm_str[10], dd_str[10], yy_str[10];
  char org_name[80], org_num[5], cday[20];
  char **meta_val, **arg_list;
  char *metadata, **sds_names, **tmp_sds_names;
  char tile[10], jday[10], esdt[10], loc_gid[80];
  char in_fname[MAX_PATH_LENGTH], m_str[MAX_STR_LEN];
  char name[MAX_SDS_NAME_LEN], dim_str[MAX_STR_LEN], sds_name[MAX_SDS_NAME_LEN];

  status = 0;  

  if (argc == 1)
    {
      status = 1;
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_FAILURE);
    }
  
  /* display help message when user use the -help option */

  
  if ((argc == 2) && (strcmp(argv[1], "-help") == 0))
    {
      fprintf(stderr, "%s\n", HELP);
      exit(EXIT_SUCCESS);
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
                      fprintf(stdout, "Valid SDS names, dimension and data type in file: %s\n", 
			      argv[i]);
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


  if ((arg_list = (char **)Calloc2D(3*MAX_NUM_OP+3, MAX_PARAM_LENGTH, sizeof(char))) == NULL)
    {
      status = 1;
      fprintf(stderr, "Cannot allocate memory for arg_list in mask_sds: main()\n");
    }

  if ((sds_names = (char **)Calloc2D(MAX_NSDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    {
      status = 1;
      fprintf(stderr, "Cannot allocate memory for sds_names mask_sds: main()\n");
    }

  if ((arg_list != NULL) && (sds_names != NULL))
    {
      if ((st = parse_cmd_mask_sds(argc, argv, &arg_cnt, arg_list, opt, fqa_l2g, m_str, sds_names, 
				      &sds_cnt, in_fname, &m_opt, &f_opt, &fill_val)) == -1)
	{
	  status = 1;
	  fprintf(stderr, "%s\n", USAGE);
	}
      else if (st != 0)     /* valid command line input */ 
	{ 
	  if ((meta_val = (char **)Calloc2D(5, MAX_GID_LEN, sizeof(char))) == NULL)
	    {
	      status = 1;
	      fprintf(stderr, "Cannot allocate memory for meta_val mask_sds: main()\n");
	    }
	  if ((tmp_sds_names = (char **)Calloc2D(MAX_NSDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
	    {
	      status = 1;
	      fprintf(stderr, "Cannot allocate memory for tmp_sds_names mask_sds: main()\n");
	    }
	  if ((meta_val != NULL) && (tmp_sds_names != NULL)) 
	    {
	      if ((out_sd_id = SDstart(arg_list[arg_cnt-2], DFACC_CREATE)) == FAIL)
		{
		  status = 1;
		  fprintf(stderr, "Cannot create output file %s in mask_sds: main() \n", 
			  arg_list[arg_cnt-2]);
		}
	      else
		{
		  strcpy(arg_list[0], in_fname);
		  fprintf(stderr, "\nProcessing HDF file: %s\n", in_fname);
		  fprintf(stderr, "--------------------------------------------------------------------\n");
		  
		  /* Search for the tile/granule indentifier */
		  meta_cnt = 0;
		  if ((metadata = get_attr_metadata(in_fname, "CoreMetadata.0")) != NULL)
		    {
		      strcpy(meta_name, "LOCALGRANULEID");
		      get_sel_metadata(metadata, meta_name, meta_val, &meta_cnt, 0);
		      free(metadata);
		      if(meta_cnt != 0)
			{
			  strcpy(loc_gid, meta_val[0]);
			  p1 = sd_strpos(meta_val[0], ".A", 0);
			  sd_strmid(meta_val[0], 1, p1-1, esdt);
			  sd_strmid(meta_val[0], p1+2, 7, jday);
			  sd_strmid(meta_val[0], p1+10, 6, tile);
			}
		    }
		  
		  /* is input an l2g file */
		  meta_cnt = 0;
		  if ((metadata = get_attr_metadata(in_fname, "ArchiveMetadata.0")) != NULL)
		    {
		      strcpy(meta_name, "NUMBEROFOVERLAPGRANULES");
		      get_sel_metadata(metadata, meta_name, meta_val, &meta_cnt, 0); 
		    }
		  fin_l2g = ((meta_cnt != 0) && (atoi(meta_val[0]) != 0)) ? 1 : 0;
		  if (fin_l2g == 1)
		    {
		      meta_cnt = 0;
		      strcpy(meta_name, "MAXIMUMOBSERVATIONS");
		      get_sel_metadata(metadata, meta_name, meta_val, &meta_cnt, 0);
		      if (meta_cnt != 0) nobs = (int)atoi(meta_val[0]);
		    }
		  
		  /* if the product is MOD14A1 or MYD14A1 get the actual day number for the SDS */
		  if ((strstr(loc_gid, "MOD14A1") != NULL) || (strstr(loc_gid, "MYD14A1") != NULL))
		    {
		      meta_cnt = 0;
		      strcpy(meta_name, "NUMBEROFDAYS");
		      get_sel_metadata(metadata, meta_name, meta_val, &meta_cnt, 0);
		      strcpy(meta_name, "DAYSOFYEAR");
		      get_sel_metadata(metadata, meta_name, meta_val, &meta_cnt, 0);
		      if (meta_cnt == 2) 
			{
			  nday = (int)atoi(meta_val[0]);
			  p1 = sd_charpos(sds_names[0], '.', 0);
			  len = (int)strlen(sds_names[0]);
			  sd_strmid(sds_names[0], p1+1, len-p1-1, tmp_str);
			  day_id = (int)atoi(tmp_str);
			  if (day_id > nday)
			    {
			      Free2D((void **)arg_list);
			      Free2D((void **)sds_names);
			      Free2D((void **)tmp_sds_names);
			      Free2D((void **)meta_val);
			      free(metadata);
			      fprintf(stderr, "SDS layer %d not found in %s\n", day_id, sds_names[0]);
			      status = 1;
			      exit(EXIT_FAILURE);
			    }
			  else
			    {
			      for (i=0, p2=0; i<day_id; i++)
				{
				  p1 = p2 + 1;
				  p2 = sd_charpos(meta_val[1], ',', p1+1);
				}
			      if (day_id == nday) p2 = (int)strlen(meta_val[1]) - 1;
			      sd_strmid(meta_val[1], p1, p2-p1, cday);
			      sd_strtrim(cday);
			      
			      p1 = sd_charpos(cday, '-', 1);
			      sd_strmid(cday, 1, p1-1, yy_str);
			      p1++;
			      p2 = sd_charpos(cday, '-', p1);
			      sd_strmid(cday, p1, p2-p1, mm_str);
			      len = (int)strlen(cday);
			      p2++;
			      sd_strmid(cday, p2, len-p2-1, dd_str);
			      yy = (int)atoi(yy_str);
			      mm = (int)atoi(mm_str);
			      dd = (int)atoi(dd_str);
			      if (conv_date(&mm, &dd, yy) != -1)
				{
				  if (dd < 10) sprintf(jday, "%d00%d", yy, dd);
				  else if (dd < 100) sprintf(jday, "%d0%d", yy, dd);
				  else sprintf(jday, "%d%d", yy, dd);
				}
			    }
			}
		    }
		  
		  if (metadata != NULL) free(metadata);
		  
		  /* get SDS names and update the names if it is l2g/3D/4D SDS */
		  for (isds=0; isds<sds_cnt; isds++)
		    strcpy(tmp_sds_names[isds], sds_names[isds]);
		  if ((sds_cnt == 1) && (strcmp(tmp_sds_names[0], "all") == 0))
		    nsds = get_sds_names(in_fname, tmp_sds_names);
		  else nsds = sds_cnt;
		  if (nsds > 0)
		    {
		      if (fin_l2g == 0) update_nd_sdsnames(tmp_sds_names, &nsds, in_fname);
		      else update_l2g_sdsnames(tmp_sds_names, &nsds, in_fname, nobs);
		      
		      if (fin_l2g == 1)
			{
			  for (isds=0; isds<nsds; isds++)
			    {
			      len = (int)strlen(tmp_sds_names[isds]);
			      p1 = sd_charpos(tmp_sds_names[isds], '.', 0);
			      sd_strmid(tmp_sds_names[isds], 0, p1, org_name);
			      sd_strmid(tmp_sds_names[isds], p1, len-p1, org_num);
			      sprintf(tmp_sds_names[isds], "%s_1%s", org_name, org_num);
			    }
			}
		      if (mask_nsds(fin_l2g, fqa_l2g, m_str, arg_list, arg_cnt, opt, tmp_sds_names,
				    nsds, out_sd_id, m_opt, f_opt, fill_val) != 1)
			fprintf(stderr, "Mask SDS failed . . Output may be in error \n");
		    }
		  SDend(out_sd_id);
		} 
	      Free2D((void **)meta_val);
	      Free2D((void **)tmp_sds_names);
	    } /* if ((meta_val != NULL) . . .) */
	} /* if (st != -1) */
      Free2D((void **)arg_list);
      Free2D((void **)sds_names);
    } /* if ((arg_list != NULL) . . .  ) */
  fprintf(stderr, "Processing done ! \n"); 
  return status;
}

int parse_cmd_mask_sds(int argc, char **argv, int *arg_cnt, char **arg_list, int *opt, 
			  int *fqa_l2g, char *m_str, char **sds_names, int *sds_cnt, 
			  char *in_fname, int *m_opt, int *f_opt, int *fill_val)
/*
!C*********************************************************************************

!Function: parse_cmd_mask_sds
       
!Description:
   function parse_cmd_mask_sds parses the argument string to the corresponding
   count of masking strings, list of masking string in  the input mask(arg_list), 
   operator string, l2g file indicator,mask string, SDS names, count of SDSs, 
   input filename (in_fname), meta option indicator.
 
   It print out a error message if there the input file, output file or 
   input mask is missing. 

!Input Parameters: 
  argc      input argument count
  argv      input argument vector

!Input/output Parameters:

  arg_cnt   count of the masking strings in the masking string list.
  arg_list  list of masking string in the input mask
  opt       operator string (>, >=, <, <=, ==, !=)
  fqa_l2g   indicator of whether it is l2g file or not.
  m_str     masking string
  sds_name  list of SDSs to mask.
  sds_cnt   count of SDSs
  in_fname  input L3 MODIS Land product name
  m_opt     indicator of whether the meta data is to be copied to the output file.

!Output Parameters:

  st   :    1  --- successful completion
            -1  --- error exit   

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/

{
  int i, k;
  int st = 1;
  char fill_str[20];
  char out_fname[MAX_PATH_LENGTH];
  
  *sds_cnt = *m_opt = *f_opt = 0;
  m_str[0] = out_fname[0] = '\0';
  in_fname[0] = '\0';
  
  for (i=1; i<argc; i++)
    {
      if ((is_arg_id(argv[i], "-of=") == 0) || (is_arg_id(argv[i], "-o=") == 0))
	get_arg_val(argv[i], out_fname);
      else if (is_arg_id(argv[i], "-sds=") == 0)
	get_arg_val_arr(argv[i], sds_names, sds_cnt);
      else if ((is_arg_id(argv[i], "-m=") == 0) || (is_arg_id(argv[i], "-mask=") == 0))
	get_arg_val(argv[i], m_str);
      else if (strcmp(argv[i], "-meta") == 0) *m_opt = 1;
      else if (is_arg_id(argv[i], "-fill=") == 0) 
      {
	*f_opt = 1;
	get_arg_val(argv[i], fill_str);
        *fill_val = (int)atoi(fill_str);
      }
      else
	{
	  if (argv[i][0] == '-') fprintf(stderr, "Unknown option %s\n", argv[i]);
	  else strcpy(in_fname, argv[i]);
	}
    }
  if (in_fname[0] == '\0') {
    st = -1; fprintf(stderr,"Missing input file \n");
  }
  if (out_fname[0] == '\0') {
    st = -1; fprintf(stderr,"Missing output file \n");
  }
  if (m_str[0] == '\0') {
    st = -1; fprintf(stderr,"Missing mask option \n");
  }

  if (st == 1)
    {
      if (*sds_cnt == 0)
	{
	  *sds_cnt = 1;
	  strcpy(sds_names[0], "all");
	  fprintf(stderr, "No SDS name input. Masking all SDS in the input file. . \n");
	}
      if ((k = get_mask_string(m_str, arg_list, opt, fqa_l2g)) < 0) st = -1;
      if (st == 1)
	{
	  *arg_cnt = 3 + k*3;
	  strcpy(arg_list[*arg_cnt], out_fname);
	  *arg_cnt = *arg_cnt + 1;
	  arg_list[*arg_cnt][0] = '\0';
	  *arg_cnt = *arg_cnt + 1;
	}
    }
  return st;
}

int mask_nsds(int fin_l2g, int *fqa_l2g, char *m_str, char **arg_list, int arg_cnt, int *opt_arr, 
		char **sds_names, int nsds, int32 out_sd_id, int m_opt, int f_opt, int fill_val)

/*
!C*********************************************************************************

!Function: mask_nsds
       
!Description:
   Function mask_nsds apply the mask to input SDSs.

!Input/output Parameters:

   fin_l2g    indicator of whether input file is an l2g file.
   fqa_l2g    Array holds l2g information.
   m_str      masking string
   arg_cnt    count of the masking strings in the masking
   arg_list   argument list containing parameters.
   opt_arr    operator string (>, >=, <, <=, ==, !=)
   sds_names  list of SDSs to mask.
   nsds       count of SDSs
   out_sd_id  output SDS id
   m_opt      indicator of whether the meta data is to be copied to the output file. 

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/

{
  int out_hdf_st;
  uint8 *data_mask;
  int rank, data_size;
  int len, p1, obs_num_in;
  int st, i_op, j_op, n_op;
  int max = 0, min = 0, diff_max, diff_min;
  int i, j, isds, jsds, irow, nrow;
  int mask_fill[MAX_NUM_SDS], tmp_val;
  int bsq, ndata_in, ndata_out, ndata_mask;
  int sel_qa_op[MAX_NUM_OP], rel_op[MAX_NUM_OP];
  int res_s[MAX_NUM_OP], res_l[MAX_NUM_OP], obs_num[MAX_NUM_OP];
  int st_c[MAX_NSDS], offset[MAX_NSDS], n[MAX_NSDS], m[MAX_NSDS];
  char **qa_fnames, num_str[10], org_sds_name[MAX_SDS_NAME_LEN];
  int32 *data_qa_nadd[MAX_NUM_OP], *data_in_nadd = NULL; 
  int32 **in_edge, **in_start, **out_edge, **out_start;
  sds_t sds_info, in_sds_nobs_info, *in_sds_info, *in_sdsc_info = NULL, *out_sds_info = NULL;
  sds_t qa_sdsc_info[MAX_NUM_OP], qa_sds_info[MAX_NUM_OP], qa_sds_nobs_info[MAX_NUM_OP];
  unsigned long bit_mask_arr[MAX_NUM_OP], mask_val_arr[MAX_NUM_OP];
  void **data_in, **data_out, *data_qa[MAX_NUM_OP];

  out_hdf_st = (strlen(arg_list[arg_cnt-2]) > 0) ? 1 : 0;
  if ((qa_fnames = (char **)Calloc2D(MAX_NUM_OP, MAX_PARAM_LENGTH, sizeof(char))) == NULL) {
    fprintf(stderr, "Cannot allocate memory for qa_fnames mask_nsds\n"); return -1;
  }
  n_op = (arg_cnt - 3)/3;
  st = get_parameters(arg_list, n_op, sel_qa_op, qa_fnames, qa_sds_info, bit_mask_arr, 
		      mask_val_arr, opt_arr, rel_op);
  if (st != -1)
    st = get_qa_sds_info(qa_fnames, qa_sds_info, qa_sdsc_info, fqa_l2g, n_op);
  if (st != -1)
    {
      if (out_hdf_st == 1) {
	if ((out_sds_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
	  fprintf(stderr, "Cannot allocate memory for out_sds_info in mask_nsds()\n");
      }
      if ((in_sds_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
      fprintf(stderr, "Cannot allocate memory for in_sds_info in mask_nsds()\n");
      if (fin_l2g == 1) {
	if ((in_sdsc_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
	  fprintf(stderr, "Cannot allocate memory for in_sdsc_info in mask_nsds()\n");
      }
      if (in_sds_info != NULL)
	{
	  for (isds=0; isds<nsds; isds++)
	    strcpy(in_sds_info[isds].name, sds_names[isds]);
	  st = get_in_sds_info(arg_list[0], in_sds_info, in_sdsc_info, &in_sds_nobs_info, 
			       fin_l2g, nsds);
	  for (isds=0; isds<nsds; isds++)
          {
            if (f_opt == 1) {
	      if ((fill_val < in_sds_info[isds].range[0]) || (fill_val > in_sds_info[isds].range[1]))
	      {
	        mask_fill[isds] = fill_val;
		if (fill_val == in_sds_info[isds].fill_val)
		  fprintf(stderr, "Warning: The mask fill value is equal to the input SDS fill value\n"); 
	      }
	      else 
	      {
	        fprintf(stderr, "\nInvalid mask fill value %d for SDS %s\n", fill_val, in_sds_info[isds].name);
	        fprintf(stderr, "Valid range of the SDS: (%d %d)\n", in_sds_info[isds].range[0], in_sds_info[isds].range[1]);
	        fprintf(stderr, "SDS fill value: %ld\n", in_sds_info[isds].fill_val);
	        tmp_val = fill_val;
	        if ((tmp_val >= in_sds_info[isds].range[0]) && (tmp_val <= in_sds_info[isds].range[1])) {
	          fprintf(stderr, "Input a new mask fill value for the SDS (outside of the valid range): ");
	          fscanf(stdin, "%d", &tmp_val);
	        }
	        mask_fill[isds] = tmp_val;
		if (tmp_val == in_sds_info[isds].fill_val)
                  fprintf(stderr, "Warning: The mask fill value is equal to the input SDS fill value\n");
		if ((tmp_val >= in_sds_info[isds].range[0]) && (tmp_val <= in_sds_info[isds].range[1]))
		  fprintf(stderr, "Warning: The mask fill value is within the valid range of the SDS\n"); 
	      }
            }
	    else {
	      switch(in_sds_info[isds].data_type)
	      {
		case 20: min = -128; max = 127; break; 
		case 21: min = 0; max = 255; break;
		case 22: min = -32768; max = 32767; break;
		case 23: min = 0; max = 65535; break;
		case 24: min = -2147483647 - 1; max = 2147483647; break; 
		case 25: min = 0; max = 4294967295u; break; 
	      }
	      diff_max = max - in_sds_info[isds].range[1];
	      diff_min = in_sds_info[isds].range[0] - min;
	      if (   ((diff_max == 0) && (diff_min == 0)) 
		  || ((diff_max == 1) && (diff_min == 0) && (in_sds_info[isds].fill_val == max))
		  || ((diff_min == 1) && (diff_max == 0) && (in_sds_info[isds].fill_val == min))) {
		fprintf(stderr, "Warning: Problem setting mask fill value for SDS %s\n", in_sds_info[isds].name);
		fprintf(stderr, "The SDS data type and maximum data value range of data type: ");
	        switch(in_sds_info[isds].data_type)
	        {
		  case 20: fprintf(stdout, "INT8 (%d %d)\n", min, max); break; 
		  case 21: fprintf(stdout, "UINT8 (%d %d)\n", min, max); break;
		  case 22: fprintf(stdout, "INT16 (%d %d)\n", min, max); break;
		  case 23: fprintf(stdout, "UINT16 (%d %d)\n", min, max); break;
		  case 24: fprintf(stdout, "INT32 (%d %d)\n", min, max); break; 
		  case 25: fprintf(stdout, "UINT32 (%d %d)\n", min, max); break; 
	        }
		fprintf(stderr, "The SDS valid range: (%d %d)\n", in_sds_info[isds].range[0], in_sds_info[isds].range[1]);
		fprintf(stderr, "The SDS fill value: %ld\n", in_sds_info[isds].fill_val);
		fprintf(stderr, "The mask fill value is being set to SDS fill value: %ld\n", in_sds_info[isds].fill_val);
		mask_fill[isds] = in_sds_info[isds].fill_val;
	      }
	      else {
		if (diff_max == 1) {
		  if (max != in_sds_info[isds].fill_val) 
		    mask_fill[isds] = max;
		  else {
		    mask_fill[isds] = min;
		    if (mask_fill[isds] == in_sds_info[isds].fill_val) 
		      mask_fill[isds] = mask_fill[isds] + 1;
		  }
		}
		else if (diff_min == 1) {
		  if (min != in_sds_info[isds].fill_val) 
		    mask_fill[isds] = min;
		  else {
		    mask_fill[isds] = max;
		    if (mask_fill[isds] == in_sds_info[isds].fill_val) 
		      mask_fill[isds] = mask_fill[isds] - 1;
		  }
		}
		else if (diff_max > 1) {
		  mask_fill[isds] = max;
		  if (mask_fill[isds] == in_sds_info[isds].fill_val) 
		    mask_fill[isds] = mask_fill[isds] - 1;
	        }
		else if (diff_min > 1) {
		  mask_fill[isds] = min;
		  if (mask_fill[isds] == in_sds_info[isds].fill_val) 
		    mask_fill[isds] = mask_fill[isds] + 1;
		}
		fprintf(stderr, "The mask fill value is set to: %d\n", mask_fill[isds]);
	      }
	    }
	  }
          
	  if (st != -1)
	    st = get_res_factors(&in_sds_info[0], qa_sds_info, n_op, res_l, res_s);
	  if (st != -1)
	    {
	      in_edge = (int32 **)Calloc2D(nsds, 4, sizeof(int32));
	      in_start = (int32 **)Calloc2D(nsds, 4, sizeof(int32));
	      out_edge = (int32 **)Calloc2D(nsds, 4, sizeof(int32));
	      out_start = (int32 **)Calloc2D(nsds, 4, sizeof(int32));
	      if (fin_l2g == 1)
		{
		  strcpy(sds_info.name, "nadd_obs_row"); sds_info.rank = 1; 
		  sds_info.sd_id = in_sds_info[0].sd_id; sds_info.sds_id = -1;
		  sds_info.dim_size[0] = in_sds_info[0].dim_size[0];
		  if ((data_in_nadd = (int32 *)calloc(sds_info.dim_size[0], sizeof(int32))) == NULL)
		    fprintf(stderr, "Cannot allocate memory for data_in_nadd in mask_nsds\n");
		  else
		    {
		      get_sds_data(&sds_info, data_in_nadd);
		      SDendaccess(sds_info.sds_id);
		    }
		}
	      
	      if (out_hdf_st == 1)
		for (isds=0; isds<nsds; isds++) {  
		  out_sds_info[isds].sd_id = out_sd_id; 
		  out_sds_info[isds].sds_id = -1; 
		}
	      create_out_sds(in_sds_info, out_sds_info, nsds, arg_list[arg_cnt], m_str, n, m, 
			     out_sd_id, out_hdf_st, mask_fill);

	      for (isds=0; isds<nsds; isds++)
		{
		  rank = in_sds_info[isds].rank;
		  if (rank == 2)
		    {
		      in_edge[isds][0] = out_edge[isds][0] = 1;
		      in_edge[isds][1] = out_edge[isds][1] = in_sds_info[isds].dim_size[1];
		    }
		  else if (rank > 2)
		    {
		      for (i=0; i<rank; i++)
			in_edge[isds][i] = out_edge[isds][i] = in_sds_info[isds].dim_size[i];
		      if ((n[isds] == -1) && (m[isds] == -1))
			{
			  if (in_sds_info[isds].dim_size[0] > in_sds_info[isds].dim_size[rank-1])
			    in_edge[isds][0] = out_edge[isds][0] = 1;
			  else in_edge[isds][rank-2] = out_edge[isds][rank-2] = 1;
			}
		      else
			{
			  out_edge[isds][0] = 1; out_edge[isds][1] = out_sds_info[isds].dim_size[1];
			  if (in_sds_info[isds].dim_size[0] > in_sds_info[isds].dim_size[rank-1])
			    { 
			      in_edge[isds][0] = 1;
			    }
			  else 
			    {
			      in_edge[isds][rank-2] = 1;
			    }
			}
		    }
		  compute_sds_start_offset(&in_sds_info[isds], n[isds], m[isds], &st_c[isds], &offset[isds]); 
		} /* for (isds=0; . . . ) */

	      get_ndata_vals(&in_sds_info[0], &bsq, &nrow, &ndata_in, &ndata_mask, &ndata_out, n[0], m[0]);
	      
	      for (isds=0, data_size=in_sds_info[0].data_size; isds<nsds; isds++)
		if (in_sds_info[isds].data_size > data_size) data_size = in_sds_info[isds].data_size;
	      if ((data_in = (void **)Calloc2D(nsds, ndata_in, data_size)) == NULL)
		fprintf(stderr, "Cannot allocate memory for data_in in mask_nsds()\n"); 
	      if ((data_mask = (uint8 *)calloc(ndata_mask, sizeof(uint8))) == NULL)
		fprintf(stderr, "Cannot allocate memory for data_mask in mask_nsds()\n"); 
	      if (out_hdf_st == 1) {
          if ((data_out = (void **)Calloc2D(nsds, ndata_out, data_size)) == NULL)
            fprintf(stderr, "Cannot allocate memory for data_out in mask_nsds()\n"); 
	      }
	      else data_out = NULL;
	      
	      if ((data_in != NULL) && (data_mask != NULL))
		{	
		  st = open_qa_sds_nsds(arg_list[0], in_sds_info, in_sdsc_info, &in_sds_nobs_info, nsds, qa_fnames,
					qa_sds_info, qa_sdsc_info, qa_sds_nobs_info, fqa_l2g, n_op);
		  if (st != -1)
		    st = malloc_qa_sds(qa_sds_info, n_op, fqa_l2g, data_qa, data_qa_nadd);
		  
		  if (st != -1)
		    {
		      rank = in_sds_info[0].rank;
		      for (i_op=0; i_op<=n_op; i_op++)
			if (fqa_l2g[i_op] == 1)
			  {
			    len = (int)strlen(qa_sdsc_info[i_op].name);
			    p1 = sd_charpos(qa_sdsc_info[i_op].name, '.', 0);
			    sd_strmid(qa_sdsc_info[i_op].name, p1+1, len-p1-1, num_str);
			    obs_num[i_op] = (int)atoi(num_str);
			  }
			else obs_num[i_op] = 1;
		      
		      for (irow=0; irow<nrow; irow++)
			{
			  read_qa_sds(qa_sds_info, qa_sdsc_info, qa_sds_nobs_info, n_op, data_qa, data_qa_nadd,
				      irow, res_l, fqa_l2g, obs_num);
			  
			  process_mask_data(data_qa, ndata_mask, qa_sds_info, n_op, sel_qa_op, bit_mask_arr, 
					    mask_val_arr, rel_op, res_s, data_mask, YES, NO, MASK_FILL);
			  
			  for (isds=0; isds<nsds; isds++)
			    {
			      if ((rank == 2) || (in_sds_info[isds].dim_size[0] > in_sds_info[isds].dim_size[rank-1]))
				in_start[isds][0] = irow;
			      else in_start[isds][rank-2] = irow;
			      if (fin_l2g == 1)
				{
				  len = (int)strlen(in_sdsc_info[isds].name);
				  p1 = sd_charpos(in_sdsc_info[isds].name, '.', 0);
				  sd_strmid(in_sdsc_info[isds].name, p1+1, len-p1-1, num_str);
				  obs_num_in = (int)atoi(num_str);
				  if (obs_num_in > 1)
				    read_sdsc_data(&in_sdsc_info[isds], &in_sds_nobs_info, data_in[isds], 
						   data_in_nadd, irow, obs_num_in);
				  else
				    if (SDreaddata(in_sds_info[isds].sds_id, in_start[isds], NULL, in_edge[isds], 
						   (VOIDP)data_in[isds]) == FAIL)
				      fprintf(stderr, "Cannot read data line from SDS %s in mask_sds()\n", 
					      in_sds_info[isds].name);
				}
			      else
				if (SDreaddata(in_sds_info[isds].sds_id, in_start[isds], NULL, in_edge[isds], 
					       (VOIDP)data_in[isds]) == FAIL)
				  fprintf(stderr, "Cannot read data line from SDS %s in mask_sds()\n", 
					  in_sds_info[isds].name);
			    } /* for (isds=0; . . .  */
			  
			  mask_nsds_data_row(data_in, data_out, data_mask, ndata_out, ndata_mask, nsds, bsq, 
					     st_c, offset, in_sds_info, mask_fill);
			  
			  if ((out_hdf_st == 1) && (out_sd_id != -1))
			    {
			      for (isds=0; isds<nsds; isds++)
				{
				  if ((rank>2) && (n[isds] == -1) && (m[isds] == -1) && 
				      (in_sds_info[isds].dim_size[rank-1] > in_sds_info[isds].dim_size[0]))
				    out_start[isds][rank-2] = irow;
				  else out_start[isds][0] = irow;
				  if (SDwritedata(out_sds_info[isds].sds_id, out_start[isds], NULL, out_edge[isds], 
						  data_out[isds]) == FAIL)
				    fprintf(stderr, "Cannot write data line to SDS %s in mask_nsds()\n", 
					    out_sds_info[isds].name);
				}
			    }
			} /* for (irow=0;  . .  ) */
		      
		      if ((m_opt == 1) && (out_hdf_st == 1))
			copy_metadata(in_sds_info[0].sd_id, out_sds_info[0].sd_id);
		      
		      /* close all SDS and HDF files */
		      if (fin_l2g == 1) SDendaccess(in_sds_nobs_info.sds_id);
		      for (i_op=0; i_op<=n_op; i_op++)
			if (fqa_l2g[i_op] == 1)
			  {
			    if ((qa_sds_nobs_info[i_op].sd_id == in_sds_nobs_info.sd_id) && 
				(qa_sds_nobs_info[i_op].sds_id == in_sds_nobs_info.sds_id)) ;
			    else
			      {
				for (j_op=0; j_op<i_op; j_op++)
                                  /* The following if statement does nothing...
				  if ((qa_sds_nobs_info[i_op].sd_id == qa_sds_nobs_info[j_op].sd_id) &&
				      (qa_sds_nobs_info[i_op].sds_id == qa_sds_nobs_info[j_op].sds_id)) ;
                                  */
				if (j_op >= i_op) SDendaccess(qa_sds_nobs_info[i_op].sds_id);
			      }
			  }
		      
		      for (i_op=0; i_op<=n_op; i_op++)
			if (fqa_l2g[i_op] == 1)
			  {
			    for (j_op=0; j_op<i_op; j_op++)
			      if (qa_sds_info[i_op].sd_id == qa_sds_info[j_op].sd_id) break;
			    if (j_op >= i_op) free(data_qa_nadd[i_op]);
			  }
		      
		      if (out_hdf_st == 1) 
			for (isds=0; isds<nsds; isds++)
			  if (out_sds_info[isds].sds_id != -1) SDendaccess(out_sds_info[isds].sds_id);
		      close_qa_hdf_nsds(arg_list[0], in_sds_info, nsds, qa_fnames, qa_sds_info, n_op);
		      for (isds=0; isds<nsds; isds++)
			if (in_sds_info[isds].sds_id != -1) 
			  {
			    SDendaccess(in_sds_info[isds].sds_id);
			    if (fin_l2g == 1)  SDendaccess(in_sdsc_info[isds].sds_id);
			    for (jsds=isds+1; jsds<nsds; jsds++)
			      if (in_sds_info[isds].sds_id == in_sds_info[jsds].sds_id)
				{
				  in_sds_info[jsds].sds_id = -1; 
				  if (fin_l2g == 1)  SDendaccess(in_sdsc_info[jsds].sds_id);
				}
			    in_sds_info[isds].sds_id = -1;
			    if (fin_l2g == 1)
			      in_sdsc_info[isds].sds_id = -1;
			  }
		      SDend(in_sds_info[0].sd_id);
		      
		      if (data_qa[0] != NULL) free(data_qa[0]);
		      for (i=1; i<=n_op; i++)
			{
			  
			  for (i=0; i<=n_op; i++)
			    {
			      p1 = sd_charpos(qa_sds_info[i].name, '.', 0);
			      if (p1 != -1)
				{
				  sd_strmid(qa_sds_info[i].name, 0, p1, org_sds_name);
				  strcpy(qa_sds_info[i].name, org_sds_name);
				}
			    }
			  
			  for (j=0; j<i; j++)
			    if ((strcmp(qa_fnames[i], qa_fnames[j]) == 0) &&
				(strcmp(qa_sds_info[i].name, qa_sds_info[j].name) == 0)) break;
			  if (j >= i) { if (data_qa[i] != NULL) free(data_qa[i]); }
			} /* for (i=1; . . .  */
		    } /* if (st != -1) . . .  */
		  free(data_mask);
		  Free2D((void **)data_in);
		  if (out_hdf_st == 1)
		    Free2D((void **)data_out);
		} /* if ((data_in != NULL) && . . .  ) */
	      Free2D((void **)in_edge);
	      Free2D((void **)in_start);
	      Free2D((void **)out_edge);
	      Free2D((void **)out_start);
	    } /* if (st != -1)  */
	  free(in_sds_info);
	  if (out_hdf_st == 1) free(out_sds_info);
	  if (fin_l2g == 1) { free(in_sdsc_info);  free(data_in_nadd); }
	} /* if (out_sds_info != . . ) */
    } /* if (st != -1) */ 
  Free2D((void **)qa_fnames);
  return st;
}

void mask_nsds_data_row(void **data_in, void **data_out, uint8 *data_mask, 
			int ndata_out, int ndata_mask, int nsds, int bsq, int *st_c, 
			int *offset, sds_t *sds_info, int *mask_fill)
/*
!C*********************************************************************************

!Function: mask_nsds_data_row
       
!Description:
   Function mask_nsds_data_row masked input sds by the mask and put the result
   into output sds.

!Input/output Parameters:
  data_in    input data to be masked
  data_out   masked data
  data_mask  mask to be used.
  ndata_out  dimension size of the output data
  ndata_mask dimension size of the mask
  nsds       count of input SDS
  bsq        indicator of whether the input SDS is band sequence or not
  st_c       starting indices of subset to be extracted.
  offset     offset of the subset to be extracted. 
  sds_info   Array of input SDS information structure.

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/

{
  int sds_val = 0;
  int i, j, k;
  int k1, jj, isds;
  int ii[MAX_NSDS];
  int ic[MAX_NSDS];

  k = ndata_out/ndata_mask;
  if (bsq == 1)
  {
    for (isds=0; isds<nsds; isds++)
      ii[isds] = st_c[isds];
    for (k1=0, j=0; k1<k; k1++)
    {
      for (i=0; i<ndata_mask; i++, j++)
      {
	for (isds=0; isds<nsds; isds++)
	{
	  jj = ii[isds];
          if (data_mask[i] == YES)
          {
            switch(sds_info[isds].data_type)
            {
              case 20: ((int8 **)data_out)[isds][j] = ((int8 **)data_in)[isds][jj]; break;
              case 21: ((uint8 **)data_out)[isds][j] = ((uint8 **)data_in)[isds][jj]; break;
              case 22: ((int16 **)data_out)[isds][j] = ((int16 **)data_in)[isds][jj]; break;
              case 23: ((uint16 **)data_out)[isds][j] = ((uint16 **)data_in)[isds][jj]; break;
              case 24: ((int32 **)data_out)[isds][j] = ((int32 **)data_in)[isds][jj]; break;
              case 25: ((uint32 **)data_out)[isds][j] = ((uint32 **)data_in)[isds][jj]; break;
	      default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
            }
          }
	  else if (data_mask[i] == MASK_FILL)
          {
	    switch(sds_info[isds].data_type)
            {
              case 20: sds_val = ((int8 **)data_in)[isds][jj]; break;
              case 21: sds_val = ((uint8 **)data_in)[isds][jj]; break;
              case 22: sds_val = ((int16 **)data_in)[isds][jj]; break;
              case 23: sds_val = ((uint16 **)data_in)[isds][jj]; break;
              case 24: sds_val = ((int32 **)data_in)[isds][jj]; break;
              case 25: sds_val = ((uint32 **)data_in)[isds][jj]; break;
              default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
            }
	    if (sds_val != sds_info[isds].fill_val)
	    {
              switch(sds_info[isds].data_type)
              {
                case 20: ((int8 **)data_out)[isds][j] = (int8)mask_fill[isds]; break;
                case 21: ((uint8 **)data_out)[isds][j] = (uint8)mask_fill[isds]; break;
                case 22: ((int16 **)data_out)[isds][j] = (int16)mask_fill[isds]; break;
                case 23: ((uint16 **)data_out)[isds][j] = (uint16)mask_fill[isds]; break;
                case 24: ((int32 **)data_out)[isds][j] = (int32)mask_fill[isds]; break;
                case 25: ((uint32 **)data_out)[isds][j] = (uint32)mask_fill[isds]; break;
                default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
              }
	    }
	    else
	    {
              switch(sds_info[isds].data_type)
              {
                case 20: ((int8 **)data_out)[isds][j] = (int8)sds_info[isds].fill_val; break;
                case 21: ((uint8 **)data_out)[isds][j] = (uint8)sds_info[isds].fill_val; break;
                case 22: ((int16 **)data_out)[isds][j] = (int16)sds_info[isds].fill_val; break;
                case 23: ((uint16 **)data_out)[isds][j] = (uint16)sds_info[isds].fill_val; break;
                case 24: ((int32 **)data_out)[isds][j] = (int32)sds_info[isds].fill_val; break;
                case 25: ((uint32 **)data_out)[isds][j] = (uint32)sds_info[isds].fill_val; break;
                default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
              }
	    }
          }
          else
	  {
            switch(sds_info[isds].data_type)
            {
              case 20: ((int8 **)data_out)[isds][j] = (int8)mask_fill[isds]; break;
              case 21: ((uint8 **)data_out)[isds][j] = (uint8)mask_fill[isds]; break;
              case 22: ((int16 **)data_out)[isds][j] = (int16)mask_fill[isds]; break;
              case 23: ((uint16 **)data_out)[isds][j] = (uint16)mask_fill[isds]; break;
              case 24: ((int32 **)data_out)[isds][j] = (int32)mask_fill[isds]; break;
              case 25: ((uint32 **)data_out)[isds][j] = (uint32)mask_fill[isds]; break;
	      default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);	
            }
	  }
	  ii[isds] += offset[isds];
	} /* for (isds=0; . . .) */
      } /* for (i=0; . . . ) */
    } /* for (k1=0; . .  .) */
  }
  else
  {
    for (isds=0; isds<nsds; isds++)
    {
      ic[isds] = 0;
      ii[isds] = st_c[isds];
    }
    for (i=0; i<ndata_mask; i++)
    {
      for (isds=0; isds<nsds; ++ic[isds], isds++)
      {
        j = ic[isds];
        for (k1=0; k1<k; k1++)
        {
	  jj = ii[isds];
          if (data_mask[i] == YES)
          {
            switch(sds_info[isds].data_type)
            {
              case 20: ((int8 **)data_out)[isds][j] = ((int8 **)data_in)[isds][jj]; break;
              case 21: ((uint8 **)data_out)[isds][j] = ((uint8 **)data_in)[isds][jj]; break;
              case 22: ((int16 **)data_out)[isds][j] = ((int16 **)data_in)[isds][jj]; break;
              case 23: ((uint16 **)data_out)[isds][j] = ((uint16 **)data_in)[isds][jj]; break;
              case 24: ((int32 **)data_out)[isds][j] = ((int32 **)data_in)[isds][jj]; break;
              case 25: ((uint32 **)data_out)[isds][j] = ((uint32 **)data_in)[isds][jj]; break;
	      default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);	
            }
          }
          else if (data_mask[i] == MASK_FILL)
          {
            switch(sds_info[isds].data_type)
            {
              case 20: sds_val = ((int8 **)data_in)[isds][jj]; break;
              case 21: sds_val = ((uint8 **)data_in)[isds][jj]; break;
              case 22: sds_val = ((int16 **)data_in)[isds][jj]; break;
              case 23: sds_val = ((uint16 **)data_in)[isds][jj]; break;
              case 24: sds_val = ((int32 **)data_in)[isds][jj]; break;
              case 25: sds_val = ((uint32 **)data_in)[isds][jj]; break;
              default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
            }
            if (sds_val != sds_info[isds].fill_val)
            {
              switch(sds_info[isds].data_type)
              {
                case 20: ((int8 **)data_out)[isds][j] = (int8)mask_fill[isds]; break;
                case 21: ((uint8 **)data_out)[isds][j] = (uint8)mask_fill[isds]; break;
                case 22: ((int16 **)data_out)[isds][j] = (int16)mask_fill[isds]; break;
                case 23: ((uint16 **)data_out)[isds][j] = (uint16)mask_fill[isds]; break;
                case 24: ((int32 **)data_out)[isds][j] = (int32)mask_fill[isds]; break;
                case 25: ((uint32 **)data_out)[isds][j] = (uint32)mask_fill[isds]; break;
                default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
              }
            }
	    else
            {
              switch(sds_info[isds].data_type)
              {
                case 20: ((int8 **)data_out)[isds][j] = (int8)sds_info[isds].fill_val; break;
                case 21: ((uint8 **)data_out)[isds][j] = (uint8)sds_info[isds].fill_val; break;
                case 22: ((int16 **)data_out)[isds][j] = (int16)sds_info[isds].fill_val; break;
                case 23: ((uint16 **)data_out)[isds][j] = (uint16)sds_info[isds].fill_val; break;
                case 24: ((int32 **)data_out)[isds][j] = (int32)sds_info[isds].fill_val; break;
                case 25: ((uint32 **)data_out)[isds][j] = (uint32)sds_info[isds].fill_val; break;
                default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
              }
            }
	  }
          else
	  {
            switch(sds_info[isds].data_type)
            {
              case 20: ((int8 **)data_out)[isds][j] = (int8)mask_fill[isds]; break;
              case 21: ((uint8 **)data_out)[isds][j] = (uint8)mask_fill[isds]; break;
              case 22: ((int16 **)data_out)[isds][j] = (int16)mask_fill[isds]; break;
              case 23: ((uint16 **)data_out)[isds][j] = (uint16)mask_fill[isds]; break;
              case 24: ((int32 **)data_out)[isds][j] = (int32)mask_fill[isds]; break;
              case 25: ((uint32 **)data_out)[isds][j] = (uint32)mask_fill[isds]; break;
	      default: fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", sds_info[isds].data_type);
            }
	  }
	  ii[isds] += offset[isds];
        } /* for (k1=0; . . . ) */
      } /* for (isds=0; . . .) */ 
    } /* for (i=0; . .  .) */
  }
}
