/****************************************************************************
!C

!File: reduce_sds_rank.c

!Description:
  Contains routines for converting a mutlidimensional SDS into 
  many SDS of two dimension

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  
  Version 1.0 September, 2002

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
#include "meta.h"

#define HELP \
"NAME \n" \
"    reduce_sds_rank - convert one or more 3D/4D SDS from a MODIS Land\n"\
"                      HDF-EOS data product to many 2D SDSs. \n" \
"SYNOPSIS \n" \
"    reduce_sds_rank -help [filename] \n" \
"    reduce_sds_rank -of=<output filename>\n"\
"                    {[-sds=<SDSname> [-dim=<dimstr> [-dim=<dimstr>. . ]]]}\n"\
"                    [-all] [-meta] filename  \n" \
" \n" \
"DESCRIPTION \n" \
"    Several MODIS Land HDF-EOS data products (e.g., MOD43, MYD43) and\n"\
"    related MODIS products (e.g., MOD35) contain multidimensional SDSs.\n"\
"    This tool converts one or more multidimensional (3D or 4D) SDS to a\n"\
"    series of 2D HDF SDSs. Specific SDS layers may be selected using the\n"\
"    -sds and -dim options for each input 3D or 4D SDS.\n" \
" \n" \
"    The output file SDS names reflect the input SDS name, the dimension\n"\
"    name and the dimension element numbers. For example: \n" \
"        BRDF_Albedo_Parameters.Num_Land_Bands_Plus3_3.Num_Parameters_1 \n" \
"        (parameter 1 for land band 3 of SDS BRDF_Albedo_Parameters in\n"\
"        MOD43B1. Note that BRDF_Albedo_Parameters is the input SDS\n"\
"        name, Num_Land_Bands_Plus3 and Num_Parameters are the 3rd and\n"\
"        the 4th dimension names) \n" \
"        Surface_Refl.Num_Obs_Max_1.Num_Land_Bands_2 (1st observation of\n"\
"        land band 2 of SDS Surface_Refl in MODAGAGG.  Note that\n"\
"        Surface_Refl is the input SDS name, Num_Obs_Max and Num_Land_Bands\n"\
"        are the 3rd and the 4th dimension names)\n"\
"        Angles.Num_Obs_Max_2.Num_Angles_3 (2nd observation of the 3rd angle\n"\
"        component of SDS Angles in MODAGAGG.  Note that Angles is the input\n"\
"        SDS name, Num_Obs_Max and Num_Angles are the 3rd and the 4th\n"\
"        dimension names) \n" \
" \n" \
"    This tool supports 3D and 4D SDS(s). \n" \
" \n" \
"    The tool command arguments can be specified in any order. \n" \
" \n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If the input\n"\
"                             filename is provided with this option, then\n"\
"                             the names of all the SDS in the file with the\n" \
"                             SDS dimension name and size are displayed. \n" \
"    -of=<filename>           Output filename. \n" \
"    -sds=<SDS name>          Name of SDS to convert. By default all 3D/4D\n"\
"                             SDS in the input file are converted to 2D SDSs.\n"\
"    -dim=<dimension name, dimension element number(s)> \n" \
"                             Name of the 3rd or higher SDS dimension and \n" \
"                             the dimension element number(s) (1-based). The\n"\
"                             dimension element numbers can be separated by\n"\
"                             comma or by '-'. The -dim option must be\n"\
"                             repeated for each of the defined SDS dimension\n"\
"                             names, e.g., -dim=Band_id,1-7 -dim=obs_id,1,2.\n"\
"                             This option should follow the -sds option for\n"\
"                             each input SDS. By default, the tool reduces\n"\
"                             all of the 3rd and 4th dimensions of a\n"\
"                             specified SDS. \n" \
"    -all                     Create an additional output SDS containing all\n"\
"                             the 2D SDSs mosaiced in a single SDS in\n"\
"                             row-column order. \n" \
"    -meta                    Copy metadata from input file to output file.\n"\
"    filename                 Input filename. \n" \
" \n" \
"EXAMPLES \n" \
"    reduce_sds_rank -sds=BRDF_Albedo_Parameters\n"\
"                    -dim=Num_Land_Bands_Plus3,1-7\n"\
"                    -dim=Num_Parameters,1-3 -of=brdf_albedo_2dsds.hdf \n" \
"                    MYD43B1.A2002177.h11v11.003.2002210233848.hdf \n" \
"        {Note: This example extracts twenty-one 2D SDSs from the 4D SDS \n" \
"               BRDF_Albedo_Parameters stored in the input MYD43B1file. In\n"\
"               this example, the values of the BRDF_Albedo_Parameters 1, 2\n"\
"               and 3 are written as separate 2D SDSs for each of the first\n"\
"               7 land bands. } \n" \
" \n" \
"    reduce_sds_rank -sds=Angles -dim=Num_Obs_Max,1-4 -dim=Num_Angles,1-4\n"\
"                    -all -of=agg_angles.hdf\n"\
"                    MODAGAGG.A2002017.h20v11.004.2002203181248.hdf \n" \
"        {Note: This example extracts sixteen 2D SDSs from the 4D SDS Angles\n"\
"               stored in the input MODAGAGG file. In this example, the\n"\
"               values of the four angles (View Zenith, Solar Zenith, Solar\n"\
"               Azimuth, View Azimuth) are written as separate 2D SDSs for\n"\
"               each of the four observations). An additional output SDS\n"\
"               containing all the 2D SDSs in a single SDS in row-column\n"\
"               order is also output. }  \n" \
" \n" \
"    reduce_sds_rank -sds=Surface_Refl -dim=Num_Obs_Max,1,2 \n" \
"                    -dim=Num_Land_Bands,1-7 -sds=Band_QC\n"\
"                    -dim=Num_Obs_Max,1,4 -dim=Num_Band_QC,1-7\n"\
"                    -of=agg_angles_bandqc.hdf \n" \
"                    MODAGAGG.A2002017.h20v11.004.2002203181248.hdf \n" \
" \n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 08/08/2002 \n" \

#define USAGE \
"usage:	reduce_sds_rank -help [filename] \n" \
"       reduce_sd_rank -of=<output filename> {[-sds=<SDSname>\n"\
"       [-dim=<dimstr> [-dim=<dimstr>. . ]]]} [-all] [-meta] filename \n" \
" \n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If the input\n"\
"                             filename is provided with this option, then\n"\
"                             the names of all the SDS in the file with the\n" \
"                             SDS dimension name and size are displayed. \n" \
"    -of=<filename>           Output filename. \n" \
"    -sds=<SDS name>          Name of SDS to convert. By default all 3D/4D\n"\
"                             SDS in the input file are converted to 2D SDSs.\n"\
"    -dim=<dimension name, dimension element number(s)> \n" \
"                             Name of the 3rd or higher SDS dimension and \n" \
"                             the dimension element number(s) (1-based). The\n"\
"                             dimension element numbers can be separated by\n"\
"                             comma or by '-'. The -dim option must be\n"\
"                             repeated for each of the defined SDS dimension\n"\
"                             names, e.g., -dim=Band_id,1-7 -dim=obs_id,1,2.\n"\
"                             This option should follow the -sds option for\n"\
"                             each input SDS. By default, the tool reduces\n"\
"                             all of the 3rd and 4th dimensions of a\n"\
"                             specified SDS. \n" \
"    -all                     Create an additional output SDS containing all\n"\
"                             the 2D SDSs mosaiced in a single SDS in\n"\
"                             row-column order. \n" \
"    -meta                    Copy metadata from input file to output file.\n"\
"    filename                 Input filename. \n" \
" \n"


/******************************************************************************************
                            Prototypes.
******************************************************************************************/
int parse_cmd_reduce_sds_rank(int argc, char **argv, int *if_cnt, char **sds_names, 
			      int *sds_cnt, char ***dim_ids_nums, int *dim_ids_cnt, 
			      char *out_fname, int *out_all, int *m_opt);
int get_out_dim_id(char **dim_names, int32 *dim_size, int rank, char **dim_ids, 
	int dim_cnt, int **out_dim_ids, int *out_dim_ids_cnt);
void cvrt_rank_m2(int32 in_sd_id, int32 out_sd_id, char *sds_name, char **dim_ids_nums,
	int dim_ids_cnt, int out_all);

/*****************************************************************************************/

int main(int argc, char *argv[])
/*
!C*****************************************************************************************

!Description:
  Main is the main program for the reduce_sds_rank routine. It converts one or more 3D/4D 
  SDS from a MODISLand HDF file to many 2D SDSs.

!Input Parameters: (none)

!Output Parameters:
  (returns)  Completion status:
             0 - successful completion
             1 - error exit

!Revision History: (none)

!Team-unique Header: (see file prolog)

!References and Credits: (see file prolog)

!Design Notes: 

   If the user doesn't know the SDS information of a file, type 
   
     reduce_sds_rand -help filename

   The program will print out the names of all the SDSs in the file with the SDS dimension
   name and dimension size.

!END
********************************************************************************************/

{
  int i;
  int if_cnt;
  int fid;
  int isds;
  int m_opt;
  int status;
  int sds_cnt;
  int out_all;
  int is_all_sds;
  int32 in_sd_id;
  int32 out_sd_id;
  char **sds_names;
  char ***dim_ids_nums;
  int dim_ids_cnt[MAX_NUM_SDS];
  char in_fname[MAX_PATH_LENGTH];
  char out_fname[MAX_PATH_LENGTH];
  sds_t *in_sds_info;
  char temp_sds_names[MAX_STR_LEN];

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
  
  if ((argc > 3) && (strcmp(argv[1], "-help")==0))
    {
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_SUCCESS);
    }                
  
  /* Display SDS names of input HDF file */ 
  
  if ((argc==3) && (strcmp(argv[1], "-help")==0))
    {
      if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
        fprintf(stderr, "Cannot allocate memory for sds_names in reduce_sds_rank (main)\n");
      else
        {
          sds_cnt = get_sds_names(argv[2], sds_names);
          if (sds_cnt == 0)
            {
              fprintf(stderr, "No SDS found in %s\n", argv[2]);
              exit(EXIT_FAILURE);
            }
          else 
            {
              fprintf(stderr, "File %s contains the following SDSs : \n", argv[2]);
              i = 1;
	      in_sds_info = (sds_t *)calloc(i, sizeof(sds_t));
              for (isds=0; isds < sds_cnt; isds++ )
                {
                  fprintf(stderr, "SDS name : %s \n", sds_names[isds]);
                  strcpy(in_sds_info->name, sds_names[isds]);
		  in_sds_info->sds_id = -1;
		  in_sds_info->sd_id = SDstart(argv[2], DFACC_READ);
		  in_sds_info->sds_index = SDnametoindex(in_sds_info->sd_id, sds_names[isds]);
		  if ((in_sds_info->sds_id = SDselect(in_sds_info->sd_id, 
						      in_sds_info->sds_index)) == FAIL)
		    {
		      fprintf(stderr, "Cannot open the SDS %s \n", sds_names[isds]);
		    }
		  else 
		    {
		      if (SDgetinfo(in_sds_info->sds_id, temp_sds_names, &in_sds_info->rank,
				    in_sds_info->dim_size, &in_sds_info->data_type, 
				    &in_sds_info->nattr) == FAIL)
			{
			  fprintf(stderr, "Cannot get information for the SDS %s \n", sds_names[isds]);
			}
		      else
			{
			  print_sds_dim_size(in_sds_info);
			}
		    }
		}
	    }
	}
      exit(EXIT_SUCCESS);
    }
  
  /* initialize dim_ids_cnt */
  
  for (i=0; i<MAX_NUM_SDS; i++)
    {
      dim_ids_cnt[i] = 0;
    }
  
  if ((dim_ids_nums = (char ***)Calloc3D(MAX_NUM_SDS, MAX_NUM_DIM, MAX_DIM_NAME_LEN, sizeof(char))) == NULL)  
    fprintf(stderr, "Cannot allocate memory for dim_ids_nums in reduce_sds_rank(main)\n");
  
  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL) 
    fprintf(stderr, "Cannot allocate memory for sds_names in reduce_sds_rank(main)\n");
  
  if ((sds_names != NULL) && (dim_ids_nums != NULL))
    {
      status = parse_cmd_reduce_sds_rank(argc, argv, &if_cnt, sds_names, &sds_cnt, dim_ids_nums,
					 dim_ids_cnt, out_fname, &out_all, &m_opt);
      if (status == -1)
	fprintf(stderr, "%s\n", USAGE);
      else if (status != 0)
	{
	  if ((out_sd_id = SDstart(out_fname, DFACC_CREATE)) == FAIL)
	    fprintf(stderr, "Cannot create output hdf file %s\n", out_fname);
	  else 
	    {
	      if ((sds_cnt == 1) && (strcmp(sds_names[0], "all") == 0)) is_all_sds = 1;
	      else is_all_sds = 0;
	      for (fid=1; fid<argc; fid++)
		{
		  strcpy(in_fname, argv[fid]);
		  if (in_fname[0] == '-') continue;
		  if (is_all_sds == 1)
		    sds_cnt = get_sds_names(in_fname, sds_names);
		  if (sds_cnt != 0)
		    {
		      if ((in_sd_id = SDstart(in_fname, DFACC_READ)) == FAIL)
			fprintf(stderr, "Cannot open input HDF file %s\n", in_fname);
		      else
			{
			  fprintf(stderr, "\nProcessing input file: %s\n", in_fname);
			  fprintf(stderr, "------------------------------------------------------------\n");
			  for (isds=0; isds<sds_cnt; isds++)
			    cvrt_rank_m2(in_sd_id, out_sd_id, sds_names[isds], dim_ids_nums[isds], 
					 dim_ids_cnt[isds], out_all);
			  if (m_opt == 1)
			    copy_metadata(in_sd_id, out_sd_id);
			  SDend(in_sd_id);
			}
		    }
		} /* for (fid=0; . . . */
	      SDend(out_sd_id);
	    }
	}
    }
  Free2D((void **)sds_names);
  Free3D((void ***)dim_ids_nums);
  fprintf(stderr, "Processing done ! \n");
  return 0;
}

int parse_cmd_reduce_sds_rank(int argc, char **argv, int *if_cnt, char **sds_names, int *sds_cnt, 
			      char ***dim_ids_nums, int *dim_ids_cnt, char *out_fname, int *out_all, 
			      int *m_opt)
/*
!C******************************************************************************
    
!Function: parse_cmd_reduce_sds_rank
       
!Description:
    
  Parse the command line arguments into corresponding fields.
    
!Input Parameters: 
  argc        input argument count
  argv        input argument vector
    
!Input/output Parameters:
  if_cnt         count of input files.
  sds_names      string array containing the list of user input sds names
  sds_cnt        count of user input sds names
  dim_ids_nums   dimension element numbers
  dim_ids_cnt    count of dimension element numbers
  out_fname      output filename
  out_all        if out_all is specified, then all SDS in the HDF file are processsed 
                 and a mosaic of all the SDS is created in the output file.
  m_opt          if "-meta" option is specified by user, m_opt is set to 1, all meta 
                 data name will be copied from input file to output file.
    
!Output Parameters:
    
  status:     1 - successful completion
              0 - error exit   
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*********************************************************************************/
{
  int i;
  int status;
  int dim_cnt;
  
  status = 1;
  out_fname[0] = '\0';
  *m_opt = *if_cnt = 0;
  *sds_cnt = *out_all = 0;
  i = 1;
  while (i<argc)
    {
      if (is_arg_id(argv[i], "-sds=") == 0)
	{
	  get_arg_val(argv[i], sds_names[*sds_cnt]);
	  i++;
	  dim_cnt = 0;
	  while ((i<argc) && (is_arg_id(argv[i], "-dim") == 0))
	    {
	      get_arg_val(argv[i], dim_ids_nums[*sds_cnt][dim_cnt]);
	      dim_cnt++;
	      i++;
	    }
	  i--;
	  dim_ids_cnt[*sds_cnt] = dim_cnt;
	  ++*sds_cnt;
	}
      else if (strcmp(argv[i], "-sds") == 0) ;
      else if (is_arg_id(argv[i], "-of=") == 0)
	get_arg_val(argv[i], out_fname);
      else if (strcmp(argv[i], "-all") == 0) *out_all = 1;
      else if (strcmp(argv[i], "-meta") == 0) *m_opt = 1;
      else if (argv[i][0] == '-')
	fprintf(stderr, "Unknown option %s\n", argv[i]);
      else ++*if_cnt;
      i++;
    }
  if ((*if_cnt == 0) || (out_fname[0] == '\0')) status = -1;
  if (*if_cnt == 0) fprintf(stderr, "Missing input filename\n");
  if (out_fname[0] == '\0') fprintf(stderr, "Missing output filename\n");
  if ((status != -1) && (*sds_cnt == 0))
    {
      fprintf(stderr, "No SDS name input. Reading all SDS . . \n");
      *sds_cnt = 1;
      strcpy(sds_names[0], "all");
    }
  return status;
}

void cvrt_rank_m2(int32 in_sd_id, int32 out_sd_id, char *sds_name, char **dim_ids_nums,
		  int dim_ids_cnt, int out_all)
/*
!C******************************************************************************
    
!Function: cvrt_rank_m2
        
!Description:
  This routines converts the input SDS to 2D SDSs. Also, if the out_all is specifed, 
  a mosaic of all the 2D SDSs from the same input SDS is created in the output file.
    
!Input Parameters: 
  in_sd_id        input file descriptor.
  out_sd_id       output file descriptor.
  sds_names       user input SDS names.
  dim_ids_nums    dimension element numbers.
  dim_ids_cnt     dimension element counts.
  out_all         if the out_all is specifed, a mosaic of all the 2D SDSs from 
                  the same input SDS is created in the output file.
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
********************************************************************************/
{
  int offset_is;
  int i, j1, j2;
  int is, iline, bsq;
  int offset = 0, p1;
  int tot_num_vals;
  int idim, isds, nsds;
  int32 num_data;
  int32 data_size;
  int32 edge1[2] = {0, 0};
  int32 start1[2] = {0, 0};
  int32 edgef[2] = {0, 0};
  int32 startf[2] = {0, 0};
  int32 edge[4] = {0, 0, 0, 0};
  int32 start[4] = {0, 0, 0, 0};
  char tmp_name[MAX_NAME_LENGTH];
  sds_t in_sds_info;
  sds_t out_sdsf_info;
  sds_t *out_sds1_info;
  char **dim_names, **short_dim_names;
  void *sds_data_in;
  void **sds_data_out;
  int **out_dim_ids;
  int out_dim_ids_cnt[MAX_NUM_DIM];
  char tmp_dname1[MAX_DIM_NAME_LEN];
  char tmp_dname2[MAX_DIM_NAME_LEN];

  memset( &out_sdsf_info, 0, sizeof(sds_t) );
  
  in_sds_info.sds_id = -1;
  in_sds_info.sd_id = in_sd_id;
  strcpy(in_sds_info.name, sds_name);
  fprintf(stderr, "	Processing SDS: %s\n", sds_name);
  
  /* initalize array */
  for (i=0; i< MAX_NUM_DIM; i++)
    {
      out_dim_ids_cnt[i] = 0;
    }
  
  if (get_sds_info((char *)NULL, &in_sds_info) != -1)
    {
      data_size = in_sds_info.data_size;
      if ((dim_names = (char **)Calloc2D(in_sds_info.rank, MAX_NAME_LENGTH, sizeof(char))) 
	  == NULL)
	fprintf(stderr, "Cannot allocate memory for dim_names\n");
      
      if ((short_dim_names = (char **)Calloc2D(in_sds_info.rank, MAX_NAME_LENGTH, sizeof(char))) 
	  == NULL)
	fprintf(stderr, "Cannot allocate memory for short_dim_names\n");
      
      if ((out_dim_ids = (int **)Calloc2D(MAX_NUM_DIM, MAX_NUM_DIM, sizeof(int))) == NULL)
	fprintf(stderr, "Cannot allocate memory for out_dim_ids\n");
      
      if ((dim_names == NULL) || (short_dim_names == NULL) || (out_dim_ids == NULL))
	{
	  Free2D((void **)dim_names);
	  Free2D((void **)short_dim_names);
	  Free2D((void **)out_dim_ids);
	  return;
	}
      else 
	{      
	  get_sds_dim_name(&in_sds_info, dim_names, short_dim_names);
	}
            
      if (get_out_dim_id(short_dim_names, in_sds_info.dim_size, in_sds_info.rank, dim_ids_nums, 
			 dim_ids_cnt, out_dim_ids, out_dim_ids_cnt) == -1)
	{
	  Free2D((void **)dim_names);
	  return;
	}

      if (in_sds_info.rank > 1) 
	{
          bsq = ((in_sds_info.rank == 2) || (in_sds_info.dim_size[0] < in_sds_info.dim_size[in_sds_info.rank - 1])) ? 1 : 0;   
	  if (bsq  == 0)
	    {	
	      for (i=2,nsds=1,num_data=in_sds_info.dim_size[1]; i<in_sds_info.rank; i++)
		{
		  nsds *= out_dim_ids_cnt[i-2];
		  num_data *= in_sds_info.dim_size[i];
		}
	    }
	  else 
	    {
	      for (i=0,nsds=1,num_data=in_sds_info.dim_size[in_sds_info.rank - 1]; i<in_sds_info.rank - 2; i++)
		{
		  nsds *= out_dim_ids_cnt[i];
		  num_data *= in_sds_info.dim_size[i];
		}
	    }	
	  
	  /* allocate memory for variables out_sds1_info, sds_data_in and sds_data_out which holds 
	     the data for input sds and output sds */
	  
	  if ((out_sds1_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL) 
	    fprintf(stderr, "Cannot allocate memory for out_sds1_info\n");
	  if ((sds_data_in = (void *)calloc(num_data, data_size)) == NULL)
	    fprintf(stderr, "Cannot allocate memory for sds_data_in\n");
	  
	  if (bsq == 0)
	    {
	      if ((sds_data_out = (void **)Calloc2D(nsds, in_sds_info.dim_size[1], data_size)) == NULL)
		fprintf(stderr, "Cannot allocate memory sds_data_out\n");
	    }
	  else
	    {
	      if ((sds_data_out = (void **)Calloc2D(nsds, in_sds_info.dim_size[in_sds_info.rank - 1], data_size)) == NULL)
		fprintf(stderr, "Cannot allocate memory sds_data_out\n");
	    }
	  
	  if ((dim_names != NULL) && (sds_data_in != NULL) && (sds_data_out != NULL) 
	      && (out_sds1_info != NULL))
	    {
	      if (bsq == 0)
		{
		  for (isds=0; isds<nsds; isds++)
		    {
		      out_sds1_info[isds].rank = 2;
		      out_sds1_info[isds].sd_id = out_sd_id;
		      out_sds1_info[isds].fill_val = in_sds_info.fill_val;
		      out_sds1_info[isds].data_type = in_sds_info.data_type;
		      
		      
		      out_sds1_info[isds].dim_size[0] = in_sds_info.dim_size[0];
		      out_sds1_info[isds].dim_size[1] = in_sds_info.dim_size[1];
		      for (idim=2; idim<in_sds_info.rank; idim++)
			{
			  p1 = sd_charpos(dim_names[idim], ':', 0);
			  if (p1 != -1)
			    {
			      if (idim == 2) sd_strmid(dim_names[idim], 0, p1, tmp_dname1);
			      else if (idim == 3) sd_strmid(dim_names[idim], 0, p1, tmp_dname2);
			    }
			  else 
			    {
			      if (idim == 2) strcpy(tmp_dname1, dim_names[idim]);
			      else if (idim == 3) strcpy(tmp_dname2, dim_names[idim]);
			    }
			} /* for (idim = 2; .... */
		      
		      if (in_sds_info.rank == 3)
			sprintf(tmp_name, "%s_%s%d", sds_name, tmp_dname1, out_dim_ids[0][isds+1] + 1); 
		      else if (in_sds_info.rank == 4)
			{
			  j1 = isds/out_dim_ids_cnt[1];
			  j2 = isds%out_dim_ids_cnt[1];
			  sprintf(tmp_name, "%s_%s%d_%s%d", sds_name, tmp_dname1, out_dim_ids[0][j1+1]+1, 
				  tmp_dname2, out_dim_ids[1][j2+1] + 1); 
			}
		      strcpy(out_sds1_info[isds].name, tmp_name);
		      if (open_sds((char *)NULL, &out_sds1_info[isds], 'W') != -1)
			write_attr_fval(out_sds1_info[isds].sds_id, out_sds1_info[isds].data_type, 1, 
					out_sds1_info[isds].fill_val, ATTR_FILL_NAME);
		    }
		}
	      else
		{  
		  for (isds=0; isds<nsds; isds++)
		    {
		      out_sds1_info[isds].rank = 2;
		      out_sds1_info[isds].sd_id = out_sd_id;
		      out_sds1_info[isds].fill_val = in_sds_info.fill_val;
		      out_sds1_info[isds].data_type = in_sds_info.data_type;
		      
		      out_sds1_info[isds].dim_size[0] = in_sds_info.dim_size[in_sds_info.rank - 2 ];
		      out_sds1_info[isds].dim_size[1] = in_sds_info.dim_size[in_sds_info.rank - 1];
		      
		      for (idim = 0; idim < in_sds_info.rank -2; idim++)
			{
			  p1 = sd_charpos(dim_names[idim], ':', 0);
			  if (p1 != -1)
			    { 
			      if (idim == 0) sd_strmid(dim_names[idim], 0, p1, tmp_dname1);
			      else if (idim == 1) sd_strmid(dim_names[idim], 0, p1, tmp_dname2);
			    }
			  else 
			    {
			      if (idim == 0) strcpy(tmp_dname1, dim_names[idim]);
			      else if (idim == 1) strcpy(tmp_dname2, dim_names[idim]);
			    }
			} /* for (idim = 0; .... */
		      
		      if (in_sds_info.rank == 3)
			sprintf(tmp_name, "%s_%s%d", sds_name, tmp_dname1, out_dim_ids[0][isds+1] + 1); 
		      else if (in_sds_info.rank == 4)
			{
			  j1 = isds/out_dim_ids_cnt[1];
			  j2 = isds%out_dim_ids_cnt[1];
			  sprintf(tmp_name, "%s_%s%d_%s%d", sds_name, tmp_dname1, out_dim_ids[0][j1+1]+1, 
				  tmp_dname2, out_dim_ids[1][j2+1] + 1); 
			}
		      strcpy(out_sds1_info[isds].name, tmp_name);
		      
		      if (in_sds_info.rank == 2)
			{
			  strcpy(out_sds1_info[isds].name, in_sds_info.name);
			}
		      if (open_sds((char *)NULL, &out_sds1_info[isds], 'W') != -1)
			write_attr_fval(out_sds1_info[isds].sds_id, out_sds1_info[isds].data_type, 1, 
					out_sds1_info[isds].fill_val, ATTR_FILL_NAME);
		    } /* for (isds=0; ....... */	    
		}
	      
	      if (out_all == 1)
		{
		  out_sdsf_info.rank = 2;
		  out_sdsf_info.sd_id = out_sd_id;
		  out_sdsf_info.fill_val = in_sds_info.fill_val;
		  out_sdsf_info.data_type = in_sds_info.data_type;
		  sprintf(out_sdsf_info.name, "%s_%s", sds_name, "all");
		  
		  if (bsq == 0)
		    {
		      if (in_sds_info.rank == 3)
			{
			  out_sdsf_info.dim_size[0] = in_sds_info.dim_size[0];
			  out_sdsf_info.dim_size[1] = in_sds_info.dim_size[1]*out_dim_ids_cnt[0];
			}
		      else if (in_sds_info.rank == 4)
			{
			  out_sdsf_info.dim_size[0] = in_sds_info.dim_size[0]*out_dim_ids_cnt[0];
			  out_sdsf_info.dim_size[1] = in_sds_info.dim_size[1]*out_dim_ids_cnt[1];
			}
		      
		      if (open_sds((char *)NULL, &out_sdsf_info, 'W') != -1)
			write_attr_fval(out_sdsf_info.sds_id, out_sdsf_info.data_type, 1, out_sdsf_info.fill_val, ATTR_FILL_NAME);
		    }
		  else
		    {
		      if (in_sds_info.rank ==2)
			{
			  out_sdsf_info.dim_size[0] = in_sds_info.dim_size[0];
			  out_sdsf_info.dim_size[1] = in_sds_info.dim_size[1];
			}
		      if (in_sds_info.rank == 3)
			{
			  out_sdsf_info.dim_size[0] = in_sds_info.dim_size[in_sds_info.rank -2];
			  out_sdsf_info.dim_size[1] = in_sds_info.dim_size[in_sds_info.rank -1]*out_dim_ids_cnt[0];
			}
		      else if (in_sds_info.rank == 4)
			{
			  out_sdsf_info.dim_size[0] = in_sds_info.dim_size[in_sds_info.rank -2]*out_dim_ids_cnt[0];
			  out_sdsf_info.dim_size[1] = in_sds_info.dim_size[in_sds_info.rank -1]*out_dim_ids_cnt[1];
			}
		      if (open_sds((char *)NULL, &out_sdsf_info, 'W') != -1)
			write_attr_fval(out_sdsf_info.sds_id, out_sdsf_info.data_type, 1, out_sdsf_info.fill_val, ATTR_FILL_NAME);
		    }
		  
		} /* if (out_all == 1) */
	      
	if (bsq == 0)
	  {
	    edge[0] = edge1[0] = edgef[0] = 1;
	    for (i=1; i<in_sds_info.rank; i++)
	      edge[i] = in_sds_info.dim_size[i];
	    edge1[1] = in_sds_info.dim_size[1];
	    edgef[1] = in_sds_info.dim_size[1];
	    
	    if (in_sds_info.rank == 3)
	      tot_num_vals = in_sds_info.dim_size[2];
	    else
	      tot_num_vals = in_sds_info.dim_size[2]*in_sds_info.dim_size[3];
	    
	    for (iline=0; iline<in_sds_info.dim_size[0]; iline++)
	      { 
		start[0] = start1[0] = (int32)iline;
		if (SDreaddata(in_sds_info.sds_id, start, NULL, edge, sds_data_in) != FAIL)
		  {
		    for (is=0; is<in_sds_info.dim_size[1]; is++)
		      {
			offset_is = is*tot_num_vals;
			for (isds=0; isds<nsds; isds++)
			  {
			    if (in_sds_info.rank == 3)
			      offset = offset_is + out_dim_ids[0][1+isds];
			    else if (in_sds_info.rank == 4)
			      {
				j1 = isds/out_dim_ids_cnt[1];
				j2 = isds%out_dim_ids_cnt[1];
				offset = offset_is + j1*in_sds_info.dim_size[3] + out_dim_ids[1][j2+1];
			      }
			    switch(in_sds_info.data_type)
			      {
			      case 5:
				((float32 *)sds_data_out[isds])[is] =  ((float32 *)sds_data_in)[offset];
				break;
			      case 6:
				((float64 *)sds_data_out[isds])[is] =  ((float64 *)sds_data_in)[offset];
				break;			      
			      case 20: ((int8 *)sds_data_out[isds])[is] = ((int8 *)sds_data_in)[offset];
				break;
			      case 21: ((uint8 *)sds_data_out[isds])[is] = ((uint8 *)sds_data_in)[offset];
				break;
			      case 22: ((int16 *)sds_data_out[isds])[is] = ((int16 *)sds_data_in)[offset];
				break;
			      case 23: ((uint16 *)sds_data_out[isds])[is] = ((uint16 *)sds_data_in)[offset];
				break;
			      case 24: ((int32 *)sds_data_out[isds])[is] = ((int32 *)sds_data_in)[offset];
				break;
			      case 25: ((uint32 *)sds_data_out[isds])[is] = ((uint32 *)sds_data_in)[offset];
				break;
			      default: fprintf(stdout,
                                         "HDF datatype " LONG_INT_FMT " not supported", 
                                         in_sds_info.data_type);
	
			      }
			  } /* for (isds=0; . . . */
		      } /* for (is=0; . . . */
		    for (isds=0; isds<nsds; isds++)
		      {
			if (SDwritedata(out_sds1_info[isds].sds_id, start1, NULL, edge1, sds_data_out[isds])
			    == FAIL)
			  fprintf(stderr, "Error writing data line for %s\n", out_sds1_info[isds].name);
			if (out_all == 1)
			  {
			    if (in_sds_info.rank == 3)
			      {
				startf[0] = iline;
				startf[1] = in_sds_info.dim_size[1]*isds;
			      }
			    else
			      {
				j1 = isds/out_dim_ids_cnt[1];
				j2 = isds%out_dim_ids_cnt[1];
				startf[0] = iline + (j1*in_sds_info.dim_size[0]);
				startf[1] = in_sds_info.dim_size[1]*j2;
			      }
			    if (SDwritedata(out_sdsf_info.sds_id, startf, NULL, edgef, sds_data_out[isds]) 
				== FAIL)
			      fprintf(stderr, "Error writing data line for %s\n", out_sdsf_info.name);
			  }
		      } /* for isds =0; . . . */
		  } /* if (SDreaddata . . . */
	      } /* for (iline=0; . . . . . */
	  } /* if (bsq == 0) */
	else   /* (bsq == 1) */ 
	  { 
	    edge1[0] = edgef[0] = 1;
	    
	    for (i=0; i<in_sds_info.rank; i++) 
	      edge[i] = in_sds_info.dim_size[i]; 
	    
	    edge[in_sds_info.rank -2] = 1;
	    
	    edge1[1] = in_sds_info.dim_size[in_sds_info.rank - 1];
	    edgef[1] = in_sds_info.dim_size[in_sds_info.rank - 1];	  
	    
	    for (iline=0; iline<in_sds_info.dim_size[in_sds_info.rank - 2]; iline++)
	      { 
		start1[0] = (int32)iline; 
		start[in_sds_info.rank -2] = (int32)iline;
		
		if (SDreaddata(in_sds_info.sds_id, start, NULL, edge, sds_data_in) != FAIL)
		  { 	
		  for (is=0; is<in_sds_info.dim_size[in_sds_info.rank - 1]; is++)
		    { 
		      offset_is = is;
		      for (isds=0; isds<nsds; isds++)
			{
			  if (in_sds_info.rank == 2)
			    offset = offset_is;
			  if (in_sds_info.rank == 3)
			    offset = offset_is + 
			             out_dim_ids[0][1+isds]*in_sds_info.dim_size[in_sds_info.rank-1];
			  if (in_sds_info.rank == 4)
			    {
			      j1 = isds/out_dim_ids_cnt[1];
			      j2 = isds%out_dim_ids_cnt[1];
			      offset = offset_is + j1*in_sds_info.dim_size[1] + 
				out_dim_ids[1][j2+1]*in_sds_info.dim_size[in_sds_info.rank - 1];
			    }
		
			  switch(in_sds_info.data_type)
			    {
			    case 5:
			      ((float32 *)sds_data_out[isds])[is] =  ((float32 *)sds_data_in)[offset];
			      break;
			    case 6:
			      ((float64 *)sds_data_out[isds])[is] =  ((float64 *)sds_data_in)[offset];
			      break;
			    case 20: 
			      ((int8 *)sds_data_out[isds])[is] = ((int8 *)sds_data_in)[offset];
			      break;
			    case 21: 
			      ((uint8 *)sds_data_out[isds])[is] = ((uint8 *)sds_data_in)[offset];
			      break;
			    case 22: 
			      ((int16 *)sds_data_out[isds])[is] = ((int16 *)sds_data_in)[offset];
			      break;
			    case 23: 
			      ((uint16 *)sds_data_out[isds])[is] = ((uint16 *)sds_data_in)[offset];
			      break;
			    case 24: 
			      ((int32 *)sds_data_out[isds])[is] = ((int32 *)sds_data_in)[offset];
			      break;
			    case 25: 
			      ((uint32 *)sds_data_out[isds])[is] = ((uint32 *)sds_data_in)[offset];
			      break;
			    default: fprintf(stdout,
                                         "HDF datatype " LONG_INT_FMT " not supported", 
                                         in_sds_info.data_type);

			    }
			} /* for (isds=0; . . . */
		    } /* for (is=0; . . . */
		  
		  for (isds=0; isds<nsds; isds++)
		    {	
		      if (SDwritedata(out_sds1_info[isds].sds_id, start1, NULL, edge1, sds_data_out[isds])
			  == FAIL)  
			  fprintf(stderr, "Error writing data line for %s\n", out_sds1_info[isds].name);
		      if (out_all == 1)
			{
			  if (in_sds_info.rank == 2)
			    {
			      startf[0] = (int32)iline;
			      startf[1] = 0;
			    }
			  else 
			    {
			      if (in_sds_info.rank == 3)
				{
				  startf[0] = (int32)iline;
				  startf[1] = in_sds_info.dim_size[in_sds_info.rank - 1]*isds; 
				}
			      else
				{
				  j1 = isds/out_dim_ids_cnt[1];
				  j2 = isds%out_dim_ids_cnt[1];
				  startf[0] = iline + (j1*in_sds_info.dim_size[in_sds_info.rank - 2]); 
				  startf[1] = in_sds_info.dim_size[in_sds_info.rank - 1]*j2;
				}
			    }
			  
			  if (SDwritedata(out_sdsf_info.sds_id, startf, NULL, edgef, sds_data_out[isds]) 
			      == FAIL)
			    {
			      fprintf(stderr, "Error writing data line for %s\n", out_sdsf_info.name); }
			}
		    } /* for (isds =0; . . . */
		} /* if (SDreaddata . . . */
	    } /* for (iline=0; . . . . . */
	} /* if (bsq == 1) */
      for (isds=0; isds<nsds; isds++)
	SDendaccess(out_sds1_info[isds].sds_id);
      if (out_all == 1) SDendaccess(out_sdsf_info.sds_id);
    }  /* if (dim_names != . . . .   */ 
    Free2D((void **)dim_names);
    Free2D((void **)sds_data_out);
    Free2D((void **)out_dim_ids);
    if (sds_data_in != NULL) free(sds_data_in);
    if (out_sds1_info != NULL) free(out_sds1_info);
  } /* if (get_sds_info . . .     */
  if (in_sds_info.sds_id != -1) SDendaccess(in_sds_info.sds_id);
    }
}

int get_out_dim_id(char **dim_names, int32 *dim_size, int rank, char **dim_ids, 
		   int dim_cnt, int **out_dim_ids, int *out_dim_ids_cnt)
/*
!C******************************************************************************
    
!Function: get_out_dim_id
        
!Description:
  Get the output dimension id.
    
!Input Parameters: 
  dim_names       Dimension names
  dim_size        Dimension size of each dimension.
  rank            The rank of the input SDS.
  dim_ids         Selected dimension name.
  dim_cnt         Count of selected dimension number.
  out_dim_ids     Output dimension name.
  out_dim_ids_cnt Count of each output dimension.
                                  
!Output Parameters: 
  1               On succeed.
  -1              On failure.
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
********************************************************************************/
{
  int ir, isz, bsq;
  int cnt, len;
  int idim, icnt;
  int p1, p2, p3;
  int **sel_dim_ids;
  int one_based_num, num1, num2;
  char num_str[10];
  char num1_str[10];
  char num2_str[10];
  char sel_dim_name[50];
  int sel_dim_cnt[MAX_NUM_DIM];
  int sel_dim_name_is_valid ;

  bsq = ((rank == 2) || (dim_size[0] < dim_size[rank - 1])) ? 1 : 0;   

  if (dim_cnt == 0)
    {
      if (bsq == 0)
	{
	  for (ir=2; ir<rank; ir++)
	    {
	      out_dim_ids[ir-2][0] = ir;
	      for (isz=0; isz<dim_size[ir]; isz++)
		{
		  out_dim_ids[ir-2][1+isz] = isz; 
		}
	      out_dim_ids_cnt[ir-2] = dim_size[ir];
	    }
	}
      else
	{
	  for (ir=0; ir<rank-2; ir++)
	    {
	      out_dim_ids[ir][0] = ir;
	      for (isz=0; isz<dim_size[ir]; isz++)
		{
		  out_dim_ids[ir][1+isz] = isz; 
		}
	      out_dim_ids_cnt[ir] = dim_size[ir];
	    }
	}
    }	  
  else
  {
    if ((sel_dim_ids = (int **)Calloc2D(MAX_NUM_DIM, MAX_DIM_SIZE, sizeof(int))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for sel_dim_ids in get_out_dim_id\n"); 
      return -1;
    }
    for (idim=0; idim < dim_cnt; idim++)
    {
      p1 = cnt = 0;
      sel_dim_cnt[idim] = -1;
      sel_dim_ids[idim][0] = -1;
      len = (int)strlen(dim_ids[idim]);
      p2 = sd_charpos(dim_ids[idim], ',', p1);
      if (p2 == -1)
        sd_strmid(dim_ids[idim], 0, len, sel_dim_name);
      else
        sd_strmid(dim_ids[idim], 0, p2, sel_dim_name);
      
      sel_dim_name_is_valid = 0;
      for (ir=0; ir<rank; ir++)
	{ 
	  if (strstr(dim_names[ir],sel_dim_name) != NULL) 
	    {	      
	      sel_dim_name_is_valid = 1;
	      break;
	    }
	} /* for (ir=0;......) */
      
      if (sel_dim_name_is_valid == 0)
	{
	  fprintf(stderr, "input dim name %s is invalid. Program exit.\n", sel_dim_name);
	  exit(EXIT_FAILURE);
	}
      else
	{
	    sel_dim_ids[idim][cnt] = ir;
	}

      if (sel_dim_ids[idim][cnt] != -1)
	{
	  if (p2 == -1)
	    {
	      for (isz=0; isz<dim_size[ir]; isz++)
		sel_dim_ids[idim][++cnt] = isz; 
	      sel_dim_cnt[idim] = dim_size[idim];
	    }
	  else
	    {
	      while (p2 != -1)
		{
		  p1 = p2+1;
		  p2 = sd_charpos(dim_ids[idim], ',', p1);
		  p3 = sd_charpos(dim_ids[idim], '-', p1);
		  if ((p3 == -1) || ((p2 != -1) && (p2 < p3)))   /* user use ',' between dimension 
								    selection */
		    {
		      if (p2 == -1)
			sd_strmid(dim_ids[idim], p1, len-p1, num_str);
		      else
			sd_strmid(dim_ids[idim], p1, p2-p1, num_str);

		      one_based_num = atoi(num_str);

		      if ((one_based_num -1 >= dim_size[ir]) || ( one_based_num -1 < 0 ))
			fprintf(stderr, "Ignoring the invalid dimension size %d\n", one_based_num);
		      else sel_dim_ids[idim][++cnt] = one_based_num - 1;
		    }
		  else if ((p2 == -1) || (p2 > p3))           /* user use '-' between dimension 
								     selection */
		    {
		      sd_strmid(dim_ids[idim], p1, p3-p1, num1_str);
		      if (p2 == -1)
			sd_strmid(dim_ids[idim], p3+1, len-p3-1, num2_str);
		      else
			sd_strmid(dim_ids[idim], p3+1, p2-p3-1, num2_str);
		      num1 = atoi(num1_str);
		      num2 = atoi(num2_str);
		      
		      for (one_based_num=num1; one_based_num<=num2; one_based_num++)
			{
			  if ((one_based_num - 1 >= dim_size[ir]) || (one_based_num - 1 < 0))
				fprintf(stderr, "Ignoring the invalid dimension size %d\n", one_based_num);
			  else sel_dim_ids[idim][++cnt] = one_based_num - 1;
			}
		    }
		}
	      sel_dim_cnt[idim] = cnt;
	    }
	} /* if sel_dim_nums. . . . */
    }   /* for (idim=0; . . . */

    if (bsq == 0)
      {
	for (ir=2, cnt=0; ir<rank; ir++, cnt++)
	  {
	    out_dim_ids[cnt][0] = ir;
	    for (idim=0; idim<dim_cnt; idim++)
	      if (sel_dim_ids[idim][0] == ir)
		{
		  for (icnt=1; icnt<=sel_dim_cnt[idim]; icnt++)
		    out_dim_ids[cnt][icnt] = sel_dim_ids[idim][icnt];
		  out_dim_ids_cnt[cnt] = sel_dim_cnt[idim];
		  break;
		}
	    if (idim == dim_cnt)
	      {
		for (icnt=0; icnt<dim_size[ir]; icnt++)
		  out_dim_ids[cnt][1+icnt] = icnt;
		out_dim_ids_cnt[cnt] = dim_size[ir];
	      }
	  }	/* for (ir=2 . . .  */
      }
    else 
      {
	for (ir=0, cnt=0; ir<rank-2; ir++, cnt++)
	  {
	    out_dim_ids[cnt][0] = ir;
	    for (idim=0; idim<dim_cnt; idim++)
	      if (sel_dim_ids[idim][0] == ir)
		{
		  for (icnt=1; icnt<=sel_dim_cnt[idim]; icnt++)
		    out_dim_ids[cnt][icnt] = sel_dim_ids[idim][icnt];
		  out_dim_ids_cnt[cnt] = sel_dim_cnt[idim];
		  break;
		}
	    if (idim == dim_cnt)
	      {
		for (icnt=0; icnt<dim_size[ir]; icnt++)
		  out_dim_ids[cnt][1+icnt] = icnt;
		out_dim_ids_cnt[cnt] = dim_size[ir];
	      }
	  }	/* for (ir=0 . . .  */	
      }

    Free2D((void **)sel_dim_ids);
  }  /* else . . . */
  return 1;
}
