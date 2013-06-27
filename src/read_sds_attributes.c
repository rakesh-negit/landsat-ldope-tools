/****************************************************************************
!C

!File: read_sds_attributes.c

!Description:
  Read attributes of one or more SDS of MODIS Land HDF-EOS data product
  
!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:

  Version 1.0    April 5, 2003

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

    Yi Zhang 
    LDOPE                             Science Systems and Applications Inc.
    sophyi@ltpmail.gsfc.nasa.gov      NASA/GSFC Code 922
    phone: 301-614-5497               Greenbelt, MD 20771

    David Roy 
    LDOPE                             University of Maryland 
    droy@kratmos.gsfc.nasa.gov	      NASA/GSFC Code 922
    phone: 301-614-5571               Greenbelt, MD 20771

!Design Notes: (none)

!END
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mfhdf.h"
#include "sds_rw.h"
#include "qa_tool.h"
#include "main_util.h"
#include "alloc_mem.h"
#include "str_op.h"
#include "l2g.h"

#define MAX_NUM_VALUES 16000

#define HELP \
"NAME \n" \
"    read_sds_attributes - read attributes of one or more SDSs of input\n"\
"    MODIS Land HDF-EOS data product \n" \
" \n" \
" \n" \
"SYNOPSIS \n" \
"    read_sds_attributes -help [filename]\n" \
"    read_sds_attributes [-sds=<sds_name>] filename \n" \
" \n" \
"DESCRIPTION \n" \
"    Read attributes of one or more SDSs of input MODIS Land HDF-EOS data \n" \
"    product and output the result to stdout. The SDS attributes include: \n" \
"    fill values, units, scaling and offsets values, SDS long name, etc. \n" \
" \n" \
"OPTIONS \n" \
"    -help                    Display this help message, If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDSs in the file are\n"\
"                             displayed. \n" \
"    -sds=<SDS_list>          List of SDS to process. SDS names are\n"\
"                             separated by commas with no space. By default\n"\
"                             attributes for all SDSs in the input file are\n"\
"                             output. \n" \
"    filename                 Input filename. \n" \
" \n" \
"Examples: \n" \
"    read_sds_attributes \n" \
"        -sds='500m Surface Reflectance Band 1,\n"\
"              500m Surface Reflectance Band 3,\n"\
"              500m Surface Reflectance Band 4'\n"\
"              MOD09.A2002123.0040.003.2002125174437.hdf  \n" \
"    read_sds_attributes MOD35_L2.A2002123.0040.003.2002124023706.hdf \n" \
"\n"\
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
" \n" \
"Version 1.0, 04/05/2004 \n"

#define USAGE \
"usage:	\n" \
"    read_sds_attributes -help [filename]\n" \
"    read_sds_attributes [-sds=<sds_name>] filename \n" \
"OPTIONS \n" \
"    -help                    Display this help message, If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDSs in the file are\n"\
"                             displayed. \n" \
"    -sds=<SDS_list>          List of SDS to process. SDS names are\n"\
"                             separated by commas with no space. By default\n"\
"                             attributes for all SDSs in the input file are\n"\
"                             output. \n" \
"    filename                 Input filename. \n" \
" \n"

void read_sds_attr(char *fname, char **sds_names, int32 nsds); 
void print_sds_values(sds_t *sds_info, int **sds_values, float **sds_fvalues, int *num_val);
int parse_cmd_read_sds_attr(int argc, char **argv, char *in_fname, char **sds_names, int* sds_cnt);

int main(int argc, char **argv)
/*
!C*********************************************************************************

!read_sds_attr
       
!Description:
   read_sds_attributes display the attribute name, type and value of SDS(s) to 
   the screen.

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/
{
   char **sds_names;
 
   int i, k, status = 0;
   int isds, sds_cnt;
   int32 sd_id, sds_id, msds, nattr, rank, dim_size[4], dt;
   char in_fname[MAX_PATH_LENGTH];
   char name[MAX_SDS_NAME_LEN];
   char sds_name[MAX_SDS_NAME_LEN];
   char dim_str[MAX_STR_LEN];
	
   if (argc == 1)
   {
      status = -1;
      fprintf(stderr, "Missing input file \n");
      fprintf(stderr, "%s\n", USAGE);
      exit(EXIT_FAILURE);
   }   
  
   /* display help */
   if ((argc == 2) && (strcmp(argv[1], "-help") == 0))
   {
     	fprintf(stderr, "%s\n", HELP);
     	exit(0);
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
		       for (isds=0; isds<msds; isds++) 
		       {
			   if ((sds_id = SDselect(sd_id, isds)) != FAIL) 
			   {
			       if (SDgetinfo(sds_id, name, &rank, dim_size, &dt, &nattr) != FAIL) 
			       {
			    	  sprintf(sds_name, "%s (" LONG_INT_FMT, name, dim_size[0]);
			    	  for (k=1; k<rank; k++) 
				  {
			      	      sprintf(dim_str, " x " LONG_INT_FMT, dim_size[k]);
			      	      strcat(sds_name, dim_str);
			          }
			    	  switch(dt) 
			       	  {
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
  
   if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
   {
    	fprintf(stderr, "Cannot allocate memory for sds_names in sds_range: main()\n");
    	exit(EXIT_FAILURE);
   }
   
   if (sds_names != NULL)
   {
        sds_cnt = 0;
    	status = parse_cmd_read_sds_attr(argc, argv, in_fname, sds_names, &sds_cnt);

	if (sds_cnt != 0)
	{
            msds = sds_cnt;
	}
        else
	{
	    msds = get_sds_names(in_fname, sds_names);
	}
	

        if (status == -1)
        {
    	    fprintf(stderr, "%s", USAGE);
	    exit(EXIT_FAILURE);
        }
        else if (status != 0)
        {
	    read_sds_attr(in_fname, sds_names, msds);   
        }
   }
   Free2D((void **)sds_names);
   fprintf(stderr, "Processing done ! \n");
   return status;
}
 
   
void read_sds_attr(char *fname, char **sds_names, int32 nsds)

/*
!C*********************************************************************************

!Function: read_sds_attr
       
!Description:
   function read_sds_attr display the attribute name, type and value to 
   the screen.

!Input Parameters:

  fname     Input HDF file name.
  sds_names User input sds_name. 
  sds_cnt   Number of input sds_name.

!Revision History: (see file prolog) 

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
*****************************************************************************/

{
  int isds, j, m;
  int32 k, attr_size;
  int32 attr_type, attr_cnt, rank, num_attrs, dimsizes[4], data_type;
  sds_t *sds_info;
  void *attr_buf;
  char sds_name[MAX_SDS_NAME_LEN];
  char attr_name[MAX_SDS_NAME_LEN];
  char* tempstr;

  /* allocate memory for tempstr (for output formating purpose) */
  if ((tempstr = (char *)calloc(50, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for tempstr in read_sds_attr()\n");

  /* allocate memory for sds_info */
  if ((sds_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
    	fprintf(stderr, "Cannot allocate memory for sds_info in read_sds_attr()\n");
  if (sds_info != NULL)
  {
    for (isds=0; isds<nsds; isds++)
    {
      sds_info[isds].sds_id = -1;
      strcpy(sds_info[isds].name, sds_names[isds]);
      sds_info[isds].sd_id = (isds == 0) ? -1:sds_info[0].sd_id;
  		
      if (get_sds_info(fname, &sds_info[isds]) != -1)
      {
	if (SDgetinfo(sds_info[isds].sds_id, sds_name, &rank, dimsizes, &data_type, &num_attrs ) != FAIL) 
	{ 
	  fprintf(stdout, "======================================================================\n");
	  fprintf(stdout, "SDS : %s \n", sds_name);
	  fprintf(stdout, "  %-17s%-10s%15s\n", "Attribute", "Data Type", "Value");
	   fprintf(stdout, "-----------------+-----------+---------------------------------------\n");
	  for(k=0; k < num_attrs; k++) /* k is the index of the attribute */
	  {
	    if (SDattrinfo(sds_info[isds].sds_id, k, attr_name, &attr_type, &attr_cnt) != FAIL)
	    { 
	      fprintf(stdout, "%-20s", attr_name);

	      attr_size = DFKNTsize(attr_type);
	      if ((attr_buf = (void *)calloc(attr_cnt, attr_size)) == NULL)
		{
		  fprintf(stderr, "Cannot allocate memory for attr_buf in read_sds_attr()");
		  exit(EXIT_FAILURE);		 
		}

	      if (SDreadattr(sds_info[isds].sds_id, k, attr_buf) == FAIL)
		{
		  fprintf(stderr, "Cannot read SDS attributes in read_sds_attr()");
		  exit(EXIT_FAILURE);
		}
		else 	    
		{  
		  switch(attr_type)
		    {
		    case  4:
		      if (attr_cnt == 1 )
			{
			  fprintf(stdout, "%-10s", "CHAR");
			  fprintf(stdout, "%-10s \n", (char *)attr_buf);
			}
		      else
			{
			  fprintf(stdout, "%-10s", "STRING");
			  if (attr_cnt < 255)
			    {
			      sd_strmid((char *)attr_buf, 0, 30, tempstr);
			      fprintf(stdout, "%-30s \n", tempstr);
			      for (m=1; m <= (attr_cnt/30); m++)
				{
				  sd_strmid((char *)attr_buf, 30*m, 30, tempstr);
				  fprintf(stdout, "%-30s%-30s\n", " ", tempstr);
				}
			    }
			  else 
			    {
			       fprintf(stdout, "\n");
			       fprintf(stdout, "Value =");
			       fprintf(stdout, "%s\n", (char *)attr_buf);
			    }
			}

		      break;		      
		    case  5: 
		      fprintf(stdout, "%-10s", "FLOAT32");
		      fprintf(stdout, "%f", ((float32 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", %f", ((float32 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;				        
		    case  6:   	
		      fprintf(stdout, "%-10s", "FLOAT64");
		      fprintf(stdout, "%f", ((float64 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", %f", ((float64 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    case 20: 
		      fprintf(stdout, "%-10s", "INT8");
		      fprintf(stdout, "%d", ((int8 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", %d", ((int8 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    case 21:
		      fprintf(stdout, "%-10s", "UINT8");
		      fprintf(stdout, "%d", ((uint8 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", %d", ((uint8 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    case 22: 
		      fprintf(stdout, "%-10s", "INT16");
		      fprintf(stdout, "%d", ((int16 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", %d", ((int16 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    case 23: 
		      fprintf(stdout, "%-10s", "UINT16");
		      fprintf(stdout, "%d", ((uint16 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", %d", ((uint16 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    case 24:
		      fprintf(stdout, "%-10s", "INT32");
		      fprintf(stdout, LONG_INT_FMT, ((int32 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", " LONG_INT_FMT, ((int32 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    case 25:
		      fprintf(stdout, "%-10s", "UINT32");
		      fprintf(stdout, LONG_INT_FMT, ((uint32 *)attr_buf)[0]);
		      for (j=1; j<attr_cnt; j++)
			{
			  fprintf(stdout, ", " LONG_INT_FMT, ((uint32 *)attr_buf)[j]);
			}
		      fprintf(stdout,"\n");
		      break;
		    default: 
		      fprintf(stdout, "%-10s\n", "Unknown");
		      break;
		    } /* switch(attr_type) */
		}   /* if (SDreadattr...) */
	    } /* SDattrinfo...) */
	  } /* for(k =0...) */
	} /* if (SDgetinfo */
      } /* if (get_sds_info...)  */
    } /* for (isds=0; . . ) */						
  } /* if (sds_info != NULL) */	        
}

int parse_cmd_read_sds_attr(int argc, char **argv, char *in_fname, char **sds_names, int* sds_cnt)
/*
!C*********************************************************************************

!Function: parse_cmd_read_sds_attr
       
!Description:
   function parse_cmd_read_sds_attr parse the argument string to the corresponding
   input filename (in_fname), output filename (out_fname) and SDS name
   (sds_name).
   It print out a error message if there the input file, output file or 
   SDS name is missing. 

!Input Parameters: 
  argc      input argument count
  argv      input argument vector

!Input/output Parameters:

  in_fname  Input HDF file name.
  sds_names User input sds_name. 
  sds_cnt   Number of input sds_name.

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
  
  in_fname[0] = '\0';
  
  for (i=1, ret=1; i<argc; i++)
    {
      if (is_arg_id(argv[i], "-sds") == 0) 
	{
	  *sds_cnt = 0;
	  get_arg_val_arr(argv[i], sds_names, sds_cnt);
	}
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
  return ret;
}
