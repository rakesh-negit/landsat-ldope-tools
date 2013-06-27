/****************************************************************************
!C

!File: read_pixvals.c

!Description:
  Read observations at a specific pixel in a list of files.  

!Input Parameters: (none)

!Input/Output Parameters: (none)

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
    LDOPE                             University of Maryland, 
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
#include "meta.h"

#define HELP \
"NAME \n" \
"    read_pixvals - read MODIS Land HDF-EOS data product values at specified\n"\
"    pixel locations.\n" \
"\n" \
"SYNOPSIS \n" \
"    read_pixvals [-help] [filename] \n" \
"\n" \
"    read_pixvals -xy=col[.cs[.cs]],row[.rs[.rs]]|<coordinates filename> \n" \
"                 [-res=qkm|hkm|1km] filename(s) \n" \
"\n" \
"DESCRIPTION \n" \
"    Read the pixel values at specified locations in one or more input MODIS\n"\
"    Land data product and output to stdout. The pixel values for each SDS\n"\
"    in each file are output as separate lines. \n" \
"\n" \
"    If more than one input MODIS Land HDF-EOS data product is specified\n"\
"    then they must be all L2 products or all L2G/L3/L4 MODIS Land data\n"\
"    products. \n" \
"\n" \
"    The MODIS Land data products may contain SDSs with different spatial\n"\
"    dimensions corresponding to the 250m, 500m and 1km MODIS Land pixel\n"\
"    resolutions. In this case, the -res option is used to specify which of\n"\
"    the 1km, 500m or 250m pixel resolutions is referenced. If -res is not\n"\
"    specified then the xy location is assumed to reference the coarsest\n"\
"    spatial resolution of the different SDSs. \n" \
"\n" \
"    Sub pixel locations may be output by specifying a sub pixel offset (0\n"\
"    or 1 in the x and/or y axes). If not specified a 0 pixel offset is\n"\
"    assumed. See examples below. \n" \
"\n" \
"    This tool supports 2D/3D/4D and L2G compact and full format SDSs. Pixel\n"\
"    values for each layer of the 3D/4D and L2G SDSs are output. \n" \
"\n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDSs in the file are\n" \
"                             displayed. \n" \
"    -res=qkm|hkm|1km         Reference SDS resolution (qkm=250m, hkm=500m,\n"\
"                             1km=1000m) of the pixel location specified in\n"\
"                             the <96>xy argument. If unspecified the\n" \
"                             coarsest resolution of all the SDSs in the\n"\
"                             input file list is used. \n" \
"    -xy=col[.cs[.cs]],row[.rs[.rs]]|<coordinates filename> \n" \
"                             Column and row pixel locations (0-based) or\n"\
"                             name of an ASCII coordinates file containing\n"\
"                             the column and row pixel locations.  Multiple\n"\
"                             locations may be specified by repeating the xy\n"\
"                             option or by specifying the  x and y\n"\
"                             coordinates on different lines in the ASCII\n"\
"                             coordinates file.  \n" \
"                             Sub pixel offsets for higher spatial\n"\
"                             resolution SDS may be specified as col.cs\n"\
"                             row.rs. The offsets refer to the top left\n"\
"                             corner pixel. For example: \n" \
"                                -res=hkm -xy=100.0,200.0 \n" \
"                                  (read values at pixel 100,200 from the\n"\
"                                   500m SDS, and at pixel 200,400 from the\n"\
"                                   250m SDS) \n" \
"                                -res=hkm -xy=100.1,200.1 \n" \
"                                  (read values at pixel 100,200 from the\n"\
"                                   500m SDS, and at pixel 201,401 from the\n"\
"                                   250m SDS) \n" \
"                                -res=1km -xy=100.0.1,200.0.1 \n" \
"                                  (read values at pixel 100,200 from the\n"\
"                                   1km SDS, at pixel 200,400 from the 500m\n"\
"                                   SDS, and at pixel 401,801 from the 250m\n"\
"                                   SDS)\n" \
"    filename(s)              Input filename(s) \n" \
"\n" \
"Examples \n" \
"    read_pixvals -xy=1000,200 -xy=0,0\n"\
"                 MOD09A1.A2003057.h29v11.004.2003069051044.hdf \n" \
"    read_pixvals -xy=10,10 -res=1km\n"\
"                 MOD09GHK.A2002225.h19v09.003.2002227235523.hdf \n" \
"                 MOD09GQK.A2002225.h19v09.003.2002227235424.hdf \n" \
"\n" \
"AUTHOR \n" \
"    Code: S. Devadiga and Yi Zhang \n" \
"    Documentation: S. Devadiga and D. Roy \n" \
"\n" \
"Version 1.0, 05/12/2003 \n" \


#define USAGE \
"usage:	\n" \
"   read_pixvals -xy=col[.cs[.cs]],row[.rs[.rs]]|<coordinates filename> \n" \
"                [-res=qkm|hkm|1km] filename(s) \n" \
"\n" \
"OPTIONS \n" \
"    -help [filename]         Display this help message. If the input\n"\
"                             filename is specified with this option, then\n"\
"                             the names of all the SDSs in the file are\n" \
"                             displayed. \n" \
"    -res=qkm|hkm|1km         Reference SDS resolution (qkm=250m, hkm=500m,\n"\
"                             1km=1000m) of the pixel location specified in\n"\
"                             the <96>xy argument. If unspecified the\n" \
"                             coarsest resolution of all the SDSs in the\n"\
"                             input file list is used. \n" \
"    -xy=col[.cs[.cs]],row[.rs[.rs]]|<coordinates filename> \n" \
"                             Column and row pixel locations (0-based) or\n"\
"                             name of an ASCII coordinates file containing\n"\
"                             the column and row pixel locations.  Multiple\n"\
"                             locations may be specified by repeating the xy\n"\
"                             option or by specifying the  x and y\n"\
"                             coordinates on different lines in the ASCII\n"\
"                             coordinates file.  \n" \
"                             Sub pixel offsets for higher spatial\n"\
"                             resolution SDS may be specified as col.cs\n"\
"                             row.rs. The offsets refer to the top left\n"\
"                             corner pixel. For example: \n" \
"                                -res=hkm -xy=100.0,200.0 \n" \
"                                  (read values at pixel 100,200 from the\n"\
"                                   500m SDS, and at pixel 200,400 from the\n"\
"                                   250m SDS) \n" \
"                                -res=hkm -xy=100.1,200.1 \n" \
"                                  (read values at pixel 100,200 from the\n"\
"                                   500m SDS, and at pixel 201,401 from the\n"\
"                                   250m SDS) \n" \
"                                -res=1km -xy=100.0.1,200.0.1 \n" \
"                                  (read values at pixel 100,200 from the\n"\
"                                   1km SDS, at pixel 200,400 from the 500m\n"\
"                                   SDS, and at pixel 401,801 from the 250m\n"\
"                                   SDS)\n" \
"    filename(s)              Input filename(s) \n" \
"\n"



/****************************************************************************************
                              Prototypes
****************************************************************************************/

int parse_cmd_read_pixvals(char **argv, int argc, char **xy_str, int *pt_cnt, int *res,
		char *xy_fname);
void get_xy_pts(char *xy_str, int *xy_pt, int *xy_sh, int *xy_sq);
void read_l2l3_obs_at_pts(char *fname, char **xy_str, char *xy_fname, int npt, 
			  int res, int gran_st);
void read_l2g_obs_at_pts(char *fname, char **xy_str, char *xy_fname, int npt, int res);
void print_sds_val(int **out_pnts, int32 **sds_val, sds_t *sds_info, int nsds, int l2g_st);

/*************************************************************************************/

int main(int argc, char **argv)

/*
!C************************************************************************************
       
!Description:
 The read_pixvals routine takes 3 arguments: 
     -xy       is the Column and row number(0-based) in the user specified resolution 
               or name of the file containing the column and row numbers.  
     -res      is the reference resolution for xy.
     filename  one or more input filenames(can contain either only L2 granules or only
               L2G/L3/L4 tiles.

 Then it reads the observations at the given pixels locations and output to stdout.

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
  int l2g_st = 0, gran_st = 0;
  int status, st, iarg, k;
  int res, npt, isds;
  int32 sd_id, sds_id, nsds;
  int32 rank, dim_size[4], dt, nattr;
  char *ameta_str, *cmeta_str, **xy_str;
  char name[MAX_SDS_NAME_LEN], xy_fname[MAX_PATH_LENGTH];
  char dim_str[MAX_STR_LEN], sds_name[MAX_SDS_NAME_LEN];

  status = 0;
  
  /* Display help */
  if ((argc==2) && (strcmp(argv[1], "-help")==0))
  {
    fprintf(stderr, "%s\n", HELP);
    exit(EXIT_SUCCESS);
  }

  /* Check to see how many arguments are inputed. */
  if (argc >= MAX_NUM_PARAM)
  {
    argc = MAX_NUM_PARAM;
    fprintf(stderr, "Too many input arguments. Using first %d arguments only . \n", argc);
  }

  /* print SDS names if filename input with help option */
  for (iarg=0, status=0; iarg<argc; iarg++)
    if (strcmp(argv[iarg], "-help") == 0) status = iarg;
  if (status != 0)
  {
    /* Print help or all the SDS names from the first file in the list of input files */
    nsds = 0;
    for (iarg=1; iarg<argc; iarg++)
    {
      if (argv[iarg][0] != '-')
      {
        if ((sd_id = SDstart(argv[iarg], DFACC_READ)) == FAIL)
          fprintf(stderr, "Cannot open the HDF file %s\n", argv[iarg]);
        else
        {
          if (SDfileinfo(sd_id, &nsds, &nattr) == FAIL)
            fprintf(stderr, "Error reading information for HDF file %s\n", argv[iarg]);
          else
          {
            fprintf(stdout, "Valid SDS names, dimension and data type in file: %s\n", argv[iarg]);
            for (isds=0; isds<nsds; isds++) {
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
        if (nsds == 0)
          fprintf(stdout, "Input files do not contain any valid SDS\n");
        break;
      }
    } /* for (i=2; . .  ) */
    exit(0);
  }
   
  if ((xy_str = (char **)Calloc2D(MAX_NUM_PNTS, 25, sizeof(char))) == NULL)
    {
      status = 1;
      fprintf(stderr, "Cannot allocate memory for xy_str in read_pixvals(main)\n");
      exit(EXIT_FAILURE);
    }

  if (xy_str != NULL)
    {
      if ((st = parse_cmd_read_pixvals(argv, argc, xy_str, &npt, &res, xy_fname)) == -1)
	{
	  status = 1;
	  fprintf(stderr, "%s\n", USAGE);
	  exit(EXIT_FAILURE);
	}	  
      else if (st != 0)
	{
	  iarg = 1;
	  while (iarg < argc)
	    {
	      if (argv[iarg][0] != '-')
		{
		  if ((ameta_str = get_attr_metadata(argv[iarg], "ArchiveMetadata.0")) 
		      != NULL)
		    { /* test to see if it's a L2G file or not */
		      l2g_st = (strstr(ameta_str, "NUMBEROFOVERLAPGRANULES") != NULL) ? 1 : 0;
		    }
		  
		  if ((cmeta_str = get_attr_metadata(argv[iarg], "CoreMetadata.0")) 
		      != NULL)
		    {
		      gran_st = (strstr(cmeta_str, "HORIZONTALTILENUMBER") == NULL) ? 1 : 0;
		    }  
		  
		  if (l2g_st == 1)
		    read_l2g_obs_at_pts(argv[iarg], xy_str, xy_fname, npt, res);
		  else 
		    read_l2l3_obs_at_pts(argv[iarg], xy_str, xy_fname, npt, res, gran_st);
		  
		  if (ameta_str != NULL) free(ameta_str);
		  if (cmeta_str != NULL) free(cmeta_str);
		}
	      iarg++;
	    }
	}
    }
  Free2D((void **)xy_str);
  fprintf(stderr, "Processing done ! \n");
  return status;
}

int parse_cmd_read_pixvals(char **argv, int argc, char **xy_str, int *pt_cnt, int *res, 
		char *xy_fname)

/*
!C************************************************************************************
       
!Description:
  Parse the command line arguments.

!Input Parameters:
  argc: number of arguments in the command line.
  argv: string array containing the input arguments.

!Input/Output Paramenters: (none)

!Output Parameters: 
  xy_str:    String array containing the input column and row number.
  pt_cnt:    Count of input pixels.
  res:       Reference resolution for xy. Default is the coarsest resolution of all the SDSs
             in the input files.
  xy_fname:  input file name of the file containing the input column and row number(optional)

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
  int st = 0;
  int i, fcnt, cnt, npt;
  int len_fstr;
  char res_str[10];
  
  st = 1;
  *res = npt = cnt = fcnt = 0;
  res_str[0] = xy_fname[0] = '\0';
  for (i=1; i<argc; i++)
    {
      if (is_arg_id(argv[i], "-res") == 0) get_arg_val(argv[i], res_str);
      else if (is_arg_id(argv[i], "-xy") == 0)
	{
	  if ((cnt == 0) && (sd_charpos(argv[i], ',', 0) == -1))
	    get_arg_val(argv[i], xy_fname);
	  else
	    {
	      get_arg_val_arr(argv[i], xy_str, &npt);
	      if (npt == cnt+2) cnt = npt;
	      else 
		{
		  fprintf(stderr, "Ignoring invalid entry %s\n", argv[i]);
		  npt = cnt;
		}
	    }
	}
      else if (argv[i][0] == '-') fprintf(stderr, "Ignoring invalid option %s\n", argv[i]);
      else fcnt++;
    } /* for (i=1, . . ) */
  
  if (fcnt == 0)
    {
      st = -1;
      fprintf(stderr, "Missing input files \n");
    }
  else
    {
      if (strcmp(res_str, "qkm") == 0) *res=250;
      else if (strcmp(res_str, "hkm") == 0) *res=500;
      else if (strcmp(res_str, "1km") == 0) *res=1000;
      else fprintf(stderr, "Reference resolution will be set to coarse SDS resolution\n");
      len_fstr = (int)strlen(xy_fname);
      if ((npt == 0) && (len_fstr <= 0))
	{
	  st = -1;
	  fprintf(stderr, "Input point sets or a point set file is required . . \n");
	}
    }
  *pt_cnt = npt/2;
  return st;
}

void get_xy_pts(char *xy_str, int *pt, int *xy_sh, int *xy_sq)
/*
!C************************************************************************************
       
!Description:
  Get the column and row number (0-based) from user specified resolution. If subsample 
  position is specified as additional subscript to row and column input, then set 
  the xy_sh and xy_sq.

!Input Parameters:
  xy_str:  String contains user specifed column and row number.

!Input/Output Paramenters: (none)

!Output Parameters: 

  xy_str:    string array containing the input column and row number.
  pt_cnt:
  res:       Reference resolution for xy. Default is the coarsest resolution of all the SDSs
             in the input files.
  xy_fname:  input file name of the file containing the input column and row number(optional)
  mask_str:  mask string.

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
  int len;
  int p1, p2;
  char tmp_str[10];
  
  *xy_sh = *xy_sq = 0;
  if ((p1 = sd_charpos(xy_str, '.', 0)) == -1)
    *pt = (int)atoi(xy_str);
  else
    {
      len = (int)strlen(xy_str);
      sd_strmid(xy_str, 0, p1, tmp_str);
      *pt = (int)atoi(tmp_str);
      p2 = sd_charpos(xy_str, '.', p1+1);
      if (p2 == -1)
	{
	  sd_strmid(xy_str, p1+1, len-p1-1, tmp_str);
	  *xy_sh = (int)atoi(tmp_str);
	}
      else
	{
	  sd_strmid(xy_str, p1+1, p2-p1-1, tmp_str);
	  *xy_sh = (int)atoi(tmp_str);
	  sd_strmid(xy_str, p2+1, len-p2-1, tmp_str);
	  *xy_sq = (int)atoi(tmp_str);
	}
    }
  if ((*xy_sh > 1) || (*xy_sh <0)) 
    {
      *xy_sh = 0;
      fprintf(stderr, "Invalid sub sample specified in %s. Set to default 0\n", xy_str);
    }
  if ((*xy_sq > 3) || (*xy_sq < 0)) 
    {
      *xy_sq = 0;
      fprintf(stderr, "Invalid sub sample specified in %s. Set to default 0\n", xy_str);
    }
}

void read_l2l3_obs_at_pts(char *fname, char **xy_str, char *xy_fname, int npt, int res, int gran_st)

/*
!C************************************************************************************
       
!Description:
  Read the value of pixels from MODLand L2 or L3 products and print the result to STOUT. 

!Input Parameters:
  fname:     Input filename.
  xy_str:    String contains user specifed column and row number.
  xy_fname:  Name of the file containing the column and row number of each pixel.
  npt:       Number of pixels.
  res:       Reference resolution for xy.
  gran_st:   Flag for tile information available or not.

!Revision History: (none)

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
**************************************************************************************/     

{
  FILE *fp = NULL;
  float *rf;
  int done = 1;
  int min_rows = 0, min_cols = 0;
  int max_ndata, max_isds;
  int x_sh = 0, y_sh = 0, x_sq = 0, y_sq = 0;
  int iobs, nobs, irow, icol;
  int pt_x = 0, pt_y = 0, f_val = 0, cres = 0;
  int id, idim, ipt, isds, nsds;
  int ndata, offset, **out_pnts;
  sds_t *sds_info;
  void *attr_buf, *data_in;
  int32 edge[4] = {1, 0, 0, 0};
  int32 start[4] = {0, 0, 0, 0};
  int32 attr_type, attr_cnt, **sds_val;
  char **sds_names, x_str[15], y_str[15];
  
  fprintf(stdout, "\nReading input file %s\n", fname);
  fprintf(stdout, "-----------------------------------------------------------------------\n");
  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for sds_names in read_l2l3_obs_at_pts\n");
      return;
    }
  if ((nsds = get_sds_names(fname, sds_names)) > 0)
    {
      if ((sds_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
	fprintf(stderr, "Cannot allocate memory for sds_info in read_l2l3_obs_at_pts\n");
      if ((rf = (float *)calloc(nsds, sizeof(float))) == NULL)
	fprintf(stderr, "Cannot allocate memory for rf in read_l2l3_obs_at_pts\n");
      if ((sds_info != NULL) && (rf != NULL))
	{
	  fprintf(stdout, "SDS in file and Fill values\n");
	  for (isds=0; isds<nsds; isds++)
	    {
	      sds_info[isds].sds_id = -1;
	      strcpy(sds_info[isds].name, sds_names[isds]);
	      sds_info[isds].sd_id = (isds == 0) ? -1 : sds_info[0].sd_id;
	      if (get_sds_info(fname, &sds_info[isds]) != -1)
		{
		  if ((attr_buf = get_sds_attr(sds_info[isds].sds_id, "_FillValue", &attr_type,
					       &attr_cnt)) != NULL)
		    {
		      switch(sds_info[isds].data_type)
			{
			case 20: f_val = ((int8 *)attr_buf)[0]; break;
			case 21: f_val = ((uint8 *)attr_buf)[0]; break;
			case 22: f_val = ((int16 *)attr_buf)[0]; break;
			case 23: f_val = ((uint16 *)attr_buf)[0]; break;
			case 24: f_val = ((int32 *)attr_buf)[0]; break;
			case 25: f_val = ((uint32 *)attr_buf)[0]; break;
                        default: 
			  fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported\n", 
				  sds_info[isds].data_type);
			}
		      free(attr_buf);
		      fprintf(stdout, "\t%s (Fill_Value = %d)\n", sds_info[isds].name, f_val);
		    }
		}
	    } /* for (isds=0; . . ) */
	  
	  if (res == 0)
	    {
	      min_rows = cres = sds_info[0].dim_size[0];
	      min_cols = sds_info[0].dim_size[1];
	      for (isds=1; isds<nsds; isds++)
		if (cres > sds_info[isds].dim_size[0]) 
		  {
		    min_rows = cres = sds_info[isds].dim_size[0];
		    min_cols =  sds_info[isds].dim_size[1];
		  }
	    }
	  else
	    {
	      switch(res)
		{
		case 250: cres = (gran_st == 1) ? 8120 : 4800; break;
		case 500: cres = (gran_st == 1) ? 4060 : 2400; break;
		case 1000: cres = (gran_st == 1) ? 2030 : 1200; break;
		default: fprintf(stderr, "Wrong resolution %d encountered\n", res); break;
		}
	    }
	  for (isds=0; isds<nsds; isds++)
	    rf[isds] = (float)sds_info[isds].dim_size[0]/(float)cres;
	  
	  if ((sds_val = (int32 **)Calloc2D(nsds, MAX_NUM_OBS, sizeof(int32))) == NULL)
	    fprintf(stderr, "Cannot allocate memory for sds_val in read_l2l3_obs_at_pts\n");
	  if ((out_pnts = (int **)Calloc2D(nsds, 2, sizeof(int))) == NULL)
	    fprintf(stderr, "Cannot allocate memory for out_pnts in read_l2l3_obs_at_pts\n");
	  if ((sds_val != NULL) && (out_pnts != NULL))
	    {
	      max_ndata = max_isds = 0;
	      for (isds=0; isds<nsds; isds++)
		if (sds_info[isds].sds_id != -1)
		  {
		    for (idim=1, ndata=1; idim<sds_info[isds].rank; idim++)
		      ndata *= sds_info[isds].dim_size[idim];
		    ndata *= sds_info[isds].data_size;
		    if (ndata > max_ndata) { max_ndata = ndata; max_isds = isds; }
		  }
	      ndata = max_ndata/sds_info[max_isds].data_size;
	      if ((data_in = (void *)calloc(ndata, sds_info[max_isds].data_size)) == NULL)
		fprintf(stderr, "Cannot allocate memory for data_in in read_l2l3_obs_at_pts\n");
	else
	  {
	    if (npt == 0)
	      {
		if ((fp = fopen(xy_fname, "r")) == NULL)
		  fprintf(stderr, "Cannot open file %s\n", xy_fname);
	      }
	    ipt = 0;
	    while (done)
	      {
		if (npt > 0)
		  { 
		    id = ipt*2;
		    get_xy_pts(xy_str[id], &pt_x, &x_sh, &x_sq);
		    get_xy_pts(xy_str[id+1], &pt_y, &y_sh, &y_sq);
		    ++ipt;
		    if (ipt > npt) break;
		  }
		else 
		  {
		    if (fscanf(fp, "%s %s", x_str, y_str) == 2)
		      {
			get_xy_pts(x_str, &pt_x, &x_sh, &x_sq);
			get_xy_pts(y_str, &pt_y, &y_sh, &y_sq);
		      }
		    if (feof(fp) != 0) { fclose(fp); break; }
		  }
		for (isds=0; isds<nsds; isds++)
		  {
		    irow = (int)((float)pt_y*rf[isds]);
		    icol = (int)((float)pt_x*rf[isds]);
		    if (rf[isds] == 2.0) { irow += y_sh; icol += x_sh; }
		    else if (rf[isds] == 4.0) { irow += y_sq; icol += x_sq; }
		    
		    start[0] = irow;
		    for (idim=1; idim<sds_info[isds].rank; idim++)
		      edge[idim] = sds_info[isds].dim_size[idim];
		    if (SDreaddata(sds_info[isds].sds_id, start, NULL, edge, data_in) == FAIL)
		      fprintf(stderr, "Cannot read data line for sds %s\n", sds_info[isds].name);
		    
		    for (idim=2, nobs=1; idim<sds_info[isds].rank; idim++)
		      nobs *= sds_info[isds].dim_size[idim];
		    offset = icol*nobs;
		    out_pnts[isds][0] = icol;
		    out_pnts[isds][1] = irow;
		    for (iobs=0; iobs<nobs; iobs++, offset++) 
		      {
			switch(sds_info[isds].data_type)
			  {
			  case 20: sds_val[isds][iobs] = ((int8 *)data_in)[offset]; break;
			  case 21: sds_val[isds][iobs] = ((uint8 *)data_in)[offset]; break;
			  case 22: sds_val[isds][iobs] = ((int16 *)data_in)[offset]; break;
			  case 23: sds_val[isds][iobs] = ((uint16 *)data_in)[offset]; break;
			  case 24: sds_val[isds][iobs] = ((int32 *)data_in)[offset]; break;
			  case 25: sds_val[isds][iobs] = ((uint32 *)data_in)[offset]; break;
			  default:
			    fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", 
				    sds_info[isds].data_type);
			  }
		      } /* for (iobs=0; . . ) */
		  } /* for (isds=0; . . . ) */
		if (res == 0)
		  fprintf(stdout, "Observations: pixel (%d %d) at ref resolution (%d x %d)\n", 
			  pt_x, pt_y, min_cols, min_rows);
		else
		  fprintf(stdout, "Observations: pixel (%d %d) at ref resolution %dm\n", pt_x, pt_y, res);
		print_sds_val(out_pnts, sds_val, sds_info, nsds, 0);
	      } /* for (ipt=0; . . ) */
	    free(data_in);
	  } /* else . . . */ 
	      for (isds=0; isds<nsds; isds++)
		if (sds_info[isds].sds_id != -1) 
		  SDendaccess(sds_info[isds].sds_id);
	      Free2D((void **)out_pnts);
	      Free2D((void **)sds_val);
	    }
	  
	  if (sds_info[0].sd_id != -1) 
	    SDend(sds_info[0].sd_id);
	}
      free(rf);
      free(sds_info); 
    } /* if ((nsds =  . . ) */
  Free2D((void **)sds_names);
}

void read_l2g_obs_at_pts(char *fname, char **xy_str, char *xy_fname, int npt, int res)

/*
!C************************************************************************************
       
!Description:
  Read the value of pixels from MODLand L2G products and print the result to STOUT. 

!Input Parameters:
  fname:     Input filename.
  xy_str:    String contains user specifed column and row number.
  xy_fname:  Name of the file containing the column and row number of each pixel.
  npt:       Number of pixels.
  res:       Reference resolution for xy.

!Revision History: (none)

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
**************************************************************************************/     

{
  FILE *fp = NULL;
  float *rf;
  int done = 1;
  int ic, ir, nobs_ic;
  int min_rows = 0, min_cols = 0;
  int id, ipt, iobs, nobs = 0;
  int st1, st2, f_val = 0, cres = 0;
  int pt_x = 0, pt_y = 0, **out_pnts;
  int ndata, max_size, offset;
  int x_sh = 0, y_sh = 0, x_sq = 0, y_sq = 0;
  int irow = 0, icol = 0, isds, nsds;
  int32 edge[2], start[2];
  void *attr_buf, *data_in;
  int8 *data_nobs;
  int32 *data_nadd_obs;
  int32 attr_type, attr_cnt, **sds_val;
  char **sds_names, x_str[15], y_str[15];
  sds_t *sds1_info, *sdsc_info, sds_nadd_obs_info;
  
  fprintf(stdout, "\nReading input file %s\n", fname);
  fprintf(stdout, "-----------------------------------------------------------------------\n");
  if ((sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for sds_names in read_l2g_obs_at_pts\n");
      return;
    }
  if ((nsds = get_l2g_sds_names(fname, sds_names)) > 0)
    {
      nsds++;
      if ((sds1_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
	fprintf(stderr, "Cannot allocate memory for sds1_info in read_l2g_obs_at_pts\n");
      if ((sdsc_info = (sds_t *)calloc(nsds, sizeof(sds_t))) == NULL)
	fprintf(stderr, "Cannot allocate memory for sdsc_info in read_l2g_obs_at_pts\n");
      if ((rf = (float *)calloc(nsds, sizeof(float))) == NULL)
	fprintf(stderr, "Cannot allocate memory for rf in read_l2g_obs_at_pts\n");
      if ((sds1_info != NULL) && (sdsc_info != NULL) && (rf != NULL))
	{
	  fprintf(stdout, "SDS in file and Fill values\n");
	  for (isds=0; isds<nsds; isds++)
	    {
	      sds1_info[isds].sds_id = sdsc_info[isds].sds_id = -1;
	      if (isds == 0)
		{
		  sds1_info[isds].sd_id = -1;
		  strcpy(sds1_info[isds].name, sds_name_nobs);
		  st1 = st2 = get_sds_info(fname, &sds1_info[isds]);
		}
	      else
		{
		  sprintf(sds1_info[isds].name, "%s_1", sds_names[isds-1]);
		  sprintf(sdsc_info[isds].name, "%s_c", sds_names[isds-1]);
		  sds1_info[isds].sd_id = sdsc_info[isds].sd_id = sds1_info[0].sd_id;
		  st1 = get_sds_info(fname, &sds1_info[isds]);
		  st2 = get_sds_info(fname, &sdsc_info[isds]);
		}
	      if ((st1 != -1) && (st2 != -1))
		{
		  if ((attr_buf = get_sds_attr(sds1_info[isds].sds_id, "_FillValue", &attr_type,
					       &attr_cnt)) != NULL)
		    {
		      switch(sds1_info[isds].data_type)
			{
			case 20: f_val = ((int8 *)attr_buf)[0]; break;
			case 21: f_val = ((uint8 *)attr_buf)[0]; break;
			case 22: f_val = ((int16 *)attr_buf)[0]; break;
			case 23: f_val = ((uint16 *)attr_buf)[0]; break;
			case 24: f_val = ((int32 *)attr_buf)[0]; break;
			case 25: f_val = ((uint32 *)attr_buf)[0]; break;
			default:
			  fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", 
				  sds1_info[isds].data_type);  
			}
		      free(attr_buf);
		    }
		  if (isds == 0)
		    fprintf(stdout, "\t%s (Fill_Value = %d)\n", sds1_info[isds].name, f_val);
		  else
		    fprintf(stdout, "\t%s(c) (Fill_Value = %d)\n", sds1_info[isds].name, f_val);
		}
	    } /* for (isds=0; . . ) */
	  
	  sds_nadd_obs_info.sds_id = -1;
	  sds_nadd_obs_info.sd_id = sds1_info[0].sd_id;
	  strcpy(sds_nadd_obs_info.name, sds_name_nadd_obs);
	  if (get_sds_info(fname, &sds_nadd_obs_info) == -1)
	    fprintf(stderr, "Result may be in error \n");
	  if ((data_nadd_obs = (int32 *)calloc(sds_nadd_obs_info.dim_size[0], sizeof(int32))) == NULL)
	    fprintf(stderr, "Cannot allocate memory for data_nadd_obs in read_l2g_obs_at_pts\n");
	  start[0] = 0; edge[0] = sds_nadd_obs_info.dim_size[0];
	  start[1] = edge[1] = 0;
	  if (SDreaddata(sds_nadd_obs_info.sds_id, start, NULL, edge, (VOIDP)data_nadd_obs) == FAIL)
	    fprintf(stderr, "Cannot read data line for sds %s\n", sds_nadd_obs_info.name);
	  
	  if (res == 0)
	    {
	      min_rows = cres = sds1_info[0].dim_size[0];
	      min_cols = sds1_info[0].dim_size[1];
	      for (isds=1; isds<nsds; isds++)
		if (cres > sds1_info[isds].dim_size[0]) 
		  {
		    min_rows = cres = sds1_info[isds].dim_size[0];
		    min_cols = sds1_info[isds].dim_size[1];
		  }
	    }
	  else
	    {
	      if (res == 250) cres = 4800;
	      else if (res == 500) cres = 2400;
	      else if (res == 1000) cres = 1200;
	      else fprintf(stderr, "Wrong resolution %d encountered\n", res);
	    }
	  for (isds=0; isds<nsds; isds++)
	    rf[isds] = (float)sds1_info[isds].dim_size[0]/(float)cres;
	  
      if ((sds_val = (int32 **)Calloc2D(nsds, MAX_NUM_OBS, sizeof(int32))) == NULL)
	fprintf(stderr, "Cannot allocate memory for sds_val in read_l2g_obs_at_pts\n");
      if ((out_pnts = (int **)Calloc2D(nsds, 2, sizeof(int))) == NULL)
	fprintf(stderr, "Cannot allocate memory for out_pnts in read_l2l3_obs_at_pts\n");
      if ((sds_val != NULL) && (out_pnts != NULL))
	{
	  ndata = sds1_info[0].dim_size[1];
	  if ((data_nobs = (int8 *)calloc(ndata, sizeof(int8))) == NULL)
	    fprintf(stderr, "Cannot allocate memory for data_nobs in read_l2g_obs_at_pts\n");
	  max_size = sds1_info[0].data_size;
	  for (isds=1; isds<nsds; isds++)
	    if (sds1_info[isds].data_size > max_size)
	      max_size = sds1_info[isds].data_size;
	  if ((data_in = (void *)calloc(ndata, max_size)) == NULL)
	    fprintf(stderr, "Cannot allocate memory for data_in in read_l2g_obs_at_pts\n");
	  
	  if ((data_nobs != NULL) && (data_in != NULL))
	    {
	      if (npt == 0)
		{
		  if ((fp = fopen(xy_fname, "r")) == NULL)
		    fprintf(stderr, "Cannot open file %s\n", xy_fname);
		}
	      ipt = 0;
	      while (done)
		{
		  if (npt > 0)
		    {
		      id = ipt*2;
		      get_xy_pts(xy_str[id], &pt_x, &x_sh, &x_sq);
		      get_xy_pts(xy_str[id+1], &pt_y, &y_sh, &y_sq);
		      ++ipt;
		      if (ipt > npt) break;
		    }
		  else
		    {
		      if (fscanf(fp, "%s %s", x_str, y_str) == 2)
			{
			  get_xy_pts(x_str, &pt_x, &x_sh, &x_sq);
			  get_xy_pts(y_str, &pt_y, &y_sh, &y_sq);
			}
		      if (feof(fp) != 0) { fclose(fp); break; }
		    }
		  for (isds=1; isds<nsds; isds++)
		    {
		      if ((sds1_info[isds].sds_id != -1) && (sdsc_info[isds].sds_id != -1)
			  && (sds1_info[0].sds_id != -1))
			{
			  irow = (int)((float)pt_y*rf[isds]);
			  icol = (int)((float)pt_x*rf[isds]);
			  if (rf[isds] == 2.0) { irow += y_sh; icol += x_sh; }
			  else if (rf[isds] == 4.0) { irow += y_sq; icol += x_sq; }
			  start[0] = irow; edge[0] = 1;
			  start[1] = 0; edge[1] = ndata;
			  out_pnts[isds][0] = icol;
			  out_pnts[isds][1] = irow;
			  if (isds == 1)
			    {
			      if (SDreaddata(sds1_info[0].sds_id, start, NULL, edge, (VOIDP)data_nobs) == FAIL)
				fprintf(stderr, "Cannot read data line for sds %s\n", sds1_info[0].name);
			      nobs = ((int8 *)data_nobs)[icol];
			      sds_val[0][0] = nobs;
			      out_pnts[0][0] = icol;
			      out_pnts[0][1] = irow;
			    }
			  if (nobs > 0)
			    {
			      if (SDreaddata(sds1_info[isds].sds_id, start, NULL, edge, (VOIDP)data_in) == FAIL)
				fprintf(stderr, "Cannot read data line for sds %s\n", sds1_info[isds].name);
			      switch(sds1_info[isds].data_type)
				{
				case 20: sds_val[isds][0] = ((int8 *)data_in)[icol]; break;
				case 21: sds_val[isds][0] = ((uint8 *)data_in)[icol]; break;
				case 22: sds_val[isds][0] = ((int16 *)data_in)[icol]; break;
				case 23: sds_val[isds][0] = ((uint16 *)data_in)[icol]; break;
				case 24: sds_val[isds][0] = ((int32 *)data_in)[icol]; break;
				case 25: sds_val[isds][0] = ((uint32 *)data_in)[icol]; break;
				default:
				  fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", 
					  sds1_info[isds].data_type);    
				}
			    }
			  if (nobs > 1)
			    {
			      for (ir=0, offset=0; ir<irow; ir++)
				offset += ((int32 *)data_nadd_obs)[ir];
			      for (ic=0; ic<icol; ic++)
				{
				  nobs_ic = ((int8 *)data_nobs)[ic];
				  if (nobs_ic > 1) offset = offset + nobs_ic - 1;
				}
			      start[0] = offset; edge[0] = nobs-1;
			      start[1] = edge[1] = 0;
			      if (SDreaddata(sdsc_info[isds].sds_id, start, NULL, edge, (VOIDP)data_in) == FAIL)
				fprintf(stderr, "Cannot read data line for sds %s\n", sdsc_info[isds].name);
			      for (iobs=1, offset=0; iobs<nobs; iobs++, offset++)
				switch(sdsc_info[isds].data_type)
				  {
				  case 20: sds_val[isds][iobs] = ((int8 *)data_in)[offset]; break;
				  case 21: sds_val[isds][iobs] = ((uint8 *)data_in)[offset]; break;
				  case 22: sds_val[isds][iobs] = ((int16 *)data_in)[offset]; break;
				  case 23: sds_val[isds][iobs] = ((uint16 *)data_in)[offset]; break;
				  case 24: sds_val[isds][iobs] = ((int32 *)data_in)[offset]; break;
				  case 25: sds_val[isds][iobs] = ((uint32 *)data_in)[offset]; break;
				  default:
				    fprintf(stderr, "HDF datatype " LONG_INT_FMT " not supported", 
					    sdsc_info[isds].data_type);  	
				  }
			    } /* if (nobs > 1) */
			} /* if (sds1_info . . . ) */
		    } /* for (isds=0; . . . ) */
		  if (res == 0)
		    fprintf(stdout, "Observations: pixel (%d %d) at ref resolution (%d x %d)\n", 
			    icol, irow, min_cols, min_rows);
		  else
		    fprintf(stdout, "Observations: pixel (%d %d) at ref resolution %dm\n", icol, irow, res);
		  print_sds_val(out_pnts, sds_val, sds1_info, nsds, 1);
		} /* for (ipt=0; . . . ) */
	      free(data_in);
	      free(data_nobs);
	    } /* if ((data_nobs != . . .  )) */
	  Free2D((void **)out_pnts);
	  Free2D((void **)sds_val);
	} /* if (sds_val != . . . ) */
      for (isds=0; isds<nsds; isds++)
	{
	  if (sds1_info[isds].sds_id != -1) 
	    SDendaccess(sds1_info[isds].sds_id);
	  if (sdsc_info[isds].sds_id != -1) 
	    SDendaccess(sdsc_info[isds].sds_id);
	} 
      if (sds1_info[0].sd_id != -1) 
        SDend(sds1_info[0].sd_id);
      free(rf);
      free(sds1_info); 
      free(sdsc_info); 
      free(data_nadd_obs);
	} /* if ((sds1_info != . . . )) */
    } /* if ((nsds =  . . ) */
  Free2D((void **)sds_names);
}

void print_sds_val(int **out_pnts, int32 **sds_val, sds_t *sds_info, int nsds, int l2g_st)

/*
!C************************************************************************************
       
!Description:
  Print the pixel value for each SDS to the STOUT. 

!Input Parameters:

  out_pnts:  Number of output points.
  sds_val:   Value of SDSs.
  sds_info:  Structure containing the SDS information.
  nsds:      number of SDSs.
  l2g_st:    If the input file is L2G product or not.

!Revision History: (none)

!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
**************************************************************************************/     

{
  int isds;
  int jobs, kobs, iobs;
  
  if (l2g_st == 1)
    {
      kobs = sds_val[0][0];
      fprintf(stdout, "\t%s: %d\n", sds_info[0].name, kobs);
      for (isds=1; isds<nsds; isds++)
	{
	  fprintf(stdout, "\t%s(c) at (%d %d): ", sds_info[isds].name, out_pnts[isds][0], 
		  out_pnts[isds][1]);
	  for (iobs=0; iobs<kobs; iobs++)
	    fprintf(stdout, LONG_INT_FMT " ", sds_val[isds][iobs]);
	  fprintf(stdout, "\n");
	}
    }
  else for (isds=0; isds<nsds; isds++)
    {
      fprintf(stdout, "\t%s at (%d %d): ", sds_info[isds].name, out_pnts[isds][0], out_pnts[isds][1]);
      if (sds_info[isds].rank == 2)
	fprintf(stdout, LONG_INT_FMT "\n", sds_val[isds][0]);
      else if (sds_info[isds].rank == 3)
	{
	  for (iobs=0; iobs<sds_info[isds].dim_size[2]; iobs++)
	    fprintf(stdout, LONG_INT_FMT " ", sds_val[isds][iobs]);
	  fprintf(stdout, "\n");
	}
      else if (sds_info[isds].rank == 4)
	{
	  fprintf(stdout, "\n");
	  for (jobs=0, iobs=0; jobs<sds_info[isds].dim_size[2]; jobs++)
	    {
	      fprintf(stdout, "\t\t\t");
	      for (kobs=0; kobs<sds_info[isds].dim_size[3]; kobs++, iobs++)
		fprintf(stdout, LONG_INT_FMT " ", sds_val[isds][iobs]);
	      fprintf(stdout, "\n");
	    }
	} /* else */
    } /* for (isds=0; . . . ) */
}
