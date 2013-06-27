/****************************************************************************
!C

!File: sds_rw.c

!Description:
  Routines to support SDS read and write 

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:

  Version 1.0 August, 2002

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

#include <string.h>

#include "isoc.h"
#include "mfhdf.h"
#include "str_op.h"
#include "sds_rw.h"
#include "qa_tool.h"
#include "alloc_mem.h"

int get_sds_info(char *hdf_fname, sds_t *sds_info)  
/*******************************************************************************
	Read the sds information about the names sds sds_name and store
	them in the structure sds_info. If the file and the sds are not open, 
	then open the file and then the sds. Return 1 on success and return -1 
	on failure.
*******************************************************************************/
{
  int p1;
  int32 attr_type, attr_cnt;
  char sds_name[MAX_SDS_NAME_LEN];
  char attr_name[MAX_ATTR_NAME_LEN];
  void *attr_val;

  if ((sds_info->sds_id == -1) && (strlen(sds_info->name) > 0)) 
    check_and_fix_sdsname(hdf_fname, sds_info);

  if (sds_info->sd_id == -1)
    if ((sds_info->sd_id = SDstart(hdf_fname, DFACC_READ)) == FAIL)
    {
      fprintf(stderr, "Cannot open the HDF file %s\n", hdf_fname);
      return -1;
    }
  if (sds_info->sds_id == -1)
  {
    if (strlen(sds_info->name) > 0)
    {
      if ((sds_info->sds_index = SDnametoindex(sds_info->sd_id, sds_info->name)) == FAIL)
      {
        if ((p1 = sd_charpos(sds_info->name, '.', 0)) != -1)
        {
          sd_strmid(sds_info->name, 0, p1, sds_name);
          sds_info->sds_index = SDnametoindex(sds_info->sd_id, sds_name);
        }
      }
    }
    if (sds_info->sds_index == FAIL)
    {
      fprintf(stderr, "Cannot find the SDS %s \n", sds_info->name);
      return -1;
    }

    if ((sds_info->sds_id = SDselect(sds_info->sd_id, sds_info->sds_index)) == FAIL)
    {
      fprintf(stderr, "Cannot open the SDS %s \n", sds_name);
      return -1;
    }
  }
  if (SDgetinfo(sds_info->sds_id, sds_name, &sds_info->rank, sds_info->dim_size,
        &sds_info->data_type, &sds_info->nattr) == FAIL)
  {   
    fprintf(stderr, "Cannot get information for the SDS %s \n", sds_name);
    return -1;
  }
  else
    sds_info->data_size = DFKNTsize(sds_info->data_type);
  strcpy(attr_name, "_FillValue");
  if ((attr_val = (void *)get_sds_attr(sds_info->sds_id, attr_name, &attr_type, &attr_cnt)) != NULL)
  {
    switch(attr_type)
    {
      case 5 : sds_info->fill_fval = ((float32 *)attr_val)[0]; break;
      case 20: sds_info->fill_val = ((int8 *)attr_val)[0]; break;
      case 21: sds_info->fill_val = ((uint8 *)attr_val)[0]; break;
      case 22: sds_info->fill_val = ((int16 *)attr_val)[0]; break;
      case 23: sds_info->fill_val = ((uint16 *)attr_val)[0]; break;
      case 24: sds_info->fill_val = ((int32 *)attr_val)[0]; break;
      case 25: sds_info->fill_val = ((uint32 *)attr_val)[0]; break;
    }
    free(attr_val);
  }
  strcpy(attr_name, "valid_range");
  if ((attr_val = (void *)get_sds_attr(sds_info->sds_id, attr_name, &attr_type, &attr_cnt)) != NULL)
  {
    if (attr_cnt != 2)
      fprintf(stderr, "SDS range value in error: Does not contain two values\n");
    switch(attr_type)
    {
      case 5 : sds_info->frange[0] = ((float32 *)attr_val)[0]; 
	       sds_info->frange[1] = ((float32 *)attr_val)[1];
	       break;
      case 20: sds_info->range[0] = ((int8 *)attr_val)[0]; 
	       sds_info->range[1] = ((int8 *)attr_val)[1];
	       break;
      case 21: sds_info->range[0] = ((uint8 *)attr_val)[0]; 
	       sds_info->range[1] = ((uint8 *)attr_val)[1];
	       break;
      case 22: sds_info->range[0] = ((int16 *)attr_val)[0]; 
	       sds_info->range[1] = ((int16 *)attr_val)[1];
	       break;
      case 23: sds_info->range[0] = ((uint16 *)attr_val)[0]; 
	       sds_info->range[1] = ((uint16 *)attr_val)[1];
	       break;
      case 24: sds_info->range[0] = ((int32 *)attr_val)[0]; 
	       sds_info->range[1] = ((int32 *)attr_val)[1];
	       break;
      case 25: sds_info->range[0] = ((uint32 *)attr_val)[0]; 
	       sds_info->range[1] = ((uint32 *)attr_val)[1];
	       break;
    }
    free(attr_val);
  }                            
  return 1;
}


int get_sds_names(char *fname, char **sds_names)
/*
 * Retrieve the SDS names from the input HDF file.
 */
{
  int i;
  int32 dt;
  int32 rank;
  int32 nsds;
  int32 nattr;
  int32 sd_id;
  int32 sds_id;
  int32 dim_size[3];
  int sds_cnt;
  char name[80];

  sds_cnt = 0;
  if ((sd_id = SDstart(fname, DFACC_READ)) == FAIL)
  {
    fprintf(stderr, "Cannot open the HDF file %s\n", fname);
    return sds_cnt;
  }
  if (SDfileinfo(sd_id, &nsds, &nattr) == FAIL)
    fprintf(stderr, "Cannot read information for HDF file %s\n", fname);
  else for (i=0; i<nsds; i++)
  {
    if ((sds_id = SDselect(sd_id, i)) != FAIL)
    {
      if (SDgetinfo(sds_id, name, &rank, dim_size, &dt, &nattr) != FAIL)
        strcpy(sds_names[sds_cnt++], name);
      SDendaccess(sds_id);
    }
  }
  SDend(sd_id);
  return sds_cnt;
}

void *get_sds_attr(int32 sds_id, char *attr_name, int32 *attr_type, int32 *attr_cnt)
/* retrieve the sds attribute like attribute name, attribute type, attribute count */
{
  int status;
  int32 attr_size;
  int32 attr_index;
  void *attr_buf = NULL;
  char name[80];

  status = -1;
  if ((attr_index = SDfindattr(sds_id, attr_name)) != -1)
  {
    if (SDattrinfo(sds_id, attr_index, name, attr_type, attr_cnt) == FAIL)
      fprintf(stderr, "Cannot read the attribute information for %s\n", attr_name);
    else
    {
      attr_size = DFKNTsize(*attr_type);
      if ((attr_buf = (void *)calloc(*attr_cnt, attr_size)) == NULL)
	fprintf(stderr, "Cannot allocate memory for attr_buf in get_sds_attr\n");
      else
      {
        if (SDreadattr(sds_id, attr_index, attr_buf) == FAIL)
          fprintf(stderr, "Cannot read sds_attr %s in get_sds_attr\n", attr_name);
	else 
	  status = 1;
      }
    }
  }
  if (status == -1) return (void *)NULL;
  else return attr_buf;
}

void update_nd_sdsnames(char **sds_names, int *sds_cnt, char *fname)
/* 
   Create new sds_names.
        e.g. input sds_name is Surf_refl.1.2-3
             The resulting sds_names will be Surf_refl.1.2
	                                     Surf_refl.1.3
	     The sds_cnt will be set to 2 since there are 2 resulting sds names
   Updated to fix the problem with SDS names containing "." character 
   	by Xiaoping Zhang, Sigma Space Corp. on Jun 2012
*/
{
  int nsds, isds;
  int p1, p2, p11, len;
  char **new_sds_names;
  char nd_ext[20], md_ext[20];
  char sds_name[MAX_SDS_NAME_LEN];
  sds_t sds_info;

  memset( &sds_info, 0, sizeof( sds_t ) );
  if ((new_sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for new_sds_name in update_nd_sdsnames\n");
  else
    {
      nd_ext[0] = md_ext[0] = '\0';
      for (isds=0, nsds=0; isds < *sds_cnt; isds++)
	{
	  if (sd_charpos(sds_names[isds], '(', 0) == -1)
	{      
	  if ((p1 = sd_charpos(sds_names[isds], '.', 0)) != -1)
	  {		  
	    sds_info.sd_id = sds_info.sds_id = -1;		  
	    sd_strmid(sds_names[isds], 0, p1, sds_name);
	    strcpy(sds_info.name, sds_name);
	    if (get_sds_info(fname, &sds_info) != -1)
	    {
	      p1++;
	      len = (int)strlen(sds_names[isds]);
	      if ((p2 = sd_charpos(sds_names[isds], '.', p1)) != -1) 
	      {
		sd_strmid(sds_names[isds], p1, p2-p1, nd_ext); 
		p2++;
		sd_strmid(sds_names[isds], p2, len-p2, md_ext); 
		create_names(fname, sds_name, nd_ext, md_ext, &nsds, new_sds_names);	
	      }
	      else
	      {
		sd_strmid(sds_names[isds], p1, len-p1, nd_ext);
		create_names(fname, sds_name, nd_ext, md_ext, &nsds, new_sds_names);	
	      }
	    } 
	    else
	    {
	      p1++;
	      if ((p11 = sd_charpos(sds_names[isds], '.', p1)) != -1) 
	      {		      		      		      
		sd_strmid(sds_names[isds], 0, p11, sds_name);
		strcpy(sds_info.name, sds_name);		
		len = (int)strlen(sds_names[isds]);
		if (get_sds_info((char *)NULL, &sds_info) != -1)
		{		
		  p11++;
		  if ((p2 = sd_charpos(sds_names[isds], '.', p11)) != -1) 
		  {		    
		    sd_strmid(sds_names[isds], p11, p2-p11, nd_ext); 
		    p2++;
		    sd_strmid(sds_names[isds], p2, len-p2, md_ext); 
		    create_names(fname, sds_name, nd_ext, md_ext, &nsds, new_sds_names);	
		  }
		  else
		  {
		    sd_strmid(sds_names[isds], p11, len-p11, nd_ext);
		    create_names(fname, sds_name, nd_ext, md_ext, &nsds, new_sds_names);	
		  }
	        }
		else
		  fprintf(stderr, "Cannot find the SDS %s \n", sds_info.name);			      
	      }	
	      else
	      {
		strcpy(new_sds_names[nsds], sds_names[isds]);
		nsds++;	      
	      }	    		    		    
	    }		   
	  }
	  else 
	  {
	    strcpy(new_sds_names[nsds], sds_names[isds]);
	    nsds++;
	  }
	}
	else
	{
	  strcpy(new_sds_names[nsds], sds_names[isds]);
	  nsds++;
	}
	if (sds_info.sds_id != -1) SDendaccess(sds_info.sds_id);
	if (sds_info.sd_id != -1) SDendaccess(sds_info.sd_id);	
      }
      for (isds=0; isds<nsds; isds++)
	strcpy(sds_names[isds], new_sds_names[isds]);
      *sds_cnt = nsds;
      
      Free2D((void **)new_sds_names);
    }
}


void check_and_fix_sdsname(char *hdf_fname, sds_t *sds_info)
/* check the user input sds name is valid or not, if there is a typo in 
   the input name, use the current name as the output sds name */
{
  int i;
  int32 nsds, nattr;
  int32 sd_id, sds_id;
  int32 rank, dt, dim_size[4];
  char name[MAX_SDS_NAME_LEN];

  sd_id = sds_info->sd_id;
  if (sd_id == -1) 
    sd_id = SDstart(hdf_fname, DFACC_READ);
  if (sd_id != -1)
  {
    if (SDfileinfo(sd_id, &nsds, &nattr) != FAIL)
    {
      for (i=0; i<nsds; i++)
      {
        if ((sds_id = SDselect(sd_id, i)) != FAIL)
        {
          if (SDgetinfo(sds_id, name, &rank, dim_size, &dt, &nattr) != FAIL)
          {
	        if (strcasecmp(name, sds_info->name) == 0)
	        strcpy(sds_info->name, name);
          }
          SDendaccess(sds_id);
        }
      }
    }
    if (sds_info->sd_id == -1)
      SDendaccess(sd_id);
/*    if (sds_info->sd_id == -1)
      SDend(sd_id); */
  }
}

void compute_sds_start_offset(sds_t *sds_info, int n, int m, int *st_c, int *offset)
/* Given a n and m(which represent the user selected number for dimension 3 and 4, 
   compute where to get the SDS. */
{
  int rank;

  rank = sds_info->rank;
  if ((rank == 2) || ((n == -1) && (m == -1)))
  {
    *st_c = 0;
    *offset = 1;
  }
  else
  {
    if ((n != -1) && (m == -1))
    {
      if (sds_info->dim_size[0] > sds_info->dim_size[rank-1])
      {
        *st_c = n;
        *offset = sds_info->dim_size[rank-1];
      }
      else
      {
        *st_c = n*sds_info->dim_size[rank-1];
        *offset = 1;
      }
    }
    if ((n != -1) && (m != -1))
    {
      if (sds_info->dim_size[0] > sds_info->dim_size[rank-1])
      {
        *st_c = n*sds_info->dim_size[rank-1] + m;
        *offset = sds_info->dim_size[rank-1]*sds_info->dim_size[rank-2];
      }
      else
      {
        *st_c = n*sds_info->dim_size[rank-1]*sds_info->dim_size[1] +  m*sds_info->dim_size[rank-1];
        *offset = 1;
      }
    }
  }
}


int open_sds(char *fname, sds_t *sds_info, char open_t)
/* open an SDS from a HDF file */
{
  int status;
  int p1, len;
  int32 sds_index;
  char sds_name[MAX_SDS_NAME_LEN];

  status = 1;
  if (open_t == 'R')
  {
    if (sds_info->sd_id == -1)
    {
      if ((sds_info->sd_id = SDstart(fname, DFACC_READ)) == FAIL)
      {
	fprintf(stderr, "Cannot open the hdf file %s\n", fname);
	status = -1;
      }
    }
    if (status != -1)
    {
      if ((sds_index = SDnametoindex(sds_info->sd_id, sds_info->name)) == FAIL)
      {
        if ((p1 = sd_charpos(sds_info->name, '.', 0)) != -1)
	{
          sd_strmid(sds_info->name, 0, p1, sds_name);
          sds_index = SDnametoindex(sds_info->sd_id, sds_name);
	}
      }
      if (sds_index == FAIL)
      {
        fprintf(stderr, "Cannot find the SDS %s in file %s\n", sds_name, fname);
        status = -1;
      }
      else if ((sds_info->sds_id = SDselect(sds_info->sd_id, sds_index)) == FAIL)
      {
        fprintf(stderr, "Cannot open the SDS %s in file %s\n", sds_name, fname);
        status = -1;
      }
    }
  }
  else
  {
    if (sds_info->sd_id == -1)
    {
      if ((sds_info->sd_id = SDstart(fname, DFACC_CREATE)) == FAIL)
      {
        fprintf(stderr, "Cannot create the hdf file %s\n", fname);
        status = -1;
      }
    }
    if (status != -1)
    {
      if ((sds_info->sds_id = SDcreate(sds_info->sd_id, sds_info->name, 
		sds_info->data_type, sds_info->rank, sds_info->dim_size)) == FAIL)
      {
        fprintf(stderr, "Cannot create the SDS %s\n", sds_info->name);
        status = -1;
      }
      else
      {
        len = (int)strlen(sds_info->name);
        if (SDsetattr(sds_info->sds_id, "long_name", DFNT_CHAR8, len, 
			(VOIDP)sds_info->name) == FAIL)
          fprintf(stderr, "Could not write the attribute 'long_name' in %s\n", sds_info->name);
      }
    }
  }
  return status;
}


void create_names(char *fname, char *sds_name, char *nd_ext, char *md_ext, int *nsds, char **sds_names)
/********************************************************************************************
          create sds_names with layers surfix.
*********************************************************************************************/
{
  int n, m;
  int ndim = 0, mdim = 0;
  int n_cnt, m_cnt;
  int narr[20], marr[20];
  sds_t sds_info;

  sds_info.sd_id = sds_info.sds_id = -1;
  strcpy(sds_info.name, sds_name);
  if (get_sds_info(fname, &sds_info) != -1)
    {
      n_cnt = m_cnt = 0;
      if (sds_info.rank > 2)
	{
	  if (sds_info.dim_size[0] > sds_info.dim_size[2]) 
	    { ndim = sds_info.dim_size[2]; mdim = sds_info.dim_size[3]; }
	  else 
	    { ndim = sds_info.dim_size[0]; mdim = sds_info.dim_size[1]; }
	  get_dim_num(nd_ext, narr, ndim, &n_cnt);
	  if (sds_info.rank > 3) get_dim_num(md_ext, marr, mdim, &m_cnt);
	}
      if (sds_info.rank == 3)
	{
	  for (n=0; n<n_cnt; n++, ++*nsds)
	    {
	      if (narr[n] > ndim || narr[n] < 1)
		{ 
		  fprintf(stderr, "Invalid index for SDS %s, valid index range is 1-%d \n",
			  sds_name, ndim);
		  exit(EXIT_FAILURE);
		}
	      sprintf(sds_names[*nsds], "%s.%d", sds_name, narr[n]);
	    }
	}
      else if (sds_info.rank == 4)
	{
	  for (n=0; n<n_cnt; n++)
	    {
	      if (narr[n] > ndim || narr[n] < 1)
		{ 
		  fprintf(stderr, "Invalid index for SDS %s, valid index range is 1-%d \n",
			  sds_name, ndim);
		  exit(EXIT_FAILURE);
		}
	      for (m=0; m<m_cnt; m++, ++*nsds)
		{
		  if (marr[m] > mdim || marr[m] < 1)
		    { 
		      fprintf(stderr, "Invalid index for SDS %s, valid index range is 1-%d,1-%d \n",
			      sds_name, ndim, mdim);
		      exit(EXIT_FAILURE);
		    }
		  sprintf(sds_names[*nsds], "%s.%d.%d", sds_name, narr[n], marr[m]); 
		}
	    }
	}
    }
  if (sds_info.sds_id != -1) SDendaccess(sds_info.sds_id);
  if (sds_info.sd_id != -1) SDendaccess(sds_info.sd_id);
}

void get_dim_num(char *str_arr, int *num_arr, int dim_size, int *cnt)
/* get the dimension number */
{
  int len;
  char c;
  int i, j, k;
  int bn1, bn2;
  char nstr[5];

  if (str_arr[0] == '*')
  {
    for (i=0; i<dim_size; i++)
      num_arr[i] = i+1;
    *cnt = dim_size;
  }
  else
  {
    c  = ',';
    len = (int)strlen(str_arr);
    for (i=0, j=0; i<=len; i++)
    {
      if ((i<len) && (str_arr[i] != ',') && (str_arr[i] != '-'))
        nstr[j++] = str_arr[i];    
      else 
      {
	nstr[j] = '\0';
	j = 0;
	if (c == ',') { num_arr[*cnt] = (int)atoi(nstr); ++*cnt; c = str_arr[i]; }
        else if (c == '-') 
        {
	  bn1 = num_arr[*cnt - 1];
	  bn2 = (int)atoi(nstr);
	  for (k=bn1+1; k<=bn2; k++) { num_arr[*cnt] = k; ++*cnt; }
	}	 
      }
    }
  }
} 


void write_attr_fval(int32 sds_id, int32 fval_type, int c, int attr_val, char *attr_name)
/* wrtie the fill value for a attribute */
{
  int32 val_len;
  int status = 1;
  void *fill_val = NULL;

  val_len = 1;
  switch (fval_type)
  {
    case 5: if ((fill_val = (float32 *)calloc(1, sizeof(float32))) != NULL)
#ifdef COMP_SDS_DIFF_IMG
	       ((float32 *)fill_val)[0] = (c == 0) ? (float32)COMP_FILL_VALUE_INT32 : (float32)attr_val;
#else
	       ((float32 *)fill_val)[0] = (c == 0) ? (float32)FILL_VALUE_INT32 : (float32)attr_val;
#endif
	     else status = -1;
	     break; 
    case 6: if ((fill_val = (float64 *)calloc(1, sizeof(float64))) != NULL)
#ifdef COMP_SDS_DIFF_IMG
	       ((float64 *)fill_val)[0] = (c == 0) ? (float64)COMP_FILL_VALUE_INT32 : (float64)attr_val;
#else
	       ((float64 *)fill_val)[0] = (c == 0) ? (float64)FILL_VALUE_INT32 : (float64)attr_val;
#endif
	     else status = -1;
	     break; 
    case 20: if ((fill_val = (int8 *)calloc(1, sizeof(int8))) != NULL)
	       ((int8 *)fill_val)[0] = (c == 0) ? (int8)FILL_VALUE_INT8 : (int8)attr_val;
	     else status = -1;
	     break; 
    case 21: if ((fill_val = (uint8 *)calloc(1, sizeof(uint8))) != NULL)
	       ((uint8 *)fill_val)[0] = (c == 0) ? (uint8)FILL_VALUE_UINT8 : (uint8)attr_val;
	     else status = -1;
	     break; 
    case 22: if ((fill_val = (int16 *)calloc(1, sizeof(int16))) != NULL)
#ifdef COMP_SDS_DIFF_IMG
	       ((int16 *)fill_val)[0] = (c == 0) ? (int16)COMP_FILL_VALUE_INT16 : (int16)attr_val;
#else
	       ((int16 *)fill_val)[0] = (c == 0) ? (int16)FILL_VALUE_INT16 : (int16)attr_val;
#endif
	     else status = -1;
	     break; 
    case 23: if ((fill_val = (uint16 *)calloc(1, sizeof(uint16))) != NULL)
	       ((uint16 *)fill_val)[0] = (c == 0) ? (uint16)FILL_VALUE_UINT16 : (uint16)attr_val;
	     else status = -1;
	     break; 
    case 24: if ((fill_val = (int32 *)calloc(1, sizeof(int32))) != NULL)
#ifdef COMP_SDS_DIFF_IMG
	       ((int32 *)fill_val)[0] = (c == 0) ? (int32)COMP_FILL_VALUE_INT32 : (int32)attr_val;
#else
	       ((int32 *)fill_val)[0] = (c == 0) ? (int32)FILL_VALUE_INT32 : (int32)attr_val;
#endif
	     else status = -1;
	     break; 
    case 25: if ((fill_val = (uint32 *)calloc(1, sizeof(uint32))) != NULL)
	       ((uint32 *)fill_val)[0] = (c == 0) ? (uint32)FILL_VALUE_UINT32 : (uint32)attr_val;
	     else status = -1;
	     break; 
  }
  if (status == -1)
    fprintf(stderr, "Cannot allocate memory for fill_val\n");
  else
  {
    if (SDsetattr(sds_id, attr_name, fval_type, val_len, (VOIDP)fill_val) == FAIL)
      fprintf(stderr, "Cannot write the attribute %s in write_attr_fval\n", attr_name);
    free(fill_val);
  }
}


int get_l2g_sds_names(char *fname, char **sds_names)
/* get the sds_names from a l2g file. */
{
  int i;
  int len;
  int32 dt;
  int32 rank;
  int32 nsds;
  int32 nattr;
  int32 sd_id;
  int32 sds_id;
  int32 dim_size[3];
  int sds_cnt;
  char tmp[5];
  char name[80];

  sds_cnt = 0;
  if ((sd_id = SDstart(fname, DFACC_READ)) == FAIL)
    fprintf(stderr, "Cannot open the HDF file %s\n", fname);
  else
  {
    if (SDfileinfo(sd_id, &nsds, &nattr) == FAIL)
      fprintf(stderr, "Cannot read the information of file %s\n", fname);
    else for (i=0; i<nsds; i++)
    {
      if ((sds_id = SDselect(sd_id, i)) != FAIL)
      {
        if (SDgetinfo(sds_id, name, &rank, dim_size, &dt, &nattr) != FAIL)
        {
          len = (int)strlen(name);
	  sd_strmid(name, len-2, 2, tmp);
	  if (strcmp(tmp, "_c") == 0)
          { 
	    strcpy(sds_names[sds_cnt], name);
	    sds_names[sds_cnt][len-2] = '\0';
	    sds_cnt++; 
	  }
        }
        SDendaccess(sds_id);
      }
    }
    SDend(sd_id);
  }
  return sds_cnt;
}

int get_sds_data(sds_t *sds_info, void *data)
/*********************************************************************************
	Read the data for the named sds sds->name and return the SDS data in
	data. Return 1 on success and -1 on failure.
*********************************************************************************/
{
  int idim;
  int status;
  int32 *edge;
  int32 *start;

  status = 1;
  if ((start = (int32 *)calloc(sds_info->rank, sizeof(int32))) == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for start in get_sds_data\n");
    return -1;
  }
  if ((edge = (int32 *)calloc(sds_info->rank, sizeof(int32))) == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for edge in get_sds_data\n");
    return -1;
  }
  for (idim = 0; idim < sds_info->rank; idim++)
  {
    start[idim] = 0;
    edge[idim] = sds_info->dim_size[idim];
  }

  if ((sds_info->sds_index = SDnametoindex(sds_info->sd_id, sds_info->name)) == FAIL)
  {
    fprintf(stderr, "Cannot find the SDS %s in get_sds_data\n", sds_info->name);
    status = -1;
  }
  if ((sds_info->sds_id = SDselect(sds_info->sd_id, sds_info->sds_index)) == FAIL)
  {
    fprintf(stderr, "Cannot read the SDS %s in get_sds_data\n", sds_info->name);
    status = -1;
  }
  if (SDreaddata(sds_info->sds_id, start, NULL, edge, data) == FAIL)
  {
    fprintf(stderr, "Cannot read the SDS data %s in get_sds_data\n", sds_info->name);
    status = -1;
  }

  if (sds_info->sds_id != -1) SDendaccess(sds_info->sds_id);
  free((int32 *)start);
  free((int32 *)edge);
  return status;
}

void print_sds_dim_size(sds_t *sds_info)
/* This routine prints out dimension name and dimension size of a SDS. */ 
{
  int i, rank;
  char **dim_names, **short_dim_names;
  
  rank = sds_info->rank;
  
  dim_names = (char **)Calloc2D(rank, MAX_NAME_LENGTH, sizeof(char));
  short_dim_names = (char **)Calloc2D(rank, MAX_DIM_NAME_LEN, sizeof(char));

  fprintf(stderr, "     %d dimensions : \n", rank);

  get_sds_dim_name(sds_info, dim_names, short_dim_names);

  for (i=0; i<rank; i++) 
    {
      printf("          dimension %i    dimension name : %s    dim size : " LONG_INT_FMT "\n",
           i+1, short_dim_names[i], sds_info->dim_size[i]);    
    }
  Free2D((void **)dim_names);
  Free2D((void **)short_dim_names);
}

void get_sds_dim_name(sds_t *sds_info, char **dim_names, char **short_dim_names)

     /* this routine takes input sds_info and output the dimension name and 
	short dimension names with ':'
	eg, As in file type MYD43B1, one dimension name is YDim:MOD_Grid_BRDF,
	the short dimension name is YDim */
{   
  int i, rank, p1;
  int32 cnt, nattrs, id, dt;
  
  rank = sds_info->rank;
  
  for (i=0; i<rank; i++)
    { 
      id = SDgetdimid(sds_info->sds_id, i); 
      if (SDdiminfo(id, dim_names[i], &cnt, &dt, &nattrs) == FAIL)
	{
          fprintf(stderr, "Error reading %dth dimension info in get_sds_dim_name. Program exit.\n", i);
          exit(EXIT_FAILURE);
	}

      p1 = sd_charpos(dim_names[i], ':', 0);

      if (p1 != -1)
	{
	  sd_strmid(dim_names[i], 0, p1, short_dim_names[i]);
	}
      else
	{
	  strcpy(short_dim_names[i], dim_names[i]);
	}
    }
}  


void update_l2g_sdsnames(char **sds_names, int *sds_cnt, char *fname, int nobs)
{
  sds_t sds_info;
  int p1, len, n, n_cnt;
  int nsds, isds, narr[255];
  char sds_name[MAX_SDS_NAME_LEN];
  char **new_sds_names, nd_ext[20];

  if ((new_sds_names = (char **)Calloc2D(MAX_NUM_SDS, MAX_SDS_NAME_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for new_sds_name in sds_rw: update_l2g_sdsnames()\n");
  else
  {
    nd_ext[0] = '\0';
    for (isds=0, nsds=0; isds < *sds_cnt; isds++)
    {
      if ((p1 = sd_charpos(sds_names[isds], '.', 0)) == -1)
	fprintf(stderr, "The l2g SDS name %s is not accepted in sds_rw: update_l2g_sdsnames()\n", 
			sds_names[isds]);
      else
      {
        sd_strmid(sds_names[isds], 0, p1, sds_name);
        len = (int)strlen(sds_names[isds]);
        sd_strmid(sds_names[isds], p1+1, len-p1-1, nd_ext);
        sds_info.sd_id = sds_info.sds_id = -1;
        strcpy(sds_info.name, sds_name);
	strcat(sds_info.name, "_1");
        if (get_sds_info(fname, &sds_info) != -1)
        {
          n_cnt = 0;
          get_dim_num(nd_ext, narr, nobs, &n_cnt);
          for (n=0; n<n_cnt; n++, ++nsds)
              sprintf(new_sds_names[nsds], "%s.%d", sds_name, narr[n]);
        }
        if (sds_info.sds_id != -1) SDendaccess(sds_info.sds_id);
        if (sds_info.sd_id != -1) SDendaccess(sds_info.sd_id);
      }
    }
    for (isds=0; isds<nsds; isds++)
      strcpy(sds_names[isds], new_sds_names[isds]);
    *sds_cnt = nsds;
    Free2D((void **)new_sds_names);
  }
}

void write_all_sds_attrs(int32 in_sds_id, int32 out_sds_id, int32 nattr)
/* 
 * copy all the sds_attribute of the input SDS to the output SDS.
 */
{
  void *attr_buf;
  int32 attr_index, attr_cnt;
  int32 attr_type, attr_size;
  char attr_name[MAX_ATTR_NAME_LEN];
                                                                                                                  
  for (attr_index=0; attr_index<nattr; attr_index++)
  {
    if (SDattrinfo(in_sds_id, attr_index, attr_name, &attr_type, &attr_cnt) != FAIL)
    {
      attr_size = DFKNTsize(attr_type);
      if ((attr_buf = (void *)calloc(attr_cnt, attr_size)) == NULL)
        fprintf(stdout, "Cannot allocate memory for attr_buf in write_all_sds_attrs\n");
      else
      {
        if (SDreadattr(in_sds_id, attr_index, attr_buf) != FAIL)
          SDsetattr(out_sds_id, attr_name, attr_type, attr_cnt, attr_buf);
        free(attr_buf);
      }
    }
  }
}

void compute_sds_nrows_ncols(sds_t *sds_info, int *nrows, int *ncols)
/*
 * Compute the total number of rows and columns. For a 2D SDS, just return
 * the number of rows and number of columns. For 3D/4D SDS, it's the total 
 * number of layers * the number of rows or columns.
 */
{
  int nr, nc;
  int i, rank;

  rank = sds_info->rank;
  if (rank == 2)
  {
    nr = sds_info->dim_size[0];
    nc = sds_info->dim_size[1];
  }
  else
  {
    if (sds_info->dim_size[0] < sds_info->dim_size[rank-1])
    {
      nr = sds_info->dim_size[rank-2];
      nc = sds_info->dim_size[rank-1];
      for (i=0; i<rank-2; i++)
        nc = nc*sds_info->dim_size[i];
    }
    else
    {
      nr = sds_info->dim_size[0];
      nc = sds_info->dim_size[1];
      for (i=2; i<rank; i++)
        nc = nc*sds_info->dim_size[i];
    }                                                     
  }
  *nrows = nr;
  *ncols = nc;
}

int compute_sds_ndata(sds_t *sds_info)
/*
 * compute number of rows in a SDS. e.g. for a 3D SDS, the ndata will be
 * number_of_row * number_of_layer
 */

{
  int i, rank, ndata;

  rank = sds_info->rank;
  if (rank == 2)
    ndata = sds_info->dim_size[1];
  else
  {
    if (sds_info->dim_size[0] < sds_info->dim_size[rank-1])
    {
      ndata = sds_info->dim_size[rank-1];
      for (i=0; i<rank-2; i++)
        ndata *= sds_info->dim_size[i];
    }
    else
    {
      ndata = sds_info->dim_size[1];
      for (i=2; i<rank; i++)
        ndata *= sds_info->dim_size[i];
    }
  }
  return ndata;
}

void get_sds_edge(sds_t *sds_info, int32 *edge)
/* 
 * get the edge information of the input SDS * 
 */
{
  int i, rank;

  rank = sds_info->rank;
  if (rank == 2)
  {
    edge[0] = 1;
    edge[1] = sds_info->dim_size[1];
  }
  else if (rank > 2)
  {
    for (i=0; i<rank; i++)
      edge[i] = sds_info->dim_size[i];
    if (sds_info->dim_size[0] < sds_info->dim_size[rank-1])
      edge[rank-2] = 1;
    else edge[0] = 1;
  }
}                                                           

void get_sds_param(sds_t *sds_info, int *n, int *m, int *rank, int *dim_size)
/*
 * retrieve the information of n, m(layer number of 3D or 4D SDS, rank and 
 * dimension size of the SDS.
 */

{
  int i, rank1;
  char sds_name[MAX_SDS_NAME_LEN];

  get_sdsname_dim(sds_info->name, sds_name, n, m);
  if ((*n == -1) && (*m == -1))
  {
    *rank = sds_info->rank;
    for (i=0; i < *rank; i++)
      dim_size[i] = sds_info->dim_size[i];
  }
  else
  {
    *rank = 2;
    rank1 = sds_info->rank;
    if (sds_info->dim_size[0] < sds_info->dim_size[rank1-1])
    {
      dim_size[0] = sds_info->dim_size[rank1-2];
      dim_size[1] = sds_info->dim_size[rank1-1];
    }
    else
    {
      dim_size[0] = sds_info->dim_size[0];
      dim_size[1] = sds_info->dim_size[1];
    }
  }
}

void write_metadata(int32 in_sd_id, int32 out_sd_id)
/*
 * copy the metadata from input HDF file to the output HDF file.
 */
{
  char attr_name[80];
  void *meta_attr_buf;
  int32 attr_cnt, attr_size;
  int32 nsds, nattr, attr_len;
  int32 attr_index, attr_type;                                                                                                                  
  fprintf(stdout, "\tWriting all global metadata to output file\n");
  if (SDfileinfo(in_sd_id, &nsds, &nattr) == FAIL)
    fprintf(stderr, "Cannot read file information for the input HDF file\n");
  else
  {
    for (attr_index=0; attr_index<nattr; attr_index++)
    {
      if (SDattrinfo(in_sd_id, attr_index, attr_name, &attr_type, &attr_cnt) == FAIL)
        fprintf(stdout, "Cannot read information for attribute " LONG_INT_FMT " in input HDF file\n", attr_index);
      else
      {
        attr_size = DFKNTsize(attr_type);
        attr_len = attr_size * attr_cnt + 1;
        if ((meta_attr_buf = (void *)calloc(attr_len, attr_size)) == NULL)
          fprintf(stdout, "Cannot allocate memory for metadata in write_metadata\n");
        else
        {
          if (SDreadattr(in_sd_id, attr_index, meta_attr_buf) == FAIL)
            fprintf(stdout, "Cannot read %s from input HDF file\n", attr_name);
          else
          {
            if (SDsetattr(out_sd_id, attr_name, attr_type, attr_cnt, meta_attr_buf) == FAIL)
              fprintf(stdout, "Cannot write metadata %s to output HDF file\n", attr_name);
          }
          free(meta_attr_buf);
        }
      }
    } /* for (attr_index=0; . . .      */
  }
}

void display_sds_info_of_file(char* filename)
/* Print all the SDS names in the input file */
{
  int32 sd_id, sds_id;
  int32 msds, nattr, rank, dt, dim_size[4];
  int k, isds;
  char dim_str[MAX_STR_LEN];
  char name[MAX_SDS_NAME_LEN], sds_name[MAX_SDS_NAME_LEN];
 
  msds = 0;
  if ((sd_id = SDstart(filename, DFACC_READ)) == FAIL)
    fprintf(stderr, "Cannot open the HDF file %s\n", filename);
  else
    {
      if (SDfileinfo(sd_id, &msds, &nattr) == FAIL)
	fprintf(stderr, "Error reading information for HDF file %s\n",filename);
      else
	{
	  fprintf(stdout, "Valid SDS names, dimension and data type in file: %s\n", filename);
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
  exit(EXIT_FAILURE);
}

int open_l2g_nobs_sds(sds_t *nobs_sds_info, sds_t *nadd_obs_sds_info, sds_t *sds_info)
{
  int k = 0;

  nobs_sds_info->sd_id = sds_info->sd_id;
  nobs_sds_info->sds_id = -1;
  sprintf(nobs_sds_info->name, "num_observations");
  k += get_sds_info((char *)NULL, nobs_sds_info);

  nadd_obs_sds_info->sd_id = sds_info->sd_id;
  nadd_obs_sds_info->sds_id = -1;
  sprintf(nadd_obs_sds_info->name, "nadd_obs_row");
  k += get_sds_info((char *)NULL, nadd_obs_sds_info);
  return k;

}
