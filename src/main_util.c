/****************************************************************************
!C

!File: main_util.c

!Description:
  Contains routines for follwoing tasks
    -get the product short name
    -get the pointer filename
    -locate a given modis file
    -read and parse the input pointer to get input filenames
    -routines for parsing the argument name and values
    -check the validity of bit string for a given sds
    -parse the bit string to get bit number values
    -rm pathname in a filename
    -routine used in command line parsing
    -get qa tool environment variable
    -read the time for a granule or tile number for a tile
    -read a line from the file

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original January 1998. Version 1.0
  Modified March 1998. (Changes to use metadata to find informations about the 
			product id than the very product filename)
  Modified April 1998. (Added routines used in command line parsing)
  Modified June 1998. Version 1.1 (Modifications to all routines to support batch 
		processing. None of these routines abrupty exit. Instead they all 
		return to the calling routine with appropriate error messages.
  Modified March 2000. (Fix to find_file(), get_pnt_fname() to accomodate -ptr
			option in read_l2g)
  Modified September 2001. (Fixed bugs while importing to Linux)

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

#include "mfhdf.h"

#include "str_op.h"
#include "alloc_mem.h"
#include "qa_tool.h"
#include "meta.h"
#include "main_util.h"

int get_prod_sname(char *hdf_fname, char *prod_sname);
int find_file(char **fpath, char *fname);
void get_qa_tool_env(char *env_var, char **env_val);

int get_prod_sname(char *hdf_fname, char *prod_sname)
/******************************************************************************
	Read the product SHORTNAME of the hdf files
	On success the function returns 1 and prod_sname contains the SHORTNAME
	else function returns -1 and prod_sname is not meaningfull
*******************************************************************************/
{
  int len = 0;
  int status = 0;
  int meta_cnt = 0;
  char *cmeta_str = NULL, **meta_val = NULL;
  char meta_name[MAX_META_NAME_LEN];

  status = -1;
  meta_cnt = 0;
  if ((meta_val = (char **)Calloc2D(5, MAX_META_VAL_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for meta_val in get_prod_sname\n");
  else if ((cmeta_str = get_attr_metadata(hdf_fname, "CoreMetadata.0")) != NULL)
  {
    strcpy(meta_name, "SHORTNAME");
    get_sel_metadata(cmeta_str, meta_name, meta_val, &meta_cnt, 0);
    if (meta_cnt == 0)
      fprintf(stderr, "No metadata SHORTNAME found in file %s", hdf_fname);
    else
    {
      status = 1;
      len = (int)strlen(meta_val[0]);
      sd_strmid(meta_val[0], 1, len-2, prod_sname);
    }
    free(cmeta_str);
    Free2D((void **)meta_val);
  }
  return status;
}

int get_pnt_fname(char *hdf_fname, char **modis_fpath, char *ptr_fname, char *pnt_fname)
/******************************************************************************
	Gets pointer filename used in the production of product file hdf_fname
	modis_fpath contains a list of pathnames to search for a modis product
	On success returns 1 and pnt_fname contains the pointer filename
	else returns -1 and pnt_fname is not meaningful
*******************************************************************************/
{
  FILE *fp = NULL;
  int st = 0, status = 0;
  int pos1 = 0, pos2 = 0;
  int len = 0, meta_cnt = 0;
  char *cmeta_str = NULL, **meta_val = NULL;
  char fname[MAX_PATH_LENGTH];
  char tmp_fname[MAX_PATH_LENGTH];
  char meta_name[MAX_META_NAME_LEN];
  char pnt_sname[10], prod_sname[10];

  fprintf(stdout, "%s\n", ptr_fname);
  if (ptr_fname[0] != '\0') { 
    strcpy(pnt_fname, ptr_fname);
    return 1;
  }
  status = -1;
  if ((meta_val = (char **)Calloc2D(5, MAX_META_VAL_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for meta_val in get_pnt_fname\n");
  else if ((cmeta_str = get_attr_metadata(hdf_fname, "CoreMetadata.0")) != NULL)
  {
    meta_cnt = 0; 
    strcpy(meta_name, "INPUTPOINTER");
    get_sel_metadata(cmeta_str, meta_name, meta_val, &meta_cnt, 0);
    if (meta_cnt == 0)
      fprintf(stderr, "Metadata INPUTPOINTER not found in file %s", hdf_fname);
    else
    {
      pos1 = sd_charpos(meta_val[0], '"', 0);
      while ((pos1 != -1) && (status == -1))
      {
	st = 0; pos1++;
        pos2 = sd_charpos(meta_val[0], '"', pos1);
        sd_strmid(meta_val[0], pos1, pos2-pos1, fname);
        sd_strtrim(fname);
        if (strstr(ptr_fname, fname) != NULL)
        {
	  st = 1;
	  strcpy(fname, ptr_fname);
        }
	else
	{
	  sprintf(tmp_fname, "%s%s%s", ptr_fname, CSDI_PATH_SEP_STR, fname);
	  if ((fp = fopen(tmp_fname, "r")) != NULL)
	  {
	    st = 1;
	    strcpy(fname, tmp_fname);
	    fclose(fp);
	  }
	}
	if (st == 0)
          st = find_file(modis_fpath, fname);
	if (st > 0)
        {
          if (get_prod_sname(fname, prod_sname) != -1)
          {
            len = (int)strlen(prod_sname);
	    if (len > 5)
            {
              sd_strmid(prod_sname, 0, 5, pnt_sname); 
              if ((strcmp(pnt_sname, "MODPT") == 0) || 
	      	   strcmp(pnt_sname, "MYDPT") == 0)
		   {
		  	 status = 1;
		   }	
	    }
          }
        }
        pos1 = sd_charpos(meta_val[0], '\"', pos2+1);
      }
      if (status == -1)
        fprintf(stderr, "No pointer files found in metadata INPUTPOINTER of %s\n", hdf_fname);
      else
        strcpy(pnt_fname, fname);
    }
    free(cmeta_str);
    Free2D((void **)meta_val);
  }
  return status;
}

int find_file(char **fpath, char *fname)
/******************************************************************************
	Returns 1 if file fname is found in one of the paths listed in fpath
	else return -1
*******************************************************************************/
{
  FILE *fp = NULL;
  int i = 0;
  char tmp_fname[MAX_PATH_LENGTH];

  for (i=0; i<MAX_NUM_PATH; i++)
    if (strlen(fpath[i]) > 1)
    {
      sprintf(tmp_fname, "%s%s", fpath[i], fname);
      if ((fp = fopen(tmp_fname, "r")) != NULL) break;
    }
  if (i < MAX_NUM_PATH) 
  {
    fclose(fp);
    strcpy(fname, tmp_fname);
    return 1;
  }
  else
  {
    if ((fp = fopen(fname, "r")) != NULL) 
    {
      fclose(fp);
      return 1;
    } else return -1;
  }
}

int get_input_files(char *fname, char *iptr_id, char **iptr_fnames)
/*****************************************************************************
    Returns number of input files used in the production of fname that matches
    the string iptr_id. iptr_fname contains the resulting filenames. if iptr_id 
    == 'all', all filenames are read from metadata INPUTPOINTER. else matching 
    filenames (match to iptr_id) are returned. else filenames from INPUTPOINTER
    with matching SHORTNAME (match to iptr_id) are returned.
*****************************************************************************/ 
{
  int i, is_l2g;
  int pos1, pos2;
  int m_cnt, tot_fcnt;
  int iptr_fcnt, iptr_id_len;
  int num_op_gran, num_ip_gran;
  char **tmp_fnames, **modis_fpath;
  char *cmeta_str, *ameta_str, **meta_val;
  char tmp[MAX_STR_LEN], meta_name[MAX_META_NAME_LEN];

  if ((meta_val = (char **)Calloc2D(2, MAX_META_VAL_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for meta_val in get_input_files \n");
  if ((tmp_fnames = (char **)Calloc2D(MAX_NUM_FILE, MAX_STR_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for tmp_fnames in get_input_files \n");
  if ((modis_fpath = (char **)Calloc2D(MAX_NUM_PATH, MAX_PATH_LENGTH, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for modis_fpath in get_input_files \n");
  cmeta_str = get_attr_metadata(fname, "CoreMetadata.0");
  ameta_str = get_attr_metadata(fname, "ArchiveMetadata.0");
  if ((meta_val == NULL) || (tmp_fnames == NULL) || (modis_fpath == NULL) || (ameta_str == NULL) || (cmeta_str == NULL))
  {
    Free2D((void **)meta_val);
    Free2D((void **)tmp_fnames);
    Free2D((void **)modis_fpath);
    if (ameta_str != NULL) free(ameta_str);
    if (cmeta_str != NULL) free(cmeta_str);
    return -1;
  }

  is_l2g = 0;
  tot_fcnt = 0;
  iptr_fcnt = 0;
  num_op_gran = -1;
  num_ip_gran = -1;

  /* if l2g file make a note of number of granules used */
  if (sd_strpos(ameta_str, "NUMBEROFOVERLAPGRANULES", 0) != -1)
  {
    is_l2g = 1;
    m_cnt = 0;
    strcpy(meta_name, "NUMBEROFOVERLAPGRANULES");
    get_sel_metadata(ameta_str, meta_name, meta_val, &m_cnt, 0);
    if (m_cnt != 0) num_op_gran = atoi(meta_val[0]);
    m_cnt = 0;
    strcpy(meta_name, "NUMBEROFINPUTGRANULES");
    get_sel_metadata(ameta_str, meta_name, meta_val, &m_cnt, 0);
    if (m_cnt != 0) num_ip_gran = atoi(meta_val[0]);
  }
  /* read the input pointer */
  m_cnt = 0;
  strcpy(meta_name, "NUMBEROFINPUTGRANULES");
  get_sel_metadata(cmeta_str, meta_name, meta_val, &m_cnt, 0);
  free(ameta_str);
  free(cmeta_str);

  if (m_cnt == 0)
    fprintf(stderr, "Metadata INPUTPOINTER not found in %s\n", fname);
  else
  {
    get_qa_tool_env(MODIS_ENV, modis_fpath);
    if (strcmp(iptr_id, "all") == 0)
    {
      fprintf(stderr, "\nReading all the filenames in the INPUTPOINTER\n");
      pos1 = sd_charpos(meta_val[0], '"', 0);
      while (pos1 != -1)
      {
        pos2 = sd_charpos(meta_val[0], '"', pos1+1);
        sd_strmid(meta_val[0], pos1+1, pos2-pos1-1, tmp);
        sd_strtrim(tmp);
        if (find_file(modis_fpath, tmp) == 1)
          strcpy(iptr_fnames[iptr_fcnt++], tmp);
        pos1 = sd_charpos(meta_val[0], '"', pos2+1);
      }
    }
    else
    {
      fprintf(stderr, "\nComparing %s to filenames in the INPUTPOINTER\n", iptr_id);
      iptr_id_len = (int)strlen(iptr_id);
      pos1 = sd_charpos(meta_val[0], '"', 0);
      while (pos1 != -1)
      {
        pos2 = sd_charpos(meta_val[0], '"', pos1+1);
        sd_strmid(meta_val[0], pos1+1, pos2-pos1-1, tmp);
        sd_strtrim(tmp);
        strcpy(tmp_fnames[tot_fcnt++], tmp);
        if (strncmp(iptr_id, tmp, iptr_id_len) == 0)
          if (find_file(modis_fpath, tmp) == 1)
            strcpy(iptr_fnames[iptr_fcnt++], tmp);
        pos1 = sd_charpos(meta_val[0], '"', pos2+1);
      }
      if (iptr_fcnt == 0)
      {
        fprintf(stderr, "Comparing SHORTNAME in files from INPUTPOINTER to %s\n", iptr_id);
        for (i=0; i<tot_fcnt; i++)
        {
          strcpy(tmp, tmp_fnames[i]);
          if (find_file(modis_fpath, tmp) != 0)
          {
            if ((cmeta_str = get_attr_metadata(tmp, "CoreMetadata.0")) != NULL)
	    {
              m_cnt = 0;
	      strcpy(meta_name, "SHORTNAME");
              get_sel_metadata(cmeta_str, meta_name, meta_val, &m_cnt, 0);
              if ((m_cnt != 0) && (strcmp(meta_val[0], iptr_id) == 0))
                strcpy(iptr_fnames[iptr_fcnt++], tmp);
              free(cmeta_str);
	    }
          }
        }
      }
    }
    if (iptr_fcnt == 0)
      fprintf(stderr, "No matching input files found in the INPUTPOINTER\n");
    else if (is_l2g == 1) 
    {
      if (strcmp(iptr_id, "all") == 0)
      {
	strcpy(iptr_fnames[num_op_gran], iptr_fnames[num_ip_gran]);
	iptr_fcnt = num_op_gran+1;
      }
      else if (iptr_fcnt > num_op_gran)
        iptr_fcnt = num_op_gran; /* if l2g limit to number of valid granules*/
    }
  }
  Free2D((void **)meta_val);
  Free2D((void **)tmp_fnames);
  Free2D((void **)modis_fpath);

  return iptr_fcnt;
}


int is_arg_id(char *arg_str, char *arg_id)
/******************************************************************************
	return 0 if arg_id is initial part of the arg_str else return -1
******************************************************************************/
{
  int p1 = 0;

  p1 = sd_charpos(arg_str, '=', 0);
  if (p1 == -1) return -1;
  if (strncmp(arg_str, arg_id, p1) == 0)
   return 0;
  else return -1;
}

void get_arg_val_arr(char *arg_str, char **arg_val, int *arg_cnt)
/******************************************************************************
	read the argument values from arg_str and return the number of argument 
	values as arg_cnt and the argument values in the array arg_val.
******************************************************************************/
{
  int p1 = 0, p2 = 0, len = 0;

  p1 =  0;
  p2 = sd_charpos(arg_str, '=', p1);
  if (p2 == -1) return;
  p1 = p2 + 1;
  p2 = sd_charpos(arg_str, ',', p1);
  while ((p2 != -1) && (*arg_cnt < MAX_NUM_PARAM))
  {
    sd_strmid(arg_str, p1, p2-p1, arg_val[*arg_cnt]);
    ++*arg_cnt;
    p1 = p2 + 1;
    p2 = sd_charpos(arg_str, ',', p1);
  }
  if (*arg_cnt >= MAX_NUM_PARAM)
  {
    fprintf(stderr, "Too many parameters in option %s\n", arg_str);
    fprintf(stderr, "Considering only %d number of parameter values\n", MAX_NUM_PARAM);
  }
  else
  {
    len = (int)strlen(arg_str);
    if (p1 < len)
    {
      sd_strmid(arg_str, p1, len-p1, arg_val[*arg_cnt]);
      ++*arg_cnt;
    }
  }
} 

void get_arg_val(char *arg_str, char *arg_val)
/******************************************************************************
	Read single argument val arg_val from arg_str
******************************************************************************/
{
  int p1 = 0, len = 0;

  p1 = sd_charpos(arg_str, '=', 0);
  if (p1 == -1) return;
  p1 = p1 + 1;
  len = (int)strlen(arg_str);
  if (len == p1) return;
  else
    sd_strmid(arg_str, p1, len-p1, arg_val);
}

int check_bit_str(char *fname, char *sname, char *bn_str)
/******************************************************************************
	Returns 1 if the bit numbers in the bit number string is valid for
	the named sds sname in the hdf file fname. Else returns -1
******************************************************************************/
{
  int len, len2;
  int pos1, pos2, pos3;
  int val, val1, val2, max_val;
  int32 sd_id, sds_id, sds_index;
  int32 rank, dim_sz[3], dt, nattr;
  char num[10], t1[5], t2[5];

  if ((sd_id = SDstart(fname, DFACC_READ)) == FAIL)
  {
    fprintf(stderr, "Cannot open the file: %s\n", fname);
    return -1;
  }
  if ((sds_index = SDnametoindex(sd_id, sname)) == FAIL)
  {
    fprintf(stderr, "Cannot find the %s in file %s\n", sname, fname);
    SDend(sd_id);
    return -1;
  }
  if ((sds_id = SDselect(sd_id, sds_index)) == FAIL)
  {
     fprintf(stderr, "Cannot select the sds %s in file %s\n", sname, fname);
     SDend(sd_id);
     return -1;
  }
  if (SDgetinfo(sds_id, sname, &rank, dim_sz, &dt, &nattr) == FAIL)
  {
     fprintf(stderr, "Cannot select the sds %s in file %s\n", sname, fname);
     SDendaccess(sds_id);
     SDend(sd_id);
     return -1;
  }
  SDendaccess(sds_id);
  SDend(sd_id);

  max_val = DFKNTsize(dt)*8-1;
  sd_strtrim(bn_str);
  len = (int)strlen(bn_str);
  pos1 = 0;
  while (pos1 < len)
  {
    pos2 = sd_charpos(bn_str, ',', pos1);
    if (pos2 == -1) pos2 = len;
    sd_strmid(bn_str, pos1, pos2-pos1, num);
    pos3 =  sd_charpos(num, '-', 0);
    if (pos3 == -1)
    {
      sd_strtrim(num);
      val = atoi(num);
      if (val > max_val)
      {
        fprintf(stderr, "Bit number input excceds the SDS size\n");
        return -1;
      }
    }
    else
    {
      len2 = (int)strlen(num);
      sd_strmid(num, 0, pos3, t1);
      sd_strmid(num, pos3+1, len2-pos3, t2);
      sd_strtrim(t1);
      sd_strtrim(t2);
      val1 = atoi(t1);
      val2 = atoi(t2);
      if ((val1 > max_val) || (val2 > max_val))
      {
        fprintf(stderr, "Bit number input excceds the SDS size\n");
        return -1;
      }
    }
    pos1 = pos2+1;
  }
  return 1;
}

void get_bit_num_arr(char *bn_str, int *cnt, int *bn_arr)
/*************************************************************************************
	Read the bit numbers from the bit number string bn_str into an ineger 
	array bn_arr and number of bits into cnt
*************************************************************************************/
{
  int i = 0, j = 0, tmp = 0;
  int len = 0, len2 = 0;
  int val = 0, val1 = 0, val2 = 0;
  int pos1 = 0, pos2 = 0, pos3 = 0, bcnt = 0;
  char num[10], t1[5], t2[5];

  bcnt = 0;
  pos1 = 0;
  sd_strtrim(bn_str);
  len = (int)strlen(bn_str);
  while (pos1 < len)
  {
    pos2 = sd_charpos(bn_str, ',', pos1);
    if (pos2 == -1) pos2 = len;
    sd_strmid(bn_str, pos1, pos2-pos1, num);
    pos3 =  sd_charpos(num, '-', 0);
    if (pos3 == -1)
    {
      sd_strtrim(num);
      val = atoi(num);
      bn_arr[bcnt++] = val;
    }
    else
    {
      len2 = (int)strlen(num);
      sd_strmid(num, 0, pos3, t1);
      sd_strmid(num, pos3+1, len2-pos3, t2);
      sd_strtrim(t1);
      sd_strtrim(t2);
      val1 = atoi(t1);
      val2 = atoi(t2);
      for (i=val1; i<=val2; i++)
        bn_arr[bcnt++] = i;
    }
    pos1 = pos2+1;
  }
  *cnt = bcnt;

  for (i=0; i<bcnt-1; i++)
  for (j=bcnt-1; j>i; j--)
    if (bn_arr[j] < bn_arr[j-1])
    {
      tmp = bn_arr[j];
      bn_arr[j] = bn_arr[j-1];
      bn_arr[j-1] = tmp;
    }
}

void rm_path(char *fname)
/*****************************************************************************
	Remove the path from the file name fname and return the filename in
	fname itself
****************************************************************************/
{
  int i = 0, j = 0, k = 0;
  int len = 0;

  len = (int)strlen(fname);
  for (i=len; i>0; i--)
    if (fname[i] == '/' || fname[i] == '\\') break;
  if (i != 0)
  {
    for (j = 0, k = i+1; k < len; j++, k++)
      fname[j] = fname[k];
    fname[j] = '\0';
  }
}

int parse_stdin(int *argc, char **argv)
/*****************************************************************************
 Modified on 05/19/1998 by S. Devadiga
        Reads the argument from file input as stdin. Modified to follow exact
        unix syntax (space, line breaks, quotes and combination of quotes etc).
        Assumes that the stdin conatins only arguments.
******************************************************************************/
{
  int cnt = 0;
  int len = 0;
  int i = 0, j = 0;
  int st_quote = 0;
  char quote = 0;
  char tmp[MAX_STR_LEN];

  cnt = 1;
  j = 0;
  st_quote = 0;
  while (fscanf(stdin, "%s", tmp) == 1)
  {
    len = (int)strlen(tmp);
    if (st_quote == 1) argv[cnt][j++] = ' ';
    for (i=0; i<len; i++)
    {
      if ((tmp[i] != '"') && (tmp[i] != '\'')) argv[cnt][j++] = tmp[i];
      else
      {
        if (st_quote == 0)
        {
          st_quote = 1;
          quote = tmp[i];
          if (j > 0)
          {
            argv[cnt][j] = '\0';
            j = 0;
            cnt++;
          }
        }
        else if (tmp[i] == quote)
        {
          st_quote = 0;
          if (j > 0)
          {
            argv[cnt][j] = '\0';
            j = 0;
            cnt++;
          }
        }
        else argv[cnt][j++] = tmp[i];
      }
    }
    if ((st_quote != 1) && (j > 0))
    {
      argv[cnt][j] = '\0';
      j = 0;
      cnt++;
    }
    if (cnt == MAX_NUM_PARAM) break;
  }
  if (cnt == MAX_NUM_PARAM)
    fprintf(stderr, "Too many input arguments. Reading only first %d arguments . . .\n", cnt);
  *argc = cnt;
  if (st_quote == 1) return -1;
  else return 1;
}

void get_qa_tool_env(char *env_var, char **env_val)
/******************************************************************************
	Get the value of the environment variable env_var and return it in
	env_val. The environment variable values follow unix syntax, except
    windows, where it follows windows syntax.
******************************************************************************/
{
  int cnt = 0;
  int len = 0;
  int p1 = 0, p2 = 0;
  char *tmp_val = NULL;
  char tmp_path[MAX_PATH_LENGTH];

  tmp_val = getenv(env_var);
  if (tmp_val != NULL)
  {
    cnt = 0;
    p1 = 0;
    while ((p2 = sd_charpos(tmp_val, CSDI_PATHS_SEP_CHAR, p1)) != -1)
    {
      sd_strmid(tmp_val, p1, p2-p1, tmp_path);
      sd_strtrim(tmp_path);
      sprintf(env_val[cnt++], "%s%s", tmp_path, CSDI_PATH_SEP_STR);
      p1 = p2 + 1;
    }
    len = (int)strlen(tmp_val);
    sd_strmid(tmp_val, p1, len-p1, tmp_path);
    sd_strtrim(tmp_path);
    sprintf(env_val[cnt++], "%s%s",  tmp_path, CSDI_PATH_SEP_STR);
  }
}

int get_day_time_tile_info(char *gran_id, char *day, char *time_tile)
{
  int len;
  int p_level;
  int i, j, k;
  char fn[80];
  char tmp[10];

  len = (int)strlen(gran_id);
  sd_strmid(gran_id, 1, len-2, fn);

  i = 0;
  while (fn[i] != '.') i++;
  for (j=i+6, k=0; j<i+9; j++, k++)
    tmp[k] = fn[j];
  tmp[k] = '\0';
  strcpy(day, tmp);

  k = 0;
  i = j+1;
  p_level = 2;
  while (fn[i] != '.')
  {
    if ((fn[i] != 'v') && (fn[i] != 'h'))
      tmp[k++] = fn[i];
    else
      p_level = 3; 
    i++;
  }                                     
  tmp[k] = '\0';
  strcpy(time_tile, tmp);
  return p_level;
}               

int get_line(FILE *fp, char *s)
{
  int i, c;

  i = 0; 
  while ((c = fgetc(fp)) != EOF)
  {
    s[i++] = (char)c;
    if (c == '\n') break;
  }
  if (i != 0) s[i++] = '\0';
  return i;
}

int get_numbers(char **num_str, int cnt, int *sel_num, int max, char *msg)
{
  int i, p1, len;
  int bn, bn1, bn2, bn_cnt;
  char tmp_str[10];

  if (cnt == 0)
    for (bn=0, bn_cnt=max; bn<max; bn++) sel_num[bn] = 1;
  else
  {
    for (i=0, bn_cnt=0; i<cnt; i++)
      if ((p1 = sd_charpos(num_str[i], '-', 0)) != -1)
      {
        sd_strmid(num_str[i], 0, p1, tmp_str);
        bn1 = (int)atoi(tmp_str);
        len = (int)strlen(num_str[i]);
        sd_strmid(num_str[i], p1+1, len-p1-1, tmp_str);
        bn2 = (int)atoi(tmp_str);
        for (bn=bn1; bn<=bn2; bn++)
          if ((bn <= 0) || (bn > max))
            fprintf(stderr, "Ignoring invalid %s number %d\n", msg, bn);
          else
          {
            sel_num[bn-1] = 1;
            bn_cnt++;
          }
      }
      else
      {
        bn = (int)atoi(num_str[i]);
        if ((bn <= 0) || (bn > max))
          fprintf(stderr, "Ignoring invalid %s number %d\n", msg, bn);
        else
        {
          sel_num[bn-1] = 1;
          bn_cnt++;
        }
      }
  }
  return bn_cnt;
}

/* sort the input values(integer) to ascending order. Total input  
 * values is n */ 
void sort_values(int *values, int n)
{
  int i = 0, j = 0, tmp = 0;
 
  for (i=n-1; i>0; i--)
    for (j=i-1; j>=0; j--)
      if (values[i] < values[j])
      {
	tmp = values[i];
        values[i] = values[j];
        values[j] = tmp;
      }
}

/* sort the input values(float) to ascending order. Total input 
 * values is n */
void sort_fvalues(float *values, int n)
{
  int i = 0, j = 0;
  float tmp = 0.;
  
  for (i=n-1; i>0; i--)
    for (j=i-1; j>=0; j--)
      if (values[i] < values[j])
      {
	tmp = values[i];
        values[i] = values[j];
        values[j] = tmp;
      }
}

/* read in the color table */
int read_clr_table(char *in_fname, int *clr_table)
{
  FILE *fp = NULL;
  int st = -1;
  int i = 0, id = 0, r = 0, g = 0, b = 0;

  if ((fp = fopen(in_fname, "r")) == NULL)
    fprintf(stderr, "Cannot open input file: %s\n", in_fname);
  else
  {
    i = id = 0;
    while ((fscanf(fp, "%d %d %d", &r, &g, &b) == 3) && (i < 256))
    {
      clr_table[id++] = r;
      clr_table[id++] = g;
      clr_table[id++] = b;
      ++i;
    }
    fclose(fp);
    fprintf(stderr, "Number of colors read from color table: %d\n", i);
    st = i;
  }
  return st;
}                             

void get_esdt_tileid(char *fname, char *esdt, char *tile_id, char *jday,
                     char *ver_id)
{
  int i = 0, j = 0, k = 0, len = 0;

  esdt[0] = tile_id[0] = '\0';
  jday[0] = ver_id[0] = '\0';
  len = (int)strlen(fname);
  for (k=len-1; k>=0; k--)
    if (fname[k] == '/') break;
  k++;
  for (j=k, i=0; j<len; j++, i++)
  {
    if (fname[j] == '.') { esdt[i] = '\0'; break; }
    else esdt[i] = fname[j];
  }
  k = j+2;
  for (j=k, i=0; j<len; j++, i++)
  {
    if (fname[j] == '.') { jday[i] = '\0'; break; }
    else jday[i] = fname[j];
  }
  k = j+1;
  for (j=k, i=0; j<len; j++, i++)
  {
    if (fname[j] == '.') { tile_id[i] = '\0'; break; }
    else tile_id[i] = fname[j];
  }
  k = j+1;
  for (j=k, i=0; j<len; j++, i++)
  {
    if (fname[j] == '.') { ver_id[i] = '\0'; break; }
    else ver_id[i] = fname[j];
  }
  if ((k = sd_charpos(esdt, '_', 0)) != -1)
    esdt[k] = '\0';
}

