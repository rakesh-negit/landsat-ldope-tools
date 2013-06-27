/****************************************************************************
!C

!File: meta.c

!Description:
  Contains routines for reading metadata.

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original March 1998. Version 1.0
  Modified June 1998. Version 1.1 (Updated to support batch processing. None
	of the routines abruptly exit, instead returns to calling routines
	with error messages.

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

#include "isoc.h"
#include "meta.h"
#include "hdf.h"
#include "mfhdf.h"
#include "qa_tool.h"
#include "str_op.h"

char *get_attr_metadata(char *hdf_fname, char *meta_gname)
/**********************************************************************************     
	Get core or archive metadata from the HDF file hdf_fname.
	On success the functions returns a pointer to the character string 
	containing the metadata else returns a NULL pointer on failure

!Input Parameters: 
  hdf_name     - input HDF file name
  meta_gname   - metadata name

!Output Parameters: 
  (returns)  Completion status: 
             On success return a pointer to the character string containing the metadata
             On failure, return NULL.
    
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END
***********************************************************************************/
{
  int status = 1;
  int32 sd_id;
  int32 attr_cnt;        /* count (number of characters) of the attribute */
  int32 attr_size;       /* attribute size */
  int32 attr_type;       /* attribute data type */
  int32 attr_index;      /* attribute index */
  char attr_name[20];    /* attribute name */
  char *meta_attr_buf = NULL;   /* data buffer for holding metadta information */

  if ((sd_id = SDstart(hdf_fname, DFACC_READ)) == FAIL)
  {
    fprintf(stderr, "Cannot open the input HDF file %s\n", hdf_fname);
    status = -1;
  }
  else
  {
    /* find metaname index in input HDF file */
    if ((attr_index = SDfindattr(sd_id, meta_gname)) == FAIL)
    {
      fprintf(stderr, "Cannot find %s in %s\n", meta_gname, hdf_fname);
      status = -1;
    }
    else
    {
      /* obtain the attribute index, attribute name, attribute data type and 
         count of attribute */
      if (SDattrinfo(sd_id, attr_index, attr_name, &attr_type, &attr_cnt) == FAIL)
      {
        fprintf(stderr, "Cannot read %s in %s\n", meta_gname, hdf_fname);
        status = -1;
      }
      else
      {
      /* determine the attrbute size */
        attr_size = DFKNTsize(attr_type) * attr_cnt + 1;
        if ((meta_attr_buf = (char *)calloc(attr_size, sizeof(char))) == NULL)
        {
          fprintf(stderr, "Cannot allocate memory for metadata in get_attr_metadata\n");
          status = -1;
        }
        else if (SDreadattr(sd_id, attr_index, (VOIDP)meta_attr_buf) == FAIL)
        {
          fprintf(stderr, "Cannot read %s from %s\n", meta_gname, hdf_fname);
	  free(meta_attr_buf);
          status = -1;
        }
      }
    }
    SDend(sd_id);
  }
  if (status == -1) return (char *)NULL;
  else return meta_attr_buf;
}

void get_all_metadata(char *meta_str, char **meta_name, char **meta_val, int *meta_cnt)
/***************************************************************************************
       
!Description:

	Read all individual metadata names into the character array meta_name and the 
	corresponding values into the character array meta_val from the single character 
	string meta_str containing all the core or archive metadata. meta_cnt contains
	the number of metadata in meta_str.

!Input Parameters: 
  meta_str     - input metadata structures
  
!Input/Output Parameters:
  meta_name    - individual metadata names
  meta_val     - meta data values
  meta_cnt     - number of metadata

!Output Parameters: (none)
    
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END

****************************************************************************************/
{
  int i;
  int cnt;                   /* count of number of metadata */
  int len, len_t;
  int pos22, pos33;
  int pos1, pos2, pos3;

  cnt = *meta_cnt;
  pos3 = 0;
  len = (int)strlen(meta_str);
  while (pos3 < len)
  {
    /*  retrieve OBJECT name  */
    pos1 = sd_strpos(meta_str, "OBJECT", pos3);
    if (pos1 != -1)
    {
      pos2 = sd_charpos(meta_str, '=', pos1) + 2;
      pos3 = sd_charpos(meta_str, '\n', pos2);
      i = pos3 + 1;
      while (meta_str[i] == ' ') i++;
      /*  NUM_VAL found */
      if (meta_str[i] == 'N')
      {
        sd_strmid(meta_str, pos2, pos3-pos2, meta_name[cnt]);
        pos3 = i;
        /* VALUE found */
        pos1 = sd_strpos(meta_str, "VALUE", pos3);
        pos2 = sd_charpos(meta_str, '=', pos1) + 2;
        if (meta_str[pos2] == '(')
          pos3 = sd_charpos(meta_str, ')', pos2) + 1;
        else
          pos3 = sd_charpos(meta_str, '\n', pos2);
        sd_strmid(meta_str, pos2, pos3-pos2, meta_val[cnt]);
	if (meta_str[pos3-1] == ')')
	  sd_rm_ln_in_str(meta_val[cnt]);
	  /* get ADDITIONALATTRIBUTENAME */
        if (strcmp(meta_name[cnt], "ADDITIONALATTRIBUTENAME") == 0)
        {
	  len_t = (int)strlen(meta_val[cnt]);
	  sd_strmid(meta_val[cnt], 1, len_t-2, meta_name[cnt]);
	  cnt++;
        }
        /* get PARAMETERVALUE */
        else if (strcmp(meta_name[cnt], "PARAMETERVALUE") == 0)
	  strcpy(meta_val[cnt-1], meta_val[cnt]);
        else cnt++;
      }
      /* CLASS */
      else if (meta_str[i] == 'C')
      {
        pos22 = sd_charpos(meta_str, '=', i);
        pos33 = sd_charpos(meta_str, '\n', pos22);
        i = pos33 + 1;
        while (meta_str[i] == ' ') i++;
        /* NUM_VAL */
        if (meta_str[i] == 'N')
        {
          sd_strmid(meta_str, pos2, pos3-pos2, meta_name[cnt]);
          pos1 = sd_strpos(meta_str, "VALUE", pos33);
          pos2 = sd_charpos(meta_str, '=', pos1) + 2;
          if (meta_str[pos2] == '(')
            pos3 = sd_charpos(meta_str, ')', pos2) + 1;
          else
            pos3 = sd_charpos(meta_str, '\n', pos2);
          sd_strmid(meta_str, pos2, pos3-pos2, meta_val[cnt]);
	  if (meta_str[pos3-1] == ')')
	    sd_rm_ln_in_str(meta_val[cnt]);
          if (strcmp(meta_name[cnt], "ADDITIONALATTRIBUTENAME") == 0)
          {
	    len_t = (int)strlen(meta_val[cnt]);
	    sd_strmid(meta_val[cnt], 1, len_t-2, meta_name[cnt]);
            cnt++;
          }
          else if (strcmp(meta_name[cnt], "PARAMETERVALUE") == 0)
            strcpy(meta_val[cnt-1], meta_val[cnt]);
          else cnt++;
	}
      }
    }
    else pos3 = len + 1;
  }
  *meta_cnt = cnt;
}

void get_sel_metadata(char *meta_str, char *meta_name, char **meta_val, int *meta_cnt, int case_ch)
/***************************************************************************************
       
!Description:

	Read all instances of metadata meta_name and read the metadata values into the 
	character array meta_val from the single character string meta_str containing 
	all the core or archive metadata. meta_cnt contains the number of instances of 
	metadata meta_name found in meta_str.

    if case_ch = 0, case sensitive. otherwise case insensitive for meta_name.
        
!Input Parameters: 
  meta_str     - input metadata structures
  case_ch      - case option
  
!Input/Output Parameters:
  meta_name    - individual metadata names
  meta_val     - meta data values
  meta_cnt     - number of metadata

!Output Parameters: (none)
    
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END

****************************************************************************************/
{
  int cnt, st;
  char tmp_name[MAX_STR_LEN];
  int pos1, pos2, pos3, pos4;

  cnt = *meta_cnt;
  pos1 = (case_ch == 0) ? sd_strpos(meta_str, meta_name, 0) : sd_strcasepos(meta_str, meta_name, 0); 
  while (pos1 != -1)
  {
    pos1--;
    if ((meta_str[pos1] == '\"') || (meta_str[pos1] == ' ')) pos1++;
    pos2 = sd_charpos(meta_str, '\n', pos1);
    if (meta_str[pos2-1] == '\"') pos2--; 
    sd_strmid(meta_str, pos1, pos2-pos1, tmp_name);
    st = (case_ch == 0) ? strcmp(meta_name, tmp_name) : strcasecmp(meta_name, tmp_name);
    if (st == 0)
    {
      if (cnt == *meta_cnt) strcpy(meta_name, tmp_name);
      pos3 = sd_strpos(meta_str, "END_OBJECT", pos2);
      pos3 = sd_charpos(meta_str, '=', pos3) + 2;
      pos4 = sd_charpos(meta_str, '\n', pos3);
      sd_strmid(meta_str, pos3, pos4-pos3, tmp_name);
      if (strcmp(tmp_name, "ADDITIONALATTRIBUTENAME") == 0)
      {
	pos2 = sd_strpos(meta_str, "PARAMETERVALUE", pos4);
	pos2 = sd_charpos(meta_str, '\n', pos2);
      }
      pos2 = sd_strpos(meta_str, "VALUE", pos2);
      pos2 = sd_charpos(meta_str, '=', pos2) + 2;
      if (meta_str[pos2] == '(')
        pos3 = sd_charpos(meta_str, ')', pos2) + 1;
      else
        pos3 = sd_charpos(meta_str, '\n', pos2);
      sd_strmid(meta_str, pos2, pos3-pos2, meta_val[cnt]);
      cnt++;
      pos1 = sd_strpos(meta_str, "END_OBJECT", pos3);
      pos1 = sd_charpos(meta_str, '\n', pos1);
      pos1 = sd_strpos(meta_str, meta_name, pos1);
    } 
    else
      pos1 = sd_strpos(meta_str, meta_name, pos2);
  }
  *meta_cnt = cnt;
}

void copy_metadata(int32 in_sd_id, int32 out_sd_id)
/********************************************************************************************
       
!Description:
        Copy all global metadata from input to output 

!Input Parameters: 
  in_sd_id      - input SD interface identifer.
  out_sd_id     - output SD interface identifer.

!Output Parameters: (none)
    
!Revision History: 
  Added error testing for SDfileinfo.      Aug 13, 2002   Y. Zhang  /SSAI
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

!END

********************************************************************************************/
{
  long attr_len;
  int32 nsds;                     /* number of data sets in the file */
  int32 nattr;                    /* number of global attributes in the file */
  int32 attr_cnt;                 /* count of number of attribute */
  int32 attr_index;               /* attribute index */
  int32 attr_type, attr_size;     /* attribute type and size */
  char attr_name[40];             /* attribute name */
  char *attr_buf;                 /* data buffer for attribute */

  fprintf(stdout, "\tWriting metadata to output\n");
  SDfileinfo(in_sd_id, &nsds, &nattr);
  for (attr_index=0; attr_index<nattr; attr_index++)
  {
     /* test to see if the file or dataset do indeed contain attributes */
     if(SDattrinfo(in_sd_id,attr_index,attr_name,&attr_type,&attr_cnt) == FAIL)
        fprintf(stdout,
                "Cannot read information for attribute " LONG_INT_FMT "\n",
                attr_index);
     else
     {
        attr_size = DFKNTsize(attr_type);
        attr_len = attr_size*attr_cnt;
        if ((attr_type == DFNT_CHAR) || (attr_type == DFNT_CHAR8))
           attr_len++;
        if ((attr_buf = (void *)calloc(attr_size, attr_len)) == NULL)
           fprintf(stdout,
                   "Cannot allocate memory for attr_buf in copy_metadata\n");
        else
        {
           if (SDreadattr(in_sd_id, attr_index, (VOIDP)attr_buf) == FAIL)
              fprintf(stdout, "Cannot read %s from input file\n", attr_name);
           else
           {
              /* copy the metadata from attr_buf to out_sd_id. */
              if (SDsetattr(out_sd_id, attr_name, attr_type, attr_cnt,
                            (VOIDP)attr_buf) == FAIL)
                 fprintf(stdout,
                       "Cannot write metadata %s to output file\n", attr_name);
           }
        free(attr_buf);
        }
     }
  } /* for (attr_index=0; . . .      */
  return;
}

void write_modss_metadata(int32 in_sd_id, int32 out_sd_id, char **meta_names, char **meta_vals)
/********************************************************************************************
    Copy metadata from input to output L2G subset. Update all neccessary metadata.
********************************************************************************************/
{
  int imeta;
  int32 attr_cnt, attr_index;
  int32 attr_size, attr_type;
  char attr_name[20], *meta_attr_buf, *new_mbuf;
  char *meta_gname[] = {"CoreMetadata.0", "ArchiveMetadata.0"};

  for (imeta=0; imeta<2; imeta++)
  {
    if ((attr_index = SDfindattr(in_sd_id, meta_gname[imeta])) == FAIL)
      fprintf(stderr, "Cannot find %s in input HDF file\n", meta_gname[imeta]);
    else
    {
      if (SDattrinfo(in_sd_id, attr_index, attr_name, &attr_type, &attr_cnt) == FAIL)
        fprintf(stderr, "Cannot read information for %s from input HDF file", meta_gname[imeta]);
      else
      {
        attr_size = DFKNTsize(attr_type) * attr_cnt + 1;
        if ((meta_attr_buf = (char *)calloc(attr_size, sizeof(char))) == NULL)
          fprintf(stderr, "Cannot allocate memory for metadata in write_modss_metadata()\n");
        else
        {
          if (SDreadattr(in_sd_id, attr_index, (VOIDP)meta_attr_buf) == FAIL)
            fprintf(stderr, "Cannot read %s from input HDF file \n", meta_gname[imeta]);
          else
          {
	    new_mbuf = update_modss_metadata(meta_attr_buf, attr_cnt, imeta, meta_names, meta_vals);
	    attr_cnt = (int)strlen(new_mbuf);
            if (SDsetattr(out_sd_id, meta_gname[imeta], attr_type, attr_cnt, (VOIDP)new_mbuf) == FAIL)
              fprintf(stderr, "Cannot write metadata %s to output file \n", meta_gname[imeta]);
            free(new_mbuf);
          }
          free(meta_attr_buf);
        }
      }
    }
  } /* for (imeta=0; . . . )      */
}                                      

char *update_modss_metadata(char *attr_buf, int org_len, int imeta, char **meta_names, char **meta_vals)
/****************************************************************************************
   Create the new metadata string by modifying the values of required metadata
****************************************************************************************/
{
  char end_c, *meta_str;
  int m1, m2, m_cnt;
  int i, j, k, ii, jj;
  int p1, p2, p3, len;
  int st_id[10], st_p[10], end_p[10];

  /* get index of core metadata and esdt metadata to be updated */

  if (imeta == 0) { m1 = 0; m2 = 6; m_cnt = 6; }
  else { m1 = 6; m2 = 10; m_cnt = 4; }

  /* Find the starting and ending positions of metadata values in the string */

  for (i=m1; i<m2; i++)
  {
    if ((p1 = sd_strpos(attr_buf, meta_names[i], 0)) != -1)
    {
      p2 = sd_strpos(attr_buf, "VALUE", p1);
      p3 = sd_charpos(attr_buf, '=', p2);
      p3 += 2;
      st_p[i-m1] = p3;
      end_c = (attr_buf[p3] == '(') ? ')' : '\n';
      while (attr_buf[p3] != end_c) p3++;
      end_p[i-m1] = p3;
    }                                    
    else st_p[i-m1] = end_p[i-m1] = -1;
  }

  /* sort the metadata based on their position in the string */

  for (i=0; i<m_cnt; i++)
    st_id[i] = i;
  for (i=0; i<m_cnt; i++)
  {
    ii = st_id[i];
    for (j=i+1; j<m_cnt; j++)
    {
      jj = st_id[j];
      if (st_p[ii] > st_p[jj])
        ii = st_id[j];
    }
    if (ii != st_id[i])
    {
      k = st_id[i];
      st_id[i] = ii;
      st_id[ii] = k;
    }
  }

  /* obtain the new length */

  for (i=0, k=0, len=0; i<m_cnt; i++)
  {
    ii = st_id[i];
    if (i == 0) p1 = 0;
    else
    {
      jj = st_id[i-1];
      p1 = end_p[jj];
    }
    p2 = st_p[ii];
    len = len + p2 - p1 + 1;

    ii += m1;
    len += (int)strlen(meta_vals[ii]);
  }                                  
  ii = st_id[i-1];
  len = len + (org_len - end_p[ii]) + 1;

  /* update metadata and shift the string appropriately */ 

  if ((meta_str = (char *)calloc(len, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for meta_str in update_modss_metadata()\n");
  else
  {
    for (i=0, k=0; i<m_cnt; i++)
    {
      ii = st_id[i];
      if (i == 0) p1 = 0;
      else
      {
        jj = st_id[i-1];
        p1 = end_p[jj];
      }
      p2 = st_p[ii];
      for (j=p1; j<=p2; j++, k++)
        meta_str[k] = attr_buf[j];
      ii += m1;
      len = (int)strlen(meta_vals[ii]);
      for (j=0; j<len; j++, k++)
        meta_str[k] = meta_vals[ii][j];
    }
    jj = st_id[i-1];
    p1 = end_p[jj]; p2 = org_len;
    for (j=p1; j<=p2; j++, k++)
      meta_str[k] = attr_buf[j];
  }                                                    
  return meta_str;
}

