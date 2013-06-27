/***************************************************************************e
!C

!File: mask_sds_lib.c

!Description:
  This file contains the library routines for creating mask SDS

!Input Parameters: (none)

!Output Parameters: (none)

!Revision History:

    Version 1.0    September, 2003
    
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
#include "meta.h"
#include "mask_sds_lib.h"

static unsigned int BIT[] = {
  0x1,        0x2,        0x4,        0x8,
  0x10,       0x20,       0x40,       0x80,
  0x100,      0x200,      0x400,      0x800,
  0x1000,     0x2000,     0x4000,     0x8000,
  0x10000,    0x20000,    0x40000,    0x80000,
  0x100000,   0x200000,   0x400000,   0x800000,
  0x1000000,  0x2000000,  0x4000000,  0x8000000,
  0x10000000, 0x20000000, 0x40000000, 0x80000000
};

int get_mask_string(char *m_str, char **arg_mask_str, int *val_opt, int *l2g_st)

/*
!C******************************************************************************
    
!Function: get_mask_string
       
!Description:
    
  Retrieve the mask string.
    
!Input Parameters: 
  m_str           String contains the user specified masking logic.
  arg_mask_str    String contains the information of input masking filename, masking 
                  SDS name, bit descriptor, equality operatior(<. > , <=, >=, !=, ==) 
		  and binary bit value string to match, and logical operator (AND or OR)
  val_opt         Flag of if the equality opertor is presented in the mask. 
  l2g_st          Flag of whether the input mask file is in L2G file format.

    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*******************************************************************************/

{
  int p1, p2, p3, i, k, cnt, len;
  char m_op[10], m_bit_str[50], p_bit_str[50], **mask_str;
  char m_fname[MAX_PATH_LENGTH], m_sname[MAX_SDS_NAME_LEN];
  char p_sname[MAX_SDS_NAME_LEN], p_fname[MAX_PATH_LENGTH];

  if ((mask_str = (char **)Calloc2D(MAX_NUM_OP, MAX_STR_LEN, sizeof(char))) == NULL)
    fprintf(stderr, "Cannot allocate memory for mask_str in mask_sds_lib: get_mask_string() \n");

  cnt = 0;
  p1 = 0;
  p2 = sd_strpos(m_str, "OR", p1);
  p3 = sd_strpos(m_str, "AND", p1);

  while ((p2 != -1) || (p3 != -1))
  {
    if (p2 == -1) p2 = p3;
    else if ((p3 != -1) && (p3 < p2)) p2 = p3;
    p3 = sd_charpos(m_str, ',', p2);
    sd_strmid(m_str, p1, p3-p1, mask_str[cnt]);
    cnt++;
    p1 = p3 + 1;
    p2 = sd_strpos(m_str, "OR", p1);
    p3 = sd_strpos(m_str, "AND", p1);
  }
  len = (int)strlen(m_str);
  sd_strmid(m_str, p1, len-p1, mask_str[cnt]);
  cnt++;
  for (i=0, k=0; i<cnt; i++)
  {
    p1 = p2 = p3 = -1;
    m_fname[0] = m_sname[0] = m_bit_str[0] = m_op[0] = '\0';
    len = (int)strlen(mask_str[i]);
    if ((p1 = sd_charpos(mask_str[i], ',', 0)) != -1)
    {
      sd_strmid(mask_str[i], 0, p1, m_fname);
      p1++;
      p2 = sd_charpos(mask_str[i], ',', p1);
      if (p2 != -1)
      {
        sd_strmid(mask_str[i], p1, p2-p1, m_sname);
	p2++;
        if ((p3 = sd_strpos(mask_str[i], "OR", p2)) != -1)
	  strcpy(m_op, "OR");
	else
	{
          if ((p3 = sd_strpos(mask_str[i], "AND", p2)) != -1)
	    strcpy(m_op, "AND");
	  else strcpy(m_op, "NONE");
	}
        if (p3 != -1)
	  sd_strmid(mask_str[i], p2, p3-p2-1, m_bit_str);
	else 
	  sd_strmid(mask_str[i], p2, len-p2, m_bit_str);
      }
    }
    if ((m_fname[0] == '*') && (k > 0)) strcpy(m_fname, p_fname);
    if ((m_sname[0] == '*') && (k > 0)) strcpy(m_sname, p_sname);
    if ((m_bit_str[0] == '*') && (k > 0)) strcpy(m_bit_str, p_bit_str);
    if (check_fsds_bit_str_val(m_fname, m_sname, m_bit_str, &val_opt[k], &l2g_st[k]) == -1)
      fprintf(stderr,"Ignoring input masking option %s\n", mask_str[i]);
    else
    {
      strcpy(p_fname, m_fname);
      strcpy(p_sname, m_sname);
      strcpy(p_bit_str, m_bit_str);
      sprintf(arg_mask_str[1+3*k], "%s:%s", m_fname, m_sname);  /* create the mask_string for use of 
                                                                   get_file_sds_names(). */
      strcpy(arg_mask_str[2+3*k], m_bit_str);
      strcpy(arg_mask_str[3+3*k], m_op);
      k++;
    }
  }
  Free2D((void **)mask_str);
  return (k-1);
}

int check_fsds_bit_str_val(char *fname, char *sname, char *bit_str, int *opt, int *l2g_st)

/*
!C******************************************************************************
    
!Function: check_fsds_bit_str_val
       
!Description:
    
  Checks to see if the input mask string is valid.
    
!Input Parameters: 
  fname           Input masking filename.
  sname           Input masking SDS name.
  bit_str         String contains the information of input masking bit descriptor, 
                  equality operatior(<. > , <=, >=, !=, ==) and binary bit value string 
		  to match
  opt             Flag of if the equality opertor is presented in the mask.
  l2g_st          Flag of whether the input mask file is in L2G file format.

    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*******************************************************************************/

{
  int32 sz;
  sds_t sds_info;
  int st, meta_cnt;
  /* int len; */
  char **meta_val, *ametadata;
  char c, tmp[20], meta_name[MAX_PATH_LENGTH];
		
  if ((meta_val = (char **)Calloc2D(5, 10, sizeof(char))) == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for meta_val in mask_sds_lib: check_fsds_bit_str_val()\n");
    return -1;
  }

  meta_cnt = 0;
  /* len = (int)strlen(fname); */
  if ((ametadata = get_attr_metadata(fname, "ArchiveMetadata.0")) != NULL)
  {
    strcpy(meta_name, "NUMBEROFOVERLAPGRANULES");
    get_sel_metadata(ametadata, meta_name, meta_val, &meta_cnt, 0);
    free(ametadata);
  }
  *l2g_st = ((meta_cnt != 0) && (atoi(meta_val[0]) != 0)) ? 1 : 0;
  Free2D((void **)meta_val);

  if (*l2g_st == 0)
    strcpy(sds_info.name, sname);
  else
  {
    if (sd_strpos(sname, "_1.", 0) != -1)
      strcpy(sds_info.name, sname);
    else
    {
      /*       len = (int)strlen(sname); */
      /*       p1 = sd_charpos(sname, '.', 0); */
      /*       sd_strmid(sname, 0, p1, org_name); */
      /*       sd_strmid(sname, p1, len-p1, org_num); */
      /*       sprintf(sds_info.name, "%s_1%s", org_name, org_num); */
      sprintf(sds_info.name, "%s_1", sname);
      strcpy(sname, sds_info.name);
    }
  }
  sds_info.sd_id = sds_info.sds_id = -1;
  if ((st = get_sds_info(fname, &sds_info)) != -1)
  {
    sz = 8*sds_info.data_size;
    c = bit_str[0];
    if ((c == '=') || (c == '>') || (c == '<') || (c == '!'))
    {
      *opt = 0;
      sprintf(tmp, "%d-" LONG_INT_FMT "%s", 0, sz-1, bit_str);
      strcpy(bit_str, tmp);
    }
    else *opt = 1;
  }
  if (sds_info.sds_id != -1) SDendaccess(sds_info.sds_id);
  if (sds_info.sd_id != -1) SDend(sds_info.sd_id);
  return st;
}

int get_parameters(char **arg_list, int n_op, int *sel_qa_op, char **qa_fnames,
        sds_t *qa_sds_info, unsigned long *bit_mask_arr, unsigned long *mask_val_arr,
        int *opt_arr, int *rel_op)

/*
!C******************************************************************************
    
!Function: get_parameters
       
!Description: 
  retrieve the masking information.
    
!Input Parameters: 

  arg_list             String contains the user input masking string.     
  n_op                 Number of Logic operation in the user input masking string.
  sel_qa_op            Array holding the logic operations.
                           "AND"  1
			   "OR"   2

  qa_fname             The input masking filename.
  qa_sds_info          Struct contains the masking SDS information.
  bit_mask_arr         Array containing the bit mask.   
  mask_val_arr         Array containing the mask value.
  opt_arr              Flag of if the equality opertor is presented in the mask.
                           
  rel_op               Integer number represents the equality operator: 
                           == 0
			   <  1
			   >  2
			   <= 3
			   >= 4
			   != 5
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*******************************************************************************/


{
  int i;
  int status;
  char op_str[5];

  status = 1;
  for (i=0; i <= n_op; i++)
  {
    if (get_file_sds_names(arg_list[1+i*3], qa_fnames[i], qa_sds_info[i].name) == -1)
      status = -1;
    if (get_bit_num_val(arg_list[2+i*3], &bit_mask_arr[i], &mask_val_arr[i], opt_arr[i],
                &rel_op[i]) == -1)
      status = -1;
    else if (i != n_op)
    {
      strcpy(op_str, arg_list[3+i*3]);
      if (strcmp(op_str, "AND") == 0) sel_qa_op[i] = 1;
      else if (strcmp(op_str, "OR") == 0) sel_qa_op[i] = 2;
      else
      {
        fprintf(stderr, "Wrong logical operator type encountered in get_parameters\n");
        status = -1;
      }
    }
  }
  return status;
}

int get_file_sds_names(char *fsds_name, char *fname, char *sds_name)


/*
!C******************************************************************************
    
!Function: get_file_sds_names
       
!Description:
    
   Retrieve the user input file name and SDS name from fsds_name(which is created by 
   get_mask_string().
    
!Input Parameters: 

  fsds_name       String created by get_mask_string() which contains the user input
                  file name and SDS name pair.
  fname           User input file name.
  sds_name        User input SDS name for its associated input file.
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*******************************************************************************/

{
  int len, pos;

  pos = sd_charpos(fsds_name, ':', 0);
  if (pos != -1)
  {
    len = (int)strlen(fsds_name);
    sd_strmid(fsds_name, 0, pos, fname);
    sd_strmid(fsds_name, pos+1, len-pos-1, sds_name);
    sd_strtrim(fname);
    sd_strtrim(sds_name);
    return 1;
  }
  else
  {
    fprintf(stderr, "Error separating filename and sdsname in get_file_sds_names\n");
    fname[0] = '\0';
    sds_name[0] = '\0';
    return -1;
  }
}

int get_bit_num_val(char *in_str, unsigned long *bit_mask, unsigned long *mask_val,
   int opt, int *rop)

/*
!C******************************************************************************
    
!Function: get_file_sds_names
       
!Description:
    
   Retrieve the user input file name and SDS name from fsds_name(which is created by 
   get_mask_string().
    
!Input Parameters: 

  fsds_name       String created by get_mask_string() which contains the user input
                  file name and SDS name pair.
  fname           User input file name.
  sds_name        User input SDS name for its associated input file.
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*******************************************************************************/

{
  int status;
  int i, cnt;
  int len, len2;
  int val1, val2, bnum;
  int pos1, pos2, pos3;
  int bit_num[MAX_NUM_BITS];
  char t1[5], t2[5], num[5];
  char num_str[MAX_STR_LEN];
  char val_str[MAX_STR_LEN];

  status = 1;
  len = (int)strlen(in_str);
  pos1 = sd_charpos(in_str, '<', 0);
  if (pos1 == -1) pos1 = sd_charpos(in_str, '>', 0);
  if (pos1 == -1) pos1 = sd_charpos(in_str, '!', 0);
  if (pos1 == -1) pos1 = sd_charpos(in_str, '=', 0);
  if (pos1 != -1)
  {
    sd_strmid(in_str, 0, pos1, num_str);
    pos1++;
    if (sd_charpos(in_str, '=', pos1) != -1) pos1++;
    sd_strmid(in_str, pos1, len-pos1, val_str);
    if (opt == 1) sd_strrev(val_str);
    if (strstr(in_str, "<=") != NULL)  *rop = 3;
    else if (strstr(in_str, ">=") != NULL) *rop = 4;
    else if (strstr(in_str, "!=") != NULL) *rop = 5;
    else if (strstr(in_str, "==") != NULL) *rop = 0;
    else if (strstr(in_str, "<") != NULL) *rop = 1;
    else if (strstr(in_str, ">") != NULL) *rop = 2;
  }
  else
  {
    fprintf(stderr, "Error reading num_str and val_str in get_bit_num_val\n");
    num_str[0] = '\0';
    val_str[0] = '\0';
    status = -1;
  }

  if (status != -1)
  {
    len = (int)strlen(num_str);
    cnt = pos1 = 0;
    while (pos1 < len)
    {
      pos2 = sd_charpos(num_str, ',', pos1);
      if (pos2 == -1) pos2 = len;
      sd_strmid(num_str, pos1, pos2-pos1, num);
      pos3 =  sd_charpos(num, '-', 0);
      sd_strtrim(num);
      if (pos3 == -1) bit_num[cnt++] = atoi(num);
      else
      {
        len2 = (int)strlen(num);
        sd_strmid(num, 0, pos3, t1);
        sd_strmid(num, pos3+1, len2-pos3-1, t2);
        sd_strtrim(t1);
        sd_strtrim(t2);
        val1 = atoi(t1);
        val2 = atoi(t2);
        for (i=val1; i<=val2; i++)
          bit_num[cnt++] = i;
      }
      pos1 = pos2+1;
    }

    *bit_mask = 0;
    *mask_val = 0;
    for (i=0; i<cnt; i++)
    {
      bnum = bit_num[i];
      *bit_mask = *bit_mask | BIT[bnum];
    }
    if (opt == 1)
      for (i=0; i<cnt; i++)
      {
        bnum = bit_num[i];
        if (val_str[i] == '1')
          *mask_val = *mask_val | BIT[bnum];
      }
    else *mask_val = atoi(val_str);
  }
  return status;
}

int malloc_qa_sds(sds_t *qa_sds_info, int n_op, int *fqa_l2g, void **data_qa, int32 **data_qa_nadd)

/*
!C******************************************************************************
    
!Function: malloc_qa_sds
       
!Description:
    
   
    
!Input Parameters: 

  qa_sds_info
  n_op
  fqa_l2g
  data_qa
  data_qa_nadd   
    
!Revision History: (see file prolog) 
    
!Team-unique Header: (see file prolog)
    
!References and Credits: (see file prolog)
    
!Design Notes: (none)
    
!END
*******************************************************************************/


{
  int i, j, k;
  sds_t sds_info;
  int status = 1;
  int ndata_qa, rank;

  rank = qa_sds_info[0].rank;
  if ((rank == 2) || (qa_sds_info[0].dim_size[0] > qa_sds_info[0].dim_size[rank-1]))
  {
    ndata_qa = qa_sds_info[0].dim_size[1];
    for (k=2; k<rank; k++)
      ndata_qa *= qa_sds_info[0].dim_size[k];
  }
  else
  {
    ndata_qa = qa_sds_info[0].dim_size[rank-1];
    for (k=0; k<rank-2; k++)
      ndata_qa *= qa_sds_info[0].dim_size[k];
  }
  if ((data_qa[0] = (void *)calloc(ndata_qa, qa_sds_info[0].data_size)) == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for data_qa in malloc_qa_sds\n");
    status = -1;
  }
  for (i=1; i<=n_op; i++)
  {
    for (j=0; j<i; j++)
    {
      if ((qa_sds_info[i].sd_id == qa_sds_info[j].sd_id) && (qa_sds_info[i].sds_id == qa_sds_info[j].sds_id))
      {
        data_qa[i] = data_qa[j];
        break;
      }
    }
    if (j >= i)
    {
      rank = qa_sds_info[i].rank;
      if ((rank == 2) || (qa_sds_info[i].dim_size[0] > qa_sds_info[i].dim_size[rank-1]))
      {
        ndata_qa = qa_sds_info[i].dim_size[1];
        for (k=2; k<rank; k++)
          ndata_qa *= qa_sds_info[i].dim_size[k];
      }
      else
      {
        ndata_qa = qa_sds_info[i].dim_size[rank-1];
        for (k=0; k<rank-2; k++)
          ndata_qa *= qa_sds_info[i].dim_size[k];
      }

      if ((data_qa[i] = (void *)calloc(ndata_qa, qa_sds_info[i].data_size)) == NULL)
      {
        fprintf(stderr, "Cannot allocate memory for data_qa in malloc_qa_sds\n");
        status = -1;
      }
    }
  }
  if (status != -1)
  {
    for (i=0; i<=n_op; i++)
      if (fqa_l2g[i] == 1)
      {
	for (j=0; j<i; j++)
	  if (qa_sds_info[i].sd_id == qa_sds_info[j].sd_id) break;
	if (j >= i)
	{
          sds_info.rank = 1;
          strcpy(sds_info.name, "nadd_obs_row");
          sds_info.sd_id = qa_sds_info[i].sd_id;
          sds_info.dim_size[0] = qa_sds_info[i].dim_size[0];
          if ((data_qa_nadd[i] = (int32 *)calloc(sds_info.dim_size[0], sizeof(int32))) == NULL)
          {
            fprintf(stderr, "Cannot allocate memory for data_qa_nadd[] in malloc_qa_sds");
	    status = -1; break;
          }
          else 
	  {
	    get_sds_data(&sds_info, data_qa_nadd[i]);
	    SDendaccess(sds_info.sds_id);
	  }
	}
	else data_qa_nadd[i] = data_qa_nadd[j];
      }
  }
  return status;
}

int open_qa_sds(char *fname, sds_t *sds_info, char **qa_fnames, sds_t *qa_sds_info, int n_op)
{
  int status;
  int i, j, p1;
  char sdsi_name[MAX_SDS_NAME_LEN];
  char sdsj_name[MAX_SDS_NAME_LEN];

  status = 1;
  for (i=0; i<=n_op; i++)
    qa_sds_info[i].sd_id = qa_sds_info[i].sds_id = -1;
  for (i=0; i<=n_op; i++)
  {
    if ((p1 = sd_charpos(qa_sds_info[i].name, '.', 0)) == -1) 
      strcpy(sdsi_name, qa_sds_info[i].name);
    else sd_strmid(qa_sds_info[i].name, 0, p1, sdsi_name);
    if ((fname != NULL) && (strcmp(fname, qa_fnames[i]) == 0))
    {
      qa_sds_info[i].sd_id = sds_info->sd_id;
      if ((p1 = sd_charpos(sds_info->name, '.', 0)) == -1) 
	strcpy(sdsj_name, sds_info->name);
      else sd_strmid(sds_info->name, 0, p1, sdsj_name);
      if (strcmp(sdsi_name, sdsj_name) == 0) 
	qa_sds_info[i].sds_id = sds_info->sds_id;
      else
      {
        for (j=0; j<i; j++)
	{
          if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
            strcpy(sdsj_name, qa_sds_info[j].name);
          else sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
	  if ((strcmp(qa_fnames[i], qa_fnames[j]) == 0) && (strcmp(sdsi_name, sdsj_name) == 0))
          {
            qa_sds_info[i].sd_id = qa_sds_info[j].sd_id;
            qa_sds_info[i].sds_id = qa_sds_info[j].sds_id;
            break;
          }
	}
        if (j >= i)
        {
          if (open_sds(qa_fnames[i], &qa_sds_info[i], 'R') == -1)
          {
            status = -1;
            if (qa_sds_info[i].sds_id != -1) SDendaccess(qa_sds_info[i].sds_id);
            if (qa_sds_info[i].sd_id != -1) SDend(qa_sds_info[i].sd_id);
            qa_sds_info[i].sd_id = qa_sds_info[i].sds_id = -1;
          }
        }
      }
    }
    else
    {
      for (j=0; j<i; j++)
      {
        if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
	  strcpy(sdsj_name, qa_sds_info[j].name);
        else sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
        if ((strcmp(qa_fnames[i], qa_fnames[j]) == 0) && (strcmp(sdsi_name, sdsj_name) == 0))
        {
          qa_sds_info[i].sd_id = qa_sds_info[j].sd_id;
          qa_sds_info[i].sds_id = qa_sds_info[j].sds_id;
          break;
        }
      }
      if (j >= i)
      {
        if (open_sds(qa_fnames[i], &qa_sds_info[i], 'R') == -1)
        {
          status = -1;
          if (qa_sds_info[i].sds_id != -1) SDendaccess(qa_sds_info[i].sds_id);
          if (qa_sds_info[i].sd_id != -1) SDend(qa_sds_info[i].sd_id);
          qa_sds_info[i].sd_id = qa_sds_info[i].sds_id = -1;
        }
      }
    }
  }
  return status;
}

int open_qa_sds_nsds(char *fname, sds_t *sds_info, sds_t *sdsc_info, sds_t *sds_nobs_info, int nsds, 
	char **qa_fnames, sds_t *qa_sds_info, sds_t *qa_sdsc_info, sds_t *qa_sds_nobs_info, 
	int *qa_l2g, int n_op)
{
  int status;
  int i, j, p1, isds = 0, idim;
  char sdsi_name[MAX_SDS_NAME_LEN];
  char sdsj_name[MAX_SDS_NAME_LEN];

  status = 1;
  for (i=0; i<=n_op; i++)
  {
    qa_sds_info[i].sd_id = qa_sds_info[i].sds_id = -1;
    if (qa_l2g[i] == 1)
      qa_sdsc_info[i].sd_id = qa_sdsc_info[i].sds_id = -1;
  }

  for (i=0; i<=n_op; i++)
  {
    if ((p1 = sd_charpos(qa_sds_info[i].name, '.', 0)) == -1) 
      strcpy(sdsi_name, qa_sds_info[i].name);
    else sd_strmid(qa_sds_info[i].name, 0, p1, sdsi_name);
    if ((fname != NULL) && (strcmp(fname, qa_fnames[i]) == 0))
    {
      qa_sds_info[i].sd_id = sds_info[0].sd_id;
      for (isds=0; isds<nsds; isds++)
      {
        if ((p1 = sd_charpos(sds_info[isds].name, '.', 0)) == -1) 
	  strcpy(sdsj_name, sds_info[isds].name);
        else sd_strmid(sds_info[isds].name, 0, p1, sdsj_name);
        if (strcmp(sdsi_name, sdsj_name) == 0) 
        {
	  qa_sds_info[i].sds_id = sds_info[isds].sds_id;
          if (qa_l2g[i] == 1)
	  {
	    qa_sdsc_info[i].sd_id = sdsc_info[isds].sd_id;
	    qa_sdsc_info[i].sds_id = sdsc_info[isds].sds_id;
	    qa_sdsc_info[i].data_type = sdsc_info[isds].data_type;
	    qa_sds_nobs_info[i].sd_id = sds_nobs_info->sd_id;
	    qa_sds_nobs_info[i].sds_id = sds_nobs_info->sds_id;
	    qa_sds_nobs_info[i].rank = sds_nobs_info->rank;
	    for (idim=0; idim<sds_nobs_info->rank; idim++)
	      qa_sds_nobs_info[i].dim_size[idim] = sds_nobs_info->dim_size[idim];
	  }
        }
      }
      if (isds >= nsds)
      {
        for (j=0; j<i; j++)
	{
          if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
            strcpy(sdsj_name, qa_sds_info[j].name);
          else sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
	  if ((strcmp(qa_fnames[i], qa_fnames[j]) == 0) && (strcmp(sdsi_name, sdsj_name) == 0))
          {
            qa_sds_info[i].sd_id = qa_sds_info[j].sd_id;
            qa_sds_info[i].sds_id = qa_sds_info[j].sds_id;
            if (qa_l2g[i] == 1)
            {
              qa_sdsc_info[i].sd_id = sdsc_info[isds].sd_id;
              qa_sdsc_info[i].sds_id = sdsc_info[isds].sds_id;
	      qa_sdsc_info[i].data_type = sdsc_info[isds].data_type;
              qa_sds_nobs_info[i].sd_id = sds_nobs_info->sd_id;
              qa_sds_nobs_info[i].sds_id = sds_nobs_info->sds_id;
	      qa_sds_nobs_info[i].rank = sds_nobs_info->rank;
	      for (idim=0; idim<sds_nobs_info->rank; idim++)
	        qa_sds_nobs_info[i].dim_size[idim] = sds_nobs_info->dim_size[idim];
            }
            break;
          }
	}
        if (j >= i)
        {
          if (open_sds((char *)NULL, &qa_sds_info[i], 'R') == -1)
          {
            status = -1;
            if (qa_sds_info[i].sds_id != -1) SDendaccess(qa_sds_info[i].sds_id);
            if (qa_sds_info[i].sd_id != -1) SDend(qa_sds_info[i].sd_id);
            qa_sds_info[i].sd_id = qa_sds_info[i].sds_id = -1;
          }
	  else if (qa_l2g[i] == 1)
	  {
	    qa_sdsc_info[i].sds_id = qa_sds_nobs_info[i].sds_id = -1;
	    qa_sdsc_info[i].sd_id = qa_sds_nobs_info[i].sd_id = qa_sds_info[i].sd_id;
	    strcpy(qa_sds_nobs_info[i].name, "num_observations");
	    open_sds((char *)NULL, &qa_sdsc_info[i], 'R');
	    open_sds((char *)NULL, &qa_sds_nobs_info[i], 'R');
	    qa_sds_nobs_info[i].rank = qa_sds_info[i].rank;
	    for (idim=0; idim<qa_sds_info[i].rank; idim++)
	      qa_sds_nobs_info[i].dim_size[idim] = qa_sds_info[i].dim_size[idim];
	  }
        }
      }
    }
    else
    {
      for (j=0; j<i; j++)
      {
        if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
	  strcpy(sdsj_name, qa_sds_info[j].name);
        else sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
        if ((strcmp(qa_fnames[i], qa_fnames[j]) == 0) && (strcmp(sdsi_name, sdsj_name) == 0))
        {
          qa_sds_info[i].sd_id = qa_sds_info[j].sd_id;
          qa_sds_info[i].sds_id = qa_sds_info[j].sds_id;
          if (qa_l2g[i] == 1)
          {
            qa_sdsc_info[i].sd_id = sdsc_info[isds].sd_id;
            qa_sdsc_info[i].sds_id = sdsc_info[isds].sds_id;
            qa_sdsc_info[i].data_type = sdsc_info[isds].data_type;
            qa_sds_nobs_info[i].sd_id = sds_nobs_info->sd_id;
            qa_sds_nobs_info[i].sds_id = sds_nobs_info->sds_id;
	    qa_sds_nobs_info[i].rank = sds_nobs_info->rank;
	    for (idim=0; idim<sds_nobs_info->rank; idim++)
	      qa_sds_nobs_info[i].dim_size[idim] = sds_nobs_info->dim_size[idim];
          }
          break;
        }
      }
      if (j >= i)
      {
        if (open_sds(qa_fnames[i], &qa_sds_info[i], 'R') == -1)
        {
          status = -1;
          if (qa_sds_info[i].sds_id != -1) SDendaccess(qa_sds_info[i].sds_id);
          if (qa_sds_info[i].sd_id != -1) SDend(qa_sds_info[i].sd_id);
          qa_sds_info[i].sd_id = qa_sds_info[i].sds_id = -1;
        }
        else if (qa_l2g[i] ==1)
        {
          qa_sdsc_info[i].sds_id = qa_sds_nobs_info[i].sds_id = -1;
          qa_sdsc_info[i].sd_id = qa_sds_nobs_info[i].sd_id = qa_sds_info[i].sd_id;
	  strcpy(qa_sds_nobs_info[i].name, "num_observations");
          open_sds((char *)NULL, &qa_sdsc_info[i], 'R');
          open_sds((char *)NULL, &qa_sds_nobs_info[i], 'R');
	  qa_sds_nobs_info[i].rank = qa_sds_info[i].rank;
	  for (idim=0; idim<qa_sds_info[i].rank; idim++)
	    qa_sds_nobs_info[i].dim_size[idim] = qa_sds_info[i].dim_size[idim];
        }
      }
    }
  }
  return status;
}

void read_qa_sds(sds_t *qa_sds_info, sds_t *qa_sdsc_info, sds_t *qa_sds_nobs_info, int n_op, void **data_qa, int32 **data_qa_nadd, int irow, int *res_l, int *fqa_l2g, int *obs_num)
{
  int i, j, k, rank;
  int32 edge[4] = {0, 0, 0, 0};
  int32 start[4] = {0, 0, 0, 0};

  rank = qa_sds_info[0].rank;
  for (i=0; i<rank; i++)
    edge[i] = qa_sds_info[0].dim_size[i];

  if ((rank == 2) || (qa_sds_info[0].dim_size[0] > qa_sds_info[0].dim_size[rank-1]))
  { 
    start[0] = irow/res_l[0]; edge[0] = 1; 
  }
  else 
  { 
    start[rank-2] = irow/res_l[0]; edge[rank-2] = 1; 
  } 
  if ((fqa_l2g[0] == 0) || (obs_num[0] == 1))
  {
    if (SDreaddata(qa_sds_info[0].sds_id, start, NULL, edge, data_qa[0]) == FAIL)
      fprintf(stderr, "Cannot read data line from SDS %s in mask_sds_lib:read_qa_sds()\n", 
		qa_sds_info[0].name);
  }
  else
    read_sdsc_data(&qa_sdsc_info[0], &qa_sds_nobs_info[0], data_qa[0], data_qa_nadd[0], 
		start[0], obs_num[0]);

  for (i=1; i<=n_op; i++)
  {
    for (j=0; j<i; j++)
    {
      if ((qa_sds_info[i].sd_id == qa_sds_info[j].sd_id) &&
		(qa_sds_info[i].sds_id == qa_sds_info[j].sds_id))
      {
	if ((fqa_l2g[i] == 0) || (obs_num[i] == obs_num[j]))  data_qa[i] = data_qa[j];
	else 
        {
	  if (obs_num[i] == 1)
	  {
	    if (SDreaddata(qa_sds_info[i].sds_id, start, NULL, edge, data_qa[i]) == FAIL)
      	      fprintf(stderr, "Cannot read data line from SDS %s in mask_sds_lib:read_qa_sds()\n", 
			qa_sds_info[i].name);
	  }
	  else
	    read_sdsc_data(&qa_sdsc_info[i], &qa_sds_nobs_info[i], data_qa[i], data_qa_nadd[i], 
                start[0], obs_num[i]);
	}
	break;
      }
    }
    if (j >= i)
    {
      if ((res_l[i] == 1) || (irow%res_l[i] == 0))
      {
        rank = qa_sds_info[i].rank;
        for (k=0; k<4; k++)
          start[k] = edge[k] = 0;
        for (k=0; k<rank; k++)
          edge[k] = qa_sds_info[i].dim_size[k];
        if ((rank == 2) || (qa_sds_info[i].dim_size[0] > qa_sds_info[i].dim_size[rank-1]))
        { 
          start[0] = irow/res_l[i]; edge[0] = 1;
        }
        else
        { 
          start[rank-2] = irow/res_l[i]; edge[rank-2] = 1;
        }
        if ((fqa_l2g[i] == 0) || (obs_num[i] == 1))
        {
          if (SDreaddata(qa_sds_info[i].sds_id, start, NULL, edge, data_qa[i]) == FAIL)
            fprintf(stderr, "Cannot read data line from SDS %s in mask_sds_lib:read_qa_sds()\n", 
                      qa_sds_info[i].name);
        }
        else
          read_sdsc_data(&qa_sdsc_info[i], &qa_sds_nobs_info[i], data_qa[i], data_qa_nadd[i],
                      start[0], obs_num[i]);
      }
    }
  } /* for (i=1; . . . ) */
}

void read_sdsc_data(sds_t *sdsc_info, sds_t *sds_nobs_info, void *data, int32 *data_nadd, 
			int irow, int obs_num)
{
  int i, nadd_obs;
  int ir, ic, ncols;
  int nobs, obs_num_c;
  int32 edge[2] = {0, 0};
  int32 start[2] = {0, 0};
  void *data_c;
  int8 *data_nobs;

  ncols = sds_nobs_info->dim_size[1];
  if (data_nadd[irow] >= 1)
  {
    for (ir=0, nadd_obs=0; ir<irow; ir++)
      if (data_nadd[ir] >= 1) nadd_obs += data_nadd[ir];
    start[0] = nadd_obs; edge[0] = data_nadd[irow];
    if ((data_c = (void *)calloc(edge[0], sdsc_info->data_size)) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_c in mask_sds_lib: read_sdsc_update_data()\n");
    if (SDreaddata(sdsc_info->sds_id, start, NULL, edge, data_c) == FAIL)
      fprintf(stderr, "Cannot read data line from SDS %s in mask_sds_lib: read_sdsc_update_data()\n", 
			  sdsc_info->name);
  
    start[0] = irow; edge[0] = 1;
    start[1] = 0; edge[1] = ncols;
    if ((data_nobs = (int8 *)calloc(edge[1], sizeof(int8))) == NULL)
      fprintf(stderr, "Cannot allocate memory for data_nobs in mask_sds_lib: read_sdsc_update_data()\n");
    if (SDreaddata(sds_nobs_info->sds_id, start, NULL, edge, (VOIDP)data_nobs) == FAIL)
      fprintf(stderr, "Cannot read data line from SDS %s in mask_sds_lib: read_sdsc_update_data()\n", 
			  sds_nobs_info->name);
    obs_num_c = obs_num - 2;
    for (i=0, ic=0; i<ncols; i++)
    {
      nobs = data_nobs[i];
      if (nobs >= obs_num)
      {
        switch(sdsc_info->data_type)
        {
          case 20: ((int8 *)data)[i] = ((int8 *)data_c)[ic+obs_num_c]; break;
          case 21: ((uint8 *)data)[i] = ((uint8 *)data_c)[ic+obs_num_c]; break;
          case 22: ((int16 *)data)[i] = ((int16 *)data_c)[ic+obs_num_c]; break;
          case 23: ((uint16 *)data)[i] = ((uint16 *)data_c)[ic+obs_num_c]; break;
          case 24: ((int32 *)data)[i] = ((int32 *)data_c)[ic+obs_num_c]; break;
          case 25: ((uint32 *)data)[i] = ((uint32 *)data_c)[ic+obs_num_c]; break;
        }
      }
      else
      {
        switch(sdsc_info->data_type)
        {
	  case 20: ((int8 *)data)[i] = (int8)sdsc_info->fill_val; break; 
	  case 21: ((uint8 *)data)[i] = (uint8)sdsc_info->fill_val; break; 
	  case 22: ((int16 *)data)[i] = (int16)sdsc_info->fill_val; break; 
	  case 23: ((uint16 *)data)[i] = (uint16)sdsc_info->fill_val; break; 
	  case 24: ((int32 *)data)[i] = sdsc_info->fill_val; break; 
	  case 25: ((uint32 *)data)[i] = sdsc_info->fill_val; break; 
        }
      }
      if (nobs > 1) ic = ic + nobs - 1;
    } /* for (i=0; . .  ) */
    free(data_c);
    free(data_nobs);
  }
  else
  {
    for (i=0; i<ncols; i++)
      switch(sdsc_info->data_type)
      {
	case 20: ((int8 *)data)[i] = (int8)sdsc_info->fill_val; break; 
	case 21: ((uint8 *)data)[i] = (uint8)sdsc_info->fill_val; break; 
	case 22: ((int16 *)data)[i] = (int16)sdsc_info->fill_val; break; 
	case 23: ((uint16 *)data)[i] = (uint16)sdsc_info->fill_val; break; 
	case 24: ((int32 *)data)[i] = sdsc_info->fill_val; break; 
	case 25: ((uint32 *)data)[i] = sdsc_info->fill_val; break; 
      }
  }
}

void close_qa_hdf(char *hdf_fname, sds_t *sds_info, char **qa_fnames, sds_t *qa_sds_info, int n_op)
{
  int i, j, p1;
  char sdsi_name[MAX_SDS_NAME_LEN];
  char sdsj_name[MAX_SDS_NAME_LEN];

  for (i=0; i<=n_op; i++)
  {
    if ((p1 = sd_charpos(qa_sds_info[i].name, '.', 0)) == -1)
      strcpy(sdsi_name, qa_sds_info[i].name);
    else sd_strmid(qa_sds_info[i].name, 0, p1, sdsi_name);
    for (j=0; j<i; j++)
    {
      if (qa_sds_info[i].sd_id == qa_sds_info[j].sd_id)
      {
        if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
	  strcpy(sdsj_name, qa_sds_info[j].name);
        else sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
        if (strcmp(sdsi_name, sdsj_name) == 0) { qa_sds_info[i].sds_id = -1; break; }
      }
    }
    if (qa_sds_info[i].sds_id != -1)
    {
      if (SDendaccess(qa_sds_info[i].sds_id) == FAIL)
        fprintf(stderr, "Cannot close SDS %s in mask_sds_lib: close_qa_hdf()\n", qa_sds_info[i].name);
      if ((hdf_fname != NULL) && (qa_sds_info[i].sds_id == sds_info->sds_id) && (qa_sds_info[i].sd_id == sds_info->sd_id))
        qa_sds_info[i].sds_id = sds_info->sds_id = -1;
      else
        qa_sds_info[i].sds_id = -1;
    }
  }

  /* close all the qa hdf files if it is not also an input file */
  for (i=0; i<=n_op; i++)
  {
    for (j=0; j<i; j++)
      if (strcmp(qa_fnames[i], qa_fnames[j]) == 0) { qa_sds_info[i].sd_id = -1; break; }
    if ((qa_sds_info[i].sd_id != -1) && (hdf_fname != NULL) && (qa_sds_info[i].sd_id != sds_info->sd_id))
    {
      if (SDend(qa_sds_info[i].sd_id) == FAIL)
        fprintf(stderr, "Cannot close the HDF file %s in mask_sds_lib: close_qa_hdf()\n", qa_fnames[i]);
      qa_sds_info[i].sd_id = -1;
    }
  }
}

void close_qa_hdf_nsds(char *hdf_fname, sds_t *sds_info, int nsds, char **qa_fnames,
        sds_t *qa_sds_info, int n_op)
{
  int i, j, p1, isds;
  char sdsi_name[MAX_SDS_NAME_LEN];
  char sdsj_name[MAX_SDS_NAME_LEN];

  for (i=0; i<=n_op; i++)
  {
    if ((p1 = sd_charpos(qa_sds_info[i].name, '.', 0)) == -1)
      strcpy(sdsi_name, qa_sds_info[i].name);
    else sd_strmid(qa_sds_info[i].name, 0, p1, sdsi_name);
    for (isds=0; isds<nsds; isds++)
    {
      if (qa_sds_info[i].sd_id == sds_info[isds].sd_id)
      {
        if ((p1 = sd_charpos(sds_info[isds].name, '.', 0)) == -1)
          strcpy(sdsj_name, sds_info[isds].name);
        else sd_strmid(sds_info[isds].name, 0, p1, sdsj_name);
	if (strcmp(sdsi_name, sdsj_name) == 0)
	{
          if (SDendaccess(qa_sds_info[i].sds_id) == FAIL)
            fprintf(stderr, "Cannot close SDS %s in mask_sds_lib: close_qa_hdf()\n", qa_sds_info[i].name);
          qa_sds_info[i].sds_id = sds_info[isds].sds_id = -1;
          break;
	}
      }
    }
    if (qa_sds_info[i].sds_id != -1)
    {
      for (j=0; j<i; j++)
      {
        if (qa_sds_info[i].sd_id == qa_sds_info[j].sd_id)
        {
          if ((p1 = sd_charpos(qa_sds_info[j].name, '.', 0)) == -1)
	    strcpy(sdsj_name, qa_sds_info[j].name);
          else sd_strmid(qa_sds_info[j].name, 0, p1, sdsj_name);
          if (strcmp(sdsi_name, sdsj_name) == 0) { qa_sds_info[i].sds_id = -1; break; }
        }
      }
      if (qa_sds_info[i].sds_id != -1)
      {
        if (SDendaccess(qa_sds_info[i].sds_id) == FAIL)
          fprintf(stderr, "Cannot close SDS %s in mask_sds_lib: close_qa_hdf()\n", qa_sds_info[i].name);
        qa_sds_info[i].sds_id = sds_info->sds_id = -1;
      }
    }
  }

  /* close all the qa hdf files if it is not also an input file */
  for (i=0; i<=n_op; i++)
  {
    for (j=0; j<i; j++)
      if (strcmp(qa_fnames[i], qa_fnames[j]) == 0) { qa_sds_info[i].sd_id = -1; break; }
    if ((qa_sds_info[i].sd_id != -1) && (hdf_fname != NULL) && (qa_sds_info[i].sd_id != sds_info->sd_id))
    {
      if (SDend(qa_sds_info[i].sd_id) == FAIL)
        fprintf(stderr, "Cannot close the HDF file %s in mask_sds_lib: close_qa_hdf()\n", qa_fnames[i]);
      qa_sds_info[i].sd_id = -1;
    }
  }
}

void process_mask_data(void **data_qa, int ncols, sds_t *qa_sds_info, int n_op, int *sel_qa_op, 
	unsigned long *bit_mask_arr, unsigned long *mask_val_arr, int *rel_op, int *res_s, 
	uint8 *mask_row, int on_val, int off_val, int mask_fill)
{
  int mask_st;
  int sel_pix_fin;
  int offset[MAX_NUM_OP], j[MAX_NUM_OP];
  int n, m, i, jj, i_op;
  int sel_pix[MAX_NUM_OP];
  long masked_val, band_qa_at_ij;
  unsigned long u_masked_val, u_band_qa_at_ij;
  char sds_name[MAX_SDS_NAME_LEN];

  for (i_op=0; i_op<=n_op; i_op++)
  {
    get_sdsname_dim(qa_sds_info[i_op].name, sds_name, &n, &m);
    compute_sds_start_offset(&qa_sds_info[i_op], n, m, &j[i_op], &offset[i_op]);
  }

  for (i=0; i<ncols; i++)
  {
    mask_st = 1;
    for (i_op=0; i_op <= n_op; i_op++)
    {
      jj = j[i_op]/res_s[i_op];
      switch(qa_sds_info[i_op].data_type)
      {
        case 20: band_qa_at_ij = ((int8 *)data_qa[i_op])[jj];
		 if (band_qa_at_ij != qa_sds_info[i_op].fill_val)
		 {
        	   masked_val = band_qa_at_ij & bit_mask_arr[i_op];
        	   switch(rel_op[i_op])
        	   {
          	     case 0: sel_pix[i_op] = (masked_val == (int)mask_val_arr[i_op]) ? YES : NO; break;
          	     case 1: sel_pix[i_op] = (masked_val < (int)mask_val_arr[i_op]) ? YES : NO; break;
          	     case 2: sel_pix[i_op] = (masked_val > (int)mask_val_arr[i_op]) ? YES : NO; break;
          	     case 3: sel_pix[i_op] = (masked_val <= (int)mask_val_arr[i_op]) ? YES : NO; break;
          	     case 4: sel_pix[i_op] = (masked_val >= (int)mask_val_arr[i_op]) ? YES : NO; break;
          	     case 5: sel_pix[i_op] = (masked_val != (int)mask_val_arr[i_op]) ? YES : NO; break;
        	   }
		 }
		 else mask_st = -1;
		 break;
        case 21: u_band_qa_at_ij = ((uint8 *)data_qa[i_op])[jj];
                 if (u_band_qa_at_ij != (uint8)qa_sds_info[i_op].fill_val)
                 {
                   u_masked_val = u_band_qa_at_ij & bit_mask_arr[i_op];
                   switch(rel_op[i_op])
                   {
                     case 0: sel_pix[i_op] = (u_masked_val == mask_val_arr[i_op]) ? YES : NO; break;
                     case 1: sel_pix[i_op] = (u_masked_val < mask_val_arr[i_op]) ? YES : NO; break;
                     case 2: sel_pix[i_op] = (u_masked_val > mask_val_arr[i_op]) ? YES : NO; break;
                     case 3: sel_pix[i_op] = (u_masked_val <= mask_val_arr[i_op]) ? YES : NO; break;
                     case 4: sel_pix[i_op] = (u_masked_val >= mask_val_arr[i_op]) ? YES : NO; break;
                     case 5: sel_pix[i_op] = (u_masked_val != mask_val_arr[i_op]) ? YES : NO; break;
                   }
                 }
                 else mask_st = -1;
                 break;
        case 22: band_qa_at_ij = ((int16 *)data_qa[i_op])[jj];
                 if (band_qa_at_ij != qa_sds_info[i_op].fill_val)
                 {
                   masked_val = band_qa_at_ij & bit_mask_arr[i_op];
                   switch(rel_op[i_op])
                   {
                     case 0: sel_pix[i_op] = (masked_val == (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 1: sel_pix[i_op] = (masked_val < (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 2: sel_pix[i_op] = (masked_val > (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 3: sel_pix[i_op] = (masked_val <= (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 4: sel_pix[i_op] = (masked_val >= (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 5: sel_pix[i_op] = (masked_val != (int)mask_val_arr[i_op]) ? YES : NO; break;
                   }
                 }
                 else mask_st = -1;
                 break;
        case 23: u_band_qa_at_ij = ((uint16 *)data_qa[i_op])[jj];
                 if (u_band_qa_at_ij != (uint16)qa_sds_info[i_op].fill_val)
                 {
                   u_masked_val = u_band_qa_at_ij & bit_mask_arr[i_op];
                   switch(rel_op[i_op])
                   {
                     case 0: sel_pix[i_op] = (u_masked_val == mask_val_arr[i_op]) ? YES : NO; break;
                     case 1: sel_pix[i_op] = (u_masked_val < mask_val_arr[i_op]) ? YES : NO; break;
                     case 2: sel_pix[i_op] = (u_masked_val > mask_val_arr[i_op]) ? YES : NO; break;
                     case 3: sel_pix[i_op] = (u_masked_val <= mask_val_arr[i_op]) ? YES : NO; break;
                     case 4: sel_pix[i_op] = (u_masked_val >= mask_val_arr[i_op]) ? YES : NO; break;
                     case 5: sel_pix[i_op] = (u_masked_val != mask_val_arr[i_op]) ? YES : NO; break;
                   }
                 }
                 else mask_st = -1;
                 break;
        case 24: band_qa_at_ij = ((int32 *)data_qa[i_op])[jj];
                 if (band_qa_at_ij != qa_sds_info[i_op].fill_val)
                 {
                   masked_val = band_qa_at_ij & bit_mask_arr[i_op];
                   switch(rel_op[i_op])
                   {
                     case 0: sel_pix[i_op] = (masked_val == (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 1: sel_pix[i_op] = (masked_val < (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 2: sel_pix[i_op] = (masked_val > (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 3: sel_pix[i_op] = (masked_val <= (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 4: sel_pix[i_op] = (masked_val >= (int)mask_val_arr[i_op]) ? YES : NO; break;
                     case 5: sel_pix[i_op] = (masked_val != (int)mask_val_arr[i_op]) ? YES : NO; break;
                   }
                 }
                 else mask_st = -1;
                 break;
        case 25: u_band_qa_at_ij = ((uint32 *)data_qa[i_op])[jj];
                 if (u_band_qa_at_ij != (uint32)qa_sds_info[i_op].fill_val)
                 {
                   u_masked_val = u_band_qa_at_ij & bit_mask_arr[i_op];
                   switch(rel_op[i_op])
                   {
                     case 0: sel_pix[i_op] = (u_masked_val == mask_val_arr[i_op]) ? YES : NO; break;
                     case 1: sel_pix[i_op] = (u_masked_val < mask_val_arr[i_op]) ? YES : NO; break;
                     case 2: sel_pix[i_op] = (u_masked_val > mask_val_arr[i_op]) ? YES : NO; break;
                     case 3: sel_pix[i_op] = (u_masked_val <= mask_val_arr[i_op]) ? YES : NO; break;
                     case 4: sel_pix[i_op] = (u_masked_val >= mask_val_arr[i_op]) ? YES : NO; break;
                     case 5: sel_pix[i_op] = (u_masked_val != mask_val_arr[i_op]) ? YES : NO; break;
                   }
                 }
                 else mask_st = -1;
                 break;
      }
      j[i_op] += offset[i_op];
    }  /* for i_op ... */
    if (mask_st == -1) 
      mask_row[i] = (uint8)mask_fill;
    else 
    {
      sel_pix_fin = sel_pix[0];
      if (n_op > 0)
        for (i_op = 0; i_op < n_op; i_op++)
        {
          if (sel_qa_op[i_op] == 1) sel_pix_fin &= sel_pix[i_op+1];
          else if (sel_qa_op[i_op] == 2) sel_pix_fin |= sel_pix[i_op+1];
        }
      mask_row[i] = (sel_pix_fin == YES) ? (uint8)on_val : (uint8)off_val;
    }
  }
}

int get_qa_sds_info(char **fnames, sds_t *sds_info, sds_t *sdsc_info, int *l2g_st, int n_op)
{
  int st = 1;
  int i_op, p1;

  for (i_op=0; i_op<=n_op; i_op++)
  {
    sds_info[i_op].sd_id = sds_info[i_op].sds_id = -1;
    if (get_sds_info(fnames[i_op], &sds_info[i_op]) == -1)
    {
      if (sds_info[i_op].sds_id != -1) SDendaccess(sds_info[i_op].sds_id);
      if (sds_info[i_op].sd_id != -1) SDend(sds_info[i_op].sd_id);
      st = -1;
      break;
    }
    SDendaccess(sds_info[i_op].sds_id);
    SDend(sds_info[i_op].sd_id);
    sds_info[i_op].sd_id = sds_info[i_op].sds_id = -1;
    if (l2g_st[i_op] == 1)
    {
      strcpy(sdsc_info[i_op].name, sds_info[i_op].name);
      p1 = sd_strpos(sdsc_info[i_op].name, "_1", 0);
      sdsc_info[i_op].name[p1+1] = '_';
						sdsc_info[i_op].name[p1+2] = 'c';
      sdsc_info[i_op].data_type = sds_info[i_op].data_type;
      sdsc_info[i_op].data_size = sds_info[i_op].data_size;
      sdsc_info[i_op].fill_val = sds_info[i_op].fill_val;
    }
  } /* for (i_op=0; . . .  ) */
  return st;
}

int get_in_sds_info(char *hdf_fname, sds_t *sds_info, sds_t *sdsc_info, sds_t *sds_nobs_info,
	int l2g_st, int nsds)
{
  int st = 1;
  int p1, isds, jsds, idim;
  char sdsi_name[MAX_SDS_NAME_LEN];
  char sdsj_name[MAX_SDS_NAME_LEN];

  for (isds=0; isds<nsds; isds++)
  {
    sds_info[isds].sds_id = -1;
    sds_info[isds].sd_id = (isds == 0) ? -1 : sds_info[0].sd_id;

    if ((p1 = sd_charpos(sds_info[isds].name, '.', 0)) != -1)
      sd_strmid(sds_info[isds].name, 0, p1, sdsi_name);
    else strcpy(sdsi_name, sds_info[isds].name);
    for (jsds=0; jsds<isds; jsds++)
    {
      if ((p1 = sd_charpos(sds_info[jsds].name, '.', 0)) != -1)
	sd_strmid(sds_info[jsds].name, 0, p1, sdsj_name);
      else strcpy(sdsj_name, sds_info[jsds].name);
      if (strcmp(sdsi_name, sdsj_name) == 0)
      {
	sds_info[isds].rank = sds_info[jsds].rank;
	sds_info[isds].sds_id = sds_info[jsds].sds_id;
	sds_info[isds].data_type = sds_info[jsds].data_type;
	sds_info[isds].data_size = sds_info[jsds].data_size;
	sds_info[isds].fill_val = sds_info[jsds].fill_val;
	for (idim=0; idim<sds_info[isds].rank; idim++)
	  sds_info[isds].dim_size[idim] = sds_info[jsds].dim_size[idim];
	if (l2g_st == 1)
	{
	  sdsc_info[isds].sd_id = sdsc_info[jsds].sd_id;
	  sdsc_info[isds].sds_id = sdsc_info[jsds].sds_id;
	  sdsc_info[isds].data_type = sdsc_info[jsds].data_type;
	  sdsc_info[isds].fill_val = sdsc_info[jsds].fill_val;
          strcpy(sdsc_info[isds].name, sds_info[isds].name);
          p1 = sd_strpos(sdsc_info[isds].name, "_1.", 0);
          sdsc_info[isds].name[p1+1] = 'c';
	}
      }
    }
    if (sds_info[isds].sds_id == -1)
    {
      if (get_sds_info(hdf_fname, &sds_info[isds]) == -1)
      {
        if (sds_info[isds].sds_id != -1) SDendaccess(sds_info[isds].sds_id);
        if (sds_info[isds].sd_id != -1) SDend(sds_info[isds].sd_id);
        st = -1;
      }
      if ((st != -1) && (l2g_st == 1))
      {
        sdsc_info[isds].sd_id = sds_info[isds].sd_id;
        sdsc_info[isds].sds_id = -1;
        strcpy(sdsc_info[isds].name, sds_info[isds].name);
        p1 = sd_strpos(sdsc_info[isds].name, "_1.", 0);
        sdsc_info[isds].name[p1+1] = 'c';
        get_sds_info((char *)NULL, &sdsc_info[isds]);
      }
    }
    if (st == -1) break;
  } /* for (isds=0; . . .  ) */
  if ((st != -1) && (l2g_st == 1))
  {
    strcpy(sds_nobs_info->name, "num_observations");
    sds_nobs_info->sd_id = sds_info[0].sd_id;
    sds_nobs_info->sds_id = -1;
    get_sds_info((char *)NULL, sds_nobs_info);
  }
  return st;
}

int create_out_sds(sds_t *in_sds_info, sds_t *out_sds_info, int nsds, char *of_str, char *m_str, int *n, int *m,
		int32 out_sd_id, int out_hdf_st, int *mask_fill)
{
  int st = 1;
  char sds_name[MAX_SDS_NAME_LEN];
  int i, j, p1, isds, rank, len, len_m_str;

  len_m_str = (int)strlen(m_str);
  for (isds=0; isds<nsds; isds++)
  {
    get_sdsname_dim(in_sds_info[isds].name, sds_name, &n[isds], &m[isds]);
    if (out_hdf_st == 1)
    {
      rank = in_sds_info[isds].rank;
      if ((p1 = sd_strpos(sds_name, "_1.", 0)) != -1)
      {
        len = (int)strlen(sds_name);
        for (j=p1+2, i=p1; j<=len; j++, i++)
          sds_name[i] = sds_name[j];
      }
      else strcpy(sds_name, in_sds_info[isds].name);
      out_sds_info[isds].rank = ((n[isds] == -1) && (m[isds] == -1)) ? rank : 2;
      if (rank == out_sds_info[isds].rank)
        for (j=0; j<rank; j++)
          out_sds_info[isds].dim_size[j] = in_sds_info[isds].dim_size[j];
      else
      {
        if (in_sds_info[isds].dim_size[0] < in_sds_info[isds].dim_size[rank-1])
        {
          out_sds_info[isds].dim_size[0] = in_sds_info[isds].dim_size[rank-2];
          out_sds_info[isds].dim_size[1] = in_sds_info[isds].dim_size[rank-1];
        }
        else
        {
          out_sds_info[isds].dim_size[0] = in_sds_info[isds].dim_size[0];
          out_sds_info[isds].dim_size[1] = in_sds_info[isds].dim_size[1];
        }
      }
      out_sds_info[isds].data_type = in_sds_info[isds].data_type;
      out_sds_info[isds].data_size = in_sds_info[isds].data_size;
      if (strlen(of_str) <= 0)
        sprintf(out_sds_info[isds].name, "%s", sds_name);
      else {
	if (strstr(of_str, "VI") != NULL)
          sprintf(out_sds_info[isds].name, "%s %s", sds_name, of_str);
        else
          sprintf(out_sds_info[isds].name, "%s%s", sds_name, of_str);
      }
      out_sds_info[isds].sds_id = -1;
      if (out_sd_id != -1)
      {
        out_sds_info[isds].sd_id = out_sd_id;
        if (open_sds((char *)NULL, &out_sds_info[isds], 'W') != -1)
        {
          write_attr_fval(out_sds_info[isds].sds_id, in_sds_info[isds].data_type, 1, in_sds_info[isds].fill_val, ATTR_FILL_NAME);
          write_attr_fval(out_sds_info[isds].sds_id, in_sds_info[isds].data_type, 1, mask_fill[isds], MASK_FILL_NAME);
          if (SDsetattr(out_sds_info[isds].sds_id, "Mask_String", DFNT_CHAR8, len_m_str, (VOIDP)m_str) == FAIL)
            fprintf(stderr, "Cannot write attribute to output SDS in mask_nsds()\n");
        }
        else st = -1;
      }
    } /* if (out_hdf_st == 1)  */
  } /* for (isds=0; . . . )  */
  return st;
}

int get_res_factors(sds_t *sds_info, sds_t *qa_sds_info, int n_op, int *res_l, int *res_s)
{
  int status;
  int rank, i_op;
  int xdim, ydim;

  status = 1;
  rank = sds_info->rank;
  if ((rank==2) || (sds_info->dim_size[0]>sds_info->dim_size[rank-1]))
  {
    xdim = sds_info->dim_size[0];
    ydim = sds_info->dim_size[1];
  }
  else
  {
    xdim=sds_info->dim_size[rank-2];
    ydim=sds_info->dim_size[rank-1];
  }
  for (i_op=0; i_op<=n_op; i_op++)
  {
    rank = qa_sds_info[i_op].rank;
    if ((rank==2) || (qa_sds_info[i_op].dim_size[0]>qa_sds_info[i_op].dim_size[rank-1]))
    {
      res_l[i_op] = xdim/qa_sds_info[i_op].dim_size[0];
      res_s[i_op] = ydim/qa_sds_info[i_op].dim_size[1];
    }
    else
    {
      res_l[i_op] = xdim/qa_sds_info[i_op].dim_size[rank-2];
      res_s[i_op] = ydim/qa_sds_info[i_op].dim_size[rank-1];
    }
    if ((res_s[i_op] < 1) || (res_l[i_op] < 1))
    {
      fprintf(stderr, "Masking sds have higher resolution: mask_nsds()\n");
      status = -1; break;
    }
  }  /* for (i_op = 0; . .  ) */
  return status;
}

int compute_res_factors(sds_t sds_info, sds_t *qa_sds_info, int n_op, int *res_l, int *res_s)
{
  int status;
  int rank, i_op;
  int xdim, ydim;

  status = 1;
  rank = sds_info.rank;
  if ((rank==2) || (sds_info.dim_size[0]>sds_info.dim_size[rank-1]))
  {
    xdim = sds_info.dim_size[0];
    ydim = sds_info.dim_size[1];
  }
  else
  {
    xdim=sds_info.dim_size[rank-2];
    ydim=sds_info.dim_size[rank-1];
  }
  for (i_op=0; i_op<=n_op; i_op++)
  {
    rank = qa_sds_info[i_op].rank;
    if ((rank==2) || (qa_sds_info[i_op].dim_size[0]>qa_sds_info[i_op].dim_size[rank-1]))
    {
      res_l[i_op] = xdim/qa_sds_info[i_op].dim_size[0];
      res_s[i_op] = ydim/qa_sds_info[i_op].dim_size[1];
    }
    else
    {
      res_l[i_op] = xdim/qa_sds_info[i_op].dim_size[rank-2];
      res_s[i_op] = ydim/qa_sds_info[i_op].dim_size[rank-1];
    }
    if ((res_s[i_op] < 1) || (res_l[i_op] < 1))
    {
      fprintf(stderr, "Masking sds have higher resolution: mask_nsds()\n");
      status = -1; break;
    }
  }  /* for (i_op = 0; . .  ) */
  return status;
}

void get_ndata_vals(sds_t *sds_info, int *bsq, int *nrow, int *ndata_in, int *ndata_mask, int *ndata_out,
		int n0, int m0)
{
  int i, rank;

  rank = sds_info->rank;
  *bsq = ((rank == 2) || (sds_info->dim_size[0] < sds_info->dim_size[rank-1])) ? 1 : 0;
  if (rank == 2)
  {
    *ndata_mask = *ndata_in = *ndata_out = sds_info->dim_size[1];
    *nrow = sds_info->dim_size[0];
  }
  else
  {
    if (sds_info->dim_size[0] > sds_info->dim_size[rank-1])
    {
      *nrow = sds_info->dim_size[0];
      *ndata_mask = *ndata_in = sds_info->dim_size[1];
      for (i=2; i<rank; i++)
        *ndata_in *= sds_info->dim_size[i];
      if ((n0 == -1) && (m0 == -1)) *ndata_out = *ndata_in;
        else *ndata_out = sds_info->dim_size[1];
    }
    else
    {
      *nrow = sds_info->dim_size[rank-2];
      *ndata_mask = *ndata_in = sds_info->dim_size[rank-1];
      for (i=0; i<rank-2; i++)
        *ndata_in *= sds_info->dim_size[i];
      if ((n0 == -1) && (m0 == -1)) *ndata_out = *ndata_in;
        else *ndata_out = sds_info->dim_size[rank-1];
    }
  }
}

int conv_date(int *mm, int *dd, int yyyy)
/*
   Note: *mm and *dd are input/output varialbes
   Input julian day number:
        input *mm = 0;
        input *dd = julian day
        output in mm-dd-yyyy format
   Input in mm-dd-yyyy format
        output *mm = 0;
        output *dd = julian day number;
*/
{
  int nm, im;
  int st = -1;
  int ndays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if ((yyyy%400 == 0) || ((yyyy%4 == 0) && (yyyy%100 != 0))) ndays[1] = 29;
  if (*mm == 0)
  {
    if (*dd <= 0)
      fprintf(stderr, "Error in input julian date: %d %d\n", *dd, yyyy);
    else
    {
      for (im=0; ((im<12)&&(*dd>0)); im++)
        *dd -= ndays[im];
      if ((im > 12) || ((im == 12) && (*dd > 0)))
        fprintf(stderr, "Error in input julian date: %d %d\n", *dd, yyyy);
      else
      {
        *mm = im;
        *dd += ndays[*mm - 1];
        st = 1;
      }
    }
  }
  else
  {
    if ((*mm <= 0) || (*dd <= 0))
      fprintf(stderr, "Error in input date: %d %d %d\n", *mm, *dd, yyyy);
    else
    {
      nm = *mm - 1;
      for (im=0; im<nm; im++)
        *dd += ndays[im];
      *mm = 0;
      st = 1;
    }
  }
  return st;
}                                            
