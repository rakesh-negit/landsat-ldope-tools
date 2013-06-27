/****************************************************************************
!C

!File: str_op.c

!Description:
  Contains routines for various string manipulation. 

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original January/February/March 1998. Version 1.0
  Modified September 1999. Fixed a bug in strcmp_wc 
  Modified April 2000. get_sdsname_dim: (if sds name contains '(' character consider this SDS as 2D)

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
#include <ctype.h>

#include "isoc.h"
#include "qa_tool.h"
#include "str_op.h"

int sd_charpos(char *s, char c, int p)
/****************************************************************************

!Description:
   The charpos functions finds the first occurance of a character with an 
   object string starting from position p, it returns the character position
   of the match, otherwise it returns -1 

!Input Parameters:
   char c    - character to be matched
   char *s   - target object string 
   int p     - search starting position

!Output Parameters: 
   returns the character position of the match. -l if character can't be
   found in the target object string.
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

*****************************************************************************/
{
  int len, i;

  len = (int)strlen(s);
  for (i=p; i<len; i++)
    if (s[i] == c) return i;
  return -1;
}

int sd_strpos(char *s1, char *s2, int p)

/**************************************************************************** 
!Description:
   The strpos functions finds the first occurance of a substring s2 within 
   an object string s1 starting from position p, it returns the character 
   postion of the match, otherwise it returns -1 

!Input Parameters:
   char *s1   - target object string
   char *s2   - string to be matched
   int p      - search starting position

!Output Parameters: 
   returns the substring position of the match. -l if substring can't be
   found in the target object string.
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

*****************************************************************************/
{
  int i = 0, j, k, fn;
  int len1, len2;
  char tmp[MAX_STR_LEN];

  fn = 0;
  len1 = (int)strlen(s1);
  len2 = (int)strlen(s2);
  if (len1 >= len2)
    for (i = p; i <= len1-len2; i++)
    {
      if (s1[i] == s2[0]) 
      {
        for (j = i, k = 0; k < len2; k++, j++)
         tmp[k] = s1[j];
        tmp[k] = '\0'; 
        if (strcmp(tmp, s2) == 0) 
        {
          fn = 1;
          break;
        }
      }
    }
  if (fn == 1) return i;
    else return -1; 
}

int sd_strcasepos(char *s1, char *s2, int p)

/**************************************************************************** 
!Description:
   The strcasepos functions finds the first case insensitive occurance of 
   a substring s2 within an object string s1 starting from position p, it 
   returns the character postion of the match, otherwise it returns -1 

!Input Parameters:
   char *s1   - target object string
   char *s2   - string to be matched
   int p      - search starting position

!Output Parameters: 
   returns the substring position of the match. -l if substring can't be
   found in the target object string.
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

*****************************************************************************/

{
  int i = 0, j, k, fn;
  int len1, len2;
  char tmp[MAX_STR_LEN];

  fn = 0;
  len1 = (int)strlen(s1);
  len2 = (int)strlen(s2);
  if (len1 >= len2)
    for (i = p; i <= len1-len2; i++)
    {
    if( tolower(s1[i]) == tolower(s2[0] ) )
/*
      if ((s1[i] == s2[0]) || (s1[i] == tolower(s2[0])) || (s2[0] == tolower(s1[i]))) 
*/
      {
        for (j = i, k = 0; k < len2; k++, j++)
         tmp[k] = s1[j];
        tmp[k] = '\0'; 
        if (strcasecmp(tmp, s2) == 0) 
        {
          fn = 1;
          break;
        }
      }
    }
  if (fn == 1) return i;
    else return -1; 
}

void sd_strmid(char *s1, int p1, int cnt, char *s2)

/**************************************************************************** 
!Description:
   The strmid function extracts a substring from a string expression. The 
   result of the function is a string s2 of Length cnt taken from the 
   object string s1, starting at character position p1.

!Input Parameters:
   char *s1   - target object string
   int p1     - String subsetting start point.
   int cnt    - substring length.
!Output Parameters: 
   char *s2   - resulting substring of length cnt.
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)


****************************************************************************/

{
  int i, j, p2;

  p2 = p1 + cnt;
  for (i=p1, j=0; i<p2; i++, j++)
    s2[j] = s1[i];
  s2[j] = '\0';
}

void sd_strtrim(char *s)

/**************************************************************************** 
!Description:
    The strtrim function returns a copy of String s with leading and/or 
     trailing blanks removed. 

!Input Parameters:
   char *s   - target object string
 
!Output Parameters:  (none)
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

****************************************************************************/

{
  int i, j1, j2, k, len;

  len = (int)strlen(s);
  /* remove leading blanks */
  for (j1=0; j1<len; j1++)
    if ((s[j1] != ' ') && (s[j1] != '\t') && (s[j1] != '\n')) break;
  /* remove tailing blanks */
  for (j2=len-1; j2>=0; j2--)
    if ((s[j2] != ' ') && (s[j2] != '\t') && (s[j1] != '\n')) break;
  /* repack string */
  for (i=j1, k=0; i<=j2; i++, k++)
    s[k] = s[i];
  s[k] = '\0';
}

void sd_strrev(char *s)
     /* reverse the input string s. */
{
  int i, k, len;
  char s1[MAX_STR_LEN];

  len = (int)strlen(s);
  strcpy(s1, s);
  for (i=len-1, k=0; i>=0; i--, k++)
    s[k] = s1[i];
  s[k] = '\0';
}

void sd_rm_ln_in_str(char *s)

/**************************************************************************** 
 
!Description:
   The rm_ln_in_str functions returns a copy of the string s with newline 
   character removed.   

!Input Parameters:
   char *s   - target object s
   
!Output Parameters: 
   char *s   - string s with newline character removed.
    
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

****************************************************************************/

{
  int i, j;
  int len;

  len = (int)strlen(s);
  for (i=0; i<len; i++)
    /* remove newline character */
    if (s[i] == '\n')
    {
      for (j=i+1; j<len; j++)
        s[j-1] = s[j];
      s[j-1] = '\0';
      len--; 
    }
}

int sd_strcmp_wc(char *s1, char *s2)
/* searches for s2 in s1. s1 must be complete */
{
  char s_tmp[40];
  int p_s1, p_s2;
  int pt1, pt2, pt;
  int len_s1, len_s2;
  int end_len, l_tmp;
  char end_str_s1[25];
  char end_str_s2[25];

  len_s1 = (int)strlen(s1);
  len_s2 = (int)strlen(s2);
  for (p_s1=0, p_s2=0; ((p_s1 < len_s1) && (p_s2 < len_s2)); )
  {
    if ((s2[p_s2] == s1[p_s1]) || (s2[p_s2] == '?'))
    {
      p_s1++; p_s2++;
    }
    else if (s2[p_s2] == '*')
    {
      p_s2++;
      pt1 = sd_charpos(s2, '*', p_s2);
      pt2 = sd_charpos(s2, '?', p_s2);
      if ((pt1 == -1) && (pt2 == -1))
      {
	p_s1++;
	end_len = len_s2 - p_s2;
	sd_strmid(s2, p_s2, end_len, end_str_s2);
	sd_strmid(s1, len_s1-end_len, end_len, end_str_s1);
	if (strcmp(end_str_s1, end_str_s2) == 0)
	{
	  p_s1 = len_s1;
	  p_s2 = len_s2;
	}
	break;
      }
      else
      {
        if (pt1 == -1) pt = pt2;
        else if (pt2 == -1) pt = pt1;
        else pt = (pt1 > pt2) ? pt2 : pt1;
        sd_strmid(s2, p_s2, pt-p_s2, s_tmp);
        p_s1 = sd_strpos(s1, s_tmp, p_s1);
        if (p_s1 == -1) break;
        else
        {
          l_tmp = (int)strlen(s_tmp);
          p_s1 += l_tmp;
          p_s2 += l_tmp;
        }
      }
    }
    else break;
  }
  if ((p_s1 == len_s1) && (p_s2 == len_s2)) return 0;
  else return -1;
}

void sd_sort_strings(char **s, int n)
{
  int i,j;
  char tmp[MAX_PATH_LENGTH];

  for (i=0; i<n-1; i++)
    for (j=n-1; j>i; j--)
      if (strcmp(s[j], s[j-1]) < 0)
      {
	strcpy(tmp, s[j]);
	strcpy(s[j], s[j-1]);
	strcpy(s[j-1], tmp);
      }    
  return;
}

int get_sdsname_dim(char *sdsname_str, char *sds_name, int *n, int *m)

/**************************************************************************** 
!Description:
   The get_sdsname_dim functions returns the information of sds name and 
   dimentsion. 

!Input Parameters:
   sdsname_str   - input sds data structure

!Input/output Parameter:
   sds_name    - sds_name abstract from inpu sds data structure.
   n           - layer number for 3rd dimension
   m           - layer number for 4th dimenstion

!Output Parameters: 
   
   1           if SDS dimension information is retrieved.
   -l          if invalid layer number in SDS name.
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

****************************************************************************/

{
  int st = 1;
  int p1, p2, len;
  char nstr[10];

  *n = *m = -1;
  if (sd_charpos(sdsname_str, '(', 0) == -1)
  {
      /* 2D sds */
      if ((p1 = sd_charpos(sdsname_str, '.', 0)) == -1)
      strcpy(sds_name, sdsname_str);
    else
    {
      len = (int)strlen(sdsname_str);
      sd_strmid(sdsname_str, 0, p1, sds_name);
      p1++;
      /* 3D sds */
      if ((p2 = sd_charpos(sdsname_str, '.', p1)) == -1)
      {
        sd_strmid(sdsname_str, p1, len-p1, nstr);
        *n = (int)atoi(nstr);
      }
      else
      /* 4D sds */
      {
        sd_strmid(sdsname_str, p1, p2-p1, nstr);
        *n = (int)atoi(nstr);
        p2++;
        sd_strmid(sdsname_str, p2, len-p2, nstr);
        *m = (int)atoi(nstr);
      }
    }
    /* invalid number in SDS name */
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
  }
  return st;
}

int sd_getline(FILE *fp, char *s)
     /* read in a line from the input file. i is the length of the line */
{
  int c, i;

  c = fgetc(fp);
  for (i=0; i<MAX_LINE_LENGTH && c != EOF && c != '\n'; ++i)
  {
    s[i] = (char)c;
    c = fgetc(fp);
  }
  s[i] = '\0';
  return i;
}

void sd_split_string(char *tmp1_str, char **tmp_str, int *num)
{
  int p1, p2, len;

  p1 =  0;
  p2 = sd_charpos(tmp1_str, ',', p1);
  while ((p2 != -1) && (*num < MAX_NUM_PARAM))
  {
    sd_strmid(tmp1_str, p1, p2-p1, tmp_str[*num]);
    ++*num;
    p1 = p2 + 1;
    p2 = sd_charpos(tmp1_str, ',', p1);
  }
  if (*num >= MAX_NUM_PARAM)
  {
    fprintf(stderr, "Too many parameters in option %s\n", tmp1_str);
    fprintf(stderr, "Considering only %d number of parameter values\n", MAX_NUM_PARAM);
  }
  else
  {
    len = (int)strlen(tmp1_str);
    if (p1 < len)
    {
      sd_strmid(tmp1_str, p1, len-p1, tmp_str[*num]);
      ++*num;
    }
  }
}

char *sd_concat(const char *stra, const char *strb) 
   
/**************************************************************************** 
!Description:
   This function append string b to string a. 

!Input Parameters:
   char *stra   - String to be appended.
   char *strb   - target string.

!Return: 
   char *s   - resulting string as "ab".
   
!Revision History: (see file prolog)
 
!Team-unique Header: (see file prolog)
 
!References and Credits: (see file prolog)

!Design Notes: (none)

****************************************************************************/
{
  size_t len = 0;
  char *ret = NULL;
  
  if(!stra) 
    {
      if(!strb) return NULL;
      ret = (char *)malloc(sizeof(char)*(strlen(strb) + 1));
      strcpy(ret,strb);
      return ret;
    }
  
  len = strlen(stra);
  
  if(strb) 
    {
      len += strlen(strb) + 1;
    }

ret = (char*)malloc(sizeof(char) * len);
strcpy(ret,stra);

if(strb) 
  {
    strcat(ret,strb); 
  }
 
 return ret;
}

char *sd_remove_chars(const char *in, const char *char_list)
{
	char *ret = NULL;
	size_t idx = 0;
	size_t i = 0;
	size_t slen = strlen(in);

	ret = (char*)malloc(slen+1);
	memset(ret, 0, slen+1);
	if(ret == NULL)
	{
		fprintf(stderr, "Error: allocating memory in remove_chars");
		exit(-1);
	}
	for(idx = 0; idx < slen; idx++)
	{
		if(strchr(char_list, in[idx]) == NULL)
		{
			ret[i] = in[idx];
			i++;
		}
	}	
	return ret;
}
