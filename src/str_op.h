#include <stdio.h>
#ifndef _STR_OP_H_
#define _STR_OP_H_

/* Some of the functions below are named the same as what may be found on
   some operating systems.  An "sd_" prefix has been added to most functions
   defined in this module just to be consistent, although not all of them
   needed to be.  The "sd_" prefix eliminates collisions with functions that
   were named the same.
*/

int sd_charpos(char *s, char c, int p);
int sd_strpos(char *s1, char *s2, int p);
int sd_strcasepos(char *s1, char *s2, int p);
void sd_strmid(char *s1, int p1, int cnt, char *s2);
void sd_strtrim(char *s);
void sd_strrev(char *s);
void sd_rm_ln_in_str(char *s);
int sd_strcmp_wc(char *s1, char *s2);
void sd_sort_strings(char **s, int n);
int get_sdsname_dim(char *sdsname_str, char *sds_name, int *n, int *m);
int sd_getline(FILE *fp, char *s);
void sd_split_string(char *tmp1_str, char **tmp_str, int *tmp_num);

#ifdef WIN32
# define CSDI_PATH_SEP_STR   "\\"
# define CSDI_PATHS_SEP_CHAR ';'
#else
# define CSDI_PATH_SEP_STR   "/"
# define CSDI_PATHS_SEP_CHAR ':'
#endif

char *sd_concat(const char *stra, const char *strb);

/* 
 * remove_chars will remove any character listed by char_list from the in string
**/
char * sd_remove_chars(const char *in, const char * char_list);

/* Used to denote and quiet compiler warnings for unused parameters. */

#define LDOPE_UNUSED_PARAMETER(x) ((x) = (x))

#endif
