#include "hdf.h"

#ifndef _META_HDFEOS_T_
#define _META_HDFEOS_T_

typedef struct {
  int32 ns, nl;
  float64 ulCornerGrid[2], lrCornerGrid[2];
  int32 projCode;
  int32 sphereCode;
  float64 projParameters[13];
  int32 zoneCode;
  int32 originCode;
} meta_hdfeos_t;

char *get_attr_metadata(char *hdf_fname, char *meta_gname);
void get_all_metadata(char *meta_str, char **meta_name, char **meta_val, int *meta_cnt);
void get_sel_metadata(char *meta_str, char *meta_name, char **meta_val, int *meta_cnt, int case_ch);
void copy_metadata(int32 in_sd_id, int32 out_sd_id);
void update_gid_dt(char *loc_gid, char *prod_dt);
void write_modss_metadata(int32 in_sd_id, int32 out_sd_id, char **meta_names, char **meta_vals);
char *update_modss_metadata(char *attr_buf, int org_len, int imeta, char **meta_names, char **meta_vals);

#endif
