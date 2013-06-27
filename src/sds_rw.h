#include "sds_types.h"

#ifndef _SDS_RW_H_
#define _SDS_RW_H_

char *get_attr_metadata(char *in_fname, char *meta_str);
int get_sds_info(char *hdf_fname, sds_t *sds_info);
int get_sds_data(sds_t *sds, void *data);
void *get_sel_sds_data(sds_t *sds, int32 *start, int32 *edge);
int open_sds(char *fname, sds_t *sds_info, char open_t);
void close_hdf(sds_t *sds_info);
int get_l2g_sds_names(char *fname, char **sds_names);
int get_sds_names(char *fname, char **sds_names);
void *get_sds_attr(int32 sds_id, char *attr_name, int32 *attr_type, int32 *attr_cnt);
void write_attr_fval(int32 sds_id, int32 fval_type, int c, int attr_val, char *attr_name);
void write_sds_attrs(int32 in_sds_id, int32 out_sds_id, int bn);
void write_all_sds_attrs(int32 in_sds_id, int32 out_sds_id, int32 nattr);
void write_metadata(int32 in_sd_id, int32 out_sd_id);
void compute_sds_start_offset(sds_t *sds_info, int n, int m, int *st_c, int *offset);
void update_nd_sdsnames(char **sds_names, int *sds_cnt, char *fname);
void update_l2g_sdsnames(char **sds_names, int *nsds, char *fname, int nobs);
void create_names(char *fname, char *sds_name, char *nd_ext, char *md_ext, int *nsds, char **sds_names);
void get_dim_num(char *str_arr, int *num_arr, int dim_size, int *cnt);
void copy_sds(sds_t *sds_info, int32 out_sd_id, int compress);
void check_and_fix_sdsname(char *hdf_fname, sds_t *sds_info);
int compute_sds_ndata(sds_t *sds_info);
void get_sds_edge(sds_t *sds_info, int32 *edge);
void compute_sds_nrows_ncols(sds_t *sds_info, int *nrows, int *ncols);
void get_sds_param(sds_t *sds_info, int *n, int *m, int *rank, int *dim_size);
void print_sds_dim_size(sds_t *sds_info);
void get_sds_dim_name(sds_t *sds_info, char **dim_names, char **short_dim_names);
void display_sds_info_of_file(char* filename);

int open_l2g_nobs_sds(sds_t *nobs_sds_info, sds_t *nadd_obs_sds_info, sds_t *sds_info);

#endif
