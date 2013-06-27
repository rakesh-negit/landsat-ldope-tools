#include <stdio.h>

#ifndef _MAIN_UTIL_H_
#define _MAIN_UTIL_H_

int get_prod_sname(char *hdf_fname, char *prod_sname);
int get_pnt_fname(char *hdf_fname, char **modis_fpath, char *ptr_fname, char *pnt_fname);
int find_file(char **fpath, char *fname);
int get_input_files(char *fname, char *iptr_id, char **iptr_fnames);
int is_arg_id(char *arg_str, char *arg_id);
void get_arg_val_arr(char *arg_str, char **arg_val, int *arg_cnt);
void get_arg_val(char *arg_str, char *arg_val);
int check_bit_str(char *fname, char *sname, char *bn_str);
void get_bit_num_arr(char *bn_str, int *cnt, int *bn_arr);
void rm_path(char *fname);
int parse_stdin(int *argc, char **argv);
void get_qa_tool_env(char *env_var, char **env_val);
int get_day_time_tile_info(char *gran_id, char *day, char *time_tile);
int get_line(FILE *fp, char *s);
int get_numbers(char **num_str, int cnt, int *sel_num, int max, char *msg);
void sort_values(int *values, int n);
void sort_fvalues(float *values, int n);
int read_clr_table(char *in_fname, int *clr_table);
void get_esdt_tileid(char *fname, char *esdt, char *tile_id, char *jday, char *ver_id);

#endif
