/****************************************************************************
!C

!File: mask_l3_sds.h

!Description:

  This file contains header file of routines for masking SDS.

!Input Parameters: (none)

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
    LDOPE                             University of Maryland
                                      Department of Geography
    droy@kratmos.gsfc.nasa.gov        NASA/GSFC Code 922 (B32)
    phone: 301-614-5571               Greenbelt, MD 20771   

!Design Notes: (none)

    Uses library routines defined in mask_sds_lib.h

!END
*****************************************************************************/

int get_mask_string(char *m_str, char **arg_mask_str, int *val_opt, int *l2g_st);
int check_fsds_bit_str_val(char *fname, char *sname, char *bit_str, int *opt, 
			   int *l2g_st);
int get_parameters(char **arg_list, int n_op, int *sel_qa_op, char **qa_fnames,
		   sds_t *qa_sds_info, unsigned long *bit_mask_arr, 
		   unsigned long *mask_val_arr, int *opt_arr, int *rel_op);
int get_file_sds_names(char *fsds_name, char *fname, char *sds_name);
int get_bit_num_val(char *in_str, unsigned long *bit_mask, unsigned long *mask_val,
		    int opt, int *rop);
int open_qa_sds(char *fname, sds_t *sds_info, char **qa_fnames, sds_t *qa_sds_info, 
		int n_op);
int open_qa_sds_nsds(char *fname, sds_t *sds_info, sds_t *sdsc_info, 
		     sds_t *sds_nobs_info, int nsds, char **qa_fnames, 
		     sds_t *qa_sds_info, sds_t *qa_sdsc_info, 
		     sds_t *qa_sds_nobs_info, int *qa_l2g, int n_op);
int malloc_qa_sds(sds_t *qa_sds_info, int n_op, int *fqa_l2g, void **data_qa,
		  int32 **data_qa_nadd);
void read_qa_sds(sds_t *qa_sds_info, sds_t *qa_sdsc_info, sds_t *qa_sds_nobs_info, 
		 int n_op, void **data_qa, int32 **data_qa_nadd, int irow, int *res_l, 
		 int *fqa_l2g, int *obs_num);
void read_sdsc_data(sds_t *sdsc_info, sds_t *sds_nobs_info, void *data, 
		    int32 *data_nadd, int irow, int nobs);
void close_qa_hdf(char *hdf_fname, sds_t *sds_info, char **qa_fnames,
		  sds_t *qa_sds_info, int n_op);
void close_qa_hdf_nsds(char *hdf_fname, sds_t *sds_info, int nsds, char **qa_fnames,
		       sds_t *qa_sds_info, int n_op);
void process_mask_data(void **data_qa, int ncols, sds_t *qa_sds_info, int n_op,
		       int *sel_qa_op, unsigned long *bit_mask_arr, 
		       unsigned long *mask_val_arr, int *rel_op, int *res_s, 
		       uint8 *mask_row, int on_val, int off_val, int mask_fill);
int get_qa_sds_info(char **fnames, sds_t *sds_info, sds_t *sdsc_info, int *l2g_st, 
		    int n_op);
int get_in_sds_info(char *hdf_fname, sds_t *sds_info, sds_t *sdsc_info, 
		    sds_t *sds_nobs_info,
        int l2g_st, int nsds);
int create_out_sds(sds_t *in_sds_info, sds_t *out_sds_info, int nsds, char *of_str, 
	char *m_str, int *n, int *m, int32 out_sd_id, int out_hdf_st, int *mask_fill);
int get_res_factors(sds_t *sds_info, sds_t *qa_sds_info, int n_op, int *res_l, 
		    int *res_s);
void get_ndata_vals(sds_t *sds_info, int *bsq, int *nrow, int *ndata_in, 
		    int *ndata_mask, int *ndata_out, int n0, int m0);
int conv_date(int *mm, int *dd, int yyyy);
int compute_res_factors(sds_t sds_info, sds_t *qa_sds_info, int n_op, int *res_l, int *res_s);
