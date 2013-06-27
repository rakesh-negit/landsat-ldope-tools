#include "hdf.h"
#include "qa_tool.h"

#ifndef _SDS_TYPES_H_
#define _SDS_TYPES_H_
typedef struct
{
  char name[MAX_SDS_NAME_LEN];
  int32 data_type;
  int32 data_size;
  int32 rank; 
  int32 dim_size[4];
  int32 nattr;
  long fill_val;
  float32 fill_fval;
  int range[2];
  float32 frange[2];
  int32 sd_id;
  int32 sds_id;
  int32 sds_index;
} sds_t;

typedef struct
{
  char name[80];
  int32 type;
  int32 cnt;
  void *buf;
} sds_attr_t;

#endif
