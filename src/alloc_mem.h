#ifndef _ALLOC_MEM_H_
#define _ALLOC_MEM_H_
void **Calloc2D(size_t nobj1, size_t nobj2, size_t size);
void ***Calloc3D(size_t nobj1, size_t nobj2, size_t nobj3, size_t size);
void Free2D(void **p1);
void Free3D(void ***p1);
void *alloc_whole_sds(int32 data_type, int ndata_in, char *name);

#endif
