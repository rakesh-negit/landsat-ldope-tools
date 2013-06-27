/****************************************************************************
!C

!File: alloc_mem.c

!Description:
  Contains routines for allocating and freeing 2D and 3D memory.

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original January 1998. Version 1.0 Copied from MODLAND product source code
  library
  Modified January 1999
  Modified March 2000:  Added function "alloc_whole_sds" which allocates a
                        bulk storage block for an sds.  This block can be
                        used to copy and entire sds from one file to another.
                        (Curt Crandall   QSS   03/30/2000)

!Team-unique Header:
  This software was developed by:
    Land Data Operational Product Evaluation (LDOPE) Team for the
    National Aeronautics and Space Administration, Goddard Space Flight
    Center.

!References and Credits:

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)

!END
*****************************************************************************/
#include <stdio.h>
#include "mfhdf.h"
#include "alloc_mem.h"

void ***Calloc3D(size_t nobj1, size_t nobj2, size_t nobj3, size_t size)
{
  void ***p1;
  void **p2;
  void *p3;
  char *c3;
  size_t iobj1, iobj2;

  p3 = calloc((nobj1 *nobj2 * nobj3), size);
  if (p3 == NULL) return NULL;

  p1 = (void ***) Calloc2D(nobj1, nobj2, sizeof(void *));
  if (p1 == NULL) {
    free(p3);
    return NULL;
  }

  c3 = (char *) p3;
  for (iobj1 = 0; iobj1 < nobj1; iobj1++) {
    p2 = p1[iobj1];
    for (iobj2=0; iobj2<nobj2; iobj2++) {
      p2[iobj2] = (void *) c3;
      c3 += (nobj3 * size);
    }
  }

  return p1;
}

void **Calloc2D(size_t nobj1, size_t nobj2, size_t size)
{
  void **p1;
  void *p2;
  size_t iobj1;
  char *c2;

  p2 = calloc((nobj1 * nobj2), size);
  if (p2 == NULL) return NULL;

  p1 = (void **) calloc(nobj1, sizeof(void *));
  if (p1 == NULL) {
    free(p2);
    return NULL;
  }

  c2 = (char *) p2;
  for (iobj1 = 0; iobj1 < nobj1; iobj1++) {
    p1[iobj1] = (void *) c2;
    c2 += (nobj2 * size);
  }

  return p1;
}

void Free2D(void **p1)
{
  if (*p1 != NULL) free(*p1);
  if (p1 != NULL) free(p1);
  return;
}

void Free3D(void ***p1)
{
  if (**p1 != NULL) free(**p1);
  if (*p1 != NULL) free(*p1);
  free(p1);
  return;
}


/******************************************************************************
 *
 *    Function:  alloc_whole_sds
 *
 *    Description:  This function allocates memory for any valid HDF data
 *                  type.  This includes 128-bit numeric and 16-bit character
 *                  data types that are not currently supported.
 *
 *    Programmer:  Curt Crandall (QSS)
 *
 *    Prototype:  void *alloc_whole_sds(int32, void *, int, char *)
 *
 *    Arguments:  data_type -> HDF data type from sds structure
 *                data_in   -> array that will contain an SDS
 *                ndata_in  -> number of elements in data_in
 *                name      -> SDS name from sds structure
 *
 *    Development History:
 *       Version 1.0     Curt Crandall (QSS)     03/16/2000
 *       -Original Implementation.
 *
 *    Notes:
 *
 *****************************************************************************/

void *alloc_whole_sds(int32 data_type, int ndata_in, char *name)
{
   void *data_in = NULL;

   switch (data_type)
   {
     case 3:  if((data_in = (uchar8 *)calloc(ndata_in, sizeof(uchar8))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 4:  if((data_in = (char8 *)calloc(ndata_in, sizeof(char8))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 5:  if((data_in = (float32 *)calloc(ndata_in, sizeof(float32))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 6:  if((data_in = (float64 *)calloc(ndata_in, sizeof(float64))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 7:  fprintf(stderr, "The float128 data type is not supported\n");
              /*if((data_in = (float128 *)calloc(ndata_in, sizeof(float128))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     case 20: if((data_in = (int8 *)calloc(ndata_in, sizeof(int8))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 21: if((data_in = (uint8 *)calloc(ndata_in, sizeof(uint8))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 22: if((data_in = (int16 *)calloc(ndata_in, sizeof(int16))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 23: if((data_in = (uint16 *)calloc(ndata_in, sizeof(uint16))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 24: if((data_in = (int32 *)calloc(ndata_in, sizeof(int32))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 25: if((data_in = (uint32 *)calloc(ndata_in, sizeof(uint32))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }
              break;

     case 26: fprintf(stderr, "The int64 data type is not supported\n");
              /*if((data_in = (int64 *)calloc(ndata_in, sizeof(int64))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     case 27: fprintf(stderr, "The uint64 data type is not supported\n");
              /*if((data_in = (uint64 *)calloc(ndata_in, sizeof(uint64))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     case 28: fprintf(stderr, "The int128 data type is not supported\n");
              /*if((data_in = (int128 *)calloc(ndata_in, sizeof(int128))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     case 30: fprintf(stderr, "The uint128 data type is not supported\n");
              /*if((data_in = (uint128 *)calloc(ndata_in, sizeof(uint128))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     case 42: fprintf(stderr, "The char16 data type is not supported\n");
              /*if((data_in = (char16 *)calloc(ndata_in, sizeof(char16))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     case 43: fprintf(stderr, "The uchar16 data type is not supported\n");
              /*if((data_in = (uchar16 *)calloc(ndata_in, sizeof(uchar16))) == NULL)
              {
                 fprintf(stderr, "Cannot allocate memory for %s\n", name);
              }*/
              break;

     default: fprintf(stderr, "Data type is invalid\n");
              break;
   }

   return data_in;
}

