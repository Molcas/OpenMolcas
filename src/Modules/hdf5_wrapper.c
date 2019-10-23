/*
 * hdf5_qcm: HDF5 Fortran/C interface for QCMaquis driver and CD-NEVPT2
 * This interface replaces the previous Fortran module-based HDF5 F2003 interface and does not require the
 * compilation of the HDF5 library with the same compiler version as the QCMaquis driver/NEVPT2.
 *
 * Fortran-C interoperability inspired by http://fortranwiki.org/fortran/show/Generating+C+Interfaces
 * and Steven Vancoillie's hdf5_util in OpenMOLCAS
 * (C) 2017 Leon Freitag and Stefan Knecht
 *        Laboratory for Physical Chemistry, ETH Zurich
 *
 * hdf5_qcm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hdf5_qcm is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with hdf5_qcm. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef _HDF5_
#include "hdf5.h"


/* Type definitons */
#define MYH5T_INT 1
#define MYH5T_DOUBLE 2

#define H5T_READ_INT H5T_NATIVE_LONG
#define H5T_READ_DOUBLE H5T_NATIVE_DOUBLE

#define H5T_STORAGE_INT H5T_STD_I64LE
#define H5T_STORAGE_DOUBLE H5T_IEEE_F64LE

/* Dimension transpose routines from Fortran to C
 * taken from MOLCAS HDF5 interface by Steven Vancoillie (copy_cast_XXX)
 */

void transpose_f2c(int rank, const hsize_t *src, hsize_t *dst)
{
  int i;
  for (i=0; i<rank; i++)
  {
    dst[rank-(i+1)] = src[i];
  }
}

#define transpose_c2f transpose_f2c

const char* check_h5type(hid_t h5type)
{
    if (h5type == H5T_IEEE_F64LE)
      return "H5T_IEEE_F64LE";
    if (h5type == H5T_NATIVE_DOUBLE)
      return "H5T_NATIVE_DOUBLE";
    if (h5type == H5T_NATIVE_LONG)
      return "H5T_NATIVE_LONG";
    if (h5type == H5T_STD_I64LE)
      return "H5T_STD_I64LE";
    return "something else";
}

hid_t hdf5_create_c(const char *myfile)
{
  return H5Fcreate(myfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t hdf5_open_c(const char *myfile)
{
  return H5Fopen(myfile, H5F_ACC_RDWR, H5P_DEFAULT);
}

herr_t hdf5_close_c(hid_t file_id)
{
  return H5Fclose(file_id);
}

hid_t h5dopen_c(hid_t dset, const char *dname)
{
  return H5Dopen(dset, dname, H5P_DEFAULT);
}

hid_t h5screate_simple_c(int rank, const hsize_t *dims)
{
  hsize_t ndims[rank];
  transpose_f2c(rank,dims,ndims);
  return H5Screate_simple(rank, ndims, NULL);
}

hid_t h5dcreate_c(hid_t myfile_id, const char *dname,hid_t h5type,hid_t space_id)
{
  hid_t native_h5type;
  switch(h5type)
  {
    case MYH5T_INT:
      native_h5type = H5T_STORAGE_INT;
      break;
    case MYH5T_DOUBLE:
      native_h5type = H5T_STORAGE_DOUBLE;
      break;
    default:
      // not supported! TODO: handle this
      native_h5type = 0;
      fprintf(stderr, "Error in h5dcreate_c wrapper: type not supported yet!\n");
      break;
  }
  return H5Dcreate(myfile_id, dname, native_h5type, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

// hid_t h5dopen_c(hid_t myfile_id, const char *dname)
// {
//   return H5Dopen(myfileid, dname);
// }

int h5sget_simple_extent_dims_c(hid_t space_id, hsize_t *dims)
{
  int rank = H5Sget_simple_extent_ndims(space_id);
  hsize_t ndims[rank]; // is this working? one should probably use malloc
  int extentdims = H5Sget_simple_extent_dims(space_id, ndims, NULL);
  transpose_c2f(rank, ndims, dims);
  return extentdims;
}

herr_t h5dread_c(hid_t dset_id, hid_t h5type, void *buf)
{
  hid_t native_h5type;
  switch(h5type)
  {
    case MYH5T_INT:
      native_h5type = H5T_READ_INT;
      break;
    case MYH5T_DOUBLE:
      native_h5type = H5T_READ_DOUBLE;
      break;
    default:
      // not supported! TODO: handle this
      native_h5type = 0;
      fprintf(stderr, "Error in h5dread_c wrapper: type not supported yet!\n");
      break;
  }
  return H5Dread(dset_id, native_h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
}


herr_t hdf5_write_cholesky_c(hid_t dset_id, hid_t space_id, const hsize_t *start, const hsize_t *count, void *buf)
{
  int rank = 2; // for Cholesky datasets the rank should be always 2.
  // Otherwise: call H5Sget_simple_extent_ndims(space_id)
  // and allocate the arrays dynamically
  hsize_t c_start[rank], c_count[rank]; // Does this work? one should probably use malloc
  hid_t memspace_id;
  herr_t error;
  transpose_f2c(rank, start, c_start);
  transpose_f2c(rank, count, c_count);
  error = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, c_start, NULL, c_count, NULL);
  if (error < 0)
  {
    fprintf(stderr, "hdf5_write_cholesky_c: Error %i in H5Sselect_hyperslab\n", error);
    // error handling
  }
  memspace_id = H5Screate_simple(rank, c_count, NULL);
  error = H5Dwrite(dset_id, H5T_IEEE_F64LE, memspace_id, space_id, H5P_DEFAULT, buf);
  if (error < 0)
  {
    fprintf(stderr, "hdf5_write_cholesky_c: Error %i in H5Dwrite\n", error);
    // error handling
  }
  error = H5Sclose(memspace_id);
  if (error < 0)
  {
    fprintf(stderr, "hdf5_write_cholesky_c: Error %i in H5Sclose\n", error);
    // error handling
  }
  return error;
}

herr_t hdf5_put_data_c(hid_t fileid, hid_t h5type, const char *dname, int rank, const hsize_t* dims, void *buf)
{
  herr_t error;
  hid_t dspace, dset;
  hsize_t ndims[rank]; // Does this work? one should probably use malloc
  transpose_f2c(rank, dims, ndims);

  hid_t native_h5type;
  switch(h5type)
  {
    case MYH5T_INT:
      native_h5type = H5T_STORAGE_INT;
      break;
    case MYH5T_DOUBLE:
      native_h5type = H5T_STORAGE_DOUBLE;
      break;
    default:
      // not supported! TODO: handle this
      native_h5type = 0;
      fprintf(stderr, "Error in hdf5_put_data_c: type not supported yet!\n");
      break;
  }

  dspace = H5Screate_simple(rank, ndims, NULL);
  dset = H5Dcreate(fileid, dname, native_h5type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  error = H5Dwrite(dset, native_h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  if (error < 0)
  {
    fprintf(stderr, "hdf5_put_data_c: Error %i in H5Dwrite\n", error);
    // error handling
  }
  error = H5Sclose(dspace);
  if (error < 0)
  {
    fprintf(stderr, "hdf5_put_data_c: Error %i in H5Sclose\n", error);
    // error handling
  }
  error = H5Dclose(dset);
  if (error < 0)
  {
    fprintf(stderr, "hdf5_put_data_c: Error %i in H5Dclose\n", error);
    // error handling
  }
  return error;
}

herr_t hdf5_get_data_c(hid_t fileid, hid_t h5type, const char *dname, int rank, const hsize_t* dims, void *buf)
{
  herr_t error;
  hid_t dset;

  hsize_t ndims[rank]; // Does this work? one should probably use malloc
  transpose_f2c(rank, dims, ndims);

  hid_t native_h5type, read_h5type;
  switch(h5type)
  {
    case MYH5T_INT:
      native_h5type = H5T_READ_INT;//warning -- maybe a mismatch since we used H5T_STD_I64LE to write
      break;
    case MYH5T_DOUBLE:
      native_h5type = H5T_READ_DOUBLE;// warning -- maybe a mismatch since we used H5T_IEEE_F64LE to write
      break;
    default:
      // not supported! TODO: handle this
      native_h5type = 0;
      fprintf(stderr, "Error in hdf5_get_data_c: type not supported yet!\n");
      break;
  }

  dset = H5Dopen(fileid, dname, H5P_DEFAULT);
  read_h5type = H5Dget_type(dset);
  if (!H5Tequal(read_h5type, native_h5type))
    fprintf(stderr, "Warning in hdf5_get_data_c: type mismatch when reading dataset %s: %s instead of %s\n", dname, check_h5type(read_h5type), check_h5type(native_h5type));
  error = H5Dread(dset, read_h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  if (error < 0)
    fprintf(stderr, "hdf5_put_data_c: Error %i in H5Dread when reading dataset %s\n", error, dname);
    // error handling
  error = H5Tclose(read_h5type);
  error = H5Dclose(dset);

  return error;
}
#endif
