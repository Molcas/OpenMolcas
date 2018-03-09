/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
/* The Molcas HDF5 C interface */

#include <sys/stat.h>
#include <string.h>

#include "molcasversion.h"
#include "molcastype.h"
#include "hdf5.h"

/* molcas buffer data types */
#ifdef _I8_
#  define H5T_MOLCAS_INT H5T_NATIVE_LONG
#else
#  define H5T_MOLCAS_INT H5T_NATIVE_INT
#endif
#define H5T_MOLCAS_REAL H5T_NATIVE_DOUBLE

/* file storage data types */
#ifdef _I8_
#  define H5T_STORAGE_INT H5T_STD_I64LE
#else
#  define H5T_STORAGE_INT H5T_STD_I32LE
#endif
#define H5T_STORAGE_REAL H5T_IEEE_F64LE

/* 12/12/2016 Leon: Compression can be enabled only if underlying HDF5 library
 * supports compression! This will lead to strange effects otherwise!
 * TODO: Check this in CMAKE and set the flag accordingly! */

/* #define _HDF5_COMPRESSION_ */

/* limits */
#define MAX_RANK 7
#define CHUNK_SIZE 125000 /* 1MB chunk size as number of doubles */

/* check a value againts the maximum rank */
#define return_on_oob_rank(rank) do {if (rank > MAX_RANK) return -1;} while(0)

/* For HDF5 object identifiers, I'm explicitly casting from
 * hid_t to INT and back, assuming the hid_t will fit in INT.
 * Looking at the H5Ipublic.h include file, hid_t is a typedef
 * of int type, so the INT from Molcas should be big enough to
 * always hold that type. */

/* external interface */

/* files */
INT mh5c_create_file (char* filename);
INT mh5c_open_file_rw (char* filename);
INT mh5c_open_file_r (char* filename);
INT mh5c_close_file (INT file_id);

/* attributes */
INT mh5c_create_attr_scalar_int(INT obj_id, char* name);
INT mh5c_create_attr_scalar_real(INT obj_id, char* name);
INT mh5c_create_attr_scalar_str(INT obj_id, char* name, INT size);
INT mh5c_create_attr_array_int(INT file_id, char* name, INT rank, INT* dims);
INT mh5c_create_attr_array_real(INT file_id, char* name, INT rank, INT* dims);
INT mh5c_create_attr_array_str(INT file_id, char* name, INT rank, INT* dims, INT size);
INT mh5c_open_attr(INT obj_id, char* name);
INT mh5c_close_attr(INT attr_id);

INT mh5c_put_attr_scalar_int(INT attr_id, void* value);
INT mh5c_put_attr_scalar_real(INT attr_id, void* value);
INT mh5c_put_attr_scalar_str(INT attr_id, void* value);
INT mh5c_put_attr_array_int(INT attr_id, void* buffer);
INT mh5c_put_attr_array_real(INT attr_id, void* buffer);
INT mh5c_put_attr_array_str(INT attr_id, void* buffer);

INT mh5c_get_attr_scalar_int(INT attr_id, void* value);
INT mh5c_get_attr_scalar_real(INT attr_id, void* value);
INT mh5c_get_attr_scalar_str(INT attr_id, void* value);
INT mh5c_get_attr_array_int(INT attr_id, void* buffer);
INT mh5c_get_attr_array_real(INT attr_id, void* buffer);
INT mh5c_get_attr_array_str(INT attr_id, void* buffer);

/* datasets */
INT mh5c_create_dset_scalar_int(INT file_id, char* name);
INT mh5c_create_dset_scalar_real(INT file_id, char* name);
INT mh5c_create_dset_scalar_str(INT file_id, char* name, INT size);
INT mh5c_create_dset_array_int(INT file_id, char* name, INT rank, INT* dims);
INT mh5c_create_dset_array_real(INT file_id, char* name, INT rank, INT* dims);
INT mh5c_create_dset_array_str(INT file_id, char* name, INT rank, INT* dims, INT size);
INT mh5c_create_dset_array_dyn_int(INT file_id, char* name, INT rank, INT* dims);
INT mh5c_create_dset_array_dyn_real(INT file_id, char* name, INT rank, INT* dims);
INT mh5c_create_dset_array_dyn_str(INT file_id, char* name, INT rank, INT* dims, INT size);
INT mh5c_open_dset(INT file_id, char* name);
INT mh5c_close_dset(INT dset_id);

INT mh5c_put_dset_scalar_int(INT dest_id, void* value);
INT mh5c_put_dset_scalar_real(INT dest_id, void* value);
INT mh5c_put_dset_scalar_str(INT dest_id, void* value);
INT mh5c_put_dset_array_int(INT dset_id, INT* extents, INT* offsets, void* buffer);
INT mh5c_put_dset_array_real(INT dset_id, INT* extents, INT* offsets, void* buffer);
INT mh5c_put_dset_array_str(INT dset_id, INT* extents, INT* offsets, void* buffer);
INT mh5c_put_dset_array_int_full(INT dset_id, void* buffer);
INT mh5c_put_dset_array_real_full(INT dset_id, void* buffer);

INT mh5c_get_dset_scalar_int(INT dest_id, void* value);
INT mh5c_get_dset_scalar_real(INT dest_id, void* value);
INT mh5c_get_dset_scalar_str(INT dest_id, void* value);
INT mh5c_get_dset_array_int(INT dset_id, INT* extents, INT* offsets, void* buffer);
INT mh5c_get_dset_array_real(INT dset_id, INT* extents, INT* offsets, void* buffer);
INT mh5c_get_dset_array_str(INT dset_id, INT* extents, INT* offsets, void* buffer);
INT mh5c_get_dset_array_int_full(INT dset_id, void* buffer);
INT mh5c_get_dset_array_real_full(INT dset_id, void* buffer);

/* array rank/dimensions */
INT mh5c_get_dset_array_rank(INT dset_id);
INT mh5c_get_dset_array_dims(INT dset_id, INT* dims);
INT mh5c_extend_dset_array(INT dset_id, INT* dims);

/* internal interface */

void copy_cast_f2c(int rank, const INT *src, hsize_t *dst);
void copy_cast_c2f(int rank, const hsize_t *src, INT *dst);

/* attributes */
hid_t mh5c_create_attr_scalar(hid_t file_id, char* name, hid_t hdf5_type);
hid_t mh5c_create_attr_array(hid_t file_id, char* name, int rank, const INT* dims, hid_t hdf5_type);
herr_t mh5c_put_attr(hid_t dset_id, void* value, hid_t value_type);
herr_t mh5c_get_attr(hid_t dset_id, void* value, hid_t value_type);

/* datasets */
hid_t mh5c_create_dset_scalar(hid_t file_id, char* name, hid_t hdf5_type);
hid_t mh5c_create_dset_array(hid_t file_id, char* name, int rank, const INT* dims, const INT mdim, hid_t hdf5_type);
herr_t mh5c_put_dset_scalar(hid_t dset_id, void* value, hid_t value_type);
herr_t mh5c_put_dset_array(hid_t dset_id, const INT* extents, const INT* offsets, void* buffer, hid_t buffer_type);
herr_t mh5c_get_dset_scalar(hid_t dset_id, void* value, hid_t value_type);
herr_t mh5c_get_dset_array(hid_t dset_id, const INT* extents, const INT* offsets, void* buffer, hid_t buffer_type);

/* The layout of data on disk in the HDF5 file uses C-style array indexing,
 * i.e. the leading dimension as slowest index. This routine convert the
 * dimensions by transposing those that come from Fortran, and convert also
 * from an INT to hsize_t */
void copy_cast_f2c(int rank, const INT *src, hsize_t *dst) {
        for (int i=0; i<rank; i++) {
                dst[rank-(i+1)] = src[i];
        }
}

void copy_cast_c2f(int rank, const hsize_t *src, INT *dst) {
        for (int i=0; i<rank; i++) {
                dst[rank-(i+1)] = src[i];
        }
}

void chunk_dimensions(int rank, const hsize_t *src, hsize_t *dst) {
        /* build up chunk dimensions as long as the chunk size is larger
        * than the total cumulative block size. At the end, cut off the
        * last dimension if needed, then any remaining chunk dimension
        * is 1. CAUTION: changing this will affect performance! */
        hsize_t block_size = 1;
        for (int i=0; i<rank; i++) block_size *= src[i];
        if (block_size == 0) {
                for (int i=0; i<rank; i++) dst[i] = 1;
        } else {
                int idim = 0;
                hsize_t block_size = src[idim];
                while (block_size < CHUNK_SIZE && idim < rank-1) {
                        dst[idim] = src[idim];
                        block_size *= src[++idim];
                }
                hsize_t tail = CHUNK_SIZE / (block_size / src[idim]);
                dst[idim] = tail > src[idim] ? src[idim] : tail;
                for (idim += 1; idim<rank; idim++) {
                        dst[idim] = 1;
                }
        }
}

/*
 * files
 */

INT mh5c_create_file (char* filename) {
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);
        hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        hid_t attr_id = mh5c_create_attr_scalar_str(file_id, "MOLCAS_VERSION", sizeof(_MOLCAS_VERSION_));
        mh5c_put_attr_scalar_str(attr_id, _MOLCAS_VERSION_);
        mh5c_close_attr(attr_id);
        return file_id;
}

INT mh5c_open_file_rw (char* filename) {
        struct stat fileinfo;
        if (stat(filename, &fileinfo)) {
                return mh5c_create_file(filename);
        } else {
                return H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
        }
}

INT mh5c_open_file_r (char* filename) {
        return H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
}

INT mh5c_close_file (INT file_id) {
        return H5Fclose((hid_t)file_id);
}

INT mh5c_is_hdf5(char* filename) {
        return H5Fis_hdf5(filename);
}

INT mh5c_exists_dset(hid_t mh5c_id, char* name) {
        return H5Lexists(mh5c_id, name, H5P_DEFAULT);
}

INT mh5c_exists_attr(hid_t obj_id, char* name) {
        return H5Aexists(obj_id, name);
}

/*
 * attributes
 */

INT mh5c_create_attr_scalar_int(INT dset_id, char* name) {
        return mh5c_create_attr_scalar((hid_t)dset_id, name, H5T_STORAGE_INT);
}

INT mh5c_create_attr_scalar_real(INT dset_id, char* name) {
        return mh5c_create_attr_scalar((hid_t)dset_id, name, H5T_STORAGE_REAL);
}

INT mh5c_create_attr_scalar_str(INT dset_id, char* name, INT size) {
        hid_t h5t_string;
        INT attr_id;
        h5t_string = H5Tcopy(H5T_C_S1);
        H5Tset_size (h5t_string, size);
        H5Tset_strpad(h5t_string, H5T_STR_NULLPAD);
        attr_id = mh5c_create_attr_scalar(dset_id, name, h5t_string);
        H5Tclose(h5t_string);
        return attr_id;
}

INT mh5c_create_attr_array_int(INT file_id, char* name, INT rank, INT* dims) {
        return mh5c_create_attr_array(file_id, name, rank, dims, H5T_STORAGE_INT);
}
INT mh5c_create_attr_array_real(INT file_id, char* name, INT rank, INT* dims) {
        return mh5c_create_attr_array(file_id, name, rank, dims, H5T_STORAGE_REAL);
}
INT mh5c_create_attr_array_str(INT file_id, char* name, INT rank, INT* dims, INT size) {
        hid_t h5t_string;
        INT attr_id;
        h5t_string = H5Tcopy(H5T_C_S1);
        H5Tset_size (h5t_string, size);
        H5Tset_strpad(h5t_string, H5T_STR_NULLPAD);
        attr_id = mh5c_create_attr_array(file_id, name, rank, dims, h5t_string);
        H5Tclose(h5t_string);
        return attr_id;
}

INT mh5c_open_attr(INT obj_id, char* name) {
        return H5Aopen(obj_id, name, H5P_DEFAULT);
}

INT mh5c_close_attr(INT attr_id) {
        return H5Aclose(attr_id);
}


INT mh5c_put_attr_scalar_int(INT attr_id, void* value) {
        return mh5c_put_attr(attr_id, value, H5T_MOLCAS_INT);
}

INT mh5c_put_attr_scalar_real(INT attr_id, void* value) {
        return mh5c_put_attr(attr_id, value, H5T_MOLCAS_REAL);
}

INT mh5c_put_attr_scalar_str(INT attr_id, void* value) {
        hid_t h5t_string;
        herr_t rc;
        h5t_string = H5Aget_type(attr_id);
        rc = mh5c_put_attr(attr_id, value, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}

INT mh5c_put_attr_array_int(INT attr_id, void* buffer) {
        return mh5c_put_attr(attr_id, buffer, H5T_MOLCAS_INT);
}
INT mh5c_put_attr_array_real(INT attr_id, void* buffer) {
        return mh5c_put_attr(attr_id, buffer, H5T_MOLCAS_REAL);
}
INT mh5c_put_attr_array_str(INT attr_id, void* buffer) {
        hid_t h5t_string;
        h5t_string = H5Aget_type(attr_id);
        return mh5c_put_attr(attr_id, buffer, h5t_string);
        H5Tclose(h5t_string);
}

INT mh5c_get_attr_scalar_int(INT attr_id, void* value) {
        return mh5c_get_attr(attr_id, value, H5T_MOLCAS_INT);
}

INT mh5c_get_attr_scalar_real(INT attr_id, void* value) {
        return mh5c_get_attr(attr_id, value, H5T_MOLCAS_REAL);
}

INT mh5c_get_attr_scalar_str(INT attr_id, void* value) {
        herr_t rc;
        hid_t h5t_string;
        h5t_string = H5Aget_type(attr_id);
        rc = mh5c_get_attr(attr_id, value, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}

INT mh5c_get_attr_array_int(INT attr_id, void* buffer) {
        return mh5c_get_attr(attr_id, buffer, H5T_MOLCAS_INT);
}
INT mh5c_get_attr_array_real(INT attr_id, void* buffer) {
        return mh5c_get_attr(attr_id, buffer, H5T_MOLCAS_REAL);
}
INT mh5c_get_attr_array_str(INT attr_id, void* buffer) {
        hid_t h5t_string;
        h5t_string = H5Aget_type(attr_id);
        return mh5c_get_attr(attr_id, buffer, h5t_string);
        H5Tclose(h5t_string);
}

/*
 * datasets
 */

INT mh5c_create_dset_scalar_int(INT file_id, char* name) {
        return mh5c_create_dset_scalar((hid_t)file_id, name, H5T_STORAGE_INT);
}

INT mh5c_create_dset_scalar_real(INT file_id, char* name) {
        return mh5c_create_dset_scalar((hid_t)file_id, name, H5T_STORAGE_REAL);
}

INT mh5c_create_dset_scalar_str(INT file_id, char* name, INT size) {
        hid_t dset_id;
        hid_t h5t_string;
        h5t_string = H5Tcopy(H5T_C_S1);
        H5Tset_size (h5t_string, size);
        H5Tset_strpad(h5t_string, H5T_STR_NULLPAD);
        dset_id = mh5c_create_dset_scalar((hid_t)file_id, name, h5t_string);
        H5Tclose(h5t_string);
        return dset_id;
}

INT mh5c_create_dset_array_int(INT file_id, char* name, INT rank, INT* dims) {
        return mh5c_create_dset_array(file_id, name, rank, dims, 0, H5T_STORAGE_INT);
}
INT mh5c_create_dset_array_real(INT file_id, char* name, INT rank, INT* dims) {
        return mh5c_create_dset_array(file_id, name, rank, dims, 0, H5T_STORAGE_REAL);
}
INT mh5c_create_dset_array_str(INT file_id, char* name, INT rank, INT* dims, INT size) {
        hid_t dset_id;
        hid_t h5t_string;
        h5t_string = H5Tcopy(H5T_C_S1);
        H5Tset_size (h5t_string, size);
        H5Tset_strpad(h5t_string, H5T_STR_NULLPAD);
        dset_id = mh5c_create_dset_array(file_id, name, rank, dims, 0, h5t_string);
        H5Tclose(h5t_string);
        return dset_id;
}
INT mh5c_create_dset_array_dyn_int(INT file_id, char* name, INT rank, INT* dims) {
        return mh5c_create_dset_array(file_id, name, rank, dims, H5S_UNLIMITED, H5T_STORAGE_INT);
}
INT mh5c_create_dset_array_dyn_real(INT file_id, char* name, INT rank, INT* dims) {
        return mh5c_create_dset_array(file_id, name, rank, dims, H5S_UNLIMITED, H5T_STORAGE_REAL);
}
INT mh5c_create_dset_array_dyn_str(INT file_id, char* name, INT rank, INT* dims, INT size) {
        hid_t dset_id;
        hid_t h5t_string;
        h5t_string = H5Tcopy(H5T_C_S1);
        H5Tset_size (h5t_string, size);
        H5Tset_strpad(h5t_string, H5T_STR_NULLPAD);
        dset_id = mh5c_create_dset_array(file_id, name, rank, dims, H5S_UNLIMITED, h5t_string);
        H5Tclose(h5t_string);
        return dset_id;
}

INT mh5c_open_dset(INT file_id, char* name) {
        return H5Dopen(file_id, name, H5P_DEFAULT);
}

INT mh5c_close_dset(INT dset_id) {
        return H5Dclose(dset_id);
}

INT mh5c_put_dset_scalar_int(INT dset_id, void* value) {
        return mh5c_put_dset_scalar((hid_t)dset_id, value, H5T_MOLCAS_INT);
}

INT mh5c_put_dset_scalar_real(INT dset_id, void* value) {
        return mh5c_put_dset_scalar((hid_t)dset_id, value, H5T_MOLCAS_REAL);
}

INT mh5c_put_dset_scalar_str(INT dset_id, void* value) {
        hid_t h5t_string;
        INT rc;
        h5t_string = H5Dget_type (dset_id);
        rc = mh5c_put_dset_scalar((hid_t)dset_id, value, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}

INT mh5c_put_dset_array_int(INT dset_id, INT* extents, INT* offsets, void* buffer) {
        return mh5c_put_dset_array(dset_id, extents, offsets, buffer, H5T_MOLCAS_INT);
}

INT mh5c_put_dset_array_real(INT dset_id, INT* extents, INT* offsets, void* buffer) {
        return mh5c_put_dset_array(dset_id, extents, offsets, buffer, H5T_MOLCAS_REAL);
}

INT mh5c_put_dset_array_str(INT dset_id, INT* extents, INT* offsets, void* buffer) {
        hid_t h5t_string;
        INT rc;
        h5t_string = H5Dget_type (dset_id);
        rc = mh5c_put_dset_array(dset_id, extents, offsets, buffer, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}

INT mh5c_put_dset_array_int_full(INT dset_id, void* buffer) {
        return mh5c_put_dset_array(dset_id, NULL, NULL, buffer, H5T_MOLCAS_INT);
}

INT mh5c_put_dset_array_real_full(INT dset_id, void* buffer) {
        return mh5c_put_dset_array(dset_id, NULL, NULL, buffer, H5T_MOLCAS_REAL);
}

INT mh5c_put_dset_array_str_full(INT dset_id, void* buffer) {
        hid_t h5t_string;
        INT rc;
        h5t_string = H5Dget_type (dset_id);
        rc = mh5c_put_dset_array(dset_id, NULL, NULL, buffer, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}


INT mh5c_get_dset_scalar_int(INT dset_id, void* value) {
        return mh5c_get_dset_scalar((hid_t)dset_id, value, H5T_MOLCAS_INT);
}

INT mh5c_get_dset_scalar_real(INT dset_id, void* value) {
        return mh5c_get_dset_scalar((hid_t)dset_id, value, H5T_MOLCAS_REAL);
}

INT mh5c_get_dset_scalar_str(INT dset_id, void* value) {
        herr_t rc;
        hid_t h5t_string;
        h5t_string = H5Dget_type (dset_id);
        rc = mh5c_get_dset_scalar((hid_t)dset_id, value, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}

INT mh5c_get_dset_array_int(INT dset_id, INT* extents, INT* offsets, void* buffer) {
        return mh5c_get_dset_array(dset_id, extents, offsets, buffer, H5T_MOLCAS_INT);
}
INT mh5c_get_dset_array_real(INT dset_id, INT* extents, INT* offsets, void* buffer) {
        return mh5c_get_dset_array(dset_id, extents, offsets, buffer, H5T_MOLCAS_REAL);
}
INT mh5c_get_dset_array_str(INT dset_id, INT* extents, INT* offsets, void* buffer) {
        hid_t h5t_string;
        INT rc;
        h5t_string = H5Dget_type (dset_id);
        rc = mh5c_get_dset_array(dset_id, extents, offsets, buffer, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}

INT mh5c_get_dset_array_int_full(INT dset_id, void* buffer) {
        return mh5c_get_dset_array(dset_id, NULL, NULL, buffer, H5T_MOLCAS_INT);
}
INT mh5c_get_dset_array_real_full(INT dset_id, void* buffer) {
        return mh5c_get_dset_array(dset_id, NULL, NULL, buffer, H5T_MOLCAS_REAL);
}
INT mh5c_get_dset_array_str_full(INT dset_id, void* buffer) {
        hid_t h5t_string;
        INT rc;
        h5t_string = H5Dget_type (dset_id);
        rc = mh5c_get_dset_array(dset_id, NULL, NULL, buffer, h5t_string);
        H5Tclose(h5t_string);
        return rc;
}


INT mh5c_get_dset_array_rank(INT dset_id) {
        hid_t space_id = H5Dget_space(dset_id);
        int rank = H5Sget_simple_extent_ndims(space_id);
        H5Sclose(space_id);
        return rank;
}

INT mh5c_get_dset_array_dims(INT dset_id, INT* dims) {
        hsize_t hdims[MAX_RANK];
        hid_t space_id = H5Dget_space(dset_id);
        int rank = H5Sget_simple_extent_ndims(space_id);
        return_on_oob_rank(rank);
        copy_cast_f2c(rank,dims,hdims);
        rank = H5Sget_simple_extent_dims(space_id, hdims, NULL);
        copy_cast_c2f(rank,hdims,dims);
        H5Sclose(space_id);
        return rank;
}

INT mh5c_extend_dset_array(INT dset_id, INT* dims) {
        int rank = mh5c_get_dset_array_rank(dset_id);
        hsize_t hdims[MAX_RANK];
        copy_cast_f2c(rank,dims,hdims);
        return H5Dset_extent(dset_id,hdims);
}

/* attribute array info */

INT mh5c_get_attr_scalar_rank(INT dset_id) {
        hid_t space_id = H5Aget_space(dset_id);
        int rank = H5Sget_simple_extent_ndims(space_id);
        H5Sclose(space_id);
        return rank;
}

INT mh5c_get_attr_scalar_dims(INT dset_id, INT* dims) {
        hsize_t hdims[MAX_RANK];
        hid_t space_id = H5Aget_space(dset_id);
        int rank = H5Sget_simple_extent_ndims(space_id);
        return_on_oob_rank(rank);
        copy_cast_f2c(rank,dims,hdims);
        rank = H5Sget_simple_extent_dims(space_id, hdims, NULL);
        copy_cast_c2f(rank,hdims,dims);
        H5Sclose(space_id);
        return rank;
}


/*************************
 * the generic functions *
 *************************/

/* attributes */

hid_t mh5c_create_attr_scalar(hid_t obj_id, char* name, hid_t hdf5_type) {
        hid_t space_id, dset_id;
        space_id = H5Screate(H5S_SCALAR);
        dset_id = H5Acreate(obj_id, name, hdf5_type, space_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space_id);
        return dset_id;
}

hid_t mh5c_create_attr_array(hid_t obj_id, char* name, int rank, const INT* dims, hid_t hdf5_type) {
        herr_t status;
        hid_t space_id, dset_id;
        hsize_t hdims[MAX_RANK];
        return_on_oob_rank(rank);
        copy_cast_f2c(rank,dims,hdims);
        space_id = H5Screate_simple(rank, hdims, NULL);
        dset_id = H5Acreate(obj_id, name, hdf5_type, space_id, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Sclose(space_id);
        /* "use" the variable to avoid warning */
        (void)status;
        return dset_id;
}

herr_t mh5c_put_attr(hid_t attr_id, void* buffer, hid_t buffer_type) {
        herr_t rc;
        rc = H5Awrite(attr_id, buffer_type, buffer);
        H5Fflush(attr_id, H5F_SCOPE_LOCAL);
        return rc;
}

herr_t mh5c_get_attr(hid_t attr_id, void* buffer, hid_t buffer_type) {
        return H5Aread(attr_id, buffer_type, buffer);
}

/* datasets */

hid_t mh5c_create_dset_scalar(hid_t file_id, char* name, hid_t hdf5_type) {
        hid_t space_id, dset_id;
        space_id = H5Screate(H5S_SCALAR);
        dset_id = H5Dcreate(file_id, name, hdf5_type, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space_id);
        return dset_id;
}

hid_t mh5c_create_dset_array(hid_t file_id, char* name, int rank, const INT* dims, const INT mdim, hid_t hdf5_type) {
        herr_t status;
        hid_t space_id, plist_id, dset_id;
        hsize_t hdims[MAX_RANK];
        hsize_t kdims[MAX_RANK];
        hsize_t cdims[MAX_RANK];
        return_on_oob_rank(rank);
        copy_cast_f2c(rank,dims,hdims);
        if (mdim == 0) {
                space_id = H5Screate_simple(rank, hdims, NULL);
        } else {
                for (int i=0; i<rank; i++) {
                        kdims[i] = mdim;
                }
                space_id = H5Screate_simple(rank, hdims, kdims);
        }
        plist_id = H5Pcreate (H5P_DATASET_CREATE);
#ifdef _HDF5_COMPRESSION_
        chunk_dimensions(rank, hdims, cdims);
        status = H5Pset_chunk (plist_id, rank, cdims);
        status = H5Pset_deflate (plist_id, 6);
#else
        if (mdim < 0) {
                chunk_dimensions(rank, hdims, cdims);
                status = H5Pset_chunk (plist_id, rank, cdims);
        }
#endif
        dset_id = H5Dcreate(file_id, name, hdf5_type, space_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        status = H5Sclose(space_id);
        /* "use" the variable to avoid warning */
        (void)status;
        return dset_id;
}

herr_t mh5c_put_dset_scalar(hid_t dset_id, void* value, hid_t value_type) {
        herr_t status;
        status = H5Dwrite(dset_id, value_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
        H5Fflush(dset_id, H5F_SCOPE_LOCAL);
        return status;
}

herr_t mh5c_put_dset_array(hid_t dset_id, const INT* extents, const INT* offsets, void* buffer, hid_t buffer_type) {
        herr_t status;
        hid_t mem_space_id, dset_space_id;
        hsize_t hextents[MAX_RANK], hoffsets[MAX_RANK];
        int rank;
        if (extents == NULL) {
                status = H5Dwrite(dset_id, buffer_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        } else {
                /* fetch dataspace and its rank */
                dset_space_id = H5Dget_space(dset_id);
                rank = H5Sget_simple_extent_ndims(dset_space_id);
                return_on_oob_rank(rank);

                /* convert extents/offsets to proper type */
                copy_cast_f2c(rank,extents,hextents);
                copy_cast_f2c(rank,offsets,hoffsets);

                /* create a simple space for the data in memory */
                mem_space_id = H5Screate_simple(rank, hextents, NULL);

                /* create a space which is a subset of the dataset */
                status = H5Sselect_hyperslab(dset_space_id, H5S_SELECT_SET, hoffsets, NULL, hextents, NULL);

                /* write the hyperslab */
                status = H5Dwrite(dset_id, buffer_type, mem_space_id, dset_space_id, H5P_DEFAULT, buffer);

                /* clean up dataspaces */
                status = H5Sclose(dset_space_id);
                status = H5Sclose(mem_space_id);
        }
        status = H5Fflush(dset_id, H5F_SCOPE_LOCAL);
        return status;
}

herr_t mh5c_get_dset_scalar(hid_t dset_id, void* value, hid_t value_type) {
        return H5Dread(dset_id, value_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
}

herr_t mh5c_get_dset_array(hid_t dset_id, const INT* extents, const INT* offsets, void* buffer, hid_t buffer_type) {
        herr_t status;
        hid_t mem_space_id, dset_space_id;
        hsize_t hextents[MAX_RANK], hoffsets[MAX_RANK];
        int rank;
        if (extents == NULL) {
                status = H5Dread(dset_id, buffer_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        } else {
                /* fetch dataspace and its rank */
                dset_space_id = H5Dget_space(dset_id);
                rank = H5Sget_simple_extent_ndims(dset_space_id);
                return_on_oob_rank(rank);

                /* convert extents/offsets to proper type */
                copy_cast_f2c(rank,extents,hextents);
                copy_cast_f2c(rank,offsets,hoffsets);

                /* create a simple space for the data in memory */
                mem_space_id = H5Screate_simple(rank, hextents, NULL);

                /* create a space which is a subset of the dataset */
                status = H5Sselect_hyperslab(dset_space_id, H5S_SELECT_SET, hoffsets, NULL, hextents, NULL);

                /* read the hyperslab */
                status = H5Dread(dset_id, buffer_type, mem_space_id, dset_space_id, H5P_DEFAULT, buffer);

                /* clean up dataspaces */
                status = H5Sclose(dset_space_id);
                status = H5Sclose(mem_space_id);
        }
        status = H5Fflush(dset_id, H5F_SCOPE_LOCAL);
        return status;
}
