!! hdf5_qcm: HDF5 Fortran/C interface for QCMaquis driver and CD-NEVPT2
!! This interface replaces the previous Fortran module-based HDF5 F2003 interface and does not require the
!! compilation of the HDF5 library with the same compiler version as the QCMaquis driver/NEVPT2.
!!
!! Fortran-C interoperability inspired by http://fortranwiki.org/fortran/show/Generating+C+Interfaces
!! and Steven Vancoillie's hdf5_util in OpenMOLCAS
!! (C) 2017 Leon Freitag and Stefan Knecht
!!          Laboratory for Physical Chemistry, ETH Zurich
!!
!!  hdf5_qcm is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  hdf5_qcm is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with hdf5_qcm. If not, see <http://www.gnu.org/licenses/>.
!!
module hdf5_utils


#ifdef _OPENMP_
  use omp_lib
#endif
  use iso_c_binding
  use iso_fortran_env

  implicit none

  public hdf5_init
  public hdf5_exit
  public hdf5_open
  public hdf5_close
  public hdf5_create
  public hdf5_put_data_dp
  public hdf5_put_data_ix
  public hdf5_put_data
  public hdf5_get_data
  public hdf5_get_data_dp
  public hdf5_get_data_ix
  public hdf5_init_wr_cholesky
  public hdf5_init_rd_cholesky
  public hdf5_close_cholesky
  public hdf5_read_all_cholesky
  public hdf5_write_cholesky

! export the data types used by hdf5 since we're no longer including the HDF5 module
  public HID_T
  public HERR_T
  public HSIZE_T

  interface hdf5_put_data
   module procedure  hdf5_put_data_dp,  hdf5_put_data_ix, hdf5_put_data_dpscalar, hdf5_put_data_iscalar
  end interface  hdf5_put_data

  interface hdf5_get_data
   module procedure  hdf5_get_data_dp,  hdf5_get_data_ix, hdf5_get_data_dpscalar, hdf5_get_data_iscalar
  end interface  hdf5_get_data

! The definitions of HID_T, HERR_T, HSIZE_T should match those of their c counterparts:
! typedef int64_t hid_t; typedef int herr_t; typedef unsigned long long hsize_t;
! The code below should map the correct Fortran types to C types

  INTEGER, PARAMETER :: HID_T = c_int64_t
  INTEGER, PARAMETER :: HERR_T = c_int
  INTEGER, PARAMETER :: HSIZE_T = c_long_long

  INTEGER, PARAMETER :: Fortran_INTEGER_8 = 8
  INTEGER, PARAMETER :: Fortran_REAL_8 = 8

  INTEGER*4, PARAMETER :: int_kind    = SELECTED_INT_KIND(Fortran_INTEGER_8) !should map to INTEGER*8
  INTEGER*4, PARAMETER :: real_kind   = SELECTED_REAL_KIND(Fortran_REAL_8) !should map to REAL*8

  ! Cheap and dirty h5kind_to_type() replacement
  ! type definitions that get passed to the C routines and checked there. Only these two are supported for now
  INTEGER(HID_T), PARAMETER :: MYH5T_INT = 1;
  INTEGER(HID_T), PARAMETER :: MYH5T_DOUBLE = 2;

  character(len= 8), parameter, public :: filename          = "nevpt.h5"
  character(len= 7), parameter, public :: ijklname          = "ijkl.h5"
  integer,           parameter, public :: taglen            = 24
  integer,           parameter, public :: taglenx           = 25
  integer,                      public :: tagx              = -1
  character(len=taglen),        public :: datatag           = " "

  INTEGER(HSIZE_T), public             :: datadim(1:8)    = -1
  INTEGER(HID_T),   public             :: file_id(1:2)    = -1
  INTEGER(HSIZE_T), public             :: datadim_bound   = -1

  private

  INTEGER(HERR_T)                   :: error
  TYPE(C_PTR) :: f_ptr

  interface
    function hdf5_put_data_c(myfile_id, h5type, mydata, rank, dims, data_in) result(err) bind(C)
      import C_char, C_int, C_ptr, HID_T, HSIZE_T, HERR_T
      integer(HID_T), VALUE :: myfile_id
      integer(HID_T), VALUE :: h5type
      character(len=1,kind=C_char), dimension(*), intent(in) :: mydata
      integer(c_int), VALUE :: rank
      INTEGER(HSIZE_T) :: dims(*)
      type(c_ptr), VALUE :: data_in
      integer(HERR_T) :: err
    end function hdf5_put_data_c

    function hdf5_get_data_c(myfile_id, h5type, mydata, rank, dims, data_out) result(err) bind(C)
      import C_char, C_int, C_ptr, HID_T, HSIZE_T, HERR_T
      integer(HID_T), VALUE :: myfile_id
      integer(HID_T), VALUE :: h5type
      character(len=1,kind=C_char), dimension(*), intent(in) :: mydata
      integer(c_int), VALUE :: rank
      INTEGER(HSIZE_T) :: dims(*)
      type(c_ptr), VALUE :: data_out
      integer(HERR_T) :: err
    end function hdf5_get_data_c
  end interface

  contains

!------------------------------------------------------------------------------
! do nothing for hdf5_init or exit since open and close_f calls are not required in the C interface
! The routines are kept for compatibility though
      subroutine hdf5_init()

      ! Maybe put a deprecated warning here?

      end subroutine hdf5_init
!------------------------------------------------------------------------

      subroutine hdf5_exit()

      end subroutine hdf5_exit
!------------------------------------------------------------------------

      subroutine hdf5_create(myfile, myfile_id)

        interface
          function hdf5_create_c(myfile) result (fh) bind(C)
            import C_char, HID_T
            character(len=1,kind=C_char), dimension(*), intent(in) :: myfile
            integer(HID_T) :: fh
          end function
        end interface

        character(len=*), intent(in)  :: myfile
        integer(HID_T)  , intent(out) :: myfile_id

        myfile_id = hdf5_create_c(trim(myfile)//C_NULL_CHAR)
        if (myfile_id.lt.0) then
          error = -1
        else
          error = 0
        end if
      end subroutine hdf5_create
!------------------------------------------------------------------------

      subroutine hdf5_open(myfile, myfile_id)

        interface
          function hdf5_open_c(myfile) result(fh) bind(C)
            import C_char, HID_T
            character(len=1,kind=C_char), dimension(*), intent(in) :: myfile
            integer(HID_T) :: fh
          end function
        end interface

        character(len=*), intent(in)  :: myfile
        integer(HID_T)  , intent(out) :: myfile_id

        myfile_id = hdf5_open_c(trim(myfile)//C_NULL_CHAR)
        if (myfile_id.lt.0) then
          error = -1
        else
          error = 0
        end if
      end subroutine hdf5_open
!------------------------------------------------------------------------

      subroutine hdf5_close(myfile_id)

       interface
         function hdf5_close_c(myfile_id) result(err) bind(C)
         import HID_T, HERR_T
         integer(HID_T), VALUE :: myfile_id
         integer(HERR_T) :: err
         end function
       end interface

        integer(HID_T)  , intent(inout) :: myfile_id
        integer(HERR_T) :: err

        err = hdf5_close_c(myfile_id)
        ! error handling is ignored?

      end subroutine hdf5_close
!------------------------------------------------------------------------

! C interfaces to HDF5 functions used by more than one subroutine

      subroutine h5screate_simple(rank, dims, space_id, err)
        interface
          function h5screate_simple_c(rank, dims) result(space_id) bind(C)
            import HID_T, c_int, HSIZE_T
            integer(c_int), VALUE :: rank
            integer(HSIZE_T) :: dims(*)
            INTEGER(HID_T) :: space_id
          end function
        end interface

        integer(c_int) ,intent(in) :: rank
        integer(HSIZE_T) :: dims(*)
        INTEGER(HID_T),intent(inout) :: space_id
        INTEGER(HERR_T),intent(inout) :: err

        err = 0
        space_id = h5screate_simple_c(rank, dims)
        if (space_id.lt.0) err = space_id
      end subroutine

      subroutine h5dcreate(myfile_id, dname, h5type, space_id, dset_id, err)
        interface
          function h5dcreate_c(myfile_id, dname, h5type, space_id) result(dset_id) bind(C)
            import C_char, HID_T
            integer(HID_T), VALUE :: myfile_id
            character(len=1,kind=C_char), dimension(*), intent(in) :: dname
            integer(HID_T), VALUE :: h5type
            integer(HID_T), VALUE :: space_id
            INTEGER(HID_T) :: dset_id
          end function
        end interface
        integer(HID_T), intent(in) :: myfile_id
        character(len=*), intent(in)  :: dname
        integer(HID_T), intent(in) :: h5type
        integer(HID_T), intent(in) :: space_id
        INTEGER(HID_T), intent(inout) :: dset_id
        INTEGER(HERR_T),intent(inout) :: err

        err = 0
        dset_id = h5dcreate_c(myfile_id, trim(dname)//C_NULL_CHAR, h5type, space_id)
        if (dset_id.lt.0) err = dset_id
      end subroutine

      subroutine h5dopen(myfile_id, dname, dset_id, err)
        interface
          function h5dopen_c(myfile_id, dname) result(dset_id) bind(C)
          import C_char, HID_T
            integer(HID_T), VALUE :: myfile_id
            character(len=1,kind=C_char), dimension(*), intent(in) :: dname
            INTEGER(HID_T) :: dset_id
          end function
        end interface
        integer(HID_T), intent(in) :: myfile_id
        character(len=*), intent(in)  :: dname
        INTEGER(HID_T), intent(inout) :: dset_id
        INTEGER(HERR_T),intent(inout) :: err

        err = 0
        dset_id = h5dopen_c(myfile_id, trim(dname)//C_NULL_CHAR)
        if (dset_id.lt.0) err = dset_id
      end subroutine

      subroutine h5dget_space(dset_id, space_id, err)
        interface
          function h5dget_space_c(dset_id) result(space_id) bind(C, name='H5Dget_space')
            import HID_T
            INTEGER(HID_T), VALUE :: dset_id
            INTEGER(HID_T) :: space_id
          end function
        end interface
        INTEGER(HID_T), intent(in)    :: dset_id
        INTEGER(HID_T), intent(inout) :: space_id
        INTEGER(HERR_T),intent(inout) :: err

        err = 0
        space_id = h5dget_space_c(dset_id)
        if (space_id.lt.0) err = space_id
      end subroutine

      subroutine h5sget_simple_extent_dims(space_id, dims, err)
        interface
          function h5sget_simple_extent_dims_c(space_id, dims) result(err) bind(C)
            import HID_T, HSIZE_T, c_int
            INTEGER(HID_T), VALUE :: space_id
            INTEGER(HSIZE_T) :: dims(*)
            INTEGER(c_int) :: err
          end function
        end interface
        INTEGER(HID_T), intent(in) :: space_id
        INTEGER(HSIZE_T) :: dims(*)
        INTEGER(c_int), intent(inout) :: err

        err = h5sget_simple_extent_dims_c(space_id, dims)
      end subroutine

      subroutine h5dclose(dset_id, err)
        interface
          function h5dclose_c(dset_id) result (err) bind(C,name='H5Dclose')
            import HID_T,HERR_T
            integer(HID_T), VALUE :: dset_id
            INTEGER(HERR_T) :: err
          end function
        end interface

        integer(HID_T), intent(in) :: dset_id
        INTEGER(HERR_T), intent(inout) :: err

        err = h5dclose_c(dset_id)
      end subroutine

      subroutine h5sclose(space_id, err)
        interface
          function h5sclose_c(space_id) result (err) bind(C,name='H5Sclose')
            import HID_T, HERR_T
            integer(HID_T), VALUE :: space_id
            INTEGER(HERR_T) :: err
          end function
        end interface

        integer(HID_T), intent(in) :: space_id
        INTEGER(HERR_T), intent(inout) :: err

        err = h5sclose_c(space_id)
      end subroutine

      ! one must pass C_LOC of the actual argument to data_out
      subroutine h5dread(dset_id, h5type, data_out, err)
        interface
          function h5dread_c(dset_id, h5type, data_out) result(err) bind(C)
            import c_ptr, HID_T, HERR_T
            integer(HID_T), VALUE :: dset_id
            integer(HID_T), VALUE :: h5type
            type(c_ptr), value :: data_out
            INTEGER(HERR_T) :: err
          end function
        end interface
        integer(HID_T), intent(in) :: dset_id
        integer(HID_T), intent(in) :: h5type
        type(c_ptr), intent(inout) :: data_out
        INTEGER(HERR_T), intent(inout) :: err

        err = h5dread_c(dset_id, h5type, data_out)
      end subroutine

! Utilities needed for Cholesky MO vectors
!
! Initialises a record with Cholesky vectors for a given symmetry and an array containing the size of Cholesky sets

      subroutine hdf5_init_wr_cholesky(myfile_id, jsym, nlvec, nchovec, choset_id, space_id)
        integer(HID_T)  , intent(in) :: myfile_id
        integer, intent(in)          :: jsym    ! symmetry index
        integer, intent(in)          :: nlvec   ! Dimension of one Cholesky vector
        integer, intent(in)          :: nchovec ! Number of Cholesky vectors in symmetry jsym (pass NumCho(JSYM))
        integer(HID_T), intent(out)  :: choset_id ! Cholesky dataset identifier
        integer(HID_T), intent(out)  :: space_id ! dataset identifier
        character(len=8)             :: dname ! Cholesky dataset name
        integer(c_int), parameter          :: data_rank = 2

        integer(HSIZE_T), dimension(1:2) :: ldims
        ldims(1:2) = (/nchovec,nlvec/)

        ! First implementation: create a full dataset for all the Cholesky vectors
        ! If this turns out too memory hungry, then an implementation of an extensible dataset would replace it

  ! create the name for the group
        write(dname,'(A,I1)') "CHOVEC_", jsym

  ! create the dataset and the dataspace associated with the number of Cholesky sets
        call h5screate_simple(data_rank, ldims, space_id, error)
        call h5dcreate(myfile_id, trim(dname), MYH5T_DOUBLE, space_id, choset_id, error)

      end subroutine hdf5_init_wr_cholesky

      subroutine hdf5_init_rd_cholesky(myfile_id,choset_id,space_id,jsym,nlvec,nchovec)
        integer(HID_T), intent(in)  :: myfile_id
        integer(HID_T), intent(out) :: choset_id ! Cholesky dataset identifier
        integer(HID_T), intent(out) :: space_id ! dataspace identifier
        integer, intent(in)         :: jsym    ! symmetry index
        character(len=8)            :: dname ! Cholesky dataset name
        integer, intent(out) :: nlvec ! number of vectors in the Cholesky set
        integer, intent(out) :: nchovec ! size of a cholesky vector (should be norb(norb+1)/2 or something
        integer(HSIZE_T), dimension(1:2) :: ds_size ! size of the cholesky dataspace
  !        integer(HSIZE_T), dimension(1:2) :: maxcount ! necessary for the hdf call below, can be discarded

        write(dname,'(A,I1)') "CHOVEC_", jsym

        call h5dopen(myfile_id, dname, choset_id, error)
        call h5dget_space(choset_id, space_id, error)

        ! obtain the dimensions of the Cholesky set and return them
        call h5sget_simple_extent_dims(space_id, ds_size, error)

        nchovec = ds_size(1) ! right now, nchovec should be equal to the value stored in "norbtt"
        nlvec = ds_size(2)
      end subroutine hdf5_init_rd_cholesky

      subroutine hdf5_close_cholesky(choset_id, space_id)
        integer(HID_T), intent(in)  :: choset_id ! Cholesky dataset identifier
        integer(HID_T), intent(in)  :: space_id ! dataspace identifier

        call h5dclose(choset_id, error)
        call h5sclose(space_id, error)

      end subroutine hdf5_close_cholesky


     ! reads the whole Cholesky set into memory
     subroutine hdf5_read_all_cholesky(choset_id,space_id,data_out)
       integer(HID_T), intent(in)  :: choset_id ! Cholesky dataset identifier
       integer(HID_T), intent(in) :: space_id ! dataspace identifier
       real(kind=real_kind),dimension(*),intent(out),target :: data_out ! output: Cholesky vectors

       f_ptr = C_LOC(data_out)
       call h5dread(choset_id, MYH5T_DOUBLE, f_ptr, error)

     end subroutine hdf5_read_all_cholesky


     subroutine hdf5_write_cholesky(choset_id,space_id,idx,startl,nl,data_in)
       interface
         function hdf5_write_cholesky_c(dset_id, space_id, start, num, data_in) result(err) bind(C)
          import c_ptr, HID_T, HSIZE_T, herr_t
          integer(HID_T), VALUE :: dset_id
          integer(HID_T), VALUE :: space_id
          INTEGER(HSIZE_T) :: start(*)
          INTEGER(HSIZE_T) :: num(*)
          type(c_ptr), VALUE  :: data_in
          INTEGER(HERR_T) :: err
         end function
       end interface
       integer(HID_T), intent(in)  :: choset_id ! Cholesky dataset identifier
       integer(HID_T), intent(in)  :: space_id ! dataspace identifier

       integer, intent(in)         :: idx    ! index of the current element to be written (WARNING: this index starts at 0)
       integer, intent(in)         :: startl ! index of the first cholesky vector to be written
       integer, intent(in)         :: nl     ! number of Cholesky vectors in this batch
       real(kind=real_kind),dimension(*),intent(in),target :: data_in ! a column of the Cholesky vector elements with length of nl

       integer(HSIZE_T), dimension(1:2) :: start ! offset for the hyperslab, always 0 for rows since
                                                                  ! cholesky vectors are never broken between batches
       integer(HSIZE_T), dimension(1:2) :: count ! size of the hyperslab, full size for rows and nvec for columns

       ! select hyperslab corresponding to nl x 1 elements starting at (startl, idx) where the new Cholesky elements should be written
       start(1:2) = (/startl,idx/)
       count(1:2) = (/nl,1/)

       f_ptr = C_LOC(data_in)
       error = hdf5_write_cholesky_c(choset_id, space_id, start, count, f_ptr)

     end subroutine hdf5_write_cholesky

     subroutine hdf5_put_data_dpscalar(myfile_id, mydata, dims, data_in)
       integer(HID_T)  , intent(in)               :: myfile_id
       character(len=*), intent(in)               :: mydata
       INTEGER(HSIZE_T), DIMENSION(*), intent(in) :: dims
       real(kind=real_kind),           intent(in),target :: data_in
!------------------------------------------------------------------------
       integer(c_int)                                  :: data_rank
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound
       !>

       f_ptr = C_LOC(data_in)
       error = hdf5_put_data_c(myfile_id, MYH5T_DOUBLE, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_put_data_dpscalar
!------------------------------------------------------------------------

      subroutine hdf5_put_data_dp(myfile_id, mydata_in, dims, data_in)

       integer(HID_T)  , intent(in)               :: myfile_id
       character(len=*), intent(in)               :: mydata_in
       INTEGER(HSIZE_T), DIMENSION(*), intent(in) :: dims
       real(kind=real_kind), DIMENSION(*), intent(in),target :: data_in
!------------------------------------------------------------------------
       integer(c_int)                                  :: data_rank
       character(len=taglenx)                     :: mydata
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       !> character string with integral pointers is encoded with XXXXXX
       if(mydata_in == "XXXXXX")then
         mydata = trim(datatag(1:tagx))
         mydata = trim(mydata)//'d'
       else
         mydata = trim(mydata_in)
       end if

       f_ptr = C_LOC(data_in)
       error = hdf5_put_data_c(myfile_id, MYH5T_DOUBLE, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_put_data_dp
!------------------------------------------------------------------------

      subroutine hdf5_get_data_dpscalar(myfile_id, mydata, dims, data_out)

       integer(HID_T)  , intent(in)                :: myfile_id
       character(len=*), intent(in)                :: mydata
       INTEGER(HSIZE_T), DIMENSION(*), intent(in)  :: dims
       real(kind=real_kind),           intent(out),target :: data_out
!------------------------------------------------------------------------
       integer(c_int)                                   :: data_rank
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       f_ptr = C_LOC(data_out)
       error = hdf5_get_data_c(myfile_id, MYH5T_DOUBLE, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_get_data_dpscalar
!------------------------------------------------------------------------

      subroutine hdf5_get_data_dp(myfile_id, mydata_in, dims, data_out)

       integer(HID_T)  , intent(in)                :: myfile_id
       character(len=*), intent(in)                :: mydata_in
       INTEGER(HSIZE_T), DIMENSION(*), intent(in)  :: dims
       real(kind=real_kind), DIMENSION(*), intent(out),target :: data_out
!------------------------------------------------------------------------
       integer(c_int)                              :: data_rank
       character(len=taglenx)                      :: mydata
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       !> character string with integral pointers is encoded with XXXXXX
       if(mydata_in == "XXXXXX")then
         mydata = trim(datatag(1:tagx))
         mydata = trim(mydata)//'d'
       else
         mydata = trim(mydata_in)
       end if

       f_ptr = C_LOC(data_out)
       error = hdf5_get_data_c(myfile_id, MYH5T_DOUBLE, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_get_data_dp
!------------------------------------------------------------------------

      subroutine hdf5_put_data_iscalar(myfile_id, mydata_in, dims, data_in)

       integer(HID_T)  , intent(in)                      :: myfile_id
       character(len=*), intent(in)                      :: mydata_in
       INTEGER(HSIZE_T), DIMENSION(*), intent(in)        :: dims
       integer*8,         intent(in),target :: data_in
!------------------------------------------------------------------------
       integer(c_int)                                    :: data_rank
       character(len=taglenx)                            :: mydata
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       !> character string with integral pointers is encoded with XXXXXX
       if(mydata_in == "XXXXXX")then
         mydata = trim(datatag(1:tagx))
         mydata = trim(mydata)//'i'
       else
         mydata = trim(mydata_in)
       end if


       f_ptr = C_LOC(data_in)
       error = hdf5_put_data_c(myfile_id, MYH5T_INT, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_put_data_iscalar
!------------------------------------------------------------------------

      subroutine hdf5_put_data_ix(myfile_id, mydata, dims, data_in)

       integer(HID_T)  , intent(in)                      :: myfile_id
       character(len=*), intent(in)                      :: mydata
       INTEGER(HSIZE_T), DIMENSION(*), intent(in)        :: dims
       integer*8             , DIMENSION(*), intent(in),target :: data_in
!------------------------------------------------------------------------
       integer(c_int)                                    :: data_rank
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       f_ptr = C_LOC(data_in)
       error = hdf5_put_data_c(myfile_id, MYH5T_INT, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_put_data_ix
!------------------------------------------------------------------------

      subroutine hdf5_get_data_iscalar(myfile_id, mydata_in, dims, data_out)

       integer(HID_T)  , intent(in)                       :: myfile_id
       character(len=*), intent(in)                       :: mydata_in
       INTEGER(HSIZE_T), DIMENSION(*), intent(in)         :: dims
       integer*8             ,         intent(out),target :: data_out
!------------------------------------------------------------------------
       integer(c_int)                                     :: data_rank
       character(len=taglenx)                             :: mydata
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       !> character string with integral pointers is encoded with XXXXXX
       if(mydata_in == "XXXXXX")then
         mydata = trim(datatag(1:tagx))
         mydata = trim(mydata)//'i'
       else
         mydata = trim(mydata_in)
       end if

       f_ptr = C_LOC(data_out)
       error = hdf5_get_data_c(myfile_id, MYH5T_INT, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_get_data_iscalar
!------------------------------------------------------------------------

      subroutine hdf5_get_data_ix(myfile_id, mydata, dims, data_out)

       integer(HID_T)  , intent(in)                       :: myfile_id
       character(len=*), intent(in)                       :: mydata
       INTEGER(HSIZE_T), DIMENSION(*), intent(in)         :: dims
       integer*8             , DIMENSION(*), intent(out),target :: data_out
!------------------------------------------------------------------------
       integer(c_int)                                     :: data_rank
!------------------------------------------------------------------------

       !> possibly use internal conversion from integer*8 to *4
       data_rank = datadim_bound

       f_ptr = C_LOC(data_out)
       error = hdf5_get_data_c(myfile_id, MYH5T_INT, trim(mydata)//C_NULL_CHAR, data_rank, dims, f_ptr)

      end subroutine hdf5_get_data_ix
!------------------------------------------------------------------------

end module hdf5_utils
