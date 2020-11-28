************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
#include "molcastypes.fh"
#define MH5_MAX_LBL_LEN 256

*     create a HDF5 file and return a handle to it
      function mh5_create_file (filename) result(lu)
      use iso_c_binding
      implicit none
      character(len=*) :: filename
      integer :: lu
      interface
        function mh5c_create_file(filename) result(lu)
     &   bind(C, name='mh5c_create_file')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT) :: lu
        character(kind=C_CHAR) :: filename(*)
        end function
      end interface
      integer :: lrealname
      character(len=4096) :: realname
      call prgmtranslate(filename, realname, lrealname)
      realname(lrealname+1:lrealname+1)=C_NULL_CHAR
      lu = mh5c_create_file(realname)
      end function

      function mh5_open_file_rw (filename) result (lu)
      use iso_c_binding
      implicit none
      character(len=*) :: filename
      integer :: lu
      interface
        function mh5c_open_rw(filename) result(lu)
     &   bind(C, name='mh5c_open_file_rw')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT) :: lu
        character(kind=C_CHAR) :: filename(*)
        end function
      end interface
      integer :: lrealname
      character(len=4096) :: realname
      call prgmtranslate(filename, realname, lrealname)
      realname(lrealname+1:lrealname+1)=C_NULL_CHAR
      lu = mh5c_open_rw(realname)
      end function

      function mh5_open_file_r (filename) result (lu)
      use iso_c_binding
      implicit none
      character(len=*) :: filename
      integer :: lu
      interface
        function mh5c_open_r(filename) result(lu)
     &   bind(C, name='mh5c_open_file_r')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT) :: lu
        character(kind=C_CHAR) :: filename(*)
        end function
      end interface
      integer :: lrealname
      character(len=4096) :: realname
      call prgmtranslate(filename, realname, lrealname)
      realname(lrealname+1:lrealname+1)=C_NULL_CHAR
      lu = mh5c_open_r(realname)
      end function

      subroutine mh5_close_file (lu)
      use iso_c_binding
      implicit none
      integer :: lu
      interface
        function mh5c_close_file(lu) result(ierr)
     &   bind(C, name='mh5c_close_file')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      integer :: ierr
      ierr = mh5c_close_file(lu)
      end subroutine

*     check if a file is in the HDF5 format
      logical function mh5_is_hdf5 (filename)
      use iso_c_binding
      implicit none
      character(len=*) :: filename
      interface
        function mh5c_is_hdf5(filename) result (rc)
     &   bind(C, name='mh5c_is_hdf5')
        use iso_c_binding
        implicit none
        character(kind=c_char) :: filename(*)
        integer(MOLCAS_C_INT) :: rc
        end function
      end interface
      integer :: rc
      integer :: lrealname
      character(len=4096) :: realname
      logical :: exists
      call prgmtranslate(filename, realname, lrealname)
      call f_inquire(realname,exists)
      if (exists) then
        realname(lrealname+1:lrealname+1)=C_NULL_CHAR
        rc = mh5c_is_hdf5(realname)
      else
        rc = 0
      end if
      if      (rc > 0) then
        mh5_is_hdf5 = .true.
      else if (rc == 0) then
        mh5_is_hdf5 = .false.
      else
        mh5_is_hdf5 = .false.
        call abend
      end if
      end function

      function mh5_open_group (lu, groupname) result (groupid)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: groupname
      integer :: groupid
      interface
        function mh5c_open_group(lu, groupname) result(groupid)
     &   bind(C, name='mh5c_open_group')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT) :: groupid
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=C_CHAR) :: groupname(*)
        end function
      end interface
      call f2c_string(groupname,mh5_lbl)
      groupid = mh5c_open_group(lu,mh5_lbl)
      end function

      subroutine mh5_close_group (id)
      use iso_c_binding
      implicit none
      integer :: id
      interface
        function mh5c_close_group(id) result(ierr)
     &   bind(C, name='mh5c_close_group')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      integer :: ierr
      ierr = mh5c_close_group(id)
      end subroutine

*     check for existence of dataset/attribute by id,
*     where id could be a file or dataset id.
      logical function mh5_exists_dset (id, name)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: id
      character(len=*) :: name
      interface
        function mh5c_exists_dset(id, name) result (rc)
     &   bind(C, name='mh5c_exists_dset')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(len=c_char) :: name(*)
        integer(MOLCAS_C_INT) :: rc
        end function
      end interface
      integer :: rc
      call f2c_string(name,mh5_lbl)
      rc = mh5c_exists_dset(id, mh5_lbl)
      if      (rc > 0) then
        mh5_exists_dset = .true.
      else if (rc == 0) then
        mh5_exists_dset = .false.
      else
        mh5_exists_dset = .false.
        call abend
      end if
      end function

      logical function mh5_exists_attr (id, name)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: id
      character(len=*) :: name
      interface
        function mh5c_exists_attr(id, name) result (rc)
     &   bind(C, name='mh5c_exists_attr')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT) :: rc
        end function
      end interface
      integer :: rc
      call f2c_string(name,mh5_lbl)
      rc = mh5c_exists_attr(id, mh5_lbl)
      if      (rc > 0) then
        mh5_exists_attr = .true.
      else if (rc == 0) then
        mh5_exists_attr = .false.
      else
        mh5_exists_attr = .false.
        call abend
      end if
      end function
*     open/close dataset

      function mh5_open_dset (lu, dsetname) result(dsetid)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: dsetname
      integer :: dsetid
      interface
        function mh5c_open_dset(lu, dsetname) result(dsetid)
     &   bind(C, name='mh5c_open_dset')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT) :: dsetid
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: dsetname(*)
        end function
      end interface

      call f2c_string(dsetname,mh5_lbl)
      dsetid = mh5c_open_dset(lu,mh5_lbl)
      end function

      subroutine mh5_close_dset (dsetid)
      use iso_c_binding
      implicit none
      integer :: dsetid
      interface
        function mh5c_close_dset(dsetid) result(ierr)
     &   bind(C, name='mh5c_close_dset')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dsetid
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (mh5c_close_dset(dsetid) < 0) call abend
      end subroutine

      function mh5_open_attr (lu, attrname) result(attrid)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: attrname
      integer :: attrid
      interface
        function mh5c_open_attr(lu, attrname) result(attrid)
     &   bind(C, name='mh5c_open_attr')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT) :: attrid
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: attrname(*)
        end function
      end interface

      call f2c_string(attrname,mh5_lbl)
      attrid = mh5c_open_attr(lu,mh5_lbl)
      end function

      subroutine mh5_close_attr (attrid)
      use iso_c_binding
      implicit none
      integer :: attrid
      interface
        function mh5c_close_attr(attrid) result(ierr)
     &   bind(C, name='mh5c_close_attr')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: attrid
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (mh5c_close_attr(attrid) < 0) call abend
      end subroutine

*======================
*     ATTRIBUTES
*======================

      function mh5_create_attr_scalar_int (lu, name) result (dset_id)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: name
      integer :: dset_id
      interface
        function mh5c_create_attr_scalar_int(lu, name) result (dset_id)
     &   bind(C, name='mh5c_create_attr_scalar_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT) :: dset_id
        end function
      end interface

      call f2c_string(name,mh5_lbl)
      dset_id = mh5c_create_attr_scalar_int(lu, mh5_lbl)
      end function

      function mh5_create_attr_scalar_real (lu, name) result (dset_id)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: name
      integer :: dset_id
      interface
        function mh5c_create_attr_scalar_real(lu, name) result (dset_id)
     &   bind(C, name='mh5c_create_attr_scalar_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT) :: dset_id
        end function
      end interface

      call f2c_string(name,mh5_lbl)
      dset_id = mh5c_create_attr_scalar_real(lu, mh5_lbl)
      end function

      function mh5_create_attr_scalar_str (lu, name, size)
     &        result (dset_id)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: name
      integer :: size
      integer :: dset_id
      interface
        function mh5c_create_attr_scalar_str(lu, name, size)
     &          result (dset_id)
     &          bind(C, name='mh5c_create_attr_scalar_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: size
        integer(MOLCAS_C_INT) :: dset_id
        end function
      end interface

      call f2c_string(name,mh5_lbl)
      dset_id = mh5c_create_attr_scalar_str(lu, mh5_lbl, size)
      end function

      subroutine mh5_put_attr_scalar_int (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      integer :: value
      integer :: ierr
      interface
        function mh5c_put_attr_scalar_int(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_put_attr_scalar_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        integer(MOLCAS_C_INT) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_attr_scalar_int(dset_id, value)
      end subroutine

      subroutine mh5_put_attr_scalar_real (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      real*8 :: value
      integer :: ierr
      interface
        function mh5c_put_attr_scalar_real(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_put_attr_scalar_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        real(MOLCAS_C_REAL) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_attr_scalar_real(dset_id, value)
      end subroutine

      subroutine mh5_put_attr_scalar_str (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      character(len=*) :: value
      integer :: ierr
      interface
        function mh5c_put_attr_scalar_str(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_put_attr_scalar_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        character(kind=c_char) :: value(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_attr_scalar_str(dset_id, value)
      end subroutine


      subroutine mh5_get_attr_scalar_int (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      integer :: value
      integer :: ierr
      interface
        function mh5c_get_attr_scalar_int(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_get_attr_scalar_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        integer(MOLCAS_C_INT) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_attr_scalar_int(dset_id, value)
      end subroutine

      subroutine mh5_get_attr_scalar_real (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      real*8 :: value
      integer :: ierr
      interface
        function mh5c_get_attr_scalar_real(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_get_attr_scalar_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        real(MOLCAS_C_REAL) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_attr_scalar_real(dset_id, value)
      end subroutine

      subroutine mh5_get_attr_scalar_str (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      character(len=*) :: value
      integer :: ierr
      interface
        function mh5c_get_attr_scalar_str(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_get_attr_scalar_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        character(kind=c_char) :: value(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_attr_scalar_str(dset_id, value)
      end subroutine

      function mh5_create_attr_array_int (lu, name, rank, dims)
     &        result(attr_id)
      use iso_c_binding
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      integer :: attr_id
      interface
        function mh5c_create_attr_array_int(lu, name, rank, dims)
     &   result(ierr) bind(C, name='mh5c_create_attr_array_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl

      call f2c_string(name,mh5_lbl)
      attr_id = mh5c_create_attr_array_int(lu, mh5_lbl, rank, dims)
      end function

      function mh5_create_attr_array_real (lu, name, rank, dims)
     &        result(attr_id)
      use iso_c_binding
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      integer :: attr_id
      interface
        function mh5c_create_attr_array_real(lu, name, rank, dims)
     &   result(ierr) bind(C, name='mh5c_create_attr_array_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl

      call f2c_string(name,mh5_lbl)
      attr_id = mh5c_create_attr_array_real(lu, mh5_lbl, rank, dims)
      end function

      function mh5_create_attr_array_str (lu, name, rank, dims, size)
     &        result(attr_id)
      use iso_c_binding
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      integer :: size
      integer :: attr_id
      interface
        function mh5c_create_attr_array_str(lu, name, rank, dims, size)
     &   result(ierr) bind(C, name='mh5c_create_attr_array_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT), VALUE :: size
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl

      call f2c_string(name,mh5_lbl)
      attr_id = mh5c_create_attr_array_str(
     &        lu, mh5_lbl, rank, dims, size)
      end function


      subroutine mh5_put_attr_array_int (attr_id, buffer)
      use iso_c_binding
      implicit none
      integer :: attr_id
      integer :: buffer(*)
      integer :: ierr
      interface
        function mh5c_put_attr_array_int(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_attr_array_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_attr_array_int(attr_id, buffer)
      if (ierr < 0) call abend
      end subroutine

      subroutine mh5_put_attr_array_real (attr_id, buffer)
      use iso_c_binding
      implicit none
      integer :: attr_id
      real*8 :: buffer(*)
      integer :: ierr
      interface
        function mh5c_put_attr_array_real(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_attr_array_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_attr_array_real(attr_id, buffer)
      if (ierr < 0) call abend
      end subroutine

      subroutine mh5_put_attr_array_str (attr_id, buffer)
      use iso_c_binding
      implicit none
      integer :: attr_id
      character :: buffer(*)
      integer :: ierr
      interface
        function mh5c_put_attr_array_str(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_attr_array_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_attr_array_str(attr_id, buffer)
      if (ierr < 0) call abend
      end subroutine


      subroutine mh5_get_attr_array_int (attr_id, buffer)
      use iso_c_binding
      implicit none
      integer :: attr_id
      integer :: buffer(*)
      integer :: ierr
      interface
        function mh5c_get_attr_array_int(id, buffer)
     &   result(ierr) bind(C, name='mh5c_get_attr_array_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_attr_array_int(attr_id, buffer)
      if (ierr < 0) call abend
      end subroutine

      subroutine mh5_get_attr_array_real (attr_id, buffer)
      use iso_c_binding
      implicit none
      integer :: attr_id
      real*8 :: buffer(*)
      integer :: ierr
      interface
        function mh5c_get_attr_array_real(id, buffer)
     &   result(ierr) bind(C, name='mh5c_get_attr_array_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_attr_array_real(attr_id, buffer)
      if (ierr < 0) call abend
      end subroutine

      subroutine mh5_get_attr_array_str (attr_id, buffer)
      use iso_c_binding
      implicit none
      integer :: attr_id
      character :: buffer(*)
      integer :: ierr
      interface
        function mh5c_get_attr_array_str(id, buffer)
     &   result(ierr) bind(C, name='mh5c_get_attr_array_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_attr_array_str(attr_id, buffer)
      if (ierr < 0) call abend
      end subroutine

* Convenience wrappers: 'init' and 'fetch'

* init: create, put, close
      subroutine mh5_init_attr_scalar_int (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: value
      integer :: attr_id
      integer :: mh5_create_attr_scalar_int
      attr_id = mh5_create_attr_scalar_int(lu, name)
      call mh5_put_attr_scalar_int(attr_id, value)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_init_attr_scalar_real (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      real*8 :: value
      integer :: attr_id
      integer :: mh5_create_attr_scalar_real
      attr_id = mh5_create_attr_scalar_real(lu, name)
      call mh5_put_attr_scalar_real(attr_id, value)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_init_attr_scalar_str (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      character(len=*) :: value
      integer :: attr_id
      integer :: mh5_create_attr_scalar_str
      attr_id = mh5_create_attr_scalar_str(lu, name, len(value))
      call mh5_put_attr_scalar_str(attr_id, value)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_init_attr_array_int (lu, name, rank, dims, buffer)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      integer :: buffer(*)
      integer :: attr_id
      integer :: mh5_create_attr_array_int
      attr_id = mh5_create_attr_array_int(lu, name, rank, dims)
      call mh5_put_attr_array_int(attr_id, buffer)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_init_attr_array_real(lu, name, rank, dims, buffer)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      real*8 :: buffer(*)
      integer :: attr_id
      integer :: mh5_create_attr_array_real
      attr_id = mh5_create_attr_array_real(lu, name, rank, dims)
      call mh5_put_attr_array_real(attr_id, buffer)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_init_attr_array_str (
     &        lu, name, rank, dims, buffer, size)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      character :: buffer(*)
      integer :: size
      integer :: attr_id
      integer :: mh5_create_attr_array_str
      attr_id = mh5_create_attr_array_str(lu, name, rank, dims, size)
      call mh5_put_attr_array_str(attr_id, buffer)
      call mh5_close_attr(attr_id)
      end subroutine

* fetch: open, get, close
      subroutine mh5_fetch_attr_scalar_int (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: value
      integer :: attr_id
      integer :: mh5_open_attr
      attr_id = mh5_open_attr(lu, name)
      call mh5_get_attr_scalar_int(attr_id, value)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_fetch_attr_scalar_real (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      real*8 :: value
      integer :: attr_id
      integer :: mh5_open_attr
      attr_id = mh5_open_attr(lu, name)
      call mh5_get_attr_scalar_real(attr_id, value)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_fetch_attr_scalar_str (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      character(len=*) :: value
      integer :: attr_id
      integer :: mh5_open_attr
      attr_id = mh5_open_attr(lu, name)
      call mh5_get_attr_scalar_str(attr_id, value)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_fetch_attr_array_int (lu, name, buffer)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: buffer(*)
      integer :: attr_id
      integer :: mh5_open_attr
      attr_id = mh5_open_attr(lu, name)
      call mh5_get_attr_array_int(attr_id, buffer)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_fetch_attr_array_real(lu, name, buffer)
      implicit none
      integer :: lu
      character(len=*) :: name
      real*8 :: buffer(*)
      integer :: attr_id
      integer :: mh5_open_attr
      attr_id = mh5_open_attr(lu, name)
      call mh5_get_attr_array_real(attr_id, buffer)
      call mh5_close_attr(attr_id)
      end subroutine

      subroutine mh5_fetch_attr_array_str (lu, name, buffer)
      implicit none
      integer :: lu
      character(len=*) :: name
      character :: buffer(*)
      integer :: attr_id
      integer :: mh5_open_attr
      attr_id = mh5_open_attr(lu, name)
      call mh5_get_attr_array_str(attr_id, buffer)
      call mh5_close_attr(attr_id)
      end subroutine

*====================
*     DATASETS
*====================

      function mh5_create_dset_scalar_int (lu, name) result (dset_id)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: name
      integer :: dset_id
      interface
        function mh5c_create_dset_scalar_int(lu, name) result (dset_id)
     &   bind(C, name='mh5c_create_dset_scalar_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT) :: dset_id
        end function
      end interface

      call f2c_string(name,mh5_lbl)
      dset_id = mh5c_create_dset_scalar_int(lu, mh5_lbl)
      end function

      function mh5_create_dset_scalar_real (lu, name) result (dset_id)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: name
      integer :: dset_id
      interface
        function mh5c_create_dset_scalar_real(lu, name) result (dset_id)
     &   bind(C, name='mh5c_create_dset_scalar_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT) :: dset_id
        end function
      end interface

      call f2c_string(name,mh5_lbl)
      dset_id = mh5c_create_dset_scalar_real(lu, mh5_lbl)
      end function

      function mh5_create_dset_scalar_str (lu, name, size)
     &        result (dset_id)
      use iso_c_binding
      implicit none
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      integer :: lu
      character(len=*) :: name
      integer :: size
      integer :: dset_id
      interface
        function mh5c_create_dset_scalar_str(lu, name, size)
     &          result (dset_id)
     &          bind(C, name='mh5c_create_dset_scalar_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: size
        integer(MOLCAS_C_INT) :: dset_id
        end function
      end interface

      call f2c_string(name,mh5_lbl)
      dset_id = mh5c_create_dset_scalar_str(lu, mh5_lbl, size)
      end function

      subroutine mh5_put_dset_scalar_int (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      integer :: value
      integer :: ierr
      interface
        function mh5c_put_dset_scalar_int(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_put_dset_scalar_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        integer(MOLCAS_C_INT) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_dset_scalar_int(dset_id, value)
      end subroutine

      subroutine mh5_put_dset_scalar_real (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      real*8 :: value
      integer :: ierr
      interface
        function mh5c_put_dset_scalar_real(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_put_dset_scalar_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        real(MOLCAS_C_REAL) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_dset_scalar_real(dset_id, value)
      end subroutine

      subroutine mh5_put_dset_scalar_str (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      character(len=*) :: value
      integer :: ierr
      interface
        function mh5c_put_dset_scalar_str(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_put_dset_scalar_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        character(kind=c_char) :: value(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_put_dset_scalar_str(dset_id, value)
      end subroutine


      subroutine mh5_get_dset_scalar_int (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      integer :: value
      integer :: ierr
      interface
        function mh5c_get_dset_scalar_int(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_get_dset_scalar_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        integer(MOLCAS_C_INT) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_dset_scalar_int(dset_id, value)
      end subroutine

      subroutine mh5_get_dset_scalar_real (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      real*8 :: value
      integer :: ierr
      interface
        function mh5c_get_dset_scalar_real(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_get_dset_scalar_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        real(MOLCAS_C_REAL) :: value
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_dset_scalar_real(dset_id, value)
      end subroutine

      subroutine mh5_get_dset_scalar_str (dset_id, value)
      use iso_c_binding
      implicit none
      integer :: dset_id
      character(len=*) :: value
      integer :: ierr
      interface
        function mh5c_get_dset_scalar_str(dset_id, value) result (ierr)
     &   bind(C, name='mh5c_get_dset_scalar_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: dset_id
        character(kind=c_char) :: value(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_dset_scalar_str(dset_id, value)
      end subroutine

      function mh5_create_dset_array_int (lu, name, rank, dims, dyn)
     &        result(dset_id)
      use iso_c_binding
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      logical, optional :: dyn
      integer :: dset_id
      interface
        function mh5c_create_dset_array_int(lu, name, rank, dims)
     &   result(ierr) bind(C, name='mh5c_create_dset_array_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_create_dset_array_dyn_int(lu, name, rank, dims)
     &   result(ierr) bind(C, name='mh5c_create_dset_array_dyn_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      logical :: isdyn

      isdyn = .false.
      if (present(dyn)) isdyn = dyn

      call f2c_string(name,mh5_lbl)
      if (isdyn) then
        dset_id = mh5c_create_dset_array_dyn_int(lu, mh5_lbl, rank,
     &                                           dims)
      else
        dset_id = mh5c_create_dset_array_int(lu, mh5_lbl, rank, dims)
      end if
      end function

      function mh5_create_dset_array_real (lu, name, rank, dims, dyn)
     &        result(dset_id)
      use iso_c_binding
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      logical, optional :: dyn
      integer :: dset_id
      interface
        function mh5c_create_dset_array_real(lu, name, rank, dims)
     &   result(ierr) bind(C, name='mh5c_create_dset_array_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_create_dset_array_dyn_real(lu, name, rank, dims)
     &   result(ierr) bind(C, name='mh5c_create_dset_array_dyn_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      logical :: isdyn

      isdyn = .false.
      if (present(dyn)) isdyn = dyn

      call f2c_string(name,mh5_lbl)
      if (isdyn) then
        dset_id = mh5c_create_dset_array_dyn_real(lu, mh5_lbl, rank,
     &                                            dims)
      else
        dset_id = mh5c_create_dset_array_real(lu, mh5_lbl, rank, dims)
      end if
      end function

      function mh5_create_dset_array_str (lu, name, rank, dims, size,
     &                                    dyn)
     &        result(dset_id)
      use iso_c_binding
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      integer :: size
      logical, optional :: dyn
      integer :: dset_id
      interface
        function mh5c_create_dset_array_str(lu, name, rank, dims, size)
     &   result(ierr) bind(C, name='mh5c_create_dset_array_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT), VALUE :: size
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_create_dset_array_dyn_str(lu, name, rank, dims,
     &                                          size)
     &   result(ierr) bind(C, name='mh5c_create_dset_array_dyn_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: lu
        character(kind=c_char) :: name(*)
        integer(MOLCAS_C_INT), VALUE :: rank
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT), VALUE :: size
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
      logical :: isdyn

      isdyn = .false.
      if (present(dyn)) isdyn = dyn

      call f2c_string(name,mh5_lbl)
      if (isdyn) then
        dset_id = mh5c_create_dset_array_dyn_str(
     &          lu, mh5_lbl, rank, dims, size)
      else
        dset_id = mh5c_create_dset_array_str(
     &          lu, mh5_lbl, rank, dims, size)
      end if
      end function

      subroutine mh5_put_dset_array_int (dset_id, buffer, exts, offs)
      use iso_c_binding
      implicit none
      integer :: dset_id
      integer :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: ierr
      interface
        function mh5c_put_dset_array_int(id, exts, offs, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: exts(*), offs(*)
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_put_dset_array_int_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_int_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (present(exts).and.present(offs)) then
        ierr = mh5c_put_dset_array_int(dset_id, exts, offs, buffer)
      else if (present(exts).or.present(offs)) then
        ierr = -1
        call abend
      else
        ierr = mh5c_put_dset_array_int_full(dset_id, buffer)
      end if
      end subroutine

      subroutine mh5_put_dset_array_real (dset_id, buffer, exts, offs)
      use iso_c_binding
      implicit none
      integer :: dset_id
      real*8 :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: ierr
      interface
        function mh5c_put_dset_array_real(id, exts, offs, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: exts(*), offs(*)
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_put_dset_array_real_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_real_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (present(exts).and.present(offs)) then
        ierr = mh5c_put_dset_array_real(dset_id, exts, offs, buffer)
      else if (present(exts).or.present(offs)) then
        ierr = -1
        call abend
      else
        ierr = mh5c_put_dset_array_real_full(dset_id, buffer)
      end if
      end subroutine

      subroutine mh5_put_dset_array_str (dset_id, buffer, exts, offs)
      use iso_c_binding
      implicit none
      integer :: dset_id
      character :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: ierr
      interface
        function mh5c_put_dset_array_str(id, exts, offs, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: exts(*), offs(*)
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_put_dset_array_str_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_str_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (present(exts).and.present(offs)) then
        ierr = mh5c_put_dset_array_str(dset_id, exts, offs, buffer)
      else if (present(exts).or.present(offs)) then
        ierr = -1
        call abend
      else
        ierr = mh5c_put_dset_array_str_full(dset_id, buffer)
      end if
      end subroutine


      subroutine mh5_get_dset_array_int (dset_id, buffer, exts, offs)
      use iso_c_binding
      implicit none
      integer :: dset_id
      integer :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: ierr
      interface
        function mh5c_get_dset_array_int(id, exts, offs, buffer)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_int')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: exts(*), offs(*)
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_get_dset_array_int_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_int_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (present(exts).and.present(offs)) then
        ierr = mh5c_get_dset_array_int(dset_id, exts, offs, buffer)
      else if (present(exts).or.present(offs)) then
        ierr = -1
        call abend
      else
        ierr = mh5c_get_dset_array_int_full(dset_id, buffer)
      end if
      end subroutine

      subroutine mh5_get_dset_array_real (dset_id, buffer, exts, offs)
      use iso_c_binding
      implicit none
      integer :: dset_id
      real*8 :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: ierr
      interface
        function mh5c_get_dset_array_real(id, exts, offs, buffer)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_real')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: exts(*), offs(*)
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_get_dset_array_real_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_real_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (present(exts).and.present(offs)) then
        ierr = mh5c_get_dset_array_real(dset_id, exts, offs, buffer)
      else if (present(exts).or.present(offs)) then
        ierr = -1
        call abend
      else
        ierr = mh5c_get_dset_array_real_full(dset_id, buffer)
      end if
      end subroutine

      subroutine mh5_get_dset_array_str (dset_id, buffer, exts, offs)
      use iso_c_binding
      implicit none
      integer :: dset_id
      character :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: ierr
      interface
        function mh5c_get_dset_array_str(id, exts, offs, buffer)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_str')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: exts(*), offs(*)
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
        function mh5c_get_dset_array_str_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_str_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      if (present(exts).and.present(offs)) then
        ierr = mh5c_get_dset_array_str(dset_id, exts, offs, buffer)
      else if (present(exts).or.present(offs)) then
        ierr = -1
        call abend
      else
        ierr = mh5c_get_dset_array_str_full(dset_id, buffer)
      end if
      end subroutine

      subroutine mh5_extend_dset_array(dset_id, dims)
      implicit none
      integer :: dset_id
      integer :: dims(*)
      integer :: ierr
      interface
        function mh5c_extend_dset_array(id, dims)
     &   result(ierr) bind(C, name='mh5c_extend_dset_array')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_extend_dset_array(dset_id, dims)
      if (ierr .lt. 0) then
        call abend
      end if
      end subroutine

      subroutine mh5_get_dset_array_dims(dset_id, dims)
      implicit none
      integer :: dset_id
      integer :: dims(*)
      integer :: ierr
      interface
        function mh5c_get_dset_array_dims(id, dims)
     &   result(ierr) bind(C, name='mh5c_get_dset_array_dims')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: dims(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface

      ierr = mh5c_get_dset_array_dims(dset_id, dims)
      if (ierr .lt. 0) then
        call abend
      end if
      end subroutine

* Convenience wrappers: 'init' and 'fetch'

* init: create, put, close
      subroutine mh5_init_dset_scalar_int (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: value
      integer :: dset_id
      integer :: mh5_create_dset_scalar_int
      dset_id = mh5_create_dset_scalar_int(lu, name)
      call mh5_put_dset_scalar_int(dset_id, value)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_init_dset_scalar_real (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      real*8 :: value
      integer :: dset_id
      integer :: mh5_create_dset_scalar_real
      dset_id = mh5_create_dset_scalar_real(lu, name)
      call mh5_put_dset_scalar_real(dset_id, value)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_init_dset_scalar_str (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      character(len=*) :: value
      integer :: dset_id
      integer :: mh5_create_dset_scalar_str
      dset_id = mh5_create_dset_scalar_str(lu, name, len(value))
      call mh5_put_dset_scalar_str(dset_id, value)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_init_dset_array_int (lu, name, rank, dims, buffer,
     &                                    dyn)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      integer :: buffer(*)
      logical, optional :: dyn
      integer :: ierr
      integer :: dset_id
      interface
        function mh5_create_dset_array_int (lu, name, rank, dims, dyn)
     &   result(dset_id)
        use iso_c_binding
        implicit none
        integer :: lu
        character(len=*) :: name
        integer :: rank
        integer :: dims(*)
        logical, optional :: dyn
        integer :: dset_id
        end function
        function mh5c_put_dset_array_int_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_int_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        integer(MOLCAS_C_INT) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      logical :: isdyn
      isdyn = .false.
      if (present(dyn)) isdyn = dyn
      dset_id = mh5_create_dset_array_int(lu, name, rank, dims, isdyn)
      ierr = mh5c_put_dset_array_int_full(dset_id, buffer)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_init_dset_array_real(lu, name, rank, dims, buffer,
     &                                    dyn)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      real*8 :: buffer(*)
      logical, optional :: dyn
      integer :: ierr
      integer :: dset_id
      interface
        function mh5_create_dset_array_real (lu, name, rank, dims, dyn)
     &   result(dset_id)
        use iso_c_binding
        implicit none
        integer :: lu
        character(len=*) :: name
        integer :: rank
        integer :: dims(*)
        logical, optional :: dyn
        integer :: dset_id
        end function
        function mh5c_put_dset_array_real_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_real_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        real(MOLCAS_C_REAL) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      logical :: isdyn
      isdyn = .false.
      if (present(dyn)) isdyn = dyn
      dset_id = mh5_create_dset_array_real(lu, name, rank, dims, dyn)
      ierr = mh5c_put_dset_array_real_full(dset_id, buffer)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_init_dset_array_str(
     &        lu, name, rank, dims, buffer, size, dyn)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: rank
      integer :: dims(*)
      character :: buffer(*)
      integer :: size
      logical, optional :: dyn
      integer :: ierr
      integer :: dset_id
      interface
        function mh5_create_dset_array_str (lu, name, rank, dims, size,
     &                                      dyn)
     &   result(dset_id)
        use iso_c_binding
        implicit none
        integer :: lu
        character(len=*) :: name
        integer :: rank
        integer :: dims(*)
        integer :: size
        logical, optional :: dyn
        integer :: dset_id
        end function
        function mh5c_put_dset_array_str_full(id, buffer)
     &   result(ierr) bind(C, name='mh5c_put_dset_array_str_full')
        use iso_c_binding
        implicit none
        integer(MOLCAS_C_INT), VALUE :: id
        character(kind=c_char) :: buffer(*)
        integer(MOLCAS_C_INT) :: ierr
        end function
      end interface
      logical :: isdyn
      isdyn = .false.
      if (present(dyn)) isdyn = dyn
      dset_id = mh5_create_dset_array_str(lu, name, rank, dims, size,
     &                                    isdyn)
      ierr = mh5c_put_dset_array_str_full(dset_id, buffer)
      call mh5_close_dset(dset_id)
      end subroutine

* fetch: open, get, close
      subroutine mh5_fetch_dset_scalar_int (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: value
      integer :: dset_id
      integer :: mh5_open_dset
      dset_id = mh5_open_dset(lu, name)
      call mh5_get_dset_scalar_int(dset_id, value)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_fetch_dset_scalar_real (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      real*8 :: value
      integer :: dset_id
      integer :: mh5_open_dset
      dset_id = mh5_open_dset(lu, name)
      call mh5_get_dset_scalar_real(dset_id, value)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_fetch_dset_scalar_str (lu, name, value)
      implicit none
      integer :: lu
      character(len=*) :: name
      character(len=*) :: value
      integer :: dset_id
      integer :: mh5_open_dset
      dset_id = mh5_open_dset(lu, name)
      call mh5_get_dset_scalar_str(dset_id, value)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_fetch_dset_array_int (lu, name, buffer, exts, offs)
      implicit none
      integer :: lu
      character(len=*) :: name
      integer :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: dset_id
      integer :: mh5_open_dset
      interface
        subroutine mh5_get_dset_array_int (dset_id, buffer, exts, offs)
        implicit none
        integer :: dset_id
        integer :: buffer(*)
        integer, optional :: exts(*), offs(*)
        end subroutine
      end interface
      dset_id = mh5_open_dset(lu, name)
      call mh5_get_dset_array_int(dset_id, buffer, exts, offs)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_fetch_dset_array_real(lu, name, buffer, exts, offs)
      implicit none
      integer :: lu
      character(len=*) :: name
      real*8 :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: dset_id
      integer :: mh5_open_dset
      interface
        subroutine mh5_get_dset_array_real (dset_id, buffer, exts, offs)
        implicit none
        integer :: dset_id
        real*8 :: buffer(*)
        integer, optional :: exts(*), offs(*)
        end subroutine
      end interface
      dset_id = mh5_open_dset(lu, name)
      call mh5_get_dset_array_real(dset_id, buffer, exts, offs)
      call mh5_close_dset(dset_id)
      end subroutine

      subroutine mh5_fetch_dset_array_str(lu, name, buffer, exts, offs)
      implicit none
      integer :: lu
      character(len=*) :: name
      character :: buffer(*)
      integer, optional :: exts(*), offs(*)
      integer :: dset_id
      integer :: mh5_open_dset
      interface
        subroutine mh5_get_dset_array_str (dset_id, buffer, exts, offs)
        implicit none
        integer :: dset_id
        character :: buffer(*)
        integer, optional :: exts(*), offs(*)
        end subroutine
      end interface
      dset_id = mh5_open_dset(lu, name)
      call mh5_get_dset_array_str(dset_id, buffer, exts, offs)
      call mh5_close_dset(dset_id)
      end subroutine

* convert Fortran string to null-terminated C string
      subroutine f2c_string(name, lbl)
      use iso_c_binding
      implicit none
      character(len=*) :: name, lbl
      if (len(name) .gt. len(lbl)-1) then
        call AbEnd
      end if
      lbl = TRIM(name)//C_NULL_CHAR
      end subroutine

      subroutine f2c_upcase(name, lbl)
      use iso_c_binding
      implicit none
      character(len=*) :: name, lbl
      call f2c_string(name, lbl)
      call upcase(lbl)
      end subroutine
