!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module mh5
! Generic interfaces for HDF5

#include "intent.fh"

use, intrinsic :: iso_c_binding, only: c_char, c_null_char
use Definitions, only: wp, iwp, u6, MOLCAS_C_INT, MOLCAS_C_REAL

implicit none
private

public :: mh5_create_file, mh5_open_file_rw, mh5_open_file_r, mh5_close_file, mh5_is_hdf5, mh5_open_group, mh5_close_group, &
          mh5_exists_attr, mh5_open_attr, mh5_close_attr, mh5_exists_dset, mh5_open_dset, mh5_close_dset, mh5_create_attr_int, &
          mh5_create_attr_real, mh5_create_attr_str, mh5_put_attr, mh5_get_attr, mh5_init_attr, mh5_fetch_attr, &
          mh5_create_dset_int, mh5_create_dset_real, mh5_create_dset_str, mh5_put_dset, mh5_get_dset, mh5_init_dset, &
          mh5_fetch_dset, mh5_resize_dset, mh5_get_dset_dims

!======================
! Overloaded interfaces
!======================

! attr

interface mh5_create_attr_int
  module procedure :: mh5_create_attr_scalar_int, &
                      mh5_create_attr_array_int
end interface mh5_create_attr_int

interface mh5_create_attr_real
  module procedure :: mh5_create_attr_scalar_real, &
                      mh5_create_attr_array_real
end interface mh5_create_attr_real

interface mh5_create_attr_str
  module procedure :: mh5_create_attr_scalar_str, &
                      mh5_create_attr_array_str
end interface mh5_create_attr_str

interface mh5_put_attr
  module procedure :: mh5_put_attr_scalar_int, &
                      mh5_put_attr_scalar_real, &
                      mh5_put_attr_scalar_str, &
                      mh5_put_attr_array_int, &
                      mh5_put_attr_array_real, &
                      mh5_put_attr_array_str
end interface mh5_put_attr

interface mh5_get_attr
  module procedure :: mh5_get_attr_scalar_int, &
                      mh5_get_attr_scalar_real, &
                      mh5_get_attr_scalar_str, &
                      mh5_get_attr_array_int, &
                      mh5_get_attr_array_real, &
                      mh5_get_attr_array_str
end interface mh5_get_attr

! init = create + put + close
interface mh5_init_attr
  module procedure :: mh5_init_attr_scalar_int, &
                      mh5_init_attr_scalar_real, &
                      mh5_init_attr_scalar_str, &
                      mh5_init_attr_array_int, &
                      mh5_init_attr_array_real, &
                      mh5_init_attr_array_str
end interface mh5_init_attr

! fetch = open + get + close
interface mh5_fetch_attr
  module procedure :: mh5_fetch_attr_scalar_int, &
                      mh5_fetch_attr_scalar_real, &
                      mh5_fetch_attr_scalar_str, &
                      mh5_fetch_attr_array_int, &
                      mh5_fetch_attr_array_real, &
                      mh5_fetch_attr_array_str
end interface mh5_fetch_attr

! dset

interface mh5_create_dset_int
  module procedure :: mh5_create_dset_scalar_int, &
                      mh5_create_dset_array_int
end interface mh5_create_dset_int

interface mh5_create_dset_real
  module procedure :: mh5_create_dset_scalar_real, &
                      mh5_create_dset_array_real
end interface mh5_create_dset_real

interface mh5_create_dset_str
  module procedure :: mh5_create_dset_scalar_str, &
                      mh5_create_dset_array_str
end interface mh5_create_dset_str

interface mh5_put_dset
  module procedure :: mh5_put_dset_scalar_int, &
                      mh5_put_dset_scalar_real, &
                      mh5_put_dset_scalar_str, &
                      mh5_put_dset_array_int, &
                      mh5_put_dset_array_int_2d, &
                      mh5_put_dset_array_real, &
                      mh5_put_dset_array_real_2d, &
                      mh5_put_dset_array_real_3d, &
                      mh5_put_dset_array_str
end interface mh5_put_dset

interface mh5_get_dset
  module procedure :: mh5_get_dset_scalar_int, &
                      mh5_get_dset_scalar_real, &
                      mh5_get_dset_scalar_str, &
                      mh5_get_dset_array_int, &
                      mh5_get_dset_array_real, &
                      mh5_get_dset_array_str
end interface mh5_get_dset

! init = create + put + close
interface mh5_init_dset
  module procedure :: mh5_init_dset_scalar_int, &
                      mh5_init_dset_scalar_real, &
                      mh5_init_dset_scalar_str, &
                      mh5_init_dset_array_int, &
                      mh5_init_dset_array_real, &
                      mh5_init_dset_array_str
end interface mh5_init_dset

! fetch = open + get + close
interface mh5_fetch_dset
  module procedure :: mh5_fetch_dset_scalar_int, &
                      mh5_fetch_dset_scalar_real, &
                      mh5_fetch_dset_scalar_str, &
                      mh5_fetch_dset_array_int, &
                      mh5_fetch_dset_array_int_2d, &
                      mh5_fetch_dset_array_real, &
                      mh5_fetch_dset_array_real_2d, &
                      mh5_fetch_dset_array_real_3d, &
                      mh5_fetch_dset_array_str
end interface mh5_fetch_dset

interface mh5_resize_dset
  module procedure :: mh5_extend_dset_array
end interface mh5_resize_dset

interface mh5_get_dset_dims
  module procedure :: mh5_get_dset_array_dims
end interface mh5_get_dset_dims

!==================================
! Private interfaces to C functions
!==================================

interface

  ! file

  function mh5_c_create_file(filename) result(lu) &
           bind(C,name='mh5c_create_file')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: lu
    character(kind=c_char) :: filename(*)
  end function mh5_c_create_file

  function mh5_c_open_rw(filename) result(lu) &
           bind(C,name='mh5c_open_file_rw')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: lu
    character(kind=c_char) :: filename(*)
  end function mh5_c_open_rw

  function mh5_c_open_r(filename) result(lu) &
           bind(C,name='mh5c_open_file_r')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: lu
    character(kind=c_char) :: filename(*)
  end function mh5_c_open_r

  function mh5_c_close_file(lu) result(rc) &
           bind(C,name='mh5c_close_file')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: lu
  end function mh5_c_close_file

  function mh5_c_is_hdf5(filename) result(rc) &
           bind(C,name='mh5c_is_hdf5')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    character(kind=c_char) :: filename(*)
  end function mh5_c_is_hdf5

  ! group

  function mh5_c_open_group(lu,groupname) result(groupid) &
           bind(C,name='mh5c_open_group')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: groupid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: groupname(*)
  end function mh5_c_open_group

  function mh5_c_close_group(id) result(rc) &
           bind(C,name='mh5c_close_group')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: id
  end function mh5_c_close_group

  ! attr

  function mh5_c_exists_attr(id,attrname) result(rc) &
           bind(C,name='mh5c_exists_attr')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: id
    character(kind=c_char) :: attrname(*)
  end function mh5_c_exists_attr

  function mh5_c_open_attr(lu,attrname) result(attrid) &
           bind(C,name='mh5c_open_attr')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
  end function mh5_c_open_attr

  function mh5_c_close_attr(attrid) result(rc) &
           bind(C,name='mh5c_close_attr')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
  end function mh5_c_close_attr

  ! dset

  function mh5_c_exists_dset(id,dsetname) result(rc) &
           bind(C,name='mh5c_exists_dset')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: id
    character(len=c_char) :: dsetname(*)
  end function mh5_c_exists_dset

  function mh5_c_open_dset(lu,dsetname) result(dsetid) &
           bind(C,name='mh5c_open_dset')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
  end function mh5_c_open_dset

  function mh5_c_close_dset(dsetid) result(rc) &
           bind(C,name='mh5c_close_dset')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
  end function mh5_c_close_dset

  ! attr types
  ! (create|put|get)_attr_(scalar|array)_(int|real|str)

  function mh5_c_create_attr_scalar_int(lu,attrname) result(attrid) &
           bind(C,name='mh5c_create_attr_scalar_int')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
  end function mh5_c_create_attr_scalar_int

  function mh5_c_create_attr_scalar_real(lu,attrname) result(attrid) &
           bind(C,name='mh5c_create_attr_scalar_real')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
  end function mh5_c_create_attr_scalar_real

  function mh5_c_create_attr_scalar_str(lu,attrname,length) result(attrid) &
           bind(C,name='mh5c_create_attr_scalar_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
    integer(kind=MOLCAS_C_INT), value :: length
  end function mh5_c_create_attr_scalar_str

  function mh5_c_put_attr_scalar_int(attrid,val) result(rc) &
           bind(C,name='mh5c_put_attr_scalar_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    integer(kind=MOLCAS_C_INT) :: val
  end function mh5_c_put_attr_scalar_int

  function mh5_c_put_attr_scalar_real(attrid,val) result(rc) &
           bind(C,name='mh5c_put_attr_scalar_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    real(kind=MOLCAS_C_REAL) :: val
  end function mh5_c_put_attr_scalar_real

  function mh5_c_put_attr_scalar_str(attrid,val) result(rc) &
           bind(C,name='mh5c_put_attr_scalar_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    character(kind=c_char) :: val(*)
  end function mh5_c_put_attr_scalar_str

  function mh5_c_get_attr_scalar_int(attrid,val) result(rc) &
           bind(C,name='mh5c_get_attr_scalar_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    integer(kind=MOLCAS_C_INT) :: val
  end function mh5_c_get_attr_scalar_int

  function mh5_c_get_attr_scalar_real(attrid,val) result(rc) &
           bind(C,name='mh5c_get_attr_scalar_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    real(kind=MOLCAS_C_REAL) :: val
  end function mh5_c_get_attr_scalar_real

  function mh5_c_get_attr_scalar_str(attrid,val) result(rc) &
           bind(C,name='mh5c_get_attr_scalar_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    character(kind=c_char) :: val(*)
  end function mh5_c_get_attr_scalar_str

  function mh5_c_create_attr_array_int(lu,attrname,rank,dims) result(attrid) &
           bind(C,name='mh5c_create_attr_array_int')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_create_attr_array_int

  function mh5_c_create_attr_array_real(lu,attrname,rank,dims) result(attrid) &
           bind(C,name='mh5c_create_attr_array_real')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_create_attr_array_real

  function mh5_c_create_attr_array_str(lu,attrname,rank,dims,length) result(attrid) &
           bind(C,name='mh5c_create_attr_array_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: attrid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: attrname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
    integer(kind=MOLCAS_C_INT), value :: length
  end function mh5_c_create_attr_array_str

  function mh5_c_put_attr_array_int(attrid,buffer) result(rc) &
           bind(C,name='mh5c_put_attr_array_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    integer(kind=MOLCAS_C_INT) :: buffer(*)
  end function mh5_c_put_attr_array_int

  function mh5_c_put_attr_array_real(attrid,buffer) result(rc) &
           bind(C,name='mh5c_put_attr_array_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    real(kind=MOLCAS_C_REAL) :: buffer(*)
  end function mh5_c_put_attr_array_real

  function mh5_c_put_attr_array_str(attrid,buffer) result(rc) &
           bind(C,name='mh5c_put_attr_array_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    character(kind=c_char) :: buffer(*)
  end function mh5_c_put_attr_array_str

  function mh5_c_get_attr_array_int(attrid,buffer) result(rc) &
           bind(C,name='mh5c_get_attr_array_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    integer(kind=MOLCAS_C_INT) :: buffer(*)
  end function mh5_c_get_attr_array_int

  function mh5_c_get_attr_array_real(attrid,buffer) result(rc) &
           bind(C,name='mh5c_get_attr_array_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    real(kind=MOLCAS_C_REAL) :: buffer(*)
  end function mh5_c_get_attr_array_real

  function mh5_c_get_attr_array_str(attrid,buffer) result(rc) &
           bind(C,name='mh5c_get_attr_array_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: attrid
    character(kind=c_char) :: buffer(*)
  end function mh5_c_get_attr_array_str

  ! dset types
  ! (create|put|get)_dset_(scalar|array[_dyn])_(int|real|str)[_full]

  function mh5_c_create_dset_scalar_int(lu,dsetname) result(dsetid) &
           bind(C,name='mh5c_create_dset_scalar_int')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
  end function mh5_c_create_dset_scalar_int

  function mh5_c_create_dset_scalar_real(lu,dsetname) result(dsetid) &
           bind(C,name='mh5c_create_dset_scalar_real')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT) :: dsetid
  end function mh5_c_create_dset_scalar_real

  function mh5_c_create_dset_scalar_str(lu,dsetname,length) result(dsetid) &
           bind(C,name='mh5c_create_dset_scalar_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: length
  end function mh5_c_create_dset_scalar_str

  function mh5_c_put_dset_scalar_int(dsetid,val) result(rc) &
           bind(C,name='mh5c_put_dset_scalar_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: val
  end function mh5_c_put_dset_scalar_int

  function mh5_c_put_dset_scalar_real(dsetid,val) result(rc) &
           bind(C,name='mh5c_put_dset_scalar_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    real(kind=MOLCAS_C_REAL) :: val
  end function mh5_c_put_dset_scalar_real

  function mh5_c_put_dset_scalar_str(dsetid,val) result(rc) &
           bind(C,name='mh5c_put_dset_scalar_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    character(kind=c_char) :: val(*)
  end function mh5_c_put_dset_scalar_str

  function mh5_c_get_dset_scalar_int(dsetid,val) result(rc) &
           bind(C,name='mh5c_get_dset_scalar_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: val
  end function mh5_c_get_dset_scalar_int

  function mh5_c_get_dset_scalar_real(dsetid,val) result(rc) &
           bind(C,name='mh5c_get_dset_scalar_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    real(kind=MOLCAS_C_REAL) :: val
  end function mh5_c_get_dset_scalar_real

  function mh5_c_get_dset_scalar_str(dsetid,val) result(rc) &
           bind(C,name='mh5c_get_dset_scalar_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    character(kind=c_char) :: val(*)
  end function mh5_c_get_dset_scalar_str

  function mh5_c_create_dset_array_int(lu,dsetname,rank,dims) result(dsetid) &
           bind(C,name='mh5c_create_dset_array_int')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_create_dset_array_int

  function mh5_c_create_dset_array_dyn_int(lu,dsetname,rank,dims) result(dsetid) &
           bind(C,name='mh5c_create_dset_array_dyn_int')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_create_dset_array_dyn_int

  function mh5_c_create_dset_array_real(lu,dsetname,rank,dims) result(dsetid) &
           bind(C,name='mh5c_create_dset_array_real')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_create_dset_array_real

  function mh5_c_create_dset_array_dyn_real(lu,dsetname,rank,dims) result(dsetid) &
           bind(C,name='mh5c_create_dset_array_dyn_real')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_create_dset_array_dyn_real

  function mh5_c_create_dset_array_str(lu,dsetname,rank,dims,length) result(dsetid) &
           bind(C,name='mh5c_create_dset_array_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
    integer(kind=MOLCAS_C_INT), value :: length
  end function mh5_c_create_dset_array_str

  function mh5_c_create_dset_array_dyn_str(lu,dsetname,rank,dims,length) result(dsetid) &
           bind(C,name='mh5c_create_dset_array_dyn_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: dsetid
    integer(kind=MOLCAS_C_INT), value :: lu
    character(kind=c_char) :: dsetname(*)
    integer(kind=MOLCAS_C_INT), value :: rank
    integer(kind=MOLCAS_C_INT) :: dims(*)
    integer(kind=MOLCAS_C_INT), value :: length
  end function mh5_c_create_dset_array_dyn_str

  function mh5_c_put_dset_array_int(dsetid,exts,offs,buffer) result(rc) &
           bind(C,name='mh5c_put_dset_array_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: exts(*), offs(*)
    integer(kind=MOLCAS_C_INT) :: buffer(*)
  end function mh5_c_put_dset_array_int

  function mh5_c_put_dset_array_int_full(dsetid,buffer) result(rc) &
           bind(C,name='mh5c_put_dset_array_int_full')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: buffer(*)
  end function mh5_c_put_dset_array_int_full

  function mh5_c_put_dset_array_real(dsetid,exts,offs,buffer) result(rc) &
           bind(C,name='mh5c_put_dset_array_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: exts(*), offs(*)
    real(kind=MOLCAS_C_REAL) :: buffer(*)
  end function mh5_c_put_dset_array_real

  function mh5_c_put_dset_array_real_full(dsetid,buffer) result(rc) &
           bind(C,name='mh5c_put_dset_array_real_full')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    real(kind=MOLCAS_C_REAL) :: buffer(*)
  end function mh5_c_put_dset_array_real_full

  function mh5_c_put_dset_array_str(dsetid,exts,offs,buffer) result(rc) &
           bind(C,name='mh5c_put_dset_array_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: exts(*), offs(*)
    character(kind=c_char) :: buffer(*)
  end function mh5_c_put_dset_array_str

  function mh5_c_put_dset_array_str_full(dsetid,buffer) result(rc) &
           bind(C,name='mh5c_put_dset_array_str_full')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    character(kind=c_char) :: buffer(*)
  end function mh5_c_put_dset_array_str_full

  function mh5_c_get_dset_array_int(dsetid,exts,offs,buffer) result(rc) &
           bind(C,name='mh5c_get_dset_array_int')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: exts(*), offs(*)
    integer(kind=MOLCAS_C_INT) :: buffer(*)
  end function mh5_c_get_dset_array_int

  function mh5_c_get_dset_array_int_full(dsetid,buffer) result(rc) &
           bind(C,name='mh5c_get_dset_array_int_full')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: buffer(*)
  end function mh5_c_get_dset_array_int_full

  function mh5_c_get_dset_array_real(dsetid,exts,offs,buffer) result(rc) &
           bind(C,name='mh5c_get_dset_array_real')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: exts(*), offs(*)
    real(kind=MOLCAS_C_REAL) :: buffer(*)
  end function mh5_c_get_dset_array_real

  function mh5_c_get_dset_array_real_full(dsetid,buffer) result(rc) &
           bind(C,name='mh5c_get_dset_array_real_full')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    real(kind=MOLCAS_C_REAL) :: buffer(*)
  end function mh5_c_get_dset_array_real_full

  function mh5_c_get_dset_array_str(dsetid,exts,offs,buffer) result(rc) &
           bind(C,name='mh5c_get_dset_array_str')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: exts(*), offs(*)
    character(kind=c_char) :: buffer(*)
  end function mh5_c_get_dset_array_str

  function mh5_c_get_dset_array_str_full(dsetid,buffer) result(rc) &
           bind(C,name='mh5c_get_dset_array_str_full')
    import :: MOLCAS_C_INT, c_char
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    character(kind=c_char) :: buffer(*)
  end function mh5_c_get_dset_array_str_full

  ! dset extend array / get dims

  function mh5_c_extend_dset_array(dsetid,dims) result(rc) &
           bind(C,name='mh5c_extend_dset_array')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_extend_dset_array

  function mh5_c_get_dset_array_dims(dsetid,dims) result(rc) &
           bind(C,name='mh5c_get_dset_array_dims')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: rc
    integer(kind=MOLCAS_C_INT), value :: dsetid
    integer(kind=MOLCAS_C_INT) :: dims(*)
  end function mh5_c_get_dset_array_dims

end interface

contains

! convert Fortran string to null-terminated C string
subroutine f2c_string(string,c_string)
  character(len=*), intent(in) :: string
  character(len=*), intent(out) :: c_string
  if (len_trim(string) > len(c_string)-1) then
    write(u6,*) 'f2c_string: input string too long'
    call abend()
  end if
  c_string = trim(string)//c_null_char
end subroutine f2c_string

!==========================
! Specific Fortran wrappers
!==========================

#define MH5_MAX_LBL_LEN 256

! file

function mh5_create_file(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
  call prgmtranslate(filename,realname,lrealname)
  call f2c_string(realname,c_name)
  lu = mh5_c_create_file(c_name)
end function mh5_create_file

function mh5_open_file_rw(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
  call prgmtranslate(filename,realname,lrealname)
  call f2c_string(realname,c_name)
  lu = mh5_c_open_rw(c_name)
end function mh5_open_file_rw

function mh5_open_file_r(filename) result(lu)
  integer(kind=iwp) :: lu
  character(len=*) :: filename
  integer(kind=iwp) :: lrealname
  character(len=4096) :: realname, c_name
  call prgmtranslate(filename,realname,lrealname)
  call f2c_string(realname,c_name)
  lu = mh5_c_open_r(c_name)
end function mh5_open_file_r

subroutine mh5_close_file(lu)
  integer(kind=iwp), intent(in) :: lu
  if (mh5_c_close_file(lu) < 0) call abend()
end subroutine mh5_close_file

function mh5_is_hdf5(filename) result(ishdf5)
  logical(kind=iwp) :: ishdf5
  character(len=*), intent(in) :: filename
  integer(kind=iwp) :: rc, lrealname
  character(len=4096) :: realname
  logical(kind=iwp) :: exists
  call prgmtranslate(filename,realname,lrealname)
  call f_inquire(realname,exists)
  if (exists) then
    realname(lrealname+1:lrealname+1) = c_null_char
    rc = mh5_c_is_hdf5(realname)
  else
    rc = 0
  end if
  if (rc > 0) then
    ishdf5 = .true.
  else if (rc == 0) then
    ishdf5 = .false.
  else
    ishdf5 = .false.
    call abend()
  end if
end function mh5_is_hdf5

! group

function mh5_open_group(lu,groupname) result(groupid)
  integer(kind=iwp) :: groupid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: groupname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(groupname,mh5_lbl)
  groupid = mh5_c_open_group(lu,mh5_lbl)
end function mh5_open_group

subroutine mh5_close_group(id)
  integer(kind=iwp), intent(in) :: id
  if (mh5_c_close_group(id) < 0) call abend()
end subroutine mh5_close_group

! attr

function mh5_exists_attr(id,attrname) result(exists)
  logical(kind=iwp) :: exists
  integer(kind=iwp), intent(in) :: id
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  integer(kind=iwp) :: rc
  call f2c_string(attrname,mh5_lbl)
  rc = mh5_c_exists_attr(id,mh5_lbl)
  if (rc > 0) then
    exists = .true.
  else if (rc == 0) then
    exists = .false.
  else
    exists = .false.
    call abend()
  end if
end function mh5_exists_attr

function mh5_open_attr(lu,attrname) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_open_attr(lu,mh5_lbl)
end function mh5_open_attr

subroutine mh5_close_attr(attrid)
  integer(kind=iwp), intent(in) :: attrid
  if (mh5_c_close_attr(attrid) < 0) call abend()
end subroutine mh5_close_attr

! dset

function mh5_exists_dset(id,dsetname) result(exists)
  logical(kind=iwp) :: exists
  integer(kind=iwp), intent(in) :: id
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  integer(kind=iwp) :: rc
  call f2c_string(dsetname,mh5_lbl)
  rc = mh5_c_exists_dset(id,mh5_lbl)
  if (rc > 0) then
    exists = .true.
  else if (rc == 0) then
    exists = .false.
  else
    exists = .false.
    call abend()
  end if
end function mh5_exists_dset

function mh5_open_dset(lu,dsetname) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(dsetname,mh5_lbl)
  dsetid = mh5_c_open_dset(lu,mh5_lbl)
end function mh5_open_dset

subroutine mh5_close_dset(dsetid)
  integer(kind=iwp), intent(in) :: dsetid
  if (mh5_c_close_dset(dsetid) < 0) call abend()
end subroutine mh5_close_dset

! attr types

function mh5_create_attr_scalar_int(lu,attrname) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_create_attr_scalar_int(lu,mh5_lbl)
end function mh5_create_attr_scalar_int

function mh5_create_attr_scalar_real(lu,attrname) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_create_attr_scalar_real(lu,mh5_lbl)
end function mh5_create_attr_scalar_real

function mh5_create_attr_scalar_str(lu,attrname,length) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, length
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_create_attr_scalar_str(lu,mh5_lbl,length)
end function mh5_create_attr_scalar_str

subroutine mh5_put_attr_scalar_int(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(in) :: val
  if (mh5_c_put_attr_scalar_int(attrid,val) < 0) call abend()
end subroutine mh5_put_attr_scalar_int

subroutine mh5_put_attr_scalar_real(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(in) :: val
  if (mh5_c_put_attr_scalar_real(attrid,val) < 0) call abend()
end subroutine mh5_put_attr_scalar_real

subroutine mh5_put_attr_scalar_str(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  character(len=*), intent(in) :: val
  if (mh5_c_put_attr_scalar_str(attrid,val) < 0) call abend()
end subroutine mh5_put_attr_scalar_str

subroutine mh5_get_attr_scalar_int(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(out) :: val
  if (mh5_c_get_attr_scalar_int(attrid,val) < 0) call abend()
end subroutine mh5_get_attr_scalar_int

subroutine mh5_get_attr_scalar_real(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(out) :: val
  if (mh5_c_get_attr_scalar_real(attrid,val) < 0) call abend()
end subroutine mh5_get_attr_scalar_real

subroutine mh5_get_attr_scalar_str(attrid,val)
  integer(kind=iwp), intent(in) :: attrid
  character(len=*), intent(out) :: val
  if (mh5_c_get_attr_scalar_str(attrid,val) < 0) call abend()
end subroutine mh5_get_attr_scalar_str

function mh5_create_attr_array_int(lu,attrname,rank,dims) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_create_attr_array_int(lu,mh5_lbl,rank,dims)
end function mh5_create_attr_array_int

function mh5_create_attr_array_real(lu,attrname,rank,dims) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_create_attr_array_real(lu,mh5_lbl,rank,dims)
end function mh5_create_attr_array_real

function mh5_create_attr_array_str(lu,attrname,rank,dims,length) result(attrid)
  integer(kind=iwp) :: attrid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: attrname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(attrname,mh5_lbl)
  attrid = mh5_c_create_attr_array_str(lu,mh5_lbl,rank,dims,length)
end function mh5_create_attr_array_str

subroutine mh5_put_attr_array_int(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(in) :: buffer(*)
  if (mh5_c_put_attr_array_int(attrid,buffer) < 0) call abend()
end subroutine mh5_put_attr_array_int

subroutine mh5_put_attr_array_real(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(in) :: buffer(*)
  if (mh5_c_put_attr_array_real(attrid,buffer) < 0) call abend()
end subroutine mh5_put_attr_array_real

subroutine mh5_put_attr_array_str(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  character, intent(in) :: buffer(*)
  if (mh5_c_put_attr_array_str(attrid,buffer) < 0) call abend()
end subroutine mh5_put_attr_array_str

subroutine mh5_get_attr_array_int(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  integer(kind=iwp), intent(_OUT_) :: buffer(*)
  if (mh5_c_get_attr_array_int(attrid,buffer) < 0) call abend()
end subroutine mh5_get_attr_array_int

subroutine mh5_get_attr_array_real(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  real(kind=wp), intent(_OUT_) :: buffer(*)
  if (mh5_c_get_attr_array_real(attrid,buffer) < 0) call abend()
end subroutine mh5_get_attr_array_real

subroutine mh5_get_attr_array_str(attrid,buffer)
  integer(kind=iwp), intent(in) :: attrid
  character, intent(_OUT_) :: buffer(*)
  if (mh5_c_get_attr_array_str(attrid,buffer) < 0) call abend()
end subroutine mh5_get_attr_array_str

subroutine mh5_init_attr_scalar_int(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(in) :: val
  integer(kind=iwp) :: attrid
  attrid = mh5_create_attr_scalar_int(lu,attrname)
  call mh5_put_attr_scalar_int(attrid,val)
  call mh5_close_attr(attrid)
end subroutine mh5_init_attr_scalar_int

subroutine mh5_init_attr_scalar_real(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(in) :: val
  integer(kind=iwp) :: attrid
  attrid = mh5_create_attr_scalar_real(lu,attrname)
  call mh5_put_attr_scalar_real(attrid,val)
  call mh5_close_attr(attrid)
end subroutine mh5_init_attr_scalar_real

subroutine mh5_init_attr_scalar_str(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=*), intent(in) :: val
  integer(kind=iwp) :: attrid
  attrid = mh5_create_attr_scalar_str(lu,attrname,len(val))
  call mh5_put_attr_scalar_str(attrid,val)
  call mh5_close_attr(attrid)
end subroutine mh5_init_attr_scalar_str

subroutine mh5_init_attr_array_int(lu,attrname,rank,dims,buffer)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(in) :: buffer(*)
  integer(kind=iwp) :: attrid
  attrid = mh5_create_attr_array_int(lu,attrname,rank,dims)
  call mh5_put_attr_array_int(attrid,buffer)
  call mh5_close_attr(attrid)
end subroutine mh5_init_attr_array_int

subroutine mh5_init_attr_array_real(lu,attrname,rank,dims,buffer)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(in) :: buffer(*)
  integer(kind=iwp) :: attrid
  attrid = mh5_create_attr_array_real(lu,attrname,rank,dims)
  call mh5_put_attr_array_real(attrid,buffer)
  call mh5_close_attr(attrid)
end subroutine mh5_init_attr_array_real

subroutine mh5_init_attr_array_str(lu,attrname,rank,dims,buffer,length)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: attrname
  character, intent(in) :: buffer(*)
  integer(kind=iwp) :: attrid
  attrid = mh5_create_attr_array_str(lu,attrname,rank,dims,length)
  call mh5_put_attr_array_str(attrid,buffer)
  call mh5_close_attr(attrid)
end subroutine mh5_init_attr_array_str

subroutine mh5_fetch_attr_scalar_int(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(out) :: val
  integer(kind=iwp) :: attrid
  attrid = mh5_open_attr(lu,attrname)
  call mh5_get_attr_scalar_int(attrid,val)
  call mh5_close_attr(attrid)
end subroutine mh5_fetch_attr_scalar_int

subroutine mh5_fetch_attr_scalar_real(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(out) :: val
  integer(kind=iwp) :: attrid
  attrid = mh5_open_attr(lu,attrname)
  call mh5_get_attr_scalar_real(attrid,val)
  call mh5_close_attr(attrid)
end subroutine mh5_fetch_attr_scalar_real

subroutine mh5_fetch_attr_scalar_str(lu,attrname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character(len=*), intent(out) :: val
  integer(kind=iwp) :: attrid
  attrid = mh5_open_attr(lu,attrname)
  call mh5_get_attr_scalar_str(attrid,val)
  call mh5_close_attr(attrid)
end subroutine mh5_fetch_attr_scalar_str

subroutine mh5_fetch_attr_array_int(lu,attrname,buffer)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  integer(kind=iwp), intent(_OUT_) :: buffer(*)
  integer(kind=iwp) :: attrid
  attrid = mh5_open_attr(lu,attrname)
  call mh5_get_attr_array_int(attrid,buffer)
  call mh5_close_attr(attrid)
end subroutine mh5_fetch_attr_array_int

subroutine mh5_fetch_attr_array_real(lu,attrname,buffer)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  real(kind=wp), intent(_OUT_) :: buffer(*)
  integer(kind=iwp) :: attrid
  attrid = mh5_open_attr(lu,attrname)
  call mh5_get_attr_array_real(attrid,buffer)
  call mh5_close_attr(attrid)
end subroutine mh5_fetch_attr_array_real

subroutine mh5_fetch_attr_array_str(lu,attrname,buffer)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: attrname
  character, intent(_OUT_) :: buffer(*)
  integer(kind=iwp) :: attrid
  attrid = mh5_open_attr(lu,attrname)
  call mh5_get_attr_array_str(attrid,buffer)
  call mh5_close_attr(attrid)
end subroutine mh5_fetch_attr_array_str

! dset types

function mh5_create_dset_scalar_int(lu,dsetname) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(dsetname,mh5_lbl)
  dsetid = mh5_c_create_dset_scalar_int(lu,mh5_lbl)
end function mh5_create_dset_scalar_int

function mh5_create_dset_scalar_real(lu,dsetname) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(dsetname,mh5_lbl)
  dsetid = mh5_c_create_dset_scalar_real(lu,mh5_lbl)
end function mh5_create_dset_scalar_real

function mh5_create_dset_scalar_str(lu,dsetname,length) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, length
  character(len=*), intent(in) :: dsetname
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  call f2c_string(dsetname,mh5_lbl)
  dsetid = mh5_c_create_dset_scalar_str(lu,mh5_lbl,length)
end function mh5_create_dset_scalar_str

subroutine mh5_put_dset_scalar_int(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: val
  if (mh5_c_put_dset_scalar_int(dsetid,val) < 0) call abend()
end subroutine mh5_put_dset_scalar_int

subroutine mh5_put_dset_scalar_real(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(in) :: val
  if (mh5_c_put_dset_scalar_real(dsetid,val) < 0) call abend()
end subroutine mh5_put_dset_scalar_real

subroutine mh5_put_dset_scalar_str(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  character(len=*), intent(in) :: val
  if (mh5_c_put_dset_scalar_str(dsetid,val) < 0) call abend()
end subroutine mh5_put_dset_scalar_str

subroutine mh5_get_dset_scalar_int(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(out) :: val
  if (mh5_c_get_dset_scalar_int(dsetid,val) < 0) call abend()
end subroutine mh5_get_dset_scalar_int

subroutine mh5_get_dset_scalar_real(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(out) :: val
  if (mh5_c_get_dset_scalar_real(dsetid,val) < 0) call abend()
end subroutine mh5_get_dset_scalar_real

subroutine mh5_get_dset_scalar_str(dsetid,val)
  integer(kind=iwp), intent(in) :: dsetid
  character(len=*), intent(out) :: val
  if (mh5_c_get_dset_scalar_str(dsetid,val) < 0) call abend()
end subroutine mh5_get_dset_scalar_str

function mh5_create_dset_array_int(lu,dsetname,rank,dims,dyn) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  logical(kind=iwp), intent(in), optional :: dyn
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  logical(kind=iwp) :: isdyn
  isdyn = .false.
  if (present(dyn)) isdyn = dyn
  call f2c_string(dsetname,mh5_lbl)
  if (isdyn) then
    dsetid = mh5_c_create_dset_array_dyn_int(lu,mh5_lbl,rank,dims)
  else
    dsetid = mh5_c_create_dset_array_int(lu,mh5_lbl,rank,dims)
  end if
end function mh5_create_dset_array_int

function mh5_create_dset_array_real(lu,dsetname,rank,dims,dyn) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  logical(kind=iwp), intent(in), optional :: dyn
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  logical(kind=iwp) :: isdyn
  isdyn = .false.
  if (present(dyn)) isdyn = dyn
  call f2c_string(dsetname,mh5_lbl)
  if (isdyn) then
    dsetid = mh5_c_create_dset_array_dyn_real(lu,mh5_lbl,rank,dims)
  else
    dsetid = mh5_c_create_dset_array_real(lu,mh5_lbl,rank,dims)
  end if
end function mh5_create_dset_array_real

function mh5_create_dset_array_str(lu,dsetname,rank,dims,length,dyn) result(dsetid)
  integer(kind=iwp) :: dsetid
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: dsetname
  logical(kind=iwp), intent(in), optional :: dyn
  character(len=MH5_MAX_LBL_LEN) :: mh5_lbl
  logical(kind=iwp) :: isdyn
  isdyn = .false.
  if (present(dyn)) isdyn = dyn
  call f2c_string(dsetname,mh5_lbl)
  if (isdyn) then
    dsetid = mh5_c_create_dset_array_dyn_str(lu,mh5_lbl,rank,dims,length)
  else
    dsetid = mh5_c_create_dset_array_str(lu,mh5_lbl,rank,dims,length)
  end if
end function mh5_create_dset_array_str

subroutine mh5_put_dset_array_int(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_put_dset_array_int(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_put_dset_array_int_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_put_dset_array_int

subroutine mh5_put_dset_array_real(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(in) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_put_dset_array_real(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_put_dset_array_real_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_put_dset_array_real

subroutine mh5_put_dset_array_str(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  character, intent(in) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_put_dset_array_str(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_put_dset_array_str_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_put_dset_array_str

subroutine mh5_get_dset_array_int(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(_OUT_) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_get_dset_array_int(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_get_dset_array_int_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_get_dset_array_int

subroutine mh5_get_dset_array_real(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(_OUT_) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_get_dset_array_real(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_get_dset_array_real_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_get_dset_array_real

subroutine mh5_get_dset_array_str(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  character, intent(_OUT_) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_get_dset_array_str(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_get_dset_array_str_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_get_dset_array_str

subroutine mh5_extend_dset_array(dsetid,dims)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: dims(*)
  if (mh5_c_extend_dset_array(dsetid,dims) < 0) call abend()
end subroutine mh5_extend_dset_array

subroutine mh5_get_dset_array_dims(dsetid,dims)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(_OUT_) :: dims(*)
  if (mh5_c_get_dset_array_dims(dsetid,dims) < 0) call abend()
end subroutine mh5_get_dset_array_dims

subroutine mh5_init_dset_scalar_int(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(in) :: val
  integer(kind=iwp) :: dsetid
  dsetid = mh5_create_dset_scalar_int(lu,dsetname)
  call mh5_put_dset_scalar_int(dsetid,val)
  call mh5_close_dset(dsetid)
end subroutine mh5_init_dset_scalar_int

subroutine mh5_init_dset_scalar_real(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(in) :: val
  integer(kind=iwp) :: dsetid
  dsetid = mh5_create_dset_scalar_real(lu,dsetname)
  call mh5_put_dset_scalar_real(dsetid,val)
  call mh5_close_dset(dsetid)
end subroutine mh5_init_dset_scalar_real

subroutine mh5_init_dset_scalar_str(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=*), intent(in) :: val
  integer(kind=iwp) :: dsetid
  dsetid = mh5_create_dset_scalar_str(lu,dsetname,len(val))
  call mh5_put_dset_scalar_str(dsetid,val)
  call mh5_close_dset(dsetid)
end subroutine mh5_init_dset_scalar_str

subroutine mh5_init_dset_array_int(lu,dsetname,rank,dims,buffer,dyn)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(in) :: buffer(*)
  logical(kind=iwp), optional :: dyn
  integer(kind=iwp) :: dsetid
  logical(kind=iwp) :: isdyn
  isdyn = .false.
  if (present(dyn)) isdyn = dyn
  dsetid = mh5_create_dset_array_int(lu,dsetname,rank,dims,isdyn)
  if (mh5_c_put_dset_array_int_full(dsetid,buffer) < 0) call abend()
  call mh5_close_dset(dsetid)
end subroutine mh5_init_dset_array_int

subroutine mh5_init_dset_array_real(lu,dsetname,rank,dims,buffer,dyn)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*)
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(in) :: buffer(*)
  logical(kind=iwp), optional :: dyn
  integer(kind=iwp) :: dsetid
  logical(kind=iwp) :: isdyn
  isdyn = .false.
  if (present(dyn)) isdyn = dyn
  dsetid = mh5_create_dset_array_real(lu,dsetname,rank,dims,isdyn)
  if (mh5_c_put_dset_array_real_full(dsetid,buffer) < 0) call abend()
  call mh5_close_dset(dsetid)
end subroutine mh5_init_dset_array_real

subroutine mh5_init_dset_array_str(lu,dsetname,rank,dims,buffer,length,dyn)
  integer(kind=iwp), intent(in) :: lu, rank, dims(*), length
  character(len=*), intent(in) :: dsetname
  character, intent(in) :: buffer(*)
  logical(kind=iwp), optional :: dyn
  integer(kind=iwp) :: dsetid
  logical(kind=iwp) :: isdyn
  isdyn = .false.
  if (present(dyn)) isdyn = dyn
  dsetid = mh5_create_dset_array_str(lu,dsetname,rank,dims,length,isdyn)
  if (mh5_c_put_dset_array_str_full(dsetid,buffer) < 0) call abend()
  call mh5_close_dset(dsetid)
end subroutine mh5_init_dset_array_str

subroutine mh5_fetch_dset_scalar_int(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(out) :: val
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  call mh5_get_dset_scalar_int(dsetid,val)
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_scalar_int

subroutine mh5_fetch_dset_scalar_real(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(out) :: val
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  call mh5_get_dset_scalar_real(dsetid,val)
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_scalar_real

subroutine mh5_fetch_dset_scalar_str(lu,dsetname,val)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character(len=*), intent(out) :: val
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  call mh5_get_dset_scalar_str(dsetid,val)
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_scalar_str

subroutine mh5_fetch_dset_array_int(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(_OUT_) :: buffer(*)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  if (present(exts) .and. present(offs)) then
    call mh5_get_dset_array_int(dsetid,buffer,exts,offs)
  else if (present(exts) .or. present(offs)) then
    call abend()
  else
    call mh5_get_dset_array_int(dsetid,buffer)
  end if
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_int

subroutine mh5_fetch_dset_array_real(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(_OUT_) :: buffer(*)
  integer(kind=iwp), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  if (present(exts) .and. present(offs)) then
    call mh5_get_dset_array_real(dsetid,buffer,exts,offs)
  else if (present(exts) .or. present(offs)) then
    call abend()
  else
    call mh5_get_dset_array_real(dsetid,buffer)
  end if
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_real

subroutine mh5_fetch_dset_array_str(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  character, intent(_OUT_) :: buffer(*)
  integer(kind=iwp), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  if (present(exts) .and. present(offs)) then
    call mh5_get_dset_array_str(dsetid,buffer,exts,offs)
  else if (present(exts) .or. present(offs)) then
    call abend()
  else
    call mh5_get_dset_array_str(dsetid,buffer)
  end if
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_str

! Workarounds for the fact that assumed-size arguments only match one-dimensional arrays
! in overloaded interfaces (would need support for assumed rank)

subroutine mh5_put_dset_array_int_2d(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  integer(kind=iwp), intent(in) :: buffer(:,:)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_put_dset_array_int(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_put_dset_array_int_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_put_dset_array_int_2d

subroutine mh5_put_dset_array_real_2d(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(in) :: buffer(:,:)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_put_dset_array_real(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_put_dset_array_real_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_put_dset_array_real_2d

subroutine mh5_put_dset_array_real_3d(dsetid,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: dsetid
  real(kind=wp), intent(in) :: buffer(:,:,:)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: rc
  if (present(exts) .and. present(offs)) then
    rc = mh5_c_put_dset_array_real(dsetid,exts,offs,buffer)
  else if (present(exts) .or. present(offs)) then
    rc = -1
  else
    rc = mh5_c_put_dset_array_real_full(dsetid,buffer)
  end if
  if (rc < 0) call abend()
end subroutine mh5_put_dset_array_real_3d

subroutine mh5_fetch_dset_array_int_2d(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  integer(kind=iwp), intent(_OUT_) :: buffer(:,:)
  integer(kind=iwp), intent(in), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  if (present(exts) .and. present(offs)) then
    call mh5_get_dset_array_int(dsetid,buffer,exts,offs)
  else if (present(exts) .or. present(offs)) then
    call abend()
  else
    call mh5_get_dset_array_int(dsetid,buffer)
  end if
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_int_2d

subroutine mh5_fetch_dset_array_real_2d(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(out) :: buffer(:,:)
  integer(kind=iwp), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  if (present(exts) .and. present(offs)) then
    call mh5_get_dset_array_real(dsetid,buffer,exts,offs)
  else if (present(exts) .or. present(offs)) then
    call abend()
  else
    call mh5_get_dset_array_real(dsetid,buffer)
  end if
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_real_2d

subroutine mh5_fetch_dset_array_real_3d(lu,dsetname,buffer,exts,offs)
  integer(kind=iwp), intent(in) :: lu
  character(len=*), intent(in) :: dsetname
  real(kind=wp), intent(out) :: buffer(:,:,:)
  integer(kind=iwp), optional :: exts(*), offs(*)
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(lu,dsetname)
  if (present(exts) .and. present(offs)) then
    call mh5_get_dset_array_real(dsetid,buffer,exts,offs)
  else if (present(exts) .or. present(offs)) then
    call abend()
  else
    call mh5_get_dset_array_real(dsetid,buffer)
  end if
  call mh5_close_dset(dsetid)
end subroutine mh5_fetch_dset_array_real_3d

end module mh5
