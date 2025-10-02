!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2014-2016, Steven Vancoillie                           *
!               2015,2020, Ignacio Fdez. Galvan                        *
!               2015, Liviu Ungur                                      *
!***********************************************************************
! stdalloc: wraps the standard allocate/deallocate Fortran intrinsics.
!
! When memory is allocated, it is registered to getmem, and when it is
! deallocated it is excluded. On allocation, a check is made that the
! available memory as reported by getmem is sufficient. Thus, all memory
! allocated here counts towards the maximum allowed by getmem.
!
! To add additional data types or dimensions, just include the file
! "mma_allo_template.fh" with appropriate values for the macros:
!   _SUBR_NAME_: The base name for the subroutines
!   _DATA_NAME_: The data type string for getmem: 'REAL' for real(kind=wp),
!                'INTE' for integer(kind=iwp), undefined for anything else
!                ['CHAR' would only be appropriate for character(len=1)]
!   _DEF_LABEL_: The default label
!   _TYPE_: The data type
!   _DIMENSIONS_: The number of dimensions
! See the end of this file.
!
! Once defined, use:
!
! * call mma_allocate(array,n1,...,nk)
!     - array:     An allocatable array of k dimensions, or an allocatable string
!     - n1,...,nk: Size of each dimension (variable number of arguments)
!                  They can be all integers, or all arrays of two integers (lower and upper bound)
!     Optional arguments:
!     - label:     Name with which to identify the array, in case of problems
!     - safe:      If present, nothing will happen if array is already allocated (otherwise it's in error)
!
! * call mma_deallocate(array)
!     Optional arguments:
!     - safe:      If present, nothing will happen if array is already deallocated (otherwise it's in error)
!
! Steven Vancoillie, December 2014
! Ignacio Fdez. Galvan, April 2015 (added label optional arg. and _lim variants)
! Liviu Ungur, May 2015 (added support for COMPLEX*16)
! Ignacio Fdez. Galvan, November 2020 (rewrote with preprocessor templates)

!#define _ENABLE_POINTERS_
module stdalloc

use, intrinsic :: iso_fortran_env, only: int32
use Definitions, only: wp, iwp, byte, u6, ItoB, RtoB

implicit none
private

integer(kind=iwp) :: MxMem

interface mma_allocate
  ! 0D allocate
  module procedure :: cmma_allo_0D
  ! 1D allocate
  module procedure :: bmma_allo_1D, bmma_allo_1D_lim, cmma_allo_1D, cmma_allo_1D_lim, dmma_allo_1D, dmma_allo_1D_lim, &
                      i4mma_allo_1D, i4mma_allo_1D_lim, imma_allo_1D, imma_allo_1D_lim, lmma_allo_1D, lmma_allo_1D_lim, &
                      zmma_allo_1D, zmma_allo_1D_lim
# ifdef _ENABLE_POINTERS_
  module procedure :: dpmma_allo_1D, dpmma_allo_1D_lim, ipmma_allo_1D, ipmma_allo_1D_lim
# endif
  ! 2D allocate
  module procedure :: bmma_allo_2D, bmma_allo_2D_lim, cmma_allo_2D, cmma_allo_2D_lim, dmma_allo_2D, dmma_allo_2D_lim, &
                      imma_allo_2D, imma_allo_2D_lim, lmma_allo_2D, lmma_allo_2D_lim, zmma_allo_2D, zmma_allo_2D_lim
  ! 3D allocate
  module procedure :: dmma_allo_3D, dmma_allo_3D_lim, imma_allo_3D, imma_allo_3D_lim, zmma_allo_3D, zmma_allo_3D_lim
  ! 4D allocate
  module procedure :: dmma_allo_4D, dmma_allo_4D_lim, imma_allo_4D, imma_allo_4D_lim, zmma_allo_4D, zmma_allo_4D_lim
  ! 5D allocate
  module procedure :: dmma_allo_5D, dmma_allo_5D_lim, imma_allo_5D, imma_allo_5D_lim, zmma_allo_5D, zmma_allo_5D_lim
  ! 6D allocate
  module procedure :: dmma_allo_6D, dmma_allo_6D_lim
  ! 7D allocate
  module procedure :: dmma_allo_7D, dmma_allo_7D_lim
end interface

interface mma_deallocate
  ! 0D deallocate
  module procedure :: cmma_free_0D
  ! 1D deallocate
  module procedure :: bmma_free_1D, cmma_free_1D, dmma_free_1D, i4mma_free_1D, imma_free_1D, lmma_free_1D, zmma_free_1D
# ifdef _ENABLE_POINTERS_
  module procedure :: dpmma_free_1D, ipmma_free_1D
# endif
  ! 2D deallocate
  module procedure :: bmma_free_2D, cmma_free_2D, dmma_free_2D, imma_free_2D, lmma_free_2D, zmma_free_2D
  ! 3D deallocate
  module procedure :: dmma_free_3D, imma_free_3D, zmma_free_3D
  ! 4D deallocate
  module procedure :: dmma_free_4D, imma_free_4D, zmma_free_4D
  ! 5D deallocate
  module procedure :: dmma_free_5D, imma_free_5D, zmma_free_5D
  ! 6D deallocate
  module procedure :: dmma_free_6D
  ! 7D deallocate
  module procedure :: dmma_free_7D
end interface

public :: mma_allocate, mma_deallocate, mma_double_allo, mma_double_free, mma_maxBYTES, mma_maxDBLE, mma_maxINT, mma_oom, mxMem

contains

#include "warnings.h"

! out-of-memory handling
subroutine mma_oom(label,bufsize,mma_avail)

  character(len=*), intent(in) :: label
  integer(kind=iwp), intent(in) :: bufsize, mma_avail

  write(u6,'(1x,a)') '?mma_allo_?D: error: out-of-memory'
  write(u6,'(1x,a,a)') 'label: ',trim(label)
  write(u6,'(1x,a,1x,i12)') ' available (kB):',nint(mma_avail*1.0e-3_wp)
  write(u6,'(1x,a,1x,i12)') ' required  (kB):',nint(bufsize*1.0e-3_wp)
  call quit(_RC_MEMORY_ERROR_)

end subroutine mma_oom

! double allocation/deallocation handling
subroutine mma_double_allo(label)

  character(len=*), intent(in) :: label

  write(u6,'(1x,a)') '?mma_allo_?D: error: double allocate'
  write(u6,'(1x,a,a)') 'label: ',trim(label)
  call quit(_RC_MEMORY_ERROR_)

end subroutine mma_double_allo

subroutine mma_double_free(label)

  character(len=*), intent(in) :: label

  write(u6,'(1x,a)') '?mma_free_?D: error: double deallocate'
  write(u6,'(1x,a,a)') 'label: ',trim(label)
  call quit(_RC_MEMORY_ERROR_)

end subroutine mma_double_free

subroutine mma_maxDBLE(mma_avail)

  integer(kind=iwp), intent(out) :: mma_avail
  integer(kind=iwp), external :: mma_avmem

  mma_avail = mma_avmem()/RtoB

end subroutine mma_maxDBLE

subroutine mma_maxINT(mma_avail)

  integer(kind=iwp), intent(out) :: mma_avail
  integer(kind=iwp), external :: mma_avmem

  mma_avail = mma_avmem()/ItoB

end subroutine mma_maxINT

subroutine mma_maxBYTES(mma_avail)

  integer(kind=iwp), intent(out) :: mma_avail
  integer(kind=iwp), external :: mma_avmem

  mma_avail = mma_avmem()

end subroutine mma_maxBYTES

#define _IN_STDALLOC_MOD_

! type-specific allocation subroutines
! each #include defines NAME_allo_xD, NAME_allo_xD_lim, and NAME_free_xD

! real(kind=wp) variants

#define _SUBR_NAME_ dmma
#define _TYPE_ real(kind=wp)
#define _DATA_NAME_ 'REAL'

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'dmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'dmma_2D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 3
#  define _DEF_LABEL_ 'dmma_3D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 4
#  define _DEF_LABEL_ 'dmma_4D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 5
#  define _DEF_LABEL_ 'dmma_5D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 6
#  define _DEF_LABEL_ 'dmma_6D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 7
#  define _DEF_LABEL_ 'dmma_7D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_
#undef _DATA_NAME_

! complex(kind=wp) variants
! (note that there is no specific _DATA_NAME_ for these)

#define _SUBR_NAME_ zmma
#define _TYPE_ complex(kind=wp)

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'zmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'zmma_2D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 3
#  define _DEF_LABEL_ 'zmma_3D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 4
#  define _DEF_LABEL_ 'zmma_4D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 5
#  define _DEF_LABEL_ 'zmma_5D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_

! integer(kind=iwp) variants

#define _SUBR_NAME_ imma
#define _TYPE_ integer(kind=iwp)
#define _DATA_NAME_ 'INTE'

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'imma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'imma_2D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 3
#  define _DEF_LABEL_ 'imma_3D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 4
#  define _DEF_LABEL_ 'imma_4D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 5
#  define _DEF_LABEL_ 'imma_5D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_
#undef _DATA_NAME_

#ifdef _I8_
! integer(kind=int32) variants

#define _SUBR_NAME_ i4mma
#define _TYPE_ integer(kind=int32)
#define _DATA_NAME_ 'INTE'

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'i4mma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_
#undef _DATA_NAME_

#endif

! byte variants

#define _SUBR_NAME_ bmma
#define _TYPE_ integer(kind=byte)

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'bmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'bmma_2D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_

! character variants
! (no _DATA_NAME_ defined to make sure size is counted in bytes)

#define _SUBR_NAME_ cmma

#define _TYPE_ character(len=:)

#  define _DIMENSIONS_ 0
#  define _DEF_LABEL_ 'cmma_0D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _TYPE_

#define _TYPE_ character(len=*)

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'cmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'cmma_2D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _TYPE_

#undef _SUBR_NAME_

! logical(kind=iwp) variants
! (note that there is no specific _DATA_NAME_ for these)

#define _SUBR_NAME_ lmma
#define _TYPE_ logical(kind=iwp)

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'lmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'lmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_

#ifdef _ENABLE_POINTERS_
! pointer variants

#define _IS_POINTER_

#  define _SUBR_NAME_ ipmma
#  define _TYPE_ integer(kind=iwp)
#  define _DATA_NAME_ 'INTE'

#    define _DIMENSIONS_ 1
#    define _DEF_LABEL_ 'ipmma_1D'
#    include "mma_allo_template.fh"
#    undef _DIMENSIONS_
#    undef _DEF_LABEL_

#  undef _SUBR_NAME_
#  undef _TYPE_
#  undef _DATA_NAME_

#  define _SUBR_NAME_ dpmma
#  define _TYPE_ real(kind=wp)
#  define _DATA_NAME_ 'REAL'

#    define _DIMENSIONS_ 1
#    define _DEF_LABEL_ 'dpmma_1D'
#    include "mma_allo_template.fh"
#    undef _DIMENSIONS_
#    undef _DEF_LABEL_

#  undef _SUBR_NAME_
#  undef _TYPE_
#  undef _DATA_NAME_

#undef _IS_POINTER_
#endif

end module stdalloc
