************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014-2016, Steven Vancoillie                           *
*               2015,2020, Ignacio Fdez. Galvan                        *
*               2015, Liviu Ungur                                      *
************************************************************************
* stdalloc: wraps the standard allocate/deallocate fortran intrinsics.
*
* When memory is allocated, it is registered to getmem, and when it is
* deallocated it is excluded. On allocation, a check is made that the
* available memory as reported by getmem is sufficient. Thus, all memory
* allocated here counts towards the maximum allowed by getmem.
*
* To add additional data types or dimensions, just include the file
* "mma_allo_template.fh" with appropriate values for the macros:
*   _SUBR_NAME_: The base name for the subroutines
*   _DATA_NAME_: The data type string for getmem: 'REAL' for real*8,
*                'INTE' for integer, undefined for anything else
*                ['CHAR' would only be appropriate for character(len=1)]
*   _DEF_LABEL_: The default label
*   _TYPE_: The data type
*   _DIMENSIONS_: The number of dimensions
* Then add the subroutines to the list in the generic interface definitions
* in src/Include/stdalloc.fh
*
* Steven Vancoillie, December 2014
* Ignacio Fdez. Galvan, April 2015 (added label optional arg. and _lim variants)
* Liviu Ungur, May 2015 (added support for COMPLEX*16)
* Ignacio Fdez. Galvan, November 2020 (rewrote with preprocessor templates)

#include "molcastypes.fh"

* out-of-memory handling
      subroutine mma_oom(bufsize,mma_avail)
        implicit none
#include "warnings.fh"
        integer :: bufsize, mma_avail
        write(6,'(1x,a)') '?mma_allo_?D: error: out-of-memory'
        write(6,'(1x,a,i12)') ' available (kB): ',
     &    nint(mma_avail * 1.0d-3)
        write(6,'(1x,a,i12)') ' required  (kB):  ',
     &    nint(bufsize  * 1.0d-3)
        call quit(_RC_MEMORY_ERROR_)
      end subroutine

* double allocation/deallocation handling
      subroutine mma_double_allo
        implicit none
#include "warnings.fh"
        write(6,'(1x,a)') '?mma_allo_?D: error: double allocate'
        call quit(_RC_MEMORY_ERROR_)
      end subroutine
      subroutine mma_double_free
        implicit none
#include "warnings.fh"
        write(6,'(1x,a)') '?mma_free_?D: error: double deallocate'
        call quit(_RC_MEMORY_ERROR_)
      end subroutine

      subroutine mma_maxDBLE(mma_avail)
        implicit none
#include "SysDef.fh"
        integer, intent(out) :: mma_avail
        integer, external :: mma_avmem
        mma_avail = mma_avmem()/RtoB
      end subroutine

      subroutine mma_maxINT(mma_avail)
        implicit none
#include "SysDef.fh"
        integer, intent(out) :: mma_avail
        integer, external :: mma_avmem
        mma_avail = mma_avmem()/ItoB
      end subroutine

      subroutine mma_maxbytes(mma_avail)
        implicit none
        integer, intent(out) :: mma_avail
        integer, external :: mma_avmem
        mma_avail = mma_avmem()
      end subroutine

* type-specific pointer-to-offset routines

#define _FUNC_NAME_ d_cptr2loff
#define _TYPE_ real*8
#define _DATA_NAME_ 'REAL'
#include "cptr2loff_template.fh"
#undef _FUNC_NAME_
#undef _TYPE_
#undef _DATA_NAME_

#define _FUNC_NAME_ z_cptr2loff
#define _TYPE_ complex*16
#include "cptr2loff_template.fh"
#undef _FUNC_NAME_
#undef _TYPE_

#define _FUNC_NAME_ i_cptr2loff
#define _TYPE_ integer
#define _DATA_NAME_ 'INTE'
#include "cptr2loff_template.fh"
#undef _FUNC_NAME_
#undef _TYPE_
#undef _DATA_NAME_

! _WITH_LEN_ enables a workaround for older gfortran
#define _FUNC_NAME_ c_cptr2loff
#define _TYPE_ character(len=*)
#define _DATA_NAME_ 'CHAR'
#define _WITH_LEN_
#include "cptr2loff_template.fh"
#undef _FUNC_NAME_
#undef _TYPE_
#undef _DATA_NAME_
#undef _WITH_LEN_

* type-specific allocation subroutines
* each #include defines NAME_allo_xD, NAME_allo_xD_lim, and NAME_free_xD

* real*8 variants

#define _SUBR_NAME_ dmma
#define _TYPE_ real*8
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

#  define _DIMENSIONS_ 7
#  define _DEF_LABEL_ 'dmma_7D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_
#undef _DATA_NAME_

* complex*16 variants
* (note that there is no specific _DATA_NAME_ for these)

#define _SUBR_NAME_ zmma
#define _TYPE_ complex*16

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

* integer variants

#define _SUBR_NAME_ imma
#define _TYPE_ integer
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

#undef _SUBR_NAME_
#undef _TYPE_
#undef _DATA_NAME_

* character variants
* (no _DATA_NAME_ defined to make sure size is counted in bytes)

#define _SUBR_NAME_ cmma
#define _TYPE_ character(len=*)

#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'cmma_1D'
#  include "mma_allo_template.fh"
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_

#undef _SUBR_NAME_
#undef _TYPE_
