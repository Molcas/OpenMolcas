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
*               2015, Ignacio Fdez. Galvan                             *
*               2015,2017, Liviu Ungur                                 *
*               2015, Yingjin Ma                                       *
************************************************************************
* stdalloc: wraps the standard allocate/deallocate fortran intrinsics.
*
* When memory is allocated, it is registered to getmem, and when it is
* deallocated it is excluded. On allocation, a check is made that the
* available memory as reported by getmem is sufficient. Thus, all memory
* allocated here counts towards the maximum allowed by getmem.
*
* To add additional data types or dimensions, just follow the existing
* reference implementations of real/integer for 1D up to 3D arrays. Then
* add the subroutines to the list in the generic interface definitions in
* src/Include/stdalloc.fh
*
* Steven Vancoillie, December 2014
* Ignacio Fdez. Galvan, April 2015 (added label optional arg. and _lim variants)
* Liviu Ungur, May 2015 (added support for  COMPLEX*16, 1D-3D)
* Liviu Ungur, May 2017 (added support for  COMPLEX*16, 4D)

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

      subroutine mma_maxbytes(mma_avail)
        implicit none
        integer, intent(out) :: mma_avail
        integer, external :: mma_avmem
        mma_avail = mma_avmem()
      end subroutine

* type-specific pointer-to-offset routines
      integer function d_cptr2loff(buffer)
        implicit none
#include "molcastypes.fh"
        real*8 :: buffer(*)
        integer, external :: cptr2woff
        integer, external :: kind2goff
        d_cptr2loff = cptr2woff('REAL',buffer(1)) + kind2goff('REAL')
      end function
      integer function i_cptr2loff(buffer)
        implicit none
#include "molcastypes.fh"
        integer :: buffer(*)
        integer, external :: cptr2woff
        integer, external :: kind2goff
        i_cptr2loff = cptr2woff('INTE',buffer(1)) + kind2goff('INTE')
      end function
      integer function c_cptr2loff(buffer)
        implicit none
#include "molcastypes.fh"
        character(*) :: buffer(*)
        integer, external :: cptr2woff
        integer, external :: kind2goff
        c_cptr2loff = cptr2woff('CHAR',buffer(1)) + kind2goff('CHAR')
      end function
      integer function dc_cptr2loff(buffer)
        implicit none
#include "molcastypes.fh"
        complex*16 :: buffer(*)
        integer, external :: cptr2woff
        integer, external :: kind2goff
        dc_cptr2loff = cptr2woff('REAL',buffer(1)) + kind2goff('REAL')
      end function

* type-specific allocation subroutines
      subroutine dmma_allo_1D(buffer,n1,label)
        implicit none
        real*8, allocatable :: buffer(:)
        integer :: n1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = rtob * n1
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1))
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1)
            else
              call getmem('dmma_1D','RGST','REAL',loffset,n1)
            end if
          end if
        end if
      end subroutine
      subroutine dmma_allo_1D_lim(buffer,l1,label)
        implicit none
        real*8, allocatable :: buffer(:)
        integer, dimension(2) :: l1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        bufsize = rtob * n1
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2)))
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(l1(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1)
            else
              call getmem('dmma_1D','RGST','REAL',loffset,n1)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_1D(buffer,n1,label)
        implicit none
        integer, allocatable :: buffer(:)
        integer :: n1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = itob * n1
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1))
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(1))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1)
            else
              call getmem('imma_1D','RGST','INTE',loffset,n1)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_1D_lim(buffer,l1,label)
        implicit none
        integer, allocatable :: buffer(:)
        integer, dimension(2) :: l1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        bufsize = itob * n1
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2)))
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(l1(1)))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1)
            else
              call getmem('imma_1D','RGST','INTE',loffset,n1)
            end if
          end if
        end if
      end subroutine
      subroutine cmma_allo_1D(buffer,n1,label)
        implicit none
        character(*), allocatable :: buffer(:)
        integer :: n1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = n1 * len(buffer)
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1))
          if (bufsize.gt.0) then
            loffset = cptr2loff(buffer(1))
            if (present(label)) then
              call getmem(label,'RGST','CHAR',loffset,bufsize)
            else
              call getmem('cmma_1D','RGST','CHAR',loffset,bufsize)
            end if
          end if
        end if
      end subroutine
      subroutine cmma_allo_1D_lim(buffer,l1,label)
        implicit none
        character(*), allocatable :: buffer(:)
        integer, dimension(2) :: l1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        bufsize = n1 * len(buffer)
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2)))
          if (bufsize.gt.0) then
            loffset = cptr2loff(buffer(l1(1)))
            if (present(label)) then
              call getmem(label,'RGST','CHAR',loffset,bufsize)
            else
              call getmem('cmma_1D','RGST','CHAR',loffset,bufsize)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_1D(buffer,n1,label)
        implicit none
        complex*16, allocatable :: buffer(:)
        integer :: n1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = CtoR * RtoB * n1
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1))
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset, CtoR*n1)
            else
              call getmem('DCmma_1D','RGST','REAL',loffset, CtoR*n1)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_1D_lim(buffer,l1,label)
        implicit none
        complex*16, allocatable :: buffer(:)
        integer, dimension(2) :: l1
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        bufsize = CtoR * rtob * n1
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2)))
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(l1(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset, CtoR*n1)
            else
              call getmem('DCmma_1D','RGST','REAL',loffset, CtoR*n1)
            end if
          end if
        end if
      end subroutine
      subroutine dmma_allo_2D(buffer,n1,n2,label)
        implicit none
        real*8, allocatable :: buffer(:,:)
        integer :: n1, n2
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = rtob * n1 * n2
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2))
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(1,1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1*n2)
            else
              call getmem('dmma_2D','RGST','REAL',loffset,n1*n2)
            end if
          end if
        end if
      end subroutine
      subroutine dmma_allo_2D_lim(buffer,l1,l2,label)
        implicit none
        real*8, allocatable :: buffer(:,:)
        integer, dimension(2) :: l1, l2
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        bufsize = rtob * n1 * n2
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2)))
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1*n2)
            else
              call getmem('dmma_2D','RGST','REAL',loffset,n1*n2)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_2D(buffer,n1,n2,label)
        implicit none
        integer, allocatable :: buffer(:,:)
        integer :: n1, n2
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = itob * n1 * n2
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2))
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(1,1))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1*n2)
            else
              call getmem('imma_2D','RGST','INTE',loffset,n1*n2)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_2D_lim(buffer,l1,l2,label)
        implicit none
        integer, allocatable :: buffer(:,:)
        integer, dimension(2) :: l1, l2
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        bufsize = itob * n1 * n2
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2)))
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1)))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1*n2)
            else
              call getmem('imma_2D','RGST','INTE',loffset,n1*n2)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_2D(buffer,n1,n2,label)
        implicit none
        complex*16, allocatable :: buffer(:,:)
        integer :: n1, n2
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = CtoR * rtob * n1 * n2
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2))
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(1,1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,CtoR*n1*n2)
            else
              call getmem('DCmma_2D','RGST','REAL',loffset,
     &                                                CtoR*n1*n2)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_2D_lim(buffer,l1,l2,label)
        implicit none
        complex*16, allocatable :: buffer(:,:)
        integer, dimension(2) :: l1, l2
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        bufsize = CtoR * rtob * n1 * n2
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2)))
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,CtoR*n1*n2)
            else
              call getmem('DCmma_2D','RGST','REAL',loffset,
     &                                                CtoR*n1*n2)
            end if
          end if
        end if
      end subroutine
      subroutine dmma_allo_3D(buffer,n1,n2,n3,label)
        implicit none
        real*8, allocatable :: buffer(:,:,:)
        integer :: n1, n2, n3
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = rtob * n1 * n2 * n3
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2,n3))
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(1,1,1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1*n2*n3)
            else
              call getmem('dmma_3D','RGST','REAL',loffset,n1*n2*n3)
            end if
          end if
        end if
      end subroutine
      subroutine dmma_allo_3D_lim(buffer,l1,l2,l3,label)
        implicit none
        real*8, allocatable :: buffer(:,:,:)
        integer, dimension(2) :: l1, l2, l3
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        n3 = l3(2)-l3(1)+1
        bufsize = rtob * n1 * n2 * n3
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2),l3(1):l3(2)))
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1),l3(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1*n2*n3)
            else
              call getmem('dmma_3D','RGST','REAL',loffset,n1*n2*n3)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_3D(buffer,n1,n2,n3,label)
        implicit none
        integer, allocatable :: buffer(:,:,:)
        integer :: n1, n2, n3
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = itob * n1 * n2 * n3
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2,n3))
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(1,1,1))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1*n2*n3)
            else
              call getmem('imma_3D','RGST','INTE',loffset,n1*n2*n3)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_3D_lim(buffer,l1,l2,l3,label)
        implicit none
        integer, allocatable :: buffer(:,:,:)
        integer, dimension(2) :: l1, l2, l3
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        n3 = l3(2)-l3(1)+1
        bufsize = itob * n1 * n2 * n3
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2),l3(1):l3(2)))
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1),l3(1)))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1*n2*n3)
            else
              call getmem('imma_3D','RGST','INTE',loffset,n1*n2*n3)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_3D(buffer,n1,n2,n3,label)
        implicit none
        complex*16, allocatable :: buffer(:,:,:)
        integer :: n1, n2, n3
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = CtoR * rtob * n1 * n2 * n3
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2,n3))
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(1,1,1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3)
            else
              call getmem('DCmma_3D','RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_3D_lim(buffer,l1,l2,l3,label)
        implicit none
        complex*16, allocatable :: buffer(:,:,:)
        integer, dimension(2) :: l1, l2, l3
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        n3 = l3(2)-l3(1)+1
        bufsize = CtoR * rtob * n1 * n2 * n3
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2),l3(1):l3(2)))
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1),l3(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3)
            else
              call getmem('DCmma_3D','RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_4D(buffer,n1,n2,n3,n4,label)
        implicit none
        complex*16, allocatable :: buffer(:,:,:,:)
        integer :: n1, n2, n3, n4
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = CtoR * rtob * n1 * n2 * n3 * n4
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2,n3,n4))
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(1,1,1,1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3*n4)
            else
              call getmem('DCmma_4D','RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3*n4)
            end if
          end if
        end if
      end subroutine
      subroutine DCmma_allo_4D_lim(buffer,l1,l2,l3,l4,label)
        implicit none
        complex*16, allocatable :: buffer(:,:,:,:)
        integer, dimension(2) :: l1, l2, l3, l4
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3, n4
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        n1 = l1(2)-l1(1)+1
        n2 = l2(2)-l2(1)+1
        n3 = l3(2)-l3(1)+1
        n4 = l4(2)-l4(1)+1
        bufsize = CtoR * rtob * n1 * n2 * n3 * n4
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(l1(1):l1(2),l2(1):l2(2),
     &                    l3(1):l3(2),l4(1):l4(2)))
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(l1(1),l2(1),l3(1),l4(1)))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3*n4)
            else
              call getmem('DCmma_4D','RGST','REAL',loffset,
     &                                        CtoR*n1*n2*n3*n4)
            end if
          end if
        end if
      end subroutine
* Please notice that the 4D part is not really used in current version -
*  need check the correctness -- Yingjin 1/2
      subroutine dmma_allo_4D(buffer,n1,n2,n3,n4,label)
        implicit none
        real*8, allocatable :: buffer(:,:,:,:)
        integer :: n1, n2, n3, n4
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = rtob * n1 * n2 * n3 * n4
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2,n3,n4))
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(1,1,1,1))
            if (present(label)) then
              call getmem(label,'RGST','REAL',loffset,n1*n2*n3*n4)
            else
              call getmem('dmma_4D','RGST','REAL',loffset,n1*n2*n3*n4)
            end if
          end if
        end if
      end subroutine
      subroutine imma_allo_4D(buffer,n1,n2,n3,n4,label)
        implicit none
        integer, allocatable :: buffer(:,:,:,:)
        integer :: n1, n2, n3, n4
        character (len=*), optional :: label
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: bufsize
        integer :: loffset
        integer :: mma_avail
        if (allocated(buffer)) then
          call mma_double_allo
        end if
        call mma_maxbytes(mma_avail)
        bufsize = itob * n1 * n2 * n3 * n4
        if (bufsize .gt. mma_avail) then
          call mma_oom(bufsize,mma_avail)
        else
          allocate(buffer(n1,n2,n3,n4))
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(1,1,1,1))
            if (present(label)) then
              call getmem(label,'RGST','INTE',loffset,n1*n2*n3*n4)
            else
              call getmem('imma_4D','RGST','INTE',loffset,n1*n2*n3*n4)
            end if
          end if
        end if
      end subroutine

* type-specific deallocation subroutines
      subroutine dmma_free_1D(buffer)
        implicit none
        real*8, allocatable :: buffer(:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: loffset
        n1 = size(buffer)
        if (allocated(buffer)) then
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1)))
            call getmem('dmma_1D','EXCL','REAL',loffset,n1)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine imma_free_1D(buffer)
        implicit none
        integer, allocatable :: buffer(:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: loffset
        n1 = size(buffer)
        if (allocated(buffer)) then
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1)))
            call getmem('imma_1D','EXCL','INTE',loffset,n1)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine cmma_free_1D(buffer)
        implicit none
        character(*), allocatable :: buffer(:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: bufsize
        integer :: loffset
        n1 = size(buffer)
        bufsize = n1 * len(buffer)
        if (allocated(buffer)) then
          if (bufsize.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1)))
            call getmem('cmma_1D','EXCL','CHAR',loffset,bufsize)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine DCmma_free_1D(buffer)
        implicit none
        complex*16, allocatable :: buffer(:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1
        integer :: loffset
        n1 = size(buffer)
        if (allocated(buffer)) then
          if (n1.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1)))
            call getmem('DCmma_1D','EXCL','REAL',loffset,CtoR*n1)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine dmma_free_2D(buffer)
        implicit none
        real*8, allocatable :: buffer(:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        if (allocated(buffer)) then
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2)))
            call getmem('dmma_2D','EXCL','REAL',loffset,n1*n2)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine imma_free_2D(buffer)
        implicit none
        integer, allocatable :: buffer(:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        if (allocated(buffer)) then
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2)))
            call getmem('imma_2D','EXCL','INTE',loffset,n1*n2)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine DCmma_free_2D(buffer)
        implicit none
        complex*16, allocatable :: buffer(:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        if (allocated(buffer)) then
          if (n1*n2.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2)))
            call getmem('DCmma_2D','EXCL','REAL',loffset,
     &                                                CtoR*n1*n2)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine dmma_free_3D(buffer)
        implicit none
        real*8, allocatable :: buffer(:,:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        n3 = size(buffer,3)
        if (allocated(buffer)) then
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2),
     &                                 lbound(buffer,3)))
            call getmem('dmma_3D','EXCL','REAL',loffset,n1*n2*n3)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine imma_free_3D(buffer)
        implicit none
        integer, allocatable :: buffer(:,:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        n3 = size(buffer,3)
        if (allocated(buffer)) then
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2),
     &                                 lbound(buffer,3)))
            call getmem('imma_3D','EXCL','INTE',loffset,n1*n2*n3)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine DCmma_free_3D(buffer)
        implicit none
        complex*16, allocatable :: buffer(:,:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        n3 = size(buffer,3)
        if (allocated(buffer)) then
          if (n1*n2*n3.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2),
     &                                 lbound(buffer,3)))
            call getmem('DCmma_3D','EXCL','REAL',loffset,
     &                                         CtoR*n1*n2*n3)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine DCmma_free_4D(buffer)
        implicit none
        complex*16, allocatable :: buffer(:,:,:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3, n4
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        n3 = size(buffer,3)
        n4 = size(buffer,4)
        if (allocated(buffer)) then
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2),
     &                                 lbound(buffer,3),
     &                                 lbound(buffer,4) ))
            call getmem('DCmma_3D','EXCL','REAL',loffset,
     &                                         CtoR*n1*n2*n3*n4)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
* Please notice that the 4D part is not really used in current version -
*  need check the correctness -- Yingjin 2/2
      subroutine dmma_free_4D(buffer)
        implicit none
        real*8, allocatable :: buffer(:,:,:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3, n4
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        n3 = size(buffer,3)
        n4 = size(buffer,4)
        if (allocated(buffer)) then
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2),
     &                                 lbound(buffer,3),
     &                                 lbound(buffer,4)))
            call getmem('dmma_4D','EXCL','REAL',loffset,n1*n2*n3*n4)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
      subroutine imma_free_4D(buffer)
        implicit none
        integer, allocatable :: buffer(:,:,:,:)
#include "SysDef.fh"
#include "molcastypes.fh"
#include "cptr2loff.fh"
        integer :: n1, n2, n3, n4
        integer :: loffset
        n1 = size(buffer,1)
        n2 = size(buffer,2)
        n3 = size(buffer,3)
        n4 = size(buffer,4)
        if (allocated(buffer)) then
          if (n1*n2*n3*n4.gt.0) then
            loffset = cptr2loff(buffer(lbound(buffer,1),
     &                                 lbound(buffer,2),
     &                                 lbound(buffer,3),
     &                                 lbound(buffer,4)))
            call getmem('imma_4D','EXCL','INTE',loffset,n1*n2*n3*n4)
          end if
          deallocate(buffer)
        else
          call mma_double_free
        end if
      end subroutine
