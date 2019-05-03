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
      subroutine copy_dkoperators(i,carray1,j,carray2)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character carray1(*),carray2(*)
      call copy_dkoperators_internal(carray1,carray2)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine copy_dkoperators_internal(carray1,carray2)
      use iso_c_binding
      character, target :: carray1(*),carray2(*)
      integer, pointer :: iarray1(:),iarray2(:)
      call c_f_pointer(c_loc(carray1(1)),iarray1,[1])
      call c_f_pointer(c_loc(carray2(1)),iarray2,[1])
      call copy_dkoperators_i(i,iarray1,j,iarray2)
      nullify(iarray1,iarray2)
      end subroutine copy_dkoperators_internal
*
      end
*
      subroutine copy_dkoperators_ic(i,iarray1,j,carray2)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension iarray1(*)
      character carray2(*)
      call copy_dkoperators_ic_internal(carray2)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine copy_dkoperators_ic_internal(carray2)
      use iso_c_binding
      character, target :: carray2(*)
      integer, pointer :: iarray2(:)
      call c_f_pointer(c_loc(carray2(1)),iarray2,[1])
      call copy_dkoperators_i(i,iarray1,j,iarray2)
      nullify(iarray2)
      end subroutine copy_dkoperators_ic_internal
*
      end
*
      subroutine copy_dkoperators_ci(i,carray1,j,iarray2)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character carray1(*)
      dimension iarray2(*)
      call copy_dkoperators_ci_internal(carray1)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine copy_dkoperators_ci_internal(carray1)
      use iso_c_binding
      character, target :: carray1(*)
      integer, pointer :: iarray1(:)
      call c_f_pointer(c_loc(carray1(1)),iarray1,[1])
      call copy_dkoperators_i(i,iarray1,j,iarray2)
      nullify(iarray1)
      end subroutine copy_dkoperators_ci_internal
*
      end
