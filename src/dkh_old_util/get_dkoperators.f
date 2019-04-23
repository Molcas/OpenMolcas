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
      subroutine get_dkoperators(i,string,carray)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*(*) string
      character carray(*)
*
      call get_dkoperators_internal(carray)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine get_dkoperators_internal(carray)
      use iso_c_binding
      character, target :: carray(*)
      integer, pointer :: iarray(:)
      call c_f_pointer(c_loc(carray(1)),iarray,[1])
      call get_dkoperators_i(i,string,iarray)
      nullify(iarray)
      return
      end subroutine get_dkoperators_internal
*
      end
