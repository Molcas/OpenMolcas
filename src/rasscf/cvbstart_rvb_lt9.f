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
      subroutine cvbstart_rvb_lt9(icode)
      use casvb_global, only: endvar, iprm, iunset, lstprm, recinp,
     &                        recinp_old, variat
      implicit real*8 (a-h,o-z)

      variat=(mod(icode,10).ne.0)
      endvar=(mod(icode,10).eq.2)
      recinp=0d0
      recinp_old=0d0
      lstprm(:) = iunset
      iprm = 0
      if(.not.variat)call casinfo1_rvb()
      return
      end
