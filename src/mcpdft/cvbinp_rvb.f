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
      subroutine cvbinp_rvb_m(icode,luinp)
      implicit real*8 (a-h,o-z)

      call cvbstart_cvb_ge9(icode+10)
      call hello_cvb()
      call parse_init_cvb(luinp)
      call input_cvb()
      call cvbfinish_cvb(icode+10)
      return
      end
