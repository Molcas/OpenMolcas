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
      subroutine wzero(n,b,idum)
      implicit real*8 (a-h,o-z)
      dimension b(1)
      do 10 i=1,n
10    b(i)=0
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(idum)
      end
