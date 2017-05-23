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
      Subroutine crprod(a,b,c)
c
c       calculates cross product:   a x b = c
c                                   -   -   -
      Implicit real*8 (a-h,o-z)
      Dimension a(3),b(3),c(3)
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      Return
      End


      Subroutine vec(threshold,u,j,k,iErr)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension r(3),u(3)
#include "g_zmatconv.fh"

      iErr = 0
      r2 = 0.0d0

      Do i = 1, 3
        r(i) = Coords(j,i) - Coords(k,i)
        r2 = r2 + r(i)*r(i)
      EndDo
      r2 = sqrt(r2)
      If (r2.LT.threshold) then
        iErr = 1
        Return
      EndIf
      Do i = 1, 3
       u(i) = r(i)/r2
      EndDo
      Return
      End
