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
      subroutine hrecur(pn,dpn,pn1,x,nn)
c
      implicit real*8 (a-h,o-z)
c
      p1 = 1.0d0
      p  = x
      dp1 = 0.0d0
      dp  = 1.0d0
      do 20 j=2,nn
        fj = (j)
        fj2 = 0.5d0*(fj-1.0d0)
        q  = x*p - fj2*p1
        dq = x*dp + p - fj2*dp1
        p1 = p
        p  = q
        dp1 = dp
        dp  = dq
   20 continue
      pn  = p
      dpn = dp
      pn1 = p1
      return
      end
