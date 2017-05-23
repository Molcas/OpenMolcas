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
      subroutine hroot(x,nn,dpn,pn1,eps)
c                  improves the approximate root  x
c                in addition we also obtain
c                    dpn = derivative of h(n) at x
c                    pn1 = value of h(n-1) at x
c
      implicit real*8 (a-h,o-z)
c
c     # iter = 5 sufficient for 8-byte accuracy up to nn = 7
      do 14 iter=1,10
        call hrecur(p,dp,pn1,x,nn)
        d  = p/dp
        x  = x - d
        if( abs(d).le.eps ) go to 16
   14 continue
   16 dpn = dp
      return
      end
