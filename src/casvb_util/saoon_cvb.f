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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine saoon_cvb(c1,c2,n2,s,n,metr)
c  Put matrix product S*C1 in C2:
      implicit real*8 (a-h,o-z)
      dimension c1(n,n2),c2(n,n2),s(*)

      if(metr.eq.0)then
        call fmove_cvb(c1,c2,n*n2)
      elseif(metr.eq.1)then
        call mxatb_cvb(s,c1,n,n,n2,c2)
      elseif(metr.eq.2)then
        call fzero(c2,n*n2)
        do 100 j=1,n2
        ik=0
        do 200 k=1,n
        do 300 i=1,k-1
        ik=ik+1
        c2(i,j)=c2(i,j)+s(ik)*c1(k,j)
        c2(k,j)=c2(k,j)+s(ik)*c1(i,j)
300     continue
        ik=ik+1
        c2(k,j)=c2(k,j)+s(ik)*c1(k,j)
200     continue
100     continue
      endif
      return
      end
c  *******************************
c  ** Vector sorting and search **
c  *******************************
