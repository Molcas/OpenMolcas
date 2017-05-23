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
C
C----------------------------------------------------------------------|
C
      subroutine dkh_woprig(n,ifodd,nw,np,wr,rw,p1,p2,q1,q2,t1,t2)
C
C Product of P(np)W(nw)=Q(np+nw)
C
      implicit none
      integer n,nw,np,i,j
      logical ifodd
      Real*8 wr(n,n),rw(n,n),p1(n,n),p2(n,n),q1(n,n),q2(n,n)
      Real*8 t1(n,n),t2(n,n)
      if(ifodd)then
        call dmxma(n,'N','N',p1,rw,t1,1.d0)
        call dmxma(n,'N','N',p2,wr,t2,1.d0)
      else
        call dmxma(n,'N','N',p1,wr,t1,1.d0)
        call dmxma(n,'N','N',p2,rw,t2,1.d0)
      end if
      do i=1,n
        do j=1,n
          q1(j,i) = t1(j,i)
          q2(j,i) = t2(j,i)
        end do
      end do
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nw)
        call Unused_integer(np)
      end if
      end
