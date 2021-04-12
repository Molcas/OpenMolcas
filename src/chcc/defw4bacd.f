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
        subroutine DefW4bacd (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (ba|cd)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:abLen,1:cdLen)
c
c        help variables
        integer a,b,c,d,ba,cd
c
        if (cSGrp.eq.dSGrp) then
c        case (b,a|c=d)
          do c=2,dimc
          cd=c*(c-1)/2
          do d=1,c-1
          cd=cd+1
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ba,cd)
              W(a,b,d,c)=W(a,b,d,c)+Wx(ba,cd)
            end do
            end do
          end do
          end do
c
          do c=1,dimc
          cd=c*(c+1)/2
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,c)=W(a,b,c,c)+Wx(ba,cd)
            end do
            end do
          end do
c
        else
c        case (b,a|c,d)
          cd=0
          do d=1,dimd
          do c=1,dimc
          cd=cd+1
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ba,cd)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(aSGrp)
        call Unused_integer(bSGrp)
      end if
        end
