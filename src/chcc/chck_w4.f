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
        subroutine Chck_W4
     c  (W4,dima,dimbe,dimb,dimga,adda,addbe,addb,addga)
c    c  (W4,dima,dimga,dimb,dimbe,adda,addga,addb,addbe)
c
c        this routine test W4
c
        implicit none
#include "chcc1.fh"
        integer dima,dimbe,dimb,dimga,adda,addbe,addb,addga
         real*8 W4(1:dima,1:dimbe,1:dimb,1:dimga)
c        real*8 W4(1:dima,1:dimga,1:dimb,1:dimbe)
c
c        help var
        integer a,b,be,ga,i,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do ga=1,dimga
        do b=1,dimb
        do be=1,dimbe
        do a=1,dima
          s=Q4(a+adda,be+addbe,b+addb,ga+addga)
          s=0.0d0
           do i=1,no
             s=s-Q3(a+adda,be+addbe,b+addb,i)*T1c(ga+addga,i)
             s=s-Q3(b+addb,ga+addga,a+adda,i)*T1c(be+addbe,i)
           end do
           if (abs(W4(a,be,b,ga)-s).gt.1.0d-10) then
c          if (abs(W4(a,ga,b,be)-s).gt.1.0d-10) then
            bad=bad+1
c            W4(a,be,b,ga)=s
          end if
          ntot=ntot+1
        end do
        end do
        end do
        end do
c
        write (6,*) ' W4 test ',bad,ntot
c
        return
        end
