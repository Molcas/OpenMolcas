!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine Chck_W3                                              &
     &  (W3,dima,dimbe,dimb,adda,addbe,addb)
!
!        this routine test W3 (a,be|b,i)
!
        implicit none
#include "chcc1.fh"
        integer dima,dimbe,dimb,adda,addbe,addb
         real*8 W3(1:dima,1:dimbe,1:dimb,1:no)
!
!        help var
        integer a,b,be,i,bad,ntot
        real*8 s
!
        bad=0
        ntot=0
!
        do i=1,no
        do b=1,dimb
        do be=1,dimbe
        do a=1,dima
          s=Q3(a+adda,be+addbe,b+addb,i)
           if (abs(W3(a,be,b,i)-s).gt.1.0d-10) then
            bad=bad+1
!            W3(a,be,b,i)=s
          end if
          ntot=ntot+1
        end do
        end do
        end do
        end do
!
        write (6,*) ' W3 test ',bad,ntot
!
        return
        end
