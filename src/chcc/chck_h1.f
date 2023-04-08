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
        subroutine Chck_H1 (H1,dim,add)
!
!        this routine test H1(i,a") = t1o(a,i)
!
        implicit none
#include "chcc1.fh"
        integer dim,add
         real*8 H1(1:no,1:dim)
!
!        help var
        integer a,i,bad,ntot
        real*8 s
!
        bad=0
        ntot=0
!
        do a=1,dim
        do i=1,no
          s=T1c(a+add,i)
           if (abs(H1(i,a)-s).gt.1.0d-10) then
            bad=bad+1
!            H1(i,a)=s
          end if
          ntot=ntot+1
        end do
        end do
!
        write (6,*) ' H1 test ',bad,ntot
!
        return
        end
