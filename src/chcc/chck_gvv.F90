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
        subroutine Chck_Gvv (Gvv)
!
!        check Gvv(be,a)
!
        implicit none
#include "chcc1.fh"
        real*8 Gvv(1:nv,1:nv)
!
!        help var
        integer i,a,b,be,bad
!       integer j
        real*8 s
!
        bad=0
!
        do be=1,nv
        do a=1,nv
!
!          s=0.0d0
!          do i=1,no
!          do j=1,no
!          do b=1,nv
!          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(be,b,i,j)
!          end do
!          end do
!          end do
!          s=-s
!
          s=Hvvc(be,a)
!
          do i=1,no
          do b=1,nv
          s=s+(2.0d0*Q3(a,be,b,i)-Q3(b,be,a,i))*T1c(b,i)
          end do
          end do
!
          Gvvc(be,a)=s
!
          if (abs(Gvv(be,a)-s).gt.1.0d-10) then
          bad=bad+1
          end if
!
        end do
        end do
!
        write (6,*) ' Gvv Chck :',bad
!
        return
        end
