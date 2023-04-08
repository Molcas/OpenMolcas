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
        subroutine Chck_Hvo (Hvo)
!
!        check Hvo
!
        implicit none
#include "chcc1.fh"
        real*8 Hvo(1:nv,1:no)
!
!        help var
        integer i,a,j,b,bad,tot
        real*8 s
!
        bad=0
        tot=0
!
        do i=1,no
        do a=1,nv
!
          s=0.0d0
!
          do j=1,no
          do b=1,nv
          s=s+(2.0d0*Q21(b,j,a,i)-Q21(b,i,a,j))*T1c(b,j)
          end do
          end do
!
          Hvoc(a,i)=s
!
          if (abs(Hvo(a,i)-s).gt.1.0d-10) then
          bad=bad+1
          end if
          tot=tot+1
!
        end do
        end do
!
        write (6,*) ' Hvo Chck :',bad
!
        return
        end
