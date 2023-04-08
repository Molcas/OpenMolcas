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
        subroutine MkE_Y3 (Va,V,dima,dimb,no)
!
!       this routine do:
!       Va(a',i,b',j) = 2 V(a',j,b',i) - V(a',i,b',j)
!
!       N.B. Kvajt odflaknute

!
        implicit none
        integer dima,dimb,no
        real*8 Va(1:dima,1:no,1:dimb,1:no)
        real*8 V(1:dima,1:no,1:dimb,1:no)
!
!       help variables
        integer a,b,i,j
!
        do j=1,no
        do b=1,dimb
        do i=1,no
        do a=1,dima
          Va(a,i,b,j)=2.0d0*V(a,j,b,i)-V(a,i,b,j)
        end do
        end do
        end do
        end do
!
        return
        end
