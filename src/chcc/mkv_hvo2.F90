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
        subroutine MkV_Hvo2 (V,V2,dimbe,dima,no)
!
!       this routine do:
!       Make AntiSymetric integrals
!        V2(a',i,be',j) <- [2 V(be',j|a'i) - V(be',i|a'j)]
!
        implicit none
        integer dimbe,dima,no
        real*8 V2(1:dima,1:no,1:dimbe,1:no)
        real*8 V(1:dimbe,1:no,1:dima,1:no)
!
!       help variables
        integer a,be,i,j
!
        do j=1,no
        do be=1,dimbe
        do i=1,no
        do a=1,dima
          V2(a,i,be,j)=2.0d0*V(be,j,a,i)-V(be,i,a,j)
        end do
        end do
        end do
        end do
!
!
        return
        end
