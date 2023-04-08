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
        subroutine AdV_G2 (G2,V,nv,dimbe,dima,no,addbe,adda,fact)
!
!       this routine do:
!       G2(be',a') <- fact sum(i) V(be',i,i,a')
!
!       N.B. Kvajt odflaknute

!
        implicit none
        integer nv,dimbe,dima,no,addbe,adda
        real*8 fact
        real*8 G2(1:nv,1:nv)
        real*8 V(1:dimbe,1:no,1:no,1:dima)
!
!       help variables
!
        integer be,a,i,afull
!
        do a=1,dima
        afull=adda+a
        do i=1,no
        do be=1,dimbe
          G2(addbe+be,afull)=G2(addbe+be,afull)+fact*V(be,i,i,a)
        end do
        end do
        end do
!
        return
        end
