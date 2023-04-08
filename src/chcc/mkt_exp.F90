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
        subroutine MkT_exp (T2red,T2full,dimbe,no)
!
!        this routine produce expanded set of amplitudes from reduced set
!        (used in parallel case)
!       T2red(be'ga',u,v) -> T2full(be',ga',u,v)
!        for beGrp=gaGrp
!
        implicit none
        integer dimbe,no
        real*8 T2full(1:dimbe,1:dimbe,1:no,1:no)
        real*8 T2red(1:dimbe*(dimbe+1)/2,1:no,1:no)
!
!        help variables
        integer i,j,be,ga,bega
!
        do j=1,no
        do i=1,no
!
        bega=0
        do be=1,dimbe
        do ga=1,be
        bega=bega+1
!
!dir          T2red(bega,i,j)=T2Full(be,ga,i,j)
!inv
          T2Full(be,ga,i,j)=T2red(bega,i,j)
          T2Full(ga,be,j,i)=T2red(bega,i,j)
!
        end do
        end do
!
        end do
        end do
!
        return
        end
