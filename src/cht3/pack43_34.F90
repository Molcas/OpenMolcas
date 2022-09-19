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
        subroutine pack43_34 (AA,BB,d1,d2,d3)
!
! this routine do :
!
! A(a,b,c,d) -> B(a,b,cd)  c>=d
!
        implicit none
        integer d1,d2,d3,a,b,c,d,cd
        real*8 AA(d1,d2,d3,d3),BB(d1,d2,(d3*(d3+1)/2))
!
        cd=0
        do c=1,d3
        do d=1,c
        cd=cd+1
        do b=1,d2
        do a=1,d1
!
        BB(a,b,cd)=AA(a,b,c,d)
        end do
        end do
        end do
        end do
!
        return
        end
