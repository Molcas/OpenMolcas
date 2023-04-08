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
        subroutine ExV_X42 (Vp,V,dimab,no)
!
!       this routine do:
!       Vp(a,b,i) <- V(a,b,i,i)
!
        implicit none
        integer dimab,no
        real*8 Vp(1:dimab,1:no)
        real*8 V(1:dimab,1:no,1:no)
!
!       help variables
        integer i,ab
!
        do i=1,no
          do ab=1,dimab
            Vp(ab,i)=V(ab,i,i)
          end do
        end do
!
!
        return
        end
