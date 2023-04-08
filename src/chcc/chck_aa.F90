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
        subroutine Chck_AA (A)
!
!        check T(a,b,i,j)
!
        implicit none
#include "chcc1.fh"
!        real*8 T(1:nv,1:nv,1:no,1:no)
         real*8 A(1:no*(no+1)/2,no,no)
!
        integer j,i,ij,u,v,bad
        real*8 s
!
        bad=0
        do v=1,no
        do u=1,no
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
!
           s=Ac(i,j,u,v)
!
          if (abs(A(ij,u,v)-s).gt.1.0d-10) then
            bad=bad+1
!    A(ij,u,v)=s
          end if
!
        end do
        end do
        end do
        end do
!
        write (6,*) ' Chck AA :',bad
!
        return
        end
