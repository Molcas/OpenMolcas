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
        subroutine Chck_V (VV)
!
!        check  V
!
        implicit none
#include "chcc1.fh"
        real*8 VV(1:nv,1:no,1:no,1:no)
!
!        help var
        integer be,v,u,j,b,bad
        real*8 s
!
        bad=0
!
        do j=1,no
        do u=1,no
        do v=1,no
        do be=1,nv
!
          s=0.0d0
          do b=1,nv
          s=s+Q22(be,b,u,j)*T1c(b,v)
          end do
!
          if (abs(VV(be,v,u,j)-s).gt.1.0d-10) then
          VV(be,v,u,j)=s
          bad=bad+1
          end if
!
        end do
        end do
        end do
        end do
!
        write (6,*) ' V  Chck :',bad
!
        return
        end
