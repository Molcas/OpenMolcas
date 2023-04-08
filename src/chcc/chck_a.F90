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
        subroutine Chck_A (AA)
!
!        check AA(ij,u,v)
!
        implicit none
#include "chcc1.fh"
        real*8 AA(1:no*(no+1)/2,1:no,1:no)
!
!        help var
        integer i,j,ij,u,v,a,b,bad
        real*8 s
!
        bad=0
!
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
!
        do u=1,no
        do v=1,no
!
          s=Q0(u,i,v,j)
!
          do a=1,nv
           s=s+Q1(a,j,u,i)*T1c(a,v)
          end do
!
          do a=1,nv
           s=s+Q1(a,i,v,j)*T1c(a,u)
          end do
!
          do a=1,nv
          do b=1,nv
           s=s+Q21(a,i,b,j)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
          end do
!
          Ac(i,j,u,v)=s
!
!          write (6,99) i,j,u,v,AA(ij,u,v),s,AA(ij,u,v)-s
!99          format (4(i2,1x),3(f15.10,1x))
          if (abs(AA(ij,u,v)-s).gt.1.0d-10) then
          bad=bad+1
          end if
!
        end do
        end do
        end do
        end do
!
!
        do u=1,no
        do v=1,no
        do i=2,no
        do j=1,i-1
          Ac(j,i,v,u)=Ac(i,j,u,v)
        end do
        end do
        end do
        end do
!
        write (6,*) ' A   Chck :',bad
!
        return
        end
