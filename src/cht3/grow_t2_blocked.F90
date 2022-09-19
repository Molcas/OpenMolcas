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
        subroutine grow_t2_blocked(t2,tmp,dima,dimb,nv,no,              &
     & lasta,lastb,length1,length2,grpa,grpb,sym,switch)
!
! this routine do :
!
!
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
!mp        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
        logical sym
        logical switch
!
!mp        write (6,*) 'grow_t2neq dima , dimb  ',dima,dimb
!mp        write (6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
!mp        write (6,*) 'grow_t2neq no           ',no
!
!mp        if (lasta.eq.lastb) then
!?        if (grpa.eq.grpb) then
!?        do j=1,no
!?        do i=1,no
!?        do a=1,dima
!?        do b=1,a
!?        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
!?        if (a.ne.b) t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
!?        end do
!?        end do
!?        end do
!?        end do
!
!?        else
!
        do j=1,no
        do i=1,no
        do b=1,dimb
        do a=1,dima
          if (.not.switch) then
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
          else
!mp!        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,j,i)
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
          end if
!mpn        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
!
        if (sym) then
        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
        end if
!
        end do
        end do
        end do
        end do
!
!?        end if
!
        return
! Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
