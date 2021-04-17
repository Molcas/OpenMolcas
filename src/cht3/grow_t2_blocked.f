************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine grow_t2_blocked(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,grpa,grpb,sym,switch)
c
c this routine do :
c
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
cmp        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
        logical sym
        logical switch
c
cmp        write (6,*) 'grow_t2neq dima , dimb  ',dima,dimb
cmp        write (6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
cmp        write (6,*) 'grow_t2neq no           ',no
c
cmp        if (lasta.eq.lastb) then
c?        if (grpa.eq.grpb) then
c?        do j=1,no
c?        do i=1,no
c?        do a=1,dima
c?        do b=1,a
c?        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
c?        if (a.ne.b) t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
c?        end do
c?        end do
c?        end do
c?        end do
c
c?        else
c
        do j=1,no
        do i=1,no
        do b=1,dimb
        do a=1,dima
          if (.not.switch) then
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
          else
cmp!        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,j,i)
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
          end if
cmpn        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
c
        if (sym) then
        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
        end if
c
        end do
        end do
        end do
        end do
c
c?        end if
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
