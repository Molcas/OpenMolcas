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
        subroutine grow_t2neq(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb)
c
c this routine do :
c
c grow amplitude file t2(a,b,i,j) by the segment in tmp
c for case sa != sb
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
        integer lasta,lastb
        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
c
cmp        write (6,*) 'grow_t2neq dima , dimb  ',dima,dimb
cmp        write (6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
cmp        write (6,*) 'grow_t2neq no           ',no
c
c?        if (lasta.eq.lastb) then
c?        do j=1,no
c?        do i=1,no
c?        do a=1,dima
c?        do b=1,a
c?        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
c?cmp        if (a.ne.b) t2(lastb+b,lasta+a,j,i)=-1.0d0*tmp(a,b,j,i)
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
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
cmp        t2(lastb+b,lasta+a,j,i)=-1.0d0*tmp(b,a,j,i)
        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
        end do
        end do
        end do
        end do
c
c?        end if
c
        return
        end
