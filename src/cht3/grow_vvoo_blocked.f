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
        subroutine grow_vvoo_blocked(AA,BB,no,nv,dima,dimb,lasta,lastb,
     & length1,length2,a,b,sym)
c
c this routine do :
c
c grow A(1324)/(vvoo) from the blocked cholesky vectors
c B(1234)/(vovo)
c
        implicit none
        integer i1,i2,i3,i4,dima,dimb,no,nv
        integer lasta,lastb
        integer length1,length2,a,b
        real*8 AA(1:length1,1:length2,1:no,1:no)
        real*8 BB(1:dima,1:no,1:dimb,1:no)
        logical sym
c
cmp        write (6,'(A,4(i3,x))') 'AA lasta, lastb, dima, dimb ',
cmp     & lasta,lastb,dima,dimb
cmp        write (6,'(A,2(i9,x))') 'chk_a ',lasta+dima,length1
cmp        write (6,'(A,2(i9,x))') 'chk_b ',lastb+dimb,length2
c
        do i4=1,no
        do i3=1,no
        do i1=1,dima
        do i2=1,dimb
        AA(lasta+i1,lastb+i2,i3,i4)=BB(i1,i3,i2,i4)
c
        if (sym) AA(lastb+i2,lasta+i1,i4,i3)=BB(i1,i3,i2,i4)
c
        end do
        end do
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(a)
        call Unused_integer(b)
      end if
        end
