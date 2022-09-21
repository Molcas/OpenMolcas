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
subroutine grow_t2anti_blocked1(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2)

implicit none
integer a, b, dima, dimb, no, i, j, ij
integer lasta, lastb
integer length1, length2
!mp real*8 t2(nv,nv,no,no)
real*8 t2(length1,length2,(no-1)*no/2)
real*8 tmp(dima,dimb,no,no)

!mp write(6,*) 'lasta+dima, length1 ',lasta+dima,length1
!mp write(6,*) 'lastb+dimb, length2 ',lastb+dimb,length2

!mp write(6,'(A,2(i10,x),i3)') 'length1, length2, no ',length1,length2,no
!mp write(6,'(A,4(i3,x))') 'dima, dimb, lasta, lastb',dima,dimb,lasta,lastb

ij = 0
do i=2,no
  do j=1,i-1
    ij = ij+1
    do b=1,dimb
      do a=1,dima

        t2(lasta+a,lastb+b,ij) = tmp(a,b,i,j)+(-1.0d0*tmp(a,b,j,i))

      end do
    end do
  end do
end do

return

end subroutine grow_t2anti_blocked1
