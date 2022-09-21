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

subroutine grow_t2_fblocked2(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2)

implicit none
integer a, b, dima, dimb, no, i, j
!integer ij
integer lasta, lastb
integer length1, length2
!mp real*8 t2(nv,nv,no,no)
real*8 t2(length1,length2,no,no)
real*8 tmp(dima,dimb,no,no)

!mp write(6,*) 'lasta+dima, length1 ',lasta+dima,length1
!mp write(6,*) 'lastb+dimb, length2 ',lastb+dimb,length2
!mp write(6,*) 'lasta, lastb ',lasta,lastb
!mp write(6,*) 'dima, dimb ',dima,dimb

!mp ij = 0
!mp do i=2,no
!mp   do j=1,i-1
!mp     ij = ij+1
do i=1,no
  do j=1,no
    do b=1,dima
      do a=1,dimb

        !mp t2(lasta+a,lastb+b,ij) = tmp(b,a,j,i)+(-1.0d0*tmp(b,a,i,j))
        t2(lasta+a,lastb+b,i,j) = tmp(b,a,j,i)

      end do
    end do
  end do
end do

return

end subroutine grow_t2_fblocked2
