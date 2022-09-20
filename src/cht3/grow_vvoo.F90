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

subroutine grow_vvoo(A,B,no,nv,dima,dimb,lasta,lastb)
! this routine does:
!
! grow A(1324)/(vvoo) from the blocked cholesky vectors B(1234)/(vovo)

implicit none
integer i1, i2, i3, i4, dima, dimb, no, nv
integer lasta, lastb
real*8 A(1:nv,1:nv,1:no,1:no), B(1:dima,1:no,1:dimb,1:no)

!!write(6,'(A,4(i3,x))') 'AA lasta, lastb, dima, dimb ',lasta,lastb,dima,dimb

do i4=1,no
  do i3=1,no
    do i1=1,dima
      do i2=1,dimb
        A(lasta+i1,lastb+i2,i3,i4) = B(i1,i3,i2,i4)
      end do
    end do
  end do
end do

return

end subroutine grow_vvoo
