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

subroutine grow_l1(l1,tmp,dima,nc,no,nv,last)
! this routine does:
!
! grow Cholesky vectors L1(m,i,a) by the segment in tmp

implicit none
integer a, dima, nc, nv, no, i, m, last
real*8 l1(nc,no,nv)
real*8 tmp(nc,no,dima)

!mp write(6,*) 'grow_l1i ',dima

do a=1,dima
  do i=1,no
    do m=1,nc
      l1(m,i,last+a) = tmp(m,i,a)
    end do
  end do
end do

return

end subroutine grow_l1
