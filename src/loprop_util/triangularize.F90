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

subroutine Triangularize(A_sq,A_tr,n,Fold)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 A_sq(n,n), A_tr(n*(n+1)/2)
logical Fold

Fact = One
if (Fold) Fact = Two
ij = 1
do i=1,n
  do j=1,i-1
    A_tr(ij) = Fact*A_sq(i,j)
    ij = ij+1
  end do
  A_tr(ij) = A_sq(i,j)
  ij = ij+1
end do

return

end subroutine Triangularize
