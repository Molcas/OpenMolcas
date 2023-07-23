!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine CD_Tester_CPPF(PDM,X,n)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: PDM(n*(n+1)/2), X(n,n)
integer(kind=iwp) :: i, j
! Statement function
integer(kind=iwp) :: iTri
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

do j=1,n
  X(j,j) = PDM(iTri(j,j))
  do i=j+1,n
    X(i,j) = PDM(iTri(i,j))
    X(j,i) = X(i,j)
  end do
end do

end subroutine CD_Tester_CPPF
