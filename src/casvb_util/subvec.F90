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

subroutine SUBVEC(A,B,C,N)

implicit real*8(A-H,O-Z)
dimension A(N), B(N), C(N)

do I=1,N
  A(I) = B(I)-C(I)
end do

return

end subroutine SUBVEC
