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

subroutine SWAPVE(VEC1,VEC2,NDIM)
! SWAP ELEMENTS OF VECTORS VEC1 AND VEC2

implicit real*8(A-H,O-Z)
dimension VEC1(*), VEC2(*)

do I=1,NDIM
  BUF = VEC1(I)
  VEC1(I) = VEC2(I)
  VEC2(I) = BUF
end do

end subroutine SWAPVE
