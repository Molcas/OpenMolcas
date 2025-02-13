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

subroutine SCLDIA(A,FACTOR,NDIM,IPACK)
! scale diagonal of square matrix A
!
! IPACK = 0 : full matrix
! IPACK /= 0 : Lower triangular packed matrix
!              assumed packed columnwise !!!!

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: A(*)
real(kind=wp), intent(in) :: FACTOR
integer(kind=iwp), intent(in) :: NDIM, IPACK
integer(kind=iwp) :: I, II

if (IPACK == 0) then
  do I=1,NDIM
    II = (I-1)*NDIM+I
    A(II) = A(II)*FACTOR
  end do
else
  II = 1
  do I=1,NDIM
    A(II) = A(II)*FACTOR
    II = II+NDIM-I+1
  end do
end if

end subroutine SCLDIA
