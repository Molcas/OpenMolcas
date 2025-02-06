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

subroutine ISTVC2(IVEC,IBASE,IFACT,NDIM)
! IVEC(I) = IBASE + IFACT * I

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: NDIM, IVEC(NDIM), IBASE, IFACT
integer(kind=iwp) :: I

do I=1,NDIM
  IVEC(I) = IBASE+IFACT*I
end do

end subroutine ISTVC2
