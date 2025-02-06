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

subroutine COPVEC(FROM,TO,NDIM)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NDIM
real(kind=wp) :: FROM(NDIM), TO(NDIM)
integer(kind=iwp) :: I

do I=1,NDIM
  TO(I) = FROM(I)
end do

end subroutine COPVEC
