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

subroutine ICOPVE(IFROM,ITO,NDIM)
! COPY INTEGER ARRAY

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IFROM(*), ITO(*), NDIM
integer(kind=iwp) :: I

do I=1,NDIM
  ITO(I) = IFROM(I)
end do

end subroutine ICOPVE
