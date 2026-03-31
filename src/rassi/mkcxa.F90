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

subroutine MKCXA(NSYM,NOSH,NCXA,TRA,CXA)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSYM, NOSH(NSYM), NCXA
real(kind=wp), intent(in) :: TRA(NCXA)
real(kind=wp), intent(out) :: CXA(NCXA)
integer(kind=iwp) :: I, ISTA, NDIMEN

ISTA = 1
do I=1,NSYM
  NDIMEN = NOSH(I)
  if (NDIMEN > 0) then
    call MKCXAL(NDIMEN,TRA(ISTA),CXA(ISTA))
    ISTA = ISTA+NDIMEN**2
  end if
end do

end subroutine MKCXA
