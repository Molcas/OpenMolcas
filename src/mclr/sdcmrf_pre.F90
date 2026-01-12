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

subroutine SDCMRF_PRE(IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB,IPACK,FACTOR,SGN)

use Index_Functions, only: nTri_Elem
use Constants, only: Two, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IATP, IBTP, IASM, IBSM, NA, NB, IDC, ISGVST
real(kind=wp), intent(in) :: PS, PL
integer(kind=iwp), intent(out) :: LDET, LCOMB, IPACK
real(kind=wp), intent(out) :: FACTOR, SGN
real(kind=wp), parameter :: SQRT2 = sqrt(Two), SQRT2I = sqrt(Half)

! Is combination array packed ?

IPACK = 0
FACTOR = One
SGN = One

if (((IDC == 2) .or. (IDC == 4)) .and. (IASM == IBSM)) then
  SGN = PS
  FACTOR = SQRT2
  if (IATP == IBTP) IPACK = 1
else if ((IDC == 4) .and. (IASM == ISGVST)) then
  if (IATP == IBTP) IPACK = 1
  SGN = PS*PL
  FACTOR = Two
end if

LDET = NA*NB
if (IPACK == 0) then
  LCOMB = LDET
else
  LCOMB = nTri_Elem(NA)
end if
if ((IDC == 4) .and. (IPACK == 0)) FACTOR = SQRT2

end subroutine SDCMRF_PRE
