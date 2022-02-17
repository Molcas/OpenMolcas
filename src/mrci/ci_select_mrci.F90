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

subroutine CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NREF, NSEL, NRROOT
real(kind=wp), intent(in) :: AREF(NREF,NREF), CISEL(NREF,*)
real(kind=wp), intent(out) :: PLEN(NREF)
integer(kind=iwp), intent(out) :: IROOT(NRROOT)
integer(kind=iwp) :: I, IR, ISEL, J, JJ, JMAX
real(kind=wp) :: PL, PMAX, SUM1, SUM2

if (NSEL == 0) return
! SELECTION BY PROJECTION ONTO SPACE SPANNED BY CISEL VECTORS. IROOT()
! IS SET TO SELECT THE NRROOT VECTORS WITH MAX PROJECTED LENGTH.
do J=1,NREF
  SUM1 = Zero
  do ISEL=1,NSEL
    SUM2 = Zero
    do I=1,NREF
      SUM2 = SUM2+AREF(I,J)*CISEL(I,ISEL)
    end do
    SUM1 = SUM1+SUM2**2
  end do
  PLEN(J) = SUM1+J*1.0e-12_wp
end do
! SELECT BY MAGNITUDE OF PLEN:
do J=1,NRROOT
  PMAX = PLEN(1)
  JMAX = 1
  do JJ=2,NREF
    if (PMAX < PLEN(JJ)) then
      PMAX = PLEN(JJ)
      JMAX = JJ
    end if
  end do
  PLEN(JMAX) = -PMAX
end do
I = 0
do IR=1,NREF
  PL = PLEN(IR)
  if (PL < Zero) then
    I = I+1
    IROOT(I) = IR
    PL = -PL
  end if
  PLEN(IR) = PL-IR*1.0e-12_wp
end do

return

end subroutine CI_SELECT_MRCI
