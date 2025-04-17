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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

subroutine GATRMT(MATIN,NROWIN,NCOLIN,MATUT,NROWUT,NCOLUT,ISCA,SSCA)
! Gather rows of transposed matrix MATIN  to  MATUT
!
! MATUT(I,J) = SSCA(I)*MATIN(J,ISCA(I)),(ISCA(I) /= 0)
!
! Jeppe Olsen, Getting LUCIA in shape, Feb1994

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NROWIN, NCOLIN, NROWUT, NCOLUT, ISCA(*)
real(kind=wp), intent(in) :: MATIN(NCOLIN,NROWIN), SSCA(*)
real(kind=wp), intent(out) :: MATUT(NROWUT,NCOLUT)
integer(kind=iwp) :: I, ICBL, ICEND, ICOFF, IROW, LBLK, NBLK
real(kind=wp) :: S

! (MATIN Transposed)

! To get rid of annoying and incorrect compiler warnings
ICOFF = 0
!LBLK = 100
LBLK = 40
NBLK = NCOLUT/LBLK
if (LBLK*NBLK < NCOLUT) NBLK = NBLK+1
do ICBL=1,NBLK
  if (ICBL == 1) then
    ICOFF = 1
  else
    ICOFF = ICOFF+LBLK
  end if
  ICEND = min(ICOFF+LBLK-1,NCOLUT)
  do I=1,NROWUT
    IROW = ISCA(I)
    if (IROW == 0) then
      MATUT(I,ICOFF:ICEND) = Zero
    else
      S = SSCA(I)
      MATUT(I,ICOFF:ICEND) = S*MATIN(ICOFF:ICEND,IROW)
    end if
  end do
end do

return

end subroutine GATRMT
