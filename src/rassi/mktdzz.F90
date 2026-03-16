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

subroutine MKTDZZ(CMOA,CMOB,TDMAB,TDMZZ,iRC)

use stdalloc, only: mma_allocate, mma_deallocate
use Cntrl, only: LSYM1, LSYM2
use Symmetry_Info, only: nSym => nIrrep, MUL
use rassi_data, only: NCMO, NTDMAB, NTDMZZ, NBASF, NOSH
use Constants, only: Zero, One

implicit none
real*8 CMOA(NCMO), CMOB(NCMO)
real*8 TDMAB(NTDMAB), TDMZZ(NTDMZZ)
integer iRC
integer ISTCMO(8)
real*8, allocatable :: SCR(:)
integer ISY12, NSCR, IST, ISY1, NO1, ISY2, ISTTA, ISTCA, ISTTZ, ISTCB, NO2, NB1, NB2

if (iRC == 0) then
  TDMZZ(:) = Zero
  return
end if
! ISTCMO()=START INDEX FOR CMO ARRAY SYMMETRY BLOCKS.
ISY12 = MUL(LSYM1,LSYM2)
! NSCR=SIZE NEEDED FOR TEMPORARY MATRIX PRODUCT.
NSCR = 0
IST = 1
do ISY1=1,NSYM
  ISTCMO(ISY1) = IST
  NO1 = NOSH(ISY1)
  IST = IST+NO1*NBASF(ISY1)
  ISY2 = MUL(ISY1,ISY12)
  NSCR = max(NSCR,NO1*NBASF(ISY2))
end do
call mma_allocate(SCR,NSCR,Label='SCR')
ISTTA = 1
ISTCA = 1
ISTTZ = 1
do ISY1=1,NSYM
  ISY2 = MUL(ISY1,ISY12)
  ISTCB = ISTCMO(ISY2)
  NO1 = NOSH(ISY1)
  NO2 = NOSH(ISY2)
  NB1 = NBASF(ISY1)
  NB2 = NBASF(ISY2)
  if (NB1*NB2 /= 0) then
    if (NO1*NO2 == 0) then
      call FZERO(TDMZZ(ISTTZ),NB1*NB2)
    else
      call DGEMM_('N','T',NO1,NB2,NO2,One,TDMAB(ISTTA),NO1,CMOB(ISTCB),NB2,Zero,SCR,NO1)
      call DGEMM_('N','N',NB1,NB2,NO1,One,CMOA(ISTCA),NB1,SCR,NO1,Zero,TDMZZ(ISTTZ),NB1)
      ISTTA = ISTTA+NO1*NO2
    end if
  end if
  ISTCA = ISTCA+NB1*NO1
  ISTTZ = ISTTZ+NB1*NB2
end do
call mma_deallocate(SCR)

end subroutine MKTDZZ
