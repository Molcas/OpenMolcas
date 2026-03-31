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

use Cntrl, only: LSYM1, LSYM2
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NBASF, NCMO, NOSH, NTDMAB, NTDMZZ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: CMOA(NCMO), CMOB(NCMO), TDMAB(NTDMAB), TDMZZ(NTDMZZ)
integer(kind=iwp) :: iRC
integer(kind=iwp) :: IST, ISTCA, ISTCB, ISTCMO(8), ISTTA, ISTTZ, ISY1, ISY12, ISY2, NB1, NB2, NO1, NO2, NSCR
real(kind=wp), allocatable :: SCR(:)

if (iRC == 0) then
  TDMZZ(:) = Zero
  return
end if
! ISTCMO()=START INDEX FOR CMO ARRAY SYMMETRY BLOCKS.
ISY12 = MUL(LSYM1,LSYM2)
! NSCR=SIZE NEEDED FOR TEMPORARY MATRIX PRODUCT.
NSCR = 0
IST = 1
do ISY1=1,nIrrep
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
do ISY1=1,nIrrep
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
