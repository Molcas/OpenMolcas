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

!#define _DEBUGPRINT_
subroutine NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSM,IBLTP,IDC,NONEW,IOCOC)
! Obtain allowed block following IATP IBTP IASM

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(inout) :: IATP, IBTP, IASM
integer(kind=iwp), intent(in) :: NOCTPA, NOCTPB, NSM, IBLTP(*), IDC, IOCOC(NOCTPA,NOCTPB)
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: IA, IB, ISM

! Initialize
ISM = IASM
IA = IATP
IB = IBTP
NONEW = 0
do
  ! Next block
  if (IB < NOCTPB) then
    IB = IB+1
  else
    IB = 1
    if (IA < NOCTPA) then
      IA = IA+1
    else
      IA = 1
      if (ISM < NSM) then
        ISM = ISM+1
      else
        NONEW = 1
      end if
    end if
  end if
  if (NONEW == 1) exit
  ! Should this block be included
  if ((IDC /= 1) .and. (IBLTP(ISM) == 0)) then
    !continue
  else if ((IDC /= 1) .and. (IBLTP(ISM) == 2) .and. (IA < IB)) then
    !continue
  else if (IOCOC(IA,IB) == 0) then
    !continue
  else
    exit
  end if
end do

IATP = IA
IBTP = IB
IASM = ISM

#ifdef _DEBUGPRINT_
write(u6,'(A,4I4)') ' NXTBLK : ISM IA IB NONEW ',IASM,IA,IB,NONEW
#endif

return

end subroutine NXTBLK_MCLR
