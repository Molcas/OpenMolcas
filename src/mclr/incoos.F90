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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine INCOOS(IDC,IBLTP,NOOS,NOCTPA,NOCTPB,ISTSM,ISTTA,ISTTB,NSM,IENSM,IENTA,IENTB,IACOOS,MXLNG,IFINI,NBLOCK,INCFST,IOCOC)
! Obtain Number of OOS blocks that can be included
! IN MXLNG word starting from block after ISTSM,ISTTA,ISTTB
! Activated blocks are given in IACOOS
! Last activated block is (IENSM,IENTA,IENTB)
! If all blocks have been accessed IFINI is returned as 1
! Diagonal blocks are expanded
!
! Jeppe Olsen, Winter of 1991

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IDC, IBLTP(*), NOCTPA, NOCTPB, NSM, NOOS(NOCTPA,NOCTPB,NSM), MXLNG, INCFST, IOCOC(NOCTPA,NOCTPB)
integer(kind=iwp), intent(inout) :: ISTSM, ISTTA, ISTTB
integer(kind=iwp), intent(out) :: IENSM, IENTA, IENTB, IACOOS(NOCTPA,NOCTPB,NSM), IFINI, NBLOCK
integer(kind=iwp) :: IA, IB, IPA, IPB, IPSM, ISM, LBLOCK, LENGTH
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ISMST
#endif
logical(kind=iwp) :: Skip

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ================'
write(u6,*) ' INCOOS in action'
write(u6,*) ' ================'
write(u6,*)
write(u6,*) ' NOOS(NOCTPA,NOCTPB,NSM) array (input)'
write(u6,*)
do ISMST=1,NSM
  write(u6,*) ' ISMST = ',ISMST
  call IWRTMA(NOOS(:,:,ISMST),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
end do
#endif

IPA = 0
IPB = 0
IPSM = 0

! Initialize
IACOOS(:,:,:) = 0
ISM = ISTSM
IA = ISTTA
IB = ISTTB
LENGTH = 0
NBLOCK = 0
IENSM = ISTSM
IENTA = ISTTA
IENTB = ISTTB
IFINI = 0
if (INCFST == 1) then
  Skip = .true.
else
  Skip = .false.
end if
do
  if (Skip) then
    Skip = .false.
  else
    ! Next block
    IPA = IA
    IPB = IB
    IPSM = ISM

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
          IFINI = 1
        end if
      end if
    end if
    if (IFINI == 1) exit
  end if
  ! Should this block be included
  if (IDC /= 1) then
    if ((IBLTP(ISM) == 0) .or. ((IBLTP(ISM) == 2) .and. (IA < IB))) cycle
  end if
  if (IOCOC(IA,IB) == 0) cycle
  !write(u6,*) ' INCOOS IDC IBLTP ',IDC,IBLTP(ISM)
  ! can this block be included
  LBLOCK = NOOS(IA,IB,ISM)
  !write(u6,*) ' IA IB ISM LBLOCK ',IA,IB,ISM,LBLOCK
  if (LENGTH+LBLOCK > MXLNG) then
    IA = IPA
    IB = IPB
    ISM = IPSM
    exit
  end if
  NBLOCK = NBLOCK+1
  LENGTH = LENGTH+LBLOCK
  IACOOS(IA,IB,ISM) = 1
  if (NBLOCK == 1) then
    ISTTA = IA
    ISTTB = IB
    ISTSM = ISM
  end if
end do

IENSM = ISM
IENTA = IA
IENTB = IB
if ((IFINI == 0) .and. (NBLOCK == 0)) then
  write(u6,*) ' Not enough scratch space to include a single Block'
  write(u6,*) ' Since I cannot proceed I will stop'
  write(u6,*) ' Insufficient buffer detected in INCOOS'
  write(u6,*) ' Alter RAS space of raise Buffer from ',MXLNG
  call SYSABENDMSG('lucia_util/incoos','Internal error',' ')
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'Output from INCOOS'
write(u6,*) '=================='
write(u6,*) ' Length and number of included blocks ',LENGTH,NBLOCK
do ISM=ISTSM,IENSM
  write(u6,*) ' Active blocks of symmetry ',ISM
  call IWRTMA(IACOOS(1,1,ISM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
end do
if (IFINI == 1) write(u6,*) ' No new blocks'
#endif

return

end subroutine INCOOS
