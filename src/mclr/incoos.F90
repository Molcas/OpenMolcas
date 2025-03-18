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

subroutine INCOOS(IDC,IBLTP,NOOS,NOCTPA,NOCTPB,ISTSM,ISTTA,ISTTB,NSMST,IENSM,IENTA,IENTB,IACOOS,MXLNG,IFINI,NBLOCK,INCFST,IOCOC)
! Obtain Number of OOS blocks that can be included
! IN MXLNG word starting from block after ISTSM,ISTTA,ISTTB
! Activated blocks are given in IACOOS
! Last activated block is (IENSM,IENTA,IENTB)
! If all blocks have been accessed IFINI is returned as 1
! Diagonal blocks are expanded
!
! Jeppe Olsen, Winter of 1991

implicit real*8(A-H,O-Z)
! Input
integer NOOS(NOCTPA,NOCTPB,NSMST)
integer IOCOC(NOCTPA,NOCTPB)
! May 7
integer IBLTP(*)
! May 7
! Output
integer IACOOS(NOCTPA,NOCTPB,NSMST)

NTEST = 00
if (NTEST >= 100) then
  write(6,*)
  write(6,*) ' ================'
  write(6,*) ' INCOOS in action'
  write(6,*) ' ================'
  write(6,*)
  write(6,*) ' NOOS(NOCTPA,NOCTPB,NSMST) array (input)'
  write(6,*)
  do ISMST=1,NSMST
    write(6,*) ' ISMST = ',ISMST
    call IWRTMA(NOOS(1,1,ISMST),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  end do
end if

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
if (INCFST == 1) goto 999
1000 continue
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
    if (ISM < NSMST) then
      ISM = ISM+1
    else
      IFINI = 1
    end if
  end if
end if
if (IFINI == 1) goto 1001
! Should this block be included
999 continue
if ((IDC /= 1) .and. (IBLTP(ISM) == 0)) goto 1000
if ((IDC /= 1) .and. (IBLTP(ISM) == 2) .and. (IA < IB)) goto 1000
if (IOCOC(IA,IB) == 0) goto 1000
!write(6,*) ' INCOOS IDC IBLTP ',IDC,IBLTP(ISM)
! can this block be included
LBLOCK = NOOS(IA,IB,ISM)
!write(6,*) ' IA IB ISM LBLOCK ',IA,IB,ISM,LBLOCK
if (LENGTH+LBLOCK <= MXLNG) then
  NBLOCK = NBLOCK+1
  LENGTH = LENGTH+LBLOCK
  IACOOS(IA,IB,ISM) = 1
  if (NBLOCK == 1) then
    ISTTA = IA
    ISTTB = IB
    ISTSM = ISM
  end if
  goto 1000
else
  IA = IPA
  IB = IPB
  ISM = IPSM
end if
1001 continue

IENSM = ISM
IENTA = IA
IENTB = IB
if ((IFINI == 0) .and. (NBLOCK == 0)) then
  write(6,*) ' Not enough scratch space to include a single Block'
  write(6,*) ' Since I cannot proceed I will stop'
  write(6,*) ' Insufficient buffer detected in INCOOS'
  write(6,*) ' Alter RAS space of raise Buffer from ',MXLNG
  call SYSABENDMSG('lucia_util/incoos','Internal error',' ')
end if

if (NTEST /= 0) then
  write(6,*) 'Output from INCOOS'
  write(6,*) '=================='
  write(6,*) ' Length and number of included blocks ',LENGTH,NBLOCK
end if
if (NTEST >= 2) then
  do ISM=ISTSM,IENSM
    write(6,*) ' Active blocks of symmetry ',ISM
    call IWRTMA(IACOOS(1,1,ISM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  end do
  if (IFINI == 1) write(6,*) ' No new blocks'
end if

return

end subroutine INCOOS
