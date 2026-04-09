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
! Copyright (C) 1995,1999, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PART_CIV2(IDC,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,IOCOC,ICSM,NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP)
! Jeppe Olsen
!
! Last update : May 1999 : ISIMSYM added
!
! Partition a CI vector into batches of blocks.
! If ISIMSYM == 1, TTS blocks that differs only in symmetry are not split.
!
! IF ICOMP == 1 the complete civector is constructed
!
! Compared to PART_CIV, the NOOS arrays have been eliminated.
! They are becoming the size defining arrays - atleast at
! the laptop
!
! Output
! ======
! NBATCH  : Number of batches
! LBATCH  : Number of blocks in a given batch
! LEBATCH : Number of elements in a given batch (packed)!
! I1BATCH : Number of first block in a given batch
! IBATCH  : TTS blocks in Start of a given TTS block with respect to start of batch
!   IBATCH(1,*) : Alpha type
!   IBATCH(2,*) : Beta sym
!   IBATCH(3,*) : Sym of alpha
!   IBATCH(4,*) : Sym of beta
!   IBATCH(5,*) : Offset of block with respect to start of block in expanded form
!   IBATCH(6,*) : Offset of block with respect to start of block in packed form
!   IBATCH(7,*) : Length of block, expandend form
!   IBATCH(8,*) : Length of block, packed form
!
! Jeppe Olsen, August 1995

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IDC, NSMST, NOCTPA, NOCTPB, NSSOA(NSMST,NOCTPA), NSSOB(NSMST,NOCTPB), IOCOC(NOCTPA,NOCTPB), ICSM, &
                                 ICOMP
integer(kind=iwp), intent(out) :: NBATCH
integer(kind=iwp), intent(_OUT_) :: LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
integer(kind=iwp) :: IA, IASM, IB, IBLOCK, IBSM, IFINI, IFRST, INC, ISM, LBLOCK, LBLOCKP, LENGTH, LENGTHP, NBLOCK, NSTA, NSTB
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: JBATCH
#endif

! Dummy initialize
INC = 0
LBLOCKP = 0

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' =================='
write(u6,*) '     PART_CIV2'
write(u6,*) ' =================='
write(u6,*) ' IDC = ',IDC
write(u6,*)
write(u6,*) ' IOCOC Array'
call IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
!write(u6,*) ' IBLTP array'
!call IWRTMA(IBLTP,1,NSMST,1,NSMST)
write(u6,*) ' NSSOA, NSSOB'
call IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
call IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
write(u6,*) 'ICOMP = ',ICOMP
#endif

! block 1

IB = 1
IA = 1
ISM = 1
IFRST = 1
NBATCH = 0
IBLOCK = 0
IFINI = 0
! Loop over batches of blocks
outer: do
  NBATCH = NBATCH+1
  LBATCH(NBATCH) = 0
  I1BATCH(NBATCH) = IBLOCK+1
  LENGTH = 0
  LENGTHP = 0
  NBLOCK = 0
  IFRST = 1
  ! Loop over blocks in batch
  do
    if (IFRST == 0) then
      ! New order : ISM,IB,IA (leftmost inner loop)
      if (ISM < NSMST) then
        ISM = ISM+1
      else
        ISM = 1
        if (IB < NOCTPB) then
          IB = IB+1
        else
          IB = 1
          if (IA < NOCTPA) then
            IA = IA+1
          else
            IFINI = 1
          end if
        end if
      end if
    end if
    IFRST = 0
    if (IFINI == 1) exit outer
    if (IOCOC(IA,IB) == 0) cycle
    ! Size of TT block (all symmetries)
    !LBLOCK_AS = 0
    !if ((ISIMSYM == 1) .and. (ISM == 1)) then
    !  do IASM=1,NSMST
    !    IBSM = Mul(IASM,ICSM)
    !    NSTA = NSSOA(IASM,IA)
    !    NSTB = NSSOB(IBSM,IB)
    !    if (IBLTP(IASM) == 0) cycle
    !    if ((IBLTP(IASM) == 2) .and. (IA < IB)) cycle
    !    LBLOCK_AS = LBLOCK_AS+NSTA*NSTB
    !  end do
    !  INC = 0
    !  !write(u6,*) ' IA IB LBLOCK_AS',IA,IB,LBLOCK_AS
    !  if ((LENGTH+LBLOCK_AS <= MXLNG) .or. (ICOMP == 1)) INC = 1
    !end if
    ! Should this block be included
    IASM = ISM
    IBSM = Mul(IASM,ICSM)
    if (IDC == 2) then
      if (IA < IB) cycle
      if ((IA == IB) .and. (IASM < IBSM)) cycle
    end if
    ! can this block be included
    NSTA = NSSOA(ISM,IA)
    NSTB = NSSOB(IBSM,IB)
    LBLOCK = NSTA*NSTB
    if ((IDC == 1) .or. (IA > IB) .or. ((IA == IB) .and. (IASM > IBSM))) then
      LBLOCKP = NSTA*NSTB
    else if ((IDC == 2) .and. (IA == IB) .and. (IASM == IBSM)) then
      LBLOCKP = nTri_Elem(NSTA)
    end if
    !write(u6,*) ' IASM IBSM IA IB LBLOCKP,LBLOCK',IASM,IBSM,IA,IB,LBLOCKP,LBLOCK

    !if (ISIMSYM == 0) then
    INC = 0
    if ((LENGTH+LBLOCK <= LBLOCK) .or. (ICOMP == 1)) INC = 1
    !end if

    if (INC == 1) then
      NBLOCK = NBLOCK+1
      IBLOCK = IBLOCK+1
      LBATCH(NBATCH) = LBATCH(NBATCH)+1
      IBATCH(1,IBLOCK) = IA
      IBATCH(2,IBLOCK) = IB
      IBATCH(3,IBLOCK) = ISM
      IBATCH(4,IBLOCK) = IBSM
      IBATCH(5,IBLOCK) = LENGTH+1
      IBATCH(6,IBLOCK) = LENGTHP+1
      IBATCH(7,IBLOCK) = LBLOCK
      IBATCH(8,IBLOCK) = LBLOCKP
      LENGTH = LENGTH+LBLOCK
      LENGTHP = LENGTHP+LBLOCKP
      LEBATCH(NBATCH) = LENGTHP
    else if ((ICOMP == 0) .and. (INC == 0) .and. (NBLOCK == 0)) then
      write(u6,*) ' Not enough space to include a single Block'
      write(u6,*) ' Since I cannot proceed I will stop'
      write(u6,*) ' Insufficient space detected in PART_CIV'
      !write(u6,*) ' Alter GAS space or raise Buffer from ',MXLNG
      !stop 'Error in PART_CIV2'
      call SYSABENDMSG('lucia_util/part_civ2','Internal error','')
    else
      ! This batch is finished, go to next batch
      exit
    end if
  end do
end do outer

#ifdef _DEBUGPRINT_
!write(u6,*) 'Output from PART_CIV'
!write(u6,*) '====================='
write(u6,*)
write(u6,*) ' Number of batches ',NBATCH
do JBATCH=1,NBATCH
  write(u6,*)
  write(u6,*) ' Info on batch ',JBATCH
  write(u6,*) ' ***********************'
  write(u6,*)
  write(u6,*) '      Number of blocks included ',LBATCH(JBATCH)
  write(u6,*) '      TTSS and offsets and lengths of each block'
  do IBLOCK=I1BATCH(JBATCH),I1BATCH(JBATCH)+LBATCH(JBATCH)-1
    write(u6,'(10X,4I3,4I8)') IBATCH(:,IBLOCK)
  end do
end do
#endif

end subroutine PART_CIV2
