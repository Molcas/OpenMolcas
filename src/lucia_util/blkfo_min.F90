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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

subroutine BLKFO_MIN(ISM,NBLK,LEN_BLK)
! Number of blocks and length of each block for CI expansion
!
! Jeppe Olsen, June 2001
!
! Input
! =====
!
! ISM : Symmetry of CI expansion
!
! Output
! ======
!
! NBLK : Number of blocks in expansion
! LEN_BLK(IBLK) : Length of block IBLK

use CandS, only: ISSPC
use lucia_data, only: Allocate_Local_Arrays, CBLTP, CI1BT, CIBT, CLBT, CLEBT, Deallocate_Local_Arrays, IDC, MXNTTS, NIRREP, &
                      NOCTYP, NSTSO
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ISM
integer(kind=iwp), intent(out) :: NBLK, LEN_BLK(MXNTTS)
integer(kind=iwp) :: I_DUMMY(1), IATP, IBTP, NBATCH, NOCTPA, NOCTPB
integer(kind=iwp), allocatable :: CIOIO(:)

I_DUMMY(1) = 0 ! jwk-cleanup
IATP = 1
IBTP = 2

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Pointers to local arrays
call Allocate_Local_Arrays(MXNTTS,NIRREP)
! Info needed for generation of block info
call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
call IAIBCM(ISSPC,CIOIO) ! Jesper
call ZBLTP(ISM,NIRREP,IDC,CBLTP,I_DUMMY)
! Allowed length of each batch (not important for final output)
!LBLOCK = max(MXSOOB,LCSBLK)
! Batches  of C vector
call PART_CIV2(IDC,NSTSO(IATP)%A,NSTSO(IBTP)%A,NOCTPA,NOCTPB,NIRREP,CIOIO,ISM,NBATCH,CLBT,CLEBT,CI1BT,CIBT,0)
! Number of BLOCKS
NBLK = CI1BT(NBATCH)+CLBT(NBATCH)-1
! Length of each block
call EXTRROW(CIBT,8,8,NBLK,LEN_BLK)

call Deallocate_Local_Arrays()
call mma_deallocate(CIOIO)

end subroutine BLKFO_MIN
