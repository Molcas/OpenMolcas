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
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************

subroutine Z_BLKFO(ISPC,ISM,IATP,IBTP,NBATCH,NBLOCK)
! Construct information about batch and block structure of CI space
! defined by ISPC,ISM,IATP,IBTP.
!
! Output is given in the form of pointers to vectors in WORK
! where the info is stored :
!
! CLBT : Length of each Batch ( in blocks)
! CLEBT : Length of each Batch ( in elements)
! CI1BT : Length of each block
! CIBT  : Info on each block
! CBLTP : BLock type for each symmetry
!
! NBATCH : Number of batches
! NBLOCK : Number of blocks
!
! Jeppe Olsen, Feb. 98

use stdalloc, only: mma_allocate, mma_deallocate
use Local_Arrays, only: CLBT, CLEBT, CI1BT, CIBT, CBLTP, Allocate_Local_Arrays
use strbas, only: NSTSO
use lucia_data, only: MXSOOB, ISMOST, MXNTTS, XISPSM
use lucia_data, only: ENVIRO, ISIMSYM, LCSBLK
use lucia_data, only: IREFSM, PSSIGN, IDC
use lucia_data, only: NOCTYP
use csm_data, only: NSMST

implicit none
integer ISPC, ISM, IATP, IBTP, NBATCH, NBLOCK
integer, allocatable :: LCIOIO(:)
integer, allocatable :: SVST(:)
integer, external :: IFRMR
integer NTEST, NOCTPA, NOCTPB, LBLOCK

! Some dummy initializations
NTEST = 0
#ifdef _DEBUGPRINT_
if (NTEST >= 100) then
  write(6,*)
  write(6,*) ' ==================='
  write(6,*) ' Output from Z_BLKFO'
  write(6,*) ' ==================='
  write(6,*)
  write(6,*) ' ISM, ISPC = ',ISM,ISPC
end if
#endif

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Allocate local arrays
call Allocate_Local_Arrays(MXNTTS,NSMST)
! These should be preserved after exit so put mark for flushing here
! Info needed for generation of block info
call mma_allocate(LCIOIO,NOCTPA*NOCTPB,Label='LCIOIO')
call IAIBCM(ISPC,LCIOIO)
call mma_allocate(SVST,1,Label='SVST')
call ZBLTP(ISMOST(1,ISM),NSMST,IDC,CBLTP,SVST)
call mma_deallocate(SVST)
! Allowed length of each batch
!if (ISIMSYM == 0) then
LBLOCK = MXSOOB
!else
!  LBLOCK = MXSOOB_AS
!end if

LBLOCK = max(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
if (ENVIRO(1:6) == 'RASSCF') then
  LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
  if (PSSIGN /= 0.0d0) LBLOCK = int(2.0d0*XISPSM(IREFSM,1))
end if

if (NTEST >= 10) write(6,*) ' LBLOCK = ',LBLOCK

! Batches of C vector
call PART_CIV2(IDC,CBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,NOCTPA,NOCTPB,NSMST,LBLOCK,LCIOIO,ISMOST(1,ISM),NBATCH,CLBT,CLEBT,CI1BT,CIBT, &
               0,ISIMSYM)
! Number of BLOCKS
NBLOCK = IFRMR(CI1BT,1,NBATCH)+IFRMR(CLBT,1,NBATCH)-1
if (NTEST >= 1) then
  write(6,*) ' Number of batches',NBATCH
  write(6,*) ' Number of blocks ',NBLOCK
end if
! Length of each block
call EXTRROW(CIBT,8,8,NBLOCK,CI1BT)

call mma_deallocate(LCIOIO)

end subroutine Z_BLKFO
