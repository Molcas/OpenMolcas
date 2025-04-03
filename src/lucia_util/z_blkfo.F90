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

!#define _DEBUGPRINT_
subroutine Z_BLKFO(ISPC,ISM,IATP,IBTP,NBATCH,NBLOCK)
! Construct information about batch and block structure of CI space
! defined by ISPC,ISM,IATP,IBTP.
!
! Output is given in the form of pointers to vectors in WORK
! where the info is stored :
!
! CLBT : Length of each Batch (in blocks)
! CLEBT : Length of each Batch (in elements)
! CI1BT : Length of each block
! CIBT  : Info on each block
! CBLTP : BLock type for each symmetry
!
! NBATCH : Number of batches
! NBLOCK : Number of blocks
!
! Jeppe Olsen, Feb. 98

use lucia_data, only: Allocate_Local_Arrays, CBLTP, CI1BT, CIBT, CLBT, CLEBT, IDC, MXNTTS, NIRREP, NSTSO, NOCTYP
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ISPC, ISM, IATP, IBTP
integer(kind=iwp), intent(out) :: NBATCH, NBLOCK
integer(kind=iwp) :: NOCTPA, NOCTPB
integer(kind=iwp), allocatable :: LCIOIO(:), SVST(:)

! Some dummy initializations
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ==================='
write(u6,*) ' Output from Z_BLKFO'
write(u6,*) ' ==================='
write(u6,*)
write(u6,*) ' ISM, ISPC = ',ISM,ISPC
#endif

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Allocate local arrays
call Allocate_Local_Arrays(MXNTTS,NIRREP)
! These should be preserved after exit so put mark for flushing here
! Info needed for generation of block info
call mma_allocate(LCIOIO,NOCTPA*NOCTPB,Label='LCIOIO')
call IAIBCM(ISPC,LCIOIO)
call mma_allocate(SVST,1,Label='SVST')
call ZBLTP(ISM,NIRREP,IDC,CBLTP,SVST)
call mma_deallocate(SVST)
! Allowed length of each batch
!if (ISIMSYM == 0) then
!  LBLOCK = MXSOOB
!else
!  LBLOCK = MXSOOB_AS
!end if

!LBLOCK = max(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
!if (ENVIRO == 'RASSCF') then
!  LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
!  if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))
!end if

!#ifdef _DEBUGPRINT_
!write(u6,*) ' LBLOCK = ',LBLOCK
!#endif

! Batches of C vector
call PART_CIV2(IDC,NSTSO(IATP)%A,NSTSO(IBTP)%A,NOCTPA,NOCTPB,NIRREP,LCIOIO,ISM,NBATCH,CLBT,CLEBT,CI1BT,CIBT,0)
! Number of BLOCKS
NBLOCK = CI1BT(NBATCH)+CLBT(NBATCH)-1
#ifdef _DEBUGPRINT_
write(u6,*) ' Number of batches',NBATCH
write(u6,*) ' Number of blocks ',NBLOCK
#endif
! Length of each block
call EXTRROW(CIBT,8,8,NBLOCK,CI1BT)

call mma_deallocate(LCIOIO)

end subroutine Z_BLKFO
