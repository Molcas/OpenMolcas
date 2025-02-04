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

subroutine T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
! Evaluate T**(NK_operator) times vector on file LUIN
! to yield vector on file LUOUT
! (NK_operator is number operator for orbital K)
!
! Note LUCIN and LUCOUT are both rewinded before read/write
! Input
! =====
!  T : Input constant
!  KORB : Orbital in symmetry order
!
!  ISM,ISPC : Symmetry and space of state on LUIN
!  C : Scratch block
!
! Jeppe Olsen, Feb. 98

use stdalloc, only: mma_allocate, mma_deallocate
use Local_Arrays, only: CIBT, CBLTP, Deallocate_Local_Arrays
use strbas, only: NSTSO
use lucia_data, only: ICISTR
use lucia_data, only: MXNSTR
use lucia_data, only: IREOST, NTOOB
use lucia_data, only: NELEC
use csm_data, only: NSMST
use Definitions, only: u6

implicit none
real*8 T
integer KORB, ISM, ISPC, LUCIN, LUCOUT

! Scratch block, must hold a batch of blocks
real*8 C(*)

integer, allocatable :: LASTR(:), LBSTR(:)
integer, allocatable :: LKAOC(:), LKBOC(:)
integer NTEST, IATP, IBTP, NAEL, NBEL, KKORB, NBATCH, NBLOCK

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' T_TO_NK_VEC speaking'
  write(u6,*) ' ISM, ISPC = ',ISM,ISPC
end if
! Set up block and batch structure of vector
IATP = 1
IBTP = 2

NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)

call Z_BLKFO(ISPC,ISM,IATP,IBTP,NBATCH,NBLOCK)
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)

call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
call mma_allocate(LKAOC,MXNSTR,Label='LKAOC')
call mma_allocate(LKBOC,MXNSTR,Label='LKBOC')
! Orbital K in type ordering
KKORB = IREOST(KORB)
call T_TO_NK_VECS(T,KKORB,C,LUCIN,LUCOUT,NSTSO(IATP)%I,NSTSO(IBTP)%I,NBLOCK,CIBT,NAEL,NBEL,LASTR,LBSTR,CBLTP,NSMST,ICISTR,NTOOB, &
                  LKAOC,LKBOC)

call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LKAOC)
call mma_deallocate(LKBOC)

call Deallocate_Local_Arrays()

end subroutine T_TO_NK_VEC
