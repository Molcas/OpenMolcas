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
subroutine TRACI_MASTER(JOBDISK,JOBIPH,CMOMO)

use CandS, only: ISSM, ISSPC
use lucia_data, only: IDISK, INT1, IREFSM, kvec3_length, LUC, LUDIA, LUHC, LUSC1, LUSC2, MXNTTS, MXSOOB, &
                      NCSF_PER_SYM, NROOT, NSMOB, NTOOBS, PSSIGN, VEC3, XISPSM, IH1FORM
use general_data, only: nCSF=>nConf, nDet, nLev
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(inout) :: JOBDISK
integer(kind=iwp), intent(in) :: JOBIPH
real(kind=wp), intent(in) :: CMOMO(*)
integer(kind=iwp) :: I, I_DUMMY(1), IADR, ICOL, IOFF, IROW, ISM, J, JDISK, JROOT, LBLK, LBLOCK, &
                     NDIM, NREC, IH1FORM_SAVE, IDISK_Save
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IREC, NUM_ELE
#endif
real(kind=wp), allocatable :: TransMat(:), LH1SAVE(:), VEC1(:), VEC2(:), VEC4(:)
integer(kind=iwp), allocatable :: lRec(:)

LBLK = -1

NDIM = nLev*nLev
nCSF = NCSF_PER_SYM(ISSM)
! JESPER: Should reduce I/O
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

! The three scratch  blocks
call mma_allocate(VEC1,LBLOCK,Label='VEC1')
call mma_allocate(VEC2,LBLOCK,Label='VEC2')
call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
call mma_allocate(VEC4,nCSF,Label='VEC4')

! MO-MO transformation matrix :
call mma_allocate(TransMat,NDIM,Label='TransMat')
! Copy of one-electron integrals
call mma_allocate(LH1SAVE,NDIM,Label='LH1SAVE')
! We are going to mess with the one-electron integrals, take a copy
LH1SAVE(:) = INT1(:)
IH1FORM_SAVE = IH1FORM

! The input transformation matrix, CMOMO is over the full basis but it
! has non-zero values just in the block-diagonal sections. Here we
! remove those blocks.

TransMat(:) = Zero
IOFF = 1
IADR = 1
ICOL = 1
do ISM=1,NSMOB
  if (NTOOBS(ISM) == 0) Cycle
  IROW = ICOL
  do I=1,NTOOBS(ISM)
    IADR = (ICOL-1)*nLev+IROW
    do J=1,NTOOBS(ISM)
      TransMat(IOFF+NTOOBS(ISM)*(J-1)+I-1) = CMOMO(IADR+J-1)
    end do
    ICOL = ICOL+1
  end do
  IOFF = IOFF+NTOOBS(ISM)**2
end do


Call mma_allocate(lRec,MXNTTS,Label='lRec')
call BLKFO_MIN(ISSM,NREC,LREC)

! Transfer the CI-vector to LUC

IDISK(LUC) = 0
IDISK(LUDIA) = 0
JDISK = JOBDISK
do JROOT=1,NROOT
  ! Read CI vector in SG format
  call DDAFILE(JOBIPH,2,VEC4,nCSF,JDISK)
  ! Convert to internal Lucia SD format
  call CSDTVC(VEC4,VEC1,1,ISSM,0)
  ! Write CI vector to the Lucia file LUC
  IDISK_Save=IDISK(LUC)
  call TODSCN(VEC1,NREC,LREC,LBLK,LUC)
  I_DUMMY(1) = -1
  call ITODS(I_DUMMY,1,LBLK,LUC)
  IDISK(LUC)=IDISK_Save
  ! Do the actual work
  IDISK(LUSC1) = 0
  ! Copy the CI vector from file LUC to LUSC1
  call COPVCD(LUC,LUSC1,VEC1,0,LBLK) ! From the current position
  call COPVCD(LUSC1,LUHC,VEC1,1,LBLK)! Rewind first

  ! Transform CI vector : Input on LUHC, output on LUDIA (!)
  IDISK_Save=IDISK(LUDIA)
  call TRACI_LUCIA(TransMat,LUHC,LUDIA,ISSPC,ISSM,VEC1,VEC2)
  IDISK(LUDIA)=IDISK_Save
  ! Copy CI-vector back to MOLCAS JOBIPH file
  call FRMDSCN(VEC1,NREC,LBLK,LUDIA)
  call CSDTVC(VEC2,VEC1,2,ISSM,0)
  call DDAFILE(JOBIPH,1,VEC2,nCSF,JOBDISK)
  call IFRMDS(I_DUMMY,1,LBLK,LUDIA)
end do

IDISK(LUDIA) = 0

! clean up time : copy 1-e integrals back in place
INT1(:) = LH1SAVE(:)
IH1FORM = IH1FORM_SAVE

call mma_deallocate(VEC1)
call mma_deallocate(VEC2)
call mma_deallocate(VEC3)
call mma_deallocate(VEC4)
call mma_deallocate(TransMat)
call mma_deallocate(LH1SAVE)
Call mma_deallocate(lRec)

end subroutine TRACI_MASTER
