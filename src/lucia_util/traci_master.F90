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

use lucia_data, only: IDISK, INT1, IREFSM, kvec3_length, LUC, LUDIA, LUHC, LUSC1, MXNTTS, MXSOOB, &
                      NCSF_PER_SYM, NROOT, NSMOB, NTOOBS, PSSIGN, VEC3, XISPSM, IH1FORM
use general_data, only: nCSF=>nConf, nDet, nLev, ISSM=>STSYM
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
real(kind=wp), allocatable :: TransMat(:), LH1SAVE(:), VEC1(:), VEC2(:), CI_Vector(:), LT(:)
integer(kind=iwp), allocatable :: lRec(:)


NDIM = nLev*nLev
nCSF = NCSF_PER_SYM(ISSM)

! JESPER: Should reduce I/O
LBLK = -1
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

! The three scratch  blocks
call mma_allocate(VEC1,LBLOCK,Label='VEC1')
call mma_allocate(VEC2,LBLOCK,Label='VEC2')
call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')

call mma_allocate(CI_Vector,nCSF,Label='CI_Vector')

! Copy of one-electron integrals
call mma_allocate(LH1SAVE,NDIM,Label='LH1SAVE')
! We are going to mess with the one-electron integrals, take a copy
LH1SAVE(:) = INT1(:)
IH1FORM_SAVE = IH1FORM

! The input transformation matrix, CMOMO is over the full basis but it
! has non-zero values just in the block-diagonal sections. Here we
! remove those blocks.

! MO-MO transformation matrix :
call mma_allocate(TransMat,NDIM,Label='TransMat')
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

! Memory allocation for a matrix T
call mma_allocate(LT,nLev**2,Label='LT')

! Obtain T matrix used for transformation, for each symmetry separately
IOFF = 0
do ISM=1,NSMOB
  if (ISM == 1) then
    IOFF = 1
  else
    IOFF = IOFF+NTOOBS(ISM-1)**2
  end if
  if (NTOOBS(ISM) > 0) call PAMTMT(TransMat(IOFF),LT(IOFF),NTOOBS(ISM))
end do

Call mma_allocate(lRec,MXNTTS,Label='lRec')
call BLKFO_MIN(ISSM,NREC,LREC)

IDISK(LUDIA) = 0

JDISK = JOBDISK
do JROOT=1,NROOT
  ! Read CI vector in SG format
  call DDAFILE(JOBIPH,2,CI_Vector,nCSF,JDISK)

  IDISK(LUC) = 0
  ! Convert to internal Lucia SD format
  call CSDTVC(CI_Vector,VEC1,1,ISSM,0)
  ! Write CI vector to the Lucia file LUC
  call TODSCN(VEC1,NREC,LREC,LBLK,LUC)
  I_DUMMY(1) = -1
  call ITODS(I_DUMMY,1,LBLK,LUC)
  IDISK(LUC)=0
  ! Do the actual work
  IDISK(LUSC1) = 0
  ! Copy the CI vector from file LUC to LUSC1
  call COPVCD(LUC,LUSC1,VEC1,0,LBLK) ! From the current position
  IDISK(LUC)=0
  IDISK(LUHC)=0
  call COPVCD(LUSC1,LUHC,VEC1,1,LBLK)! Rewind first

  ! Transform CI vector : Input on LUHC, output on LUDIA (!)
  IDISK(LUDIA)=0

  call TRACI_LUCIA(LT,LUHC,LUDIA,VEC1,VEC2)

  IDISK(LUDIA)=0
  ! Copy CI-vector back to MOLCAS JOBIPH file
  call FRMDSCN(VEC1,NREC,LBLK,LUDIA)
  call CSDTVC(VEC2,VEC1,2,ISSM,0)
  call IFRMDS(I_DUMMY,1,LBLK,LUDIA)

  call DDAFILE(JOBIPH,1,VEC2,nCSF,JOBDISK)
end do

IDISK(LUDIA) = 0
Call mma_deallocate(lRec)

! clean up time : copy 1-e integrals back in place
INT1(:) = LH1SAVE(:)
IH1FORM = IH1FORM_SAVE

call mma_deallocate(VEC1)
call mma_deallocate(VEC2)
call mma_deallocate(VEC3)

call mma_deallocate(CI_Vector)
call mma_deallocate(LT)
call mma_deallocate(TransMat)
call mma_deallocate(LH1SAVE)

end subroutine TRACI_MASTER
