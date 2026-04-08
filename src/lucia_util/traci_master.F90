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
subroutine TRACI_MASTER(JOBDISK,JOBIPH,CMOMO,lrec)

use CandS, only: ISSM, ISSPC
use lucia_data, only: Deallocate_Local_Arrays, DTOC, IDISK, INT1, IREFSM, kvec3_length, LUC, LUDIA, LUHC, LUSC1, LUSC2, MXNTTS, &
                      MXSOOB, NCSF_PER_SYM, NROOT, NSMOB, NTOOB, NTOOBS, PSSIGN, SDREO, VEC3, XISPSM
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
integer(kind=iwp), intent(out) :: LREC(MXNTTS)
integer(kind=iwp) :: I, I_DUMMY(1), IADR, IATP, IBTP, ICOL, IOFF, IROW, ISM, J, JDISK, JROOT, LBLK, LBLOCK, NBATCH, NBLOCK, NCONF, &
                     NDIM, NREC
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IREC, NUM_ELE
#endif
real(kind=wp), allocatable :: LCMOMO(:), LH1SAVE(:), VEC1(:), VEC2(:), VEC4(:)

LBLK = -1
NDIM = NTOOB*NTOOB
NCONF = NCSF_PER_SYM(ISSM)
! JESPER: Should reduce I/O
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

! The three scratch  blocks
call mma_allocate(VEC1,LBLOCK,Label='VEC1')
call mma_allocate(VEC2,LBLOCK,Label='VEC2')
call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')
call mma_allocate(VEC4,NCONF,Label='VEC4')

! Transfer the CI-vector to LUC

call BLKFO_MIN(ISSM,NREC,LREC)
IDISK(LUC) = 0
JDISK = JOBDISK
do JROOT=1,NROOT
  call DDAFILE(JOBIPH,2,VEC4,NCONF,JDISK)
  call CSDTVC(VEC4,VEC1,1,DTOC,SDREO,ISSM,0)
# ifdef _DEBUGPRINT_
  write(u6,*) 'CI-vector written to disk for root = ',JROOT
  call WRTMAT(VEC1,1,20,1,20)
  write(u6,*) 'Writing this to disk:'
  IOFF = 1
  do IREC=1,NREC
    if (LREC(IREC) >= 0) then
      call WRTMAT(VEC1(IOFF),1,LREC(IREC),1,LREC(IREC))
      IOFF = IOFF+LREC(IREC)
    end if
  end do
# endif
  call TODSCN(VEC1,NREC,LREC,LBLK,LUC)
  I_DUMMY(1) = -1
  call ITODS(I_DUMMY,1,LBLK,LUC)
end do

! MO-MO transformation matrix :
call mma_allocate(LCMOMO,NDIM,Label='LCMOMO')
! Copy of one-electron integrals
call mma_allocate(LH1SAVE,NDIM,Label='LH1SAVE')
! We are going to mess with the one-electron integrals, take a copy
LH1SAVE(:) = INT1(:)
! Set up block structure of CI space
IATP = 1
IBTP = 2
call Z_BLKFO(ISSPC,ISSM,IATP,IBTP,NBATCH,NBLOCK)

call Deallocate_Local_Arrays()

! The input transformation matrix contains a lot of zeros which
! is expected not to be there in Traci_Lucia, so remove them.

LCMOMO(:) = Zero
IOFF = 1
IADR = 1
ICOL = 1
do ISM=1,NSMOB
  if (NTOOBS(ISM) > 0) then
    IROW = ICOL
    do I=1,NTOOBS(ISM)
      IADR = (ICOL-1)*NTOOB+IROW
      do J=1,NTOOBS(ISM)
        LCMOMO(IOFF+NTOOBS(ISM)*(J-1)+I-1) = CMOMO(IADR+J-1)
      end do
      ICOL = ICOL+1
    end do
    IOFF = IOFF+NTOOBS(ISM)**2
  end if
end do

! Now the actual work

IDISK(LUC) = 0
IDISK(LUDIA) = 0
do JROOT=1,NROOT
  IDISK(LUSC1) = 0
  call COPVCD(LUC,LUSC1,VEC1,0,LBLK)
  call COPVCD(LUSC1,LUSC2,VEC1,1,LBLK)

  ! Transform CI vector : Input on LUHC, output on LUDIA (!)
  call COPVCD(LUSC1,LUHC,VEC1,1,LBLK)

  call TRACI_LUCIA(LCMOMO,LUHC,LUDIA,ISSPC,ISSM,VEC1,VEC2)
end do
! End of loop over roots
IDISK(LUDIA) = 0

! Copy CI-vector back to MOLCAS JOBIPH file

do JROOT=1,NROOT
  call FRMDSCN(VEC1,NREC,LBLK,LUDIA)
# ifdef _DEBUGPRINT_
  NUM_ELE = sum(LREC(1:NREC))
  write(u6,*) 'CI-Vector read from disk for root = ',JROOT
  call WRTMAT(VEC1,1,NUM_ELE,1,NUM_ELE)
# endif
  call CSDTVC(VEC2,VEC1,2,DTOC,SDREO,ISSM,0)
  call DDAFILE(JOBIPH,1,VEC2,NCONF,JOBDISK)
  call IFRMDS(I_DUMMY,1,LBLK,LUDIA)
end do
IDISK(LUDIA) = 0

#ifdef _DEBUGPRINT_
do JROOT=1,NROOT
  call WRTVCD(VEC1,LUDIA,0,LBLK)
end do
#endif

! clean up time : copy 1-e integrals back in place
INT1(:) = LH1SAVE(:)

call mma_deallocate(VEC1)
call mma_deallocate(VEC2)
call mma_deallocate(VEC3)
call mma_deallocate(VEC4)
call mma_deallocate(LCMOMO)
call mma_deallocate(LH1SAVE)

end subroutine TRACI_MASTER
