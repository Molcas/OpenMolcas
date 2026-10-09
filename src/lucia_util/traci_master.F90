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

use lucia_data, only: NROOT, NTOOBS
use general_data, only: nCSF=>nConf, nLev, STSYM, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(inout) :: JOBDISK
integer(kind=iwp), intent(in) :: JOBIPH
real(kind=wp), intent(in) :: CMOMO(*)
integer(kind=iwp) :: I, IADR, ICOL, IOFF, IROW, ISM, J, JDISK, JROOT
real(kind=wp), allocatable :: TransMat(:), CI_Vector(:), LT(:)

call mma_allocate(CI_Vector,nCSF,Label='CI_Vector')

! The input transformation matrix, CMOMO is over the full basis but it
! has non-zero values just in the block-diagonal sections. Here we
! remove those blocks.

! MO-MO transformation matrix :
call mma_allocate(TransMat,nLev**2,Label='TransMat')
TransMat(:) = Zero
IOFF = 1
IADR = 1
ICOL = 1
do ISM=1,nSym
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
do ISM=1,nSym
  if (ISM == 1) then
    IOFF = 1
  else
    IOFF = IOFF+NTOOBS(ISM-1)**2
  end if
  if (NTOOBS(ISM) > 0) call PAMTMT(TransMat(IOFF),LT(IOFF),NTOOBS(ISM))
end do

JDISK = JOBDISK
do JROOT=1,NROOT
  ! Read CI vector in SG format
  call DDAFILE(JOBIPH,2,CI_Vector,nCSF,JDISK)

  call TRACI_LUCIA(LT,nLev,CI_Vector,nCSF,STSYM)

  call DDAFILE(JOBIPH,1,CI_Vector,nCSF,JOBDISK)
end do

call mma_deallocate(CI_Vector)
call mma_deallocate(LT)
call mma_deallocate(TransMat)

end subroutine TRACI_MASTER
