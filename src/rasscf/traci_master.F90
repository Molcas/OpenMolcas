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

use ci_interfaces, only: Mk_Tra_Psi
use rasscf_global, only: NROOTS
use general_data, only: nSym, NASH
use rasscf_global, only: NAC, nCSF=>nConf
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
integer(kind=iwp), parameter :: iState=1

call mma_allocate(CI_Vector,nCSF,Label='CI_Vector')

! The input transformation matrix, CMOMO is over the full basis but it
! has non-zero values just in the block-diagonal sections. Here we
! remove those blocks.

! MO-MO transformation matrix :
call mma_allocate(TransMat,NAC**2,Label='TransMat')
TransMat(:) = Zero
IOFF = 1
IADR = 1
ICOL = 1
do ISM=1,nSym
  if (NASH(ISM) == 0) Cycle
  IROW = ICOL
  do I=1,NASH(ISM)
    IADR = (ICOL-1)*NAC+IROW
    do J=1,NASH(ISM)
      TransMat(IOFF+NASH(ISM)*(J-1)+I-1) = CMOMO(IADR+J-1)
    end do
    ICOL = ICOL+1
  end do
  IOFF = IOFF+NASH(ISM)**2
end do

! Memory allocation for a matrix T
call mma_allocate(LT,NAC**2,Label='LT')

! Obtain T matrix used for transformation, for each symmetry separately
IOFF = 0
do ISM=1,nSym
  if (ISM == 1) then
    IOFF = 1
  else
    IOFF = IOFF+NASH(ISM-1)**2
  end if
  if (NASH(ISM) > 0) call PAMTMT(TransMat(IOFF),LT(IOFF),NASH(ISM))
end do

JDISK = JOBDISK
do JROOT=1,NROOTS
  ! Read CI vector in SG format
  call DDAFILE(JOBIPH,2,CI_Vector,nCSF,JDISK)

!CI_Vector(:) = [ &
!   -0.9246364863243709_wp, &
!    0.3808508476526229_wp, &
!    0.0000000000000042_wp  &
!]

!LT(:) = [ &
!    1.0950887694440439_wp, &
!   -0.4463400194498248_wp, &
!    0.4075834141522822_wp, &
!    0.9131679804438881_wp  &
!]

  Call Mk_Tra_Psi(iState,LT,NAC**2,CI_Vector,nCSF)

  Call VecPhase(CI_Vector,nCSF)
  call DDAFILE(JOBIPH,1,CI_Vector,nCSF,JOBDISK)
end do

call mma_deallocate(CI_Vector)
call mma_deallocate(LT)
call mma_deallocate(TransMat)

end subroutine TRACI_MASTER
