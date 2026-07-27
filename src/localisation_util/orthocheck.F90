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
! Copyright (C) 2026, Lila Zapp                                        *
!                                                                      *
!***********************************************************************

subroutine OrthoCheck(CMO,nOrb2Loc,nBasis)
! check that the orbitals are orthonormal (if not -> the trafo matrix was not constructed correctly)

use Localisation_globals, only: Ovlp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOrb2Loc, nBasis
real(kind=wp), intent(in) :: CMO(nBasis,nOrb2Loc)
integer(kind=iwp) :: i
real(kind=wp) :: norm
real(kind=wp), allocatable :: CtS(:,:), CtSC(:,:)
real(kind=wp), parameter :: ThrOrtho = 1.0e-12_wp
real(kind=wp), external :: DDot_

call mma_allocate(CtS,nOrb2Loc,nBasis,Label='CtS')
call mma_allocate(CtSC,nOrb2Loc,nOrb2Loc,Label='CtSC')

! C^T*S
call dgemm_('T','N',nOrb2Loc,nBasis,nBasis, &
            One,CMO,nBasis, &
                Ovlp,nBasis, &
            Zero,CtS,nOrb2Loc)
! (C^T*S)*C
call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nBasis, &
            One,CtS,nOrb2Loc, &
                CMO,nBasis, &
            Zero,CtSC,nOrb2Loc)

#ifdef _DEBUGPRINT_
write(u6,'(/A)') 'Check the orthonormality of the orbitals'
write(u6,*) '========================================'
call RecPrt('C^T*S*C =',' ',CtSC,nOrb2Loc,nOrb2Loc)
#endif

! (C^T*S)*C - I
do i=1,nOrb2Loc
  CtSC(i,i) = CtSC(i,i)-One
end do

norm = sqrt(DDot_(nOrb2Loc**2,CtSC,1,CtSC,1))/real(nOrb2Loc,kind=wp)**2

if (norm > ThrOrtho) then
  write(u6,*)
  write(u6,*) 'Stopping calculation due to:'
  write(u6,*) 'norm of C^T*S*C - 1 =',norm
  write(u6,*) 'max. allowed norm   =',ThrOrtho
  write(u6,*) 'The current MOs do not seem to be sufficiently orthonormal (Bug in the transformation)'
  call Abend()
end if

call mma_Deallocate(CtS)
call mma_Deallocate(CtSC)

end subroutine OrthoCheck
