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

subroutine FTwo_Drv(nSym,nBas,nAsh,nSkipX,DI,D1A,FA,nTot1,ExFac,nBMX,CMO)

use Fock_util_interface, only: CHORAS_DRV
use Fock_util_global, only: ALGO
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8), nAsh(8), nSkipX(8), nTot1, nBMX
real(kind=wp), intent(in) :: DI(*), D1A(*), CMO(*), ExFac
real(kind=wp), intent(inout) :: FA(*)
logical(kind=iwp) :: DoCholesky
type(DSBA_Type) :: WFSQ(1)
real(kind=wp), allocatable :: Temp(:)

!                                                                      *
!***********************************************************************
!                                                                      *
call DecideOncholesky(DoCholesky)

if (DoCholesky .and. (ALGO == 2)) then

  ! Building of the Fock matrix directly from Cholesky vectors

  call Allocate_DT(WFSQ(1),nBas,nBas,nSym)
  WFSQ(1)%A0(:) = Zero

  call mma_allocate(Temp,nTot1,Label='nTot1')
  Temp(:) = Zero

  call CHORAS_DRV(nSym,nBas,nAsh,D1A,DI,Temp,ExFac,WFSQ,CMO)

  FA(1:nTot1) = FA(1:nTot1)+Temp(1:nTot1)

  call mma_deallocate(Temp)
  call Deallocate_DT(WFSQ(1))

else

  ! Standard building of the Fock matrix from Two-el integrals

  call FockTwo_Drv(nSym,nBas,nAsh,nSkipX,DI,D1A,FA,nTot1,ExFac,nBMX)

end if

return

end subroutine FTwo_Drv
