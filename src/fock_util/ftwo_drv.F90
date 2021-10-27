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

use fock_util_interface, only: CHORAS_DRV
use Data_Structures, only: Allocate_DSBA, Deallocate_DSBA, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nSym, nBas(8), nAsh(8), nSkipX(8), nTot1, nBMX
real(kind=wp) :: DI(*), D1A(*), FA(*), CMO(*), ExFac
#include "choras.fh"
logical(kind=iwp) :: DoCholesky
type(DSBA_Type) :: WFSQ
real(kind=wp), allocatable :: Temp(:)
!                                                                      *
!***********************************************************************
!                                                                      *

call DecideOncholesky(DoCholesky)

if (DoCholesky .and. (ALGO == 2)) then

  ! Building of the Fock matrix directly from Cholesky vectors

  call Allocate_DSBA(WFSQ,nBas,nBas,nSym)
  WFSQ%A0(:) = Zero

  call mma_allocate(Temp,nTot1,Label='nTot1')
  Temp(:) = Zero

  call CHORAS_DRV(nSym,nBas,nAsh,D1A,DI,Temp,ExFac,WFSQ,CMO)

  FA(1:nTot1) = FA(1:nTot1)+Temp(1:nTot1)

  call mma_deallocate(Temp)
  call Deallocate_DSBA(WFSQ)

else

  ! Standard building of the Fock matrix from Two-el integrals

  call FockTwo_Drv(nSym,nBas,nAsh,nSkipX,DI,D1A,FA,nTot1,ExFac,nBMX)

end if

return

end subroutine FTwo_Drv
