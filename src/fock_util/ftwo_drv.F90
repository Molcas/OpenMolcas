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

use Data_Structures, only: DSBA_Type, Allocate_DSBA, Deallocate_DSBA

implicit real*8(A-H,O-Z)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "real.fh"
integer nBas(8), nAsh(8), nSkipX(8)
real*8 CMO(*), D1A(*), DI(*), FA(*)
logical DoCholesky
type(DSBA_Type) :: WFSQ
real*8, allocatable :: Temp(:)
#include "choras.fh"
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine CHORAS_DRV(nSym,nBas,nOcc,W_DSQ,W_DLT,W_FLT,ExFac,FSQ,W_CMO)
    use Data_Structures, only: DSBA_Type
    integer nSym, nBas(8)
    integer, target :: nOcc(nSym)
    real*8 W_FLT(*), W_DSQ(*), W_DLT(*)
    real*8 ExFac
    type(DSBA_Type) FSQ
    real*8 W_CMO(*)
  end subroutine CHORAS_DRV
end interface
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
