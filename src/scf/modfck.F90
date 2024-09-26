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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine ModFck(Fock,Ovlp,nFO,CMO,nCMO,mynOcc)
!***********************************************************************
!                                                                      *
!     purpose: Modify Fock matrix taking into account frozen orbitals  *
!              F <- [1 - S*D(f)]*F  (symmetrized).                     *
!                                                                      *
!     input:                                                           *
!       Fock    : Fock matrix of length nFO                            *
!       Ovlp    : overlap matrix of length nFO                         *
!                                                                      *
!     output:                                                          *
!       Fock    : modified Fock matrix                                 *
!                                                                      *
!***********************************************************************

use InfSCF, only: MaxBas, nBas, nBO, nBT, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nFO, nCMO, mynOcc(*)
real(kind=wp) :: Fock(nFO), Ovlp(nFO), CMO(nCMO)
integer(kind=iwp) :: ij, iSym
real(kind=wp), allocatable :: Aux1(:), DFro(:), DFSq(:), OvSq(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for density contribution
call mma_allocate(DFro,nBT,Label='DFro')

! Allocate memory for squared density (then for squared Fock matrix)
call mma_allocate(DFSq,MaxBas**2,Label='DFSq')

! Allocate memory for squared overlap and later for 1 - S*D(f)
call mma_allocate(OvSq,MaxBas**2,Label='OvSq')

! Allocate memory for S*D(f)
call mma_allocate(Aux1,MaxBas**2,Label='Aux1')

! Generate contribution to the density mat. from frozen orbitals
call DFroz(DFro,nBT,CMO,nBO,mynOcc)

ij = 1
do iSym=1,nSym

  if (nBas(iSym) > 0) then
    ! Square density and overlap
    call DSq(DFro(ij),DFSq,1,nBas(iSym),nBas(iSym))
    call Square(Ovlp(ij),OvSq,1,nBas(iSym),nBas(iSym))

    ! Compute S*D(f)
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym), &
                One,OvSq,nBas(iSym), &
                DFSq,nBas(iSym), &
                Zero,Aux1,nBas(iSym))

    ! Put unit matrix to memory ascribed for overlap and compute 1 - S*D(f)
    call dcopy_(nBas(iSym)**2,[Zero],0,OvSq,1)
    call dcopy_(nBas(iSym),[One],0,OvSq,nBas(iSym)+1)
    call daxpy_(nBas(iSym)**2,-One,Aux1,1,OvSq,1)

    ! Square Fock matrix and store in memory ascribed for density
    call Square(Fock(ij),DFSq,1,nBas(iSym),nBas(iSym))

    ! Compute F <- [1 - S*D(f)]*F
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym), &
                One,OvSq,nBas(iSym), &
                DFSq,nBas(iSym), &
                Zero,Aux1,nBas(iSym))

    ! Symmetrize modified Fock matrix
    call Sym(Aux1,Fock(ij),nBas(iSym))
  end if

  ! Update pointers
  ij = ij+nBas(iSym)*(nBas(iSym)+1)/2

end do

! Deallocate memory
call mma_deallocate(Aux1)
call mma_deallocate(OvSq)
call mma_deallocate(DFSq)
call mma_deallocate(DFro)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine ModFck
