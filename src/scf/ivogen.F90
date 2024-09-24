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

subroutine IvoGen(OneHam,none,CMO,nCMO,EOrb,nEOrb,mynOcc)
!***********************************************************************
!                                                                      *
!     purpose: Generate improved virtual orbitals by diagonalization   *
!              of the one-electron hamiltonian in the subspace spanned *
!              by the virtual orbitals                                 *
!                                                                      *
!     input:                                                           *
!       OneHam  : one-electron hamiltonian of length nOne              *
!       CMO     : molecular orbital coefficients of length nCMO        *
!                                                                      *
!     output:                                                          *
!       CMO     : molecular orbital coefficients with virtual orbitals *
!                 modified                                             *
!       EOrb    : orbital energies (set to zero for virtual orbitals)  *
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!***********************************************************************

use InfSCF, only: MaxBas, MaxBOO, MaxOrO, nSym, nBas, nOrb
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer none, nCMO, nEOrb
real*8 OneHam(none), CMO(nCMO), EOrb(nEOrb)
integer mynOcc(*)
real*8 Dummy
integer iCMO, iDum, iEOr, iErr, ij, iSym, nFound, nOrbi
real*8, dimension(:), allocatable :: FckS, FckH, FckT, Scratch

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for squared modified Fock matrix
call mma_allocate(FckS,MaxBas**2,Label='FckS')

! Allocate memory for half transformed Fock matrix
call mma_allocate(FckH,MaxBOO,Label='FckH')

! Allocate memory for transformed Fock matrix
call mma_allocate(FckT,MaxOrO*(MaxOrO+1)/2,Label='FckT')

ij = 1
iCMO = 1
iEOr = 1
do iSym=1,nSym
  nOrbi = nOrb(iSym)-mynOcc(iSym)

  ! If nOrbi == 0 - no virtual orbitals; iCMO and iEOr must be
  ! updated anyway (occupied orbitals may exist)
  iCMO = iCMO+nBas(iSym)*mynOcc(iSym)
  iEOr = iEOr+mynOcc(iSym)

  if (nOrbi > 0) then

    ! Transform OneHam to space spanned by virtual orbitals
    call Square(OneHam(ij),FckS,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nOrbi,nBas(iSym), &
                One,FckS,nBas(iSym), &
                CMO(iCMO),nBas(iSym), &
                Zero,FckH,nBas(iSym))
    call DGEMM_Tri('T','N',nOrbi,nOrbi,nBas(iSym), &
                   One,CMO(iCMO),nBas(iSym), &
                   FckH,nBas(iSym), &
                   Zero,FckT,nOrbi)

    ! Diagonalize OneHam within virtual space and form orbital energies
    call mma_allocate(Scratch,nOrbi**2,Label='Scratch')
    Dummy = Zero
    iDum = 0
    call Diag_Driver('V','A','L',nOrbi,FckT,Scratch,nOrbi,Dummy,Dummy,iDum,iDum,EOrb(iEOr),CMO(iCMO),nBas(iSym),0,-1,'J',nFound, &
                     iErr)
    call mma_deallocate(Scratch)

    ! Orbital energies are now meaningless; set them to zero
    EOrb(iEOr:iEOr+nOrbi-1) = Zero

  end if

  ! Update pointers
  iCMO = iCMO+nOrbi*nBas(iSym)
  iEOr = iEOr+nOrbi
  ij = ij+nBas(iSym)*(nBas(iSym)+1)/2

end do

! Deallocate memory
call mma_deallocate(FckS)
call mma_deallocate(FckH)
call mma_deallocate(FckT)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine IvoGen
