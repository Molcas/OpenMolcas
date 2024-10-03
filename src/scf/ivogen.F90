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

subroutine IvoGen(OneHam,nOneH,CMO,nCMO,EOrb,nEOrb,mynOcc)
!***********************************************************************
!                                                                      *
!     purpose: Generate improved virtual orbitals by diagonalization   *
!              of the one-electron hamiltonian in the subspace spanned *
!              by the virtual orbitals                                 *
!                                                                      *
!     input:                                                           *
!       OneHam  : one-electron hamiltonian of length nOneH             *
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

use Index_Functions, only: nTri_Elem
use InfSCF, only: MaxBas, MaxBOO, MaxOrO, nBas, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOneH, nCMO, nEOrb, mynOcc(*)
real(kind=wp), intent(in) :: OneHam(nOneH)
real(kind=wp), intent(inout) :: CMO(nCMO)
real(kind=wp), intent(out) :: EOrb(nEOrb)
integer(kind=iwp) :: iCMO, iDum, i_EOr, iErr, ij, iSym, nFound, nOrbi
real(kind=wp) :: Dummy
real(kind=wp), allocatable :: FckH(:), FckS(:), FckT(:), Scratch(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for squared modified Fock matrix
call mma_allocate(FckS,MaxBas**2,Label='FckS')

! Allocate memory for half transformed Fock matrix
call mma_allocate(FckH,MaxBOO,Label='FckH')

! Allocate memory for transformed Fock matrix
call mma_allocate(FckT,nTri_Elem(MaxOrO),Label='FckT')

ij = 1
iCMO = 1
i_EOr = 1
do iSym=1,nSym
  nOrbi = nOrb(iSym)-mynOcc(iSym)

  ! If nOrbi == 0 - no virtual orbitals; iCMO and i_EOr must be
  ! updated anyway (occupied orbitals may exist)
  iCMO = iCMO+nBas(iSym)*mynOcc(iSym)
  i_EOr = i_EOr+mynOcc(iSym)

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
    call Diag_Driver('V','A','L',nOrbi,FckT,Scratch,nOrbi,Dummy,Dummy,iDum,iDum,EOrb(i_EOr),CMO(iCMO),nBas(iSym),0,-1,'J',nFound, &
                     iErr)
    call mma_deallocate(Scratch)

    ! Orbital energies are now meaningless; set them to zero
    EOrb(i_EOr:i_EOr+nOrbi-1) = Zero

  end if

  ! Update pointers
  iCMO = iCMO+nOrbi*nBas(iSym)
  i_EOr = i_EOr+nOrbi
  ij = ij+nTri_Elem(nBas(iSym))

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
