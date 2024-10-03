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

subroutine DCore(OneHam,nOH,CMO,TrMat,nCMO,E_Or,nEOr,mynOcc,Ovrlp)
!***********************************************************************
!                                                                      *
!     purpose: Diagonalize core hamiltonian to get starting orbitals.  *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: MaxBas, MaxBOF, MaxORF, nBas, nBO, nBT, nFro, nnFr, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOH, nCMO, nEOr, mynOcc(*)
real(kind=wp), intent(in) :: OneHam(nOH), TrMat(nCMO), Ovrlp(nOH)
real(kind=wp), intent(out) :: CMO(nCMO), E_Or(nEOr)
integer(kind=iwp) :: iCMO, iDum, iE_Or, iErr, iiBT, ij, iSym, nC, nFound, nOF
real(kind=wp) :: Dummy
integer(kind=iwp), allocatable :: Fermi(:)
real(kind=wp), allocatable :: EiVe(:), OHHl(:), OHSq(:), OHTr(:), OMod(:), Scratch(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for modified one-electron hamoltonian
call mma_allocate(OMod,nBT,Label='OMod')

! Allocate memory for squared one-el. Hamiltonian
call mma_allocate(OHSq,MaxBas**2,Label='OHSq')

! Allocate memory for half-transformed one-el. Hamiltonian
call mma_allocate(OHHl,MaxBOF,Label='OHHl')

! Allocate memory for transformed one-electron Hamiltonian
call mma_allocate(OHTr,nTri_Elem(MaxOrF),Label='OHTr')

! Allocate memory for eigenvectors
call mma_allocate(EiVe,MaxOrF**2,Label='EiVe')

! Allocate memory for info on electronic and muonic basis sets
call mma_allocate(Fermi,nEOr,Label='Fermi')
call Get_iArray('Fermion IDs',Fermi,nEOr)

! Modify one-electron hamiltonian
OMod(:) = OneHam(1:nBT)
if (nnFr > 0) call ModFck(OMod,Ovrlp,nBT,TrMat,nBO,mynOcc)

! Diagonalize core in non-frozen molecular basis
ij = 1
iCMO = 1
iE_Or = 1
do iSym=1,nSym

  iiBT = nTri_Elem(nBas(iSym))
  nOF = nOrb(iSym)-nFro(iSym)

  ! Copy frozen vectors to CMO array
  nC = nFro(iSym)*nBas(iSym)
  if (nC > 0) CMO(iCMO:iCMO+nC-1) = TrMat(iCMO:iCMO+nC-1)

  iCMO = iCMO+nBas(iSym)*nFro(iSym)
  iE_Or = iE_Or+nFro(iSym)

  if (nOF > 0) then

    ! Square one-el. Ham. and transform to orthonormal basis.
    !call Square(OMod(ij),nBas(iSym),OHSq,nBas(iSym))
    call Square(OMod(ij),OHSq,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nOF,nBas(iSym), &
                One,OHSq,nBas(iSym), &
                TrMat(iCMO),nBas(iSym), &
                Zero,OHHl,nBas(iSym))
    call DGEMM_Tri('T','N',nOF,nOF,nBas(iSym), &
                   One,TrMat(iCMO),nBas(iSym), &
                   OHHl,nBas(iSym), &
                   Zero,OHTr,nOF)

    ! Put a unit matrix into the eigenvector matrix
    call unitmat(EiVe,nOF)

    ! Add small random number to the one-electron Hamiltonian
    ! Not done here anymore, done in routine scram called by sorb!

    !if (Scrmbl) then
    !  ind = 1
    !  do i=1,nOF
    !    do j=1,i-1
    !      OHTr(ind) = OHTr(ind)+0.05_wp*Random_Molcas(iseed)
    !      ind = ind+1
    !    end do
    !    ind = ind+1
    !  end do
    !end if

    ! Diagonalize and form orbital energies
    call mma_allocate(Scratch,nOF**2,Label='Scratch')
    Dummy = Zero
    iDum = 0
    call Diag_Driver('V','A','L',nOF,OHTr,Scratch,nOF,Dummy,Dummy,iDum,iDum,E_Or(iE_Or),EiVe,nOF,1,-1,'J',nFound,iErr)
    call mma_deallocate(Scratch)

    ! Transform to AO basis
    call DGEMM_('N','N',nBas(iSym),nOF,nOF, &
                One,TrMat(iCMO),nBas(iSym), &
                EiVe,nOF, &
                Zero,CMO(iCMO),nBas(iSym))

  end if

  ! Update pointers
  iCMO = iCMO+nOF*nBas(iSym)
  iE_Or = iE_Or+nOF
  ij = ij+iiBT

end do

! Deallocate memory
call mma_deallocate(Fermi)
call mma_deallocate(EiVe)
call mma_deallocate(OHTr)
call mma_deallocate(OHHl)
call mma_deallocate(OHSq)
call mma_deallocate(OMod)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine DCore
