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

subroutine DCore(OneHam,nOH,CMO,TrMat,nCMO,EOr,nEOr,mynOcc,Ovrlp)
!***********************************************************************
!                                                                      *
!     purpose: Diagonalize core hamiltonian to get starting orbitals.  *
!                                                                      *
!***********************************************************************

use InfSCF, only: MaxBas, MaxBOF, MaxORF, nBO, nBT, nnFr, nSym, nBas, nFro, nOrb
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer nOH, nCMO, nEOr
real*8 OneHam(nOH), CMO(nCMO), TrMat(nCMO), EOr(nEOr), Ovrlp(nOH)
integer mynOcc(*)
real*8, dimension(:), allocatable :: OMod, OHSq, OHHl, OHTr, EiVe, Scratch
integer, dimension(:), allocatable :: Fermi
real*8 :: Dummy
integer :: iCMO, iDum, iEOr, iErr, iiBT, ij, iSym, nFound, nOF
real*8, external :: Random_Molcas

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
call mma_allocate(OHTr,MaxOrF*(MaxOrF+1)/2,Label='OHTr')

! Allocate memory for eigenvectors
call mma_allocate(EiVe,MaxOrF**2,Label='EiVe')

! Allocate memory for info on electronic and muonic basis sets
call mma_allocate(Fermi,nEOr,Label='Fermi')
call Get_iArray('Fermion IDs',Fermi,nEOr)

! Modify one-electron hamiltonian
call dcopy_(nBT,OneHam,1,OMod,1)
if (nnFr > 0) call ModFck(OMod,Ovrlp,nBT,TrMat,nBO,mynOcc)

! Diagonalize core in non-frozen molecular basis
ij = 1
iCMO = 1
iEOr = 1
do iSym=1,nSym

  iiBT = nBas(iSym)*(nBas(iSym)+1)/2
  nOF = nOrb(iSym)-nFro(iSym)

  ! Copy frozen vectors to CMO array
  if (nFro(iSym)*nBas(iSym) > 0) call dcopy_(nFro(iSym)*nBas(iSym),TrMat(iCMO),1,CMO(iCMO),1)

  iCMO = iCMO+nBas(iSym)*nFro(iSym)
  iEOr = iEOr+nFro(iSym)

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
    call dcopy_(nOF*nOF,[Zero],0,EiVe,1)
    call dcopy_(nOF,[One],0,EiVe,nOF+1)

    ! Add small random number to the one-electron Hamiltonian
    ! Not done here anymore, done in routine scram called by sorb!

    !if (Scrmbl) then
    !  ind = 1
    !  do i=1,nOF
    !    do j=1,i-1
    !      OHTr(ind) = OHTr(ind)+0.050d+00*Random_Molcas(iseed)
    !      ind = ind+1
    !    end do
    !    ind = ind+1
    !  end do
    !end if

    ! Diagonalize and form orbital energies
    call mma_allocate(Scratch,nOF**2,Label='Scratch')
    Dummy = Zero
    iDum = 0
    call Diag_Driver('V','A','L',nOF,OHTr,Scratch,nOF,Dummy,Dummy,iDum,iDum,EOr(iEOr),EiVe,nOF,1,-1,'J',nFound,iErr)
    call mma_deallocate(Scratch)

    ! Transform to AO basis
    call DGEMM_('N','N',nBas(iSym),nOF,nOF, &
                One,TrMat(iCMO),nBas(iSym), &
                EiVe,nOF, &
                Zero,CMO(iCMO),nBas(iSym))

  end if

  ! Update pointers
  iCMO = iCMO+nOF*nBas(iSym)
  iEOr = iEOr+nOF
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
