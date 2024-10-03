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

!#define _DEBUGPRINT_
subroutine FixOrb(Ovrlp,CMO,TrMat,nCMO)
!***********************************************************************
!                                                                      *
!     purpose: Diagonalize core hamiltonian to get starting orbitals.  *
!                                                                      *
!     input:                                                           *
!       OneHam  : one-electron hamiltonian of length nOH               *
!       TrMat   : matrix transforming to orthonotmal basis of          *
!                 length nCMO                                          *
!                                                                      *
!     output:                                                          *
!       CMO     : starting vectors of length nCMO                      *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: MaxBas, nBas, nFro, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCMO
real(kind=wp), intent(in) :: Ovrlp(nCMO), TrMat(nCMO)
real(kind=wp), intent(inout) :: CMO(nCMO)
integer(kind=iwp) :: iCMO, iS, iSym, iTrM, nBF, nOF
real(kind=wp), allocatable :: S(:), TT(:), TTS(:), CMO0(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory

call mma_allocate(S,MaxBas**2,Label='S')
call mma_allocate(TT,MaxBas**2,Label='TT')
call mma_allocate(TTS,MaxBas**2,Label='TTS')
call mma_allocate(CMO0,MaxBas**2,Label='CMO0')
!                                                                      *
!***********************************************************************
!                                                                      *
! Observe that TrGen has been called and that the dimensions might
! have been change. Note also that at this point work in the full
! basis.

! Loop over symmetry blocks

iS = 1
iCMO = 1
iTrM = 1
do iSym=1,nSym
# ifdef _DEBUGPRINT_
  call RecPrt('FixOrb: CMO(in)',' ',CMO(iCMO),nBas(iSym),nBas(iSym))
  call RecPrt('FixOrb: TrMat',' ',TrMat(iCMO),nBas(iSym),nOrb(iSym))
# endif

  nOF = nOrb(iSym)-nFro(iSym)
  nBF = nBas(iSym)-nFro(iSym)

  ! Skip the frozen orbitals

  iCMO = iCMO+nBas(iSym)*nFro(iSym)
  iTrM = iTrM+nBas(iSym)*nFro(iSym)

  if (nBF > 0) then

    ! Project to the reduced basis

    ! Create the project operator
    call DGEMM_('N','T',nBas(iSym),nBas(iSym),nOF, &
                One,TrMat(iTrM),nBas(iSym), &
                TrMat(iTrM),nBas(iSym), &
                Zero,TT,nBas(iSym))
    call Square(Ovrlp(iS),S,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym), &
                One,TT,nBas(iSym), &
                S,nBas(iSym), &
                Zero,TTS,nBas(iSym))
#   ifdef _DEBUGPRINT_
    call RecPrt('FixOrb: TT',' ',TT,nBas(iSym),nBas(iSym))
    call RecPrt('FixOrb: TTS',' ',TTS,nBas(iSym),nBas(iSym))
#   endif

    ! Project
    call DGEMM_('N','N',nBas(iSym),nOF,nBas(iSym), &
                One,TTS,nBas(iSym), &
                CMO(iCMO),nBas(iSym), &
                Zero,CMO0,nBas(iSym))
    CMO(iCMO:iCMO+nBas(iSym)*nOF-1) = CMO0(1:nBas(iSym)*nOF)
#   ifdef _DEBUGPRINT_
    call RecPrt('FixOrb: CMO(out)',' ',CMO(iCMO),nBas(iSym),nBas(iSym))
#   endif

  end if

  ! Update pointers
  iCMO = iCMO+nBF*nBas(iSym)
  iTrM = iTrM+nOF*nBas(iSym)
  iS = iS+nTri_Elem(nBas(iSym))

end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

call mma_deallocate(CMO0)
call mma_deallocate(TTS)
call mma_deallocate(TT)
call mma_deallocate(S)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FixOrb
