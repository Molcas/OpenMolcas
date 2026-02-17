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
!               2019, Liviu Ungur                                      *
!***********************************************************************

subroutine IvoGen_rasscf(nSym,nBas,nFro,nIsh,nAsh,nCMO,nEOrb,CMO,EOrb)
!***********************************************************************
!                                                                      *
!     purpose: Generate improved virtual orbitals by diagonalization   *
!              of the one-electron hamiltonian in the subspace spanned *
!              by the virtual orbitals                                 *
!                                                                      *
!     input:                                                           *
!       CMO     : molecular orbital coefficients of length nCMO=NTOT2  *
!       EOrb    : on input contain orbital energies as obtained in the *
!                 NATORB_RASSCF, energies  of the virtual orbitals     *
!                 are non-zero, length nEOrb=NTOT                      *
!                                                                      *
!     output:                                                          *
!       CMO     : molecular orbital coefficients with virtual orbitals *
!                 modified                                             *
!       EOrb    : orbital energies (set to zero for virtual orbitals)  *
!                                                                      *
!     called from: OutCtl, only when IVO is used in the input          *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Original MOLCAS/scf/ivogen.f written by:                         *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!                                                                      *
!     written by:                                                      *
!     L. Ungur  using the MOLCAS/scf/ivogen.f as start                 *
!     National University of Singapore,  Feb 2019                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nSym, nCMO, nEOrb, nBas(nSym), nFro(nSym), nIsh(nSym), nAsh(nSym)
real(kind=wp) :: CMO(nCMO), EOrb(nEOrb)
integer(kind=iwp) :: i_EOr, iCMO, iComp, iDum, iErr, ij, iOpt, iRc, iSyLbl, iSym, MaxBas, MaxBOO, MaxOrO, nBT, nFound, nOcc(nSym), &
                     nOrbi
real(kind=wp) :: Dummy
character(len=8) :: Label
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iBas, iOff
#endif
real(kind=wp), allocatable :: FckH(:), FckS(:), FckT(:), OneHam(:), Scratch(:)
#include "warnings.h"

nBT = 0
MaxBas = 0
MaxOrO = 0
MaxBOO = 0
nOcc(1:nSym) = 0
do iSym=1,nSym
  nOcc(iSym) = nFro(iSym)+nIsh(iSym)+nAsh(iSym)
  nBT = nBT+nTri_Elem(nBas(iSym))
  MaxBas = max(MaxBas,nBas(iSym))
  MaxOrO = max(MaxOrO,nBas(iSym)-nOcc(iSym))
  MaxBOO = max(MaxBOO,nBas(iSym)*(nBas(iSym)-nOcc(iSym)))
end do
#ifdef _DEBUGPRINT_
iCMO = 1
do iSym=1,nSym
  iCMO = iCMO+nBas(iSym)**2
  call RecPrt('IvoGen: CMO(in)',' ',CMO(iCMO),nBas(iSym),nBas(iSym))
end do
#endif
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
! Allocate memory for the core Hamiltonian
call mma_allocate(OneHam,nBT,Label='OneHam')
call dcopy_(nBT,[Zero],0,OneHam,1)
! Load bare nuclei Hamiltonian

iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'OneHam'
call RdOne(iRc,iOpt,Label,iComp,OneHam,iSyLbl)
if (iRc /= 0) then
  write(u6,*) ' RASSCF tried to construct compact virtual orbitals'
  write(u6,*) ' by diagonalization of core Hamiltonian, but ran   '
  write(u6,*) ' into a severe error: Failed to read the           '
  write(u6,*) ' Hamiltonian from the ONEINT file. Something may be'
  write(u6,*) ' wrong with the file.'
  call Quit(_RC_IO_ERROR_READ_)
end if
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' OneHam in AO basis in RASSCF'
write(u6,*) ' ---------------------'
write(u6,*)
iOff = 0
do iSym=1,nSym
  iBas = nBas(iSym)
  call TriPrt('IvoGen: OneHam:','(5G17.11)',OneHam(1+iOff),iBas)
  iOff = iOff+nTri_Elem(iBas)
end do
#endif

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
  nOrbi = nBas(iSym)-nOcc(iSym)

  ! If nOrbi == 0 - no virtual orbitals; iCMO and i_EOr must be
  ! updated anyway (occupied orbitals may exist)
  iCMO = iCMO+nBas(iSym)*nOcc(iSym)
  i_EOr = i_EOr+nOcc(iSym)

  if (nOrbi > 0) then

    ! Transform OneHam to space spanned by virtual orbitals
    call Square(OneHam(ij),FckS,1,nBas(iSym),nBas(iSym))
    ! multiply FckH = OneHam x CMO
    call DGEMM_('N','N',nBas(iSym),nOrbi,nBas(iSym),One,FckS,nBas(iSym),CMO(iCMO),nBas(iSym),Zero,FckH,nBas(iSym))
    ! multiply FckT =  CMO x FckH
    call DGEMM_Tri('T','N',nOrbi,nOrbi,nBas(iSym),One,CMO(iCMO),nBas(iSym),FckH,nBas(iSym),Zero,FckT,nOrbi)

    ! Diagonalize OneHam within virtual space and form orbital energies
    call mma_allocate(Scratch,nOrbi**2,Label='Scratch')
    Dummy = Zero
    iDum = 0
    call Diag_Driver('V','A','L',nOrbi,FckT,Scratch,nOrbi,Dummy,Dummy,iDum,iDum,EOrb(i_EOr),CMO(iCMO),nBas(iSym),0,-1,'J',nFound, &
                     iErr)
    call mma_deallocate(Scratch)

    ! Orbital energies are now meaningless; set them to zero
    call dcopy_(nOrbi,[Zero],0,EOrb(i_EOr),1)

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
call mma_deallocate(OneHam)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine IvoGen_rasscf
