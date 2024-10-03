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

subroutine OvlDel(Ovlp,nOvlp,TrMat,nTrMat)
!***********************************************************************
!                                                                      *
!     purpose: Remove near linear dependencies from basis set          *
!                                                                      *
!     input:                                                           *
!       Ovlp    : overlap in AO basis of length nOvlp                  *
!       TrMat   : unit matrix or matrix transforming from AO's to      *
!                 cartesian functions (dependently on input) of        *
!                 length nTrMat                                        *
!                                                                      *
!     output:                                                          *
!       TrMat   : input matrix modified such, that near linear de-     *
!                 pendencies (if exist) are removed                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: DelThr, MaxBas, MaxBOF, MaxOrF, MiniDn, nBas, nDel, nFro, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOvlp, nTrMat
real(kind=wp), intent(in) :: Ovlp(nOvlp)
real(kind=wp), intent(inout) :: TrMat(nTrMat)
integer(kind=iwp) :: iDum, iErr, iiBT, ij, Ind, iNew, iOld, iOrb, iSym, nFound, nOF, nOrbi
real(kind=wp) :: Dummy
real(kind=wp), allocatable :: EVal(:), EVec(:), NewB(:), OvlH(:), OvlS(:), OvlT(:), Scratch(:)

! Allocate memory for transformed overlap matrix
call mma_allocate(OvlT,nTri_Elem(MaxOrF),Label='OvlT')

! Allocate memory for half-transformed overlap matrix
call mma_allocate(OvlH,MaxBOF,Label='OvlH')

! Allocate memory for squared overlap matrix
call mma_allocate(OvlS,MaxBas**2,Label='OvlS')

! Allocate memory for eigenvectors of overlap
call mma_allocate(EVec,MaxOrF**2,Label='EVec')

! Allocate memory for eigenvalues of overlap
call mma_allocate(EVal,MaxOrF,Label='EVal')

! Allocate memory for 'basis' that diagonalizes overlap
call mma_allocate(NewB,MaxBOF,Label='NewB')

ij = 1
iOld = 1
iNew = 1
do iSym=1,nSym

  iiBT = nTri_Elem(nBas(iSym))
  nOF = nOrb(iSym)-nFro(iSym)

  ! Copy frozen vectors to the right position
  if (nFro(iSym)*nBas(iSym) > 0) then
    !write(u6,'(a,i2)') 'Copying symmetry block',iSym
    !write(u6,'(i8,a,i8)') iOld,' ->',iNew
    !write(u6,'(a,i8)') 'nTrMat =',nTrMat
    TrMat(iNew:iNew+nFro(iSym)*nBas(iSym)-1) = TrMat(iOld:iOld+nFro(iSym)*nBas(iSym)-1)
  !else
  !  write(u6,'(a,i2)') 'No copying of symmetry block',iSym
  end if

  iOld = iOld+nFro(iSym)*nBas(iSym)
  iNew = iNew+nFro(iSym)*nBas(iSym)

  if (nOF > 0) then

    ! Square overlap and transform to basis given by TrMat
    call Square(Ovlp(ij),OvlS,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nOF,nBas(iSym), &
                One,OvlS,nBas(iSym), &
                TrMat(iOld),nBas(iSym), &
                Zero,OvlH,nBas(iSym))
    call DGEMM_Tri('T','N',nOF,nOF,nBas(iSym), &
                   One,TrMat(iOld),nBas(iSym), &
                   OvlH,nBas(iSym), &
                   Zero,OvlT,nOF)

    ! Diagonalize overlap and form eigenvalues vector
    call mma_allocate(Scratch,nOF**2,Label='Scrtach')
    Dummy = Zero
    iDum = 0
    call Diag_Driver('V','A','L',nOF,OvlT,Scratch,nOF,Dummy,Dummy,iDum,iDum,EVal,EVec,nOF,1,0,'J',nFound,iErr)
    call mma_deallocate(Scratch)
    !?? Ovlt(1:nTri_Elem(nOF)) = Zero
    !?? iDiag = 0
    !?? do i=1,nOF
    !??   OvlT(i+iDiag) = EVal(i)
    !??   iDiag = iDiag+i
    !?? end do

    ! Transform to basis that diagonalizes overlap
    call DGEMM_('N','N',nBas(iSym),nOF,nOF, &
                One,TrMat(iOld),nBas(iSym), &
                EVec,nOF, &
                Zero,NewB,nBas(iSym))

    ! Remove linear dependencies
    nOrbi = nFro(iSym)
    ind = 1
    do iOrb=1,nOF
      if (EVal(iOrb) > DelThr) then
        if (EVal(iOrb) < 1.0e-5_wp) MiniDn = .false.
        TrMat(iNew:iNew+nBas(iSym)-1) = NewB(ind:ind+nBas(iSym)-1)
        iNew = iNew+nBas(iSym)
        nOrbi = nOrbi+1
      end if
      ind = ind+nBas(iSym)
    end do

    iOld = iOld+nOF*nBas(iSym)
    nDel(iSym) = nOrb(iSym)-nOrbi
    nOrb(iSym) = nOrbi

  end if

  ij = ij+iiBT

end do
call Put_iArray('nDel',nDel,nSym)

! Deallocate memory
call mma_deallocate(NewB)
call mma_deallocate(EVal)
call mma_deallocate(EVec)
call mma_deallocate(OvlS)
call mma_deallocate(OvlH)
call mma_deallocate(OvlT)

end subroutine OvlDel
