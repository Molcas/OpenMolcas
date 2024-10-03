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

subroutine Ortho(AMat,nAMat,Ovlp,nOvlp)
!***********************************************************************
!                                                                      *
!     purpose: Transform to orthonormal basis (one symmetry block      *
!              at a time)                                              *
!     input:                                                           *
!       AMat    : transformation matrix to non-orthonormal basis of    *
!                 length nAMat                                         *
!       Ovlp    : overlap in AO basis of length nOvlp                  *
!                                                                      *
!     output:                                                          *
!       AMat    : transformation matrix to orthonormal basis           *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: MaxBas, MaxBxO, MaxOrb, nSym, nOrb, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAMat, nOvlp
real(kind=wp), intent(inout) :: AMat(nAMat)
real(kind=wp), intent(in) :: Ovlp(nOvlp)
integer(kind=iwp) :: iiBO, iiBT, ij, im, iSym
real(kind=wp), allocatable :: OvlT(:), OvlH(:), OvlS(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for transformed overlap matrix
call mma_allocate(OvlT,MaxOrb**2,Label='OvlT')

! Allocate memory for half-transformed overlap matrix
call mma_allocate(OvlH,MaxBxO,Label='OvlH')

! Allocate memory for squared overlap matrix
call mma_allocate(OvlS,MaxBas**2,Label='OvlS')

ij = 1
im = 1
do iSym=1,nSym

  iiBO = nBas(iSym)*nOrb(iSym)
  iiBT = nTri_Elem(nBas(iSym))

  if (nOrb(iSym) > 0) then

    ! Square overlap and transform to the basis given by AMat
    !call Square(Ovlp(ij),nBas(iSym),OvlS,nBas(iSym))
    call Square(Ovlp(ij),OvlS,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym), &
                One,OvlS,nBas(iSym), &
                AMat(im),nBas(iSym), &
                Zero,OvlH,nBas(iSym))
    call DGEMM_('T','N',nOrb(iSym),nOrb(iSym),nBas(iSym), &
                One,AMat(im),nBas(iSym), &
                OvlH,nBas(iSym), &
                Zero,OvlT,nOrb(iSym))

    ! Orthogonalize (Gram-Schmidt)
    call Orthox(OvlT,AMat(im),nOrb(iSym),nBas(iSym))

  end if

  ! Update pointers
  im = im+iiBO
  ij = ij+iiBT

end do

! Deallocate memory
call mma_deallocate(OvlT)
call mma_deallocate(OvlH)
call mma_deallocate(OvlS)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine Ortho
