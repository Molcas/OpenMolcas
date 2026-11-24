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

subroutine ExplH_SA(iSym,DIAG)
!***********************************************************************
!                                                                      *
!     Compute and diagonalize the explicit Hamiltonian of the SA CI    *
!     preconditioner in a subspace of at most nexp_max CSFs, as        *
!     EXPLH2 does for the Davidson procedure of RASSCF. The subspace   *
!     and its Hamiltonian come from the CI library                     *
!     (MCLR_CI_PSpace -> CI_Lib_PSpace, the same selection rule as     *
!     RASSCF).                                                         *
!                                                                      *
!     calling arguments:                                               *
!     iSym    : integer                                                *
!               symmetry of the CI space                               *
!     DIAG    : array of Real                                          *
!               diagonal Hamiltonian, core energy included             *
!                                                                      *
!     output (MCLR_Data):                                              *
!     nSel    : number of CSFs in the subspace (0: diagonal only)      *
!     iSel    : the CSFs of the subspace                               *
!     ExplE   : eigenvalues of the explicit Hamiltonian, core energy   *
!               included                                               *
!     ExplV   : eigenvectors of the explicit Hamiltonian               *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use MCLR_procedures, only: MCLR_CI_PSpace
use MCLR_Data, only: ExplE, ExplV, FIMO, INT2, iSel, nexp_max, nSel
use input_mclr, only: nCSF, PotNuc, rIn_Ene
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym
real(kind=wp), intent(in) :: DIAG(nCSF(iSym))
integer(kind=iwp) :: I, MXXSEL
real(kind=wp) :: ECORE
integer(kind=iwp), allocatable :: IdxP(:)
real(kind=wp), allocatable :: EXHAM(:)

call ExplH_SA_Free()

MXXSEL = min(nexp_max,nCSF(iSym))
if (MXXSEL <= 1) return

! CONSTRUCT THE EXPLICIT HAMILTONIAN
! (selection of the subspace and <I|H|J> through the CI interface; the
! core energy is not in the block, but it is in the diagonal and in the
! energies of the states)

call mma_allocate(IdxP,MXXSEL,Label='IdxP')
call mma_allocate(EXHAM,nTri_Elem(MXXSEL),Label='EXHAM')
call MCLR_CI_PSpace(iSym,size(FIMO),FIMO,size(INT2),INT2,nCSF(iSym),DIAG,MXXSEL,nSel,IdxP,size(EXHAM),EXHAM)
! (no explicit subspace with some CI libraries: the diagonal alone)
if (nSel <= 1) then
  nSel = 0
  call mma_deallocate(IdxP)
  call mma_deallocate(EXHAM)
  return
end if
ECORE = rIn_Ene+PotNuc
do I=1,nSel
  EXHAM(nTri_Elem(I)) = EXHAM(nTri_Elem(I))+ECORE
end do

! DIAGONALIZE THE EXPLICIT HAMILTONIAN.

call mma_allocate(iSel,nSel,Label='iSel')
call mma_allocate(ExplE,nSel,Label='ExplE')
call mma_allocate(ExplV,nSel,nSel,Label='ExplV')
iSel(:) = IdxP(1:nSel)
ExplV(:,:) = Zero
do I=1,nSel
  ExplV(I,I) = One
end do
call NIdiag_new(EXHAM,ExplV,nSel,nSel)
call JACORD(EXHAM,ExplV,nSel,nSel)
do I=1,nSel
  ExplE(I) = EXHAM(nTri_Elem(I))
end do
call mma_deallocate(IdxP)
call mma_deallocate(EXHAM)

end subroutine ExplH_SA
