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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine DeDe_SCF(Dens,TwoHam,nDens,mDens)

use Index_Functions, only: nTri_Elem
use k2_arrays, only: DeDe, Dq, Fq, ipD00, ipDeDe, ipDijS, ipOffD, MaxDe, MxDij, nDeDe, pDq, pFq
use Basis_Info, only: nBas
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDens
real(kind=wp), target, intent(inout) :: Dens(nDens)
real(kind=wp), target, intent(in) :: TwoHam(nDens)
integer(kind=iwp), intent(out) :: mDens
integer(kind=iwp) :: i, ij, mDeDe, mIndij, nDeDe_Tot, nField, nIndij, nr_of_Densities
logical(kind=iwp) :: DFT_Storage, Special_NoSym

nr_of_Densities = 1  ! Hardwired option

nIndij = nTri_Elem(S%nShlls)
nField = 2+nr_of_Densities
call mma_allocate(ipOffD,nField,nIndij,label='ipOffD')

! The array with desymmetrized densities contain two additional
! fields.
! ipD00 is a null matrix, which should simplify the logic.
! ipDijS is an auxiliary memory if not the whole set of a
!  desymmetrized density could be used.

nDeDe_tot = nDeDe+MaxDe*nIrrep+MxDij
call mma_allocate(DeDe,nDeDe_tot,Label='DeDe')
ipDeDe = 1
ipD00 = ipDeDe+nDeDe
ipDijS = ipD00+MaxDe*nIrrep
DeDe(:) = Zero

Special_NoSym = .true.
DFT_Storage = .false.
call mk_DeDe(Dens,nDens,nr_of_Densities,ipOffD,nIndij,ipDeDe,ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,DeDe,nDeDe)
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of no symmetry do a temporary square copy of the
! density matrix.

! Change the folded density to pure triangular form,
! i.e. off-diagonal elements are divided by two.

if (nIrrep == 1) then
  Dens(:) = Half*Dens(:)
  ij = 0
  do i=1,nBas(0)
    ij = ij+i
    Dens(ij) = Two*Dens(ij)
  end do
  mDens = nbas(0)*nbas(0)
  call mma_allocate(Dq,mDens,Label='Dq')
  call Square(Dens,Dq,1,nbas(0),nbas(0))
  pDq => Dq(:)

  call mma_allocate(Fq,mDens,Label='Fq')
  Fq(:) = Zero
  pFq => Fq(:)
else
  mDens = nDens
  pDq => Dens(:)
  pFq => Twoham(:)
end if

end subroutine DeDe_SCF
