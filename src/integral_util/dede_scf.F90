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

#include "compiler_features.h"
#ifdef _IN_MODULE_

subroutine DeDe_SCF(Dens,TwoHam,nDens,mDens)

use k2_arrays, only: nDeDe, MaxDe, MxDij, ipDeDe, ipD00, ipDijS, DeDe, pDq, Dq, Fq, pFq, ipOffD
use Basis_Info, only: nBas
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, Half, Two
use stdalloc, only: mma_allocate

implicit none
integer nDens, mDens
real*8, target :: Dens(nDens), TwoHam(nDens)
logical Special_NoSym, DFT_Storage
integer nr_of_Densities, nIndij, nField, nDeDe_Tot, ij, i
integer mDeDe, mIndij

nr_of_Densities = 1  ! Hardwired option

nIndij = S%nShlls*(S%nShlls+1)/2
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
  call DScal_(nDens,Half,Dens,1)
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

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(DeDe_SCF)

#endif
