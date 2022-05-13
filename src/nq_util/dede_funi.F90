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

subroutine DeDe_Funi(Dens,nDens,nr_of_Densities)

use k2_arrays, only: DeDe, ipD00, ipDeDe, ipDijs, ipOffD, MaxDE, nDeDe_DFT
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDens, nr_of_Densities
real(kind=wp), intent(in) :: Dens(nDens,nr_of_Densities)
integer(kind=iwp) :: mDeDe, mIndij, nField, nIndij
logical(kind=iwp) :: DFT_Storage, Special_NoSym

nIndij = S%nShlls*(S%nShlls+1)/2
nField = 2+nr_of_Densities

call mma_allocate(ipOffD,nField,nIndij,label='ipOffD')
call mma_allocate(DeDe,nDeDe_DFT+MaxDe*nIrrep,Label='DeDe')
ipDeDe = 1
ipD00 = ipDeDe+nDeDe_DFT
ipDijs = -1  ! Dummy value
DeDe(:) = Zero

Special_NoSym = .false.
DFT_Storage = .true.
call mk_DeDe(Dens,nDens,nr_of_Densities,ipOffD,nIndij,ipDeDe,ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,DeDe,nDeDe_DFT)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine DeDe_Funi
