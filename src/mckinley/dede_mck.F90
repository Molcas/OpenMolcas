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

subroutine DeDe_Mck(FD,nFD,ipOffD,nOffD,DDen,lDDen,mDeDe,mIndij)

use k2_arrays, only: MaxDe
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nFD, nOffD, lDDen
real(kind=wp), intent(in) :: FD(nFD)
integer(kind=iwp), intent(out) :: ipOffD(3,nOffD), mDeDe, mIndij
real(kind=wp), intent(out) :: DDen(lDDen)
integer(kind=iwp) :: ipD00, ipDeDe, nr_of_Densities
logical(kind=iwp) :: DFT_Storage, Special_NoSym

Special_NoSym = .false.
DFT_Storage = .false.
nr_of_Densities = 1

ipDeDe = 1
ipD00 = 1
! ipDijS is controlled in the calling routine
call mk_DeDe(FD,nFD,nr_of_Densities,ipOffD,nOffD,ipDeDe,ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,DDen,lDDen)

return

end subroutine DeDe_Mck
