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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_Fock_CoulombOnly0_2(nD,ip_DBlocks,ip_WBlocks,AB,CD)
! Thomas Bondo Pedersen, September 2010.
!
!          OLD CODE
!
! Purpose: Compute
!
!      W(J_AB) = W(J_AB) + sum_[k_C l_D] (k_C l_D|J_AB)*D(k_C l_D)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nD, ip_DBlocks(nD), ip_WBlocks(nD), AB, CD
integer(kind=iwp) :: C, D, nkl, M, l_Int, iD, ipD, ipW
real(kind=wp), allocatable :: FuvJ2(:)
integer(kind=iwp), external :: LDF_nBas_Atom, LDF_nBasAux_Pair
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

! Get row and column dimensions of integrals
C = iWork(ip_AP_Atoms-1+2*(CD-1)+1)
D = iWork(ip_AP_Atoms-1+2*(CD-1)+2)
nkl = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
M = LDF_nBasAux_Pair(AB)

! Return if nothing to do
if ((nkl < 1) .or. (M < 1)) return

! Allocate integrals (k_C l_D | J_AB)
l_Int = nkl*M
call mma_allocate(FuvJ2,l_Int,label='LDFFuvJ2')

! Compute integrals (k_C l_D | J_AB)
call LDF_ComputeIntegrals_uvJ_2P(CD,AB,l_Int,FuvJ2)

! Compute contributions
do iD=1,nD
  ipD = iWork(ip_DBlocks(iD)-1+CD)
  ipW = iWork(ip_WBlocks(iD)-1+AB)
  call dGeMV_('T',nkl,M,One,FuvJ2,nkl,Work(ipD),1,One,Work(ipW),1)
end do

! Deallocate integrals
call mma_deallocate(FuvJ2)

end subroutine LDF_Fock_CoulombOnly0_2
