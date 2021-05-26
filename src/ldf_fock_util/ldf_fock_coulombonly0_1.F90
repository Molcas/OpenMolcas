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

subroutine LDF_Fock_CoulombOnly0_1(nD,FactC,ip_VBlocks,ip_FBlocks,AB,CD)
! Thomas Bondo Pedersen, September 2010.
!
!          OLD CODE
!
! Purpose: Compute
!
!      F(u_A v_B) = F(u_A v_B) + FactC*sum_[K_CD] (u_A v_B | K_CD)*V(K_CD)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nD, ip_VBlocks(nD), ip_FBlocks(nD), AB, CD
real(kind=wp), intent(in) :: FactC(nD)
integer(kind=iwp) :: A, B, nuv, M, l_Int, iD, ipV, ipF
real(kind=wp), allocatable :: FuvJ1(:)
integer(kind=iwp), external :: LDF_nBas_Atom, LDF_nBasAux_Pair
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

! Get row and column dimensions of integrals
A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
M = LDF_nBasAux_Pair(CD)

! Return if nothing to do
if ((nuv < 1) .or. (M < 1)) return

! Allocate integrals (u_A v_B | K_CD)
l_Int = nuv*M
call mma_allocate(FuvJ1,l_Int,label='LDFFuvJ1')

! Compute integrals (u_A v_B | K_CD)
call LDF_ComputeIntegrals_uvJ_2P(AB,CD,l_Int,FuvJ1)

! Compute contributions
do iD=1,nD
  ipV = iWork(ip_VBlocks(iD)-1+CD)
  ipF = iWork(ip_FBlocks(iD)-1+AB)
  call dGeMV_('N',nuv,M,FactC(iD),FuvJ1,nuv,Work(ipV),1,One,Work(ipF),1)
end do

! Deallocate integrals
call mma_deallocate(FuvJ1)

end subroutine LDF_Fock_CoulombOnly0_1
