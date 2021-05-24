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

implicit none
integer nD
integer ip_DBlocks(nD)
integer ip_WBlocks(nD)
integer AB, CD
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

integer LDF_nBas_Atom, LDF_nBasAux_Pair
external LDF_nBas_Atom, LDF_nBasAux_Pair

integer nkl, M
integer ip_Int, l_Int
integer iD
integer ipD, ipW

integer i, j
integer AP_Atoms
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

! Get row and column dimensions of integrals
nkl = LDF_nBas_Atom(AP_Atoms(1,CD))*LDF_nBas_Atom(AP_Atoms(2,CD))
M = LDF_nBasAux_Pair(AB)

! Return if nothing to do
if ((nkl < 1) .or. (M < 1)) return

! Allocate integrals (k_C l_D | J_AB)
l_Int = nkl*M
call GetMem('LDFFuvJ2','Allo','Real',ip_Int,l_Int)

! Compute integrals (k_C l_D | J_AB)
call LDF_ComputeIntegrals_uvJ_2P(CD,AB,l_Int,Work(ip_Int))

! Compute contributions
do iD=1,nD
  ipD = iWork(ip_DBlocks(iD)-1+CD)
  ipW = iWork(ip_WBlocks(iD)-1+AB)
  call dGeMV_('T',nkl,M,1.0d0,Work(ip_Int),nkl,Work(ipD),1,1.0d0,Work(ipW),1)
end do

! Deallocate integrals
call GetMem('LDFFuvJ2','Free','Real',ip_Int,l_Int)

end subroutine LDF_Fock_CoulombOnly0_2
