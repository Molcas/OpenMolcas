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

subroutine LDF_Fock_CoulombOnly0_3(Const,nD,ip_VBlocks,ip_WBlocks,AB,CD)
! Thomas Bondo Pedersen, September 2010.
!
!          OLD CODE
!
! Purpose: Compute
!
!      W(J_AB) = W(J_AB) + Const * sum_[K_CD] (J_AB | K_CD)*V(K_CD)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Const
integer(kind=iwp), intent(in) :: nD, ip_VBlocks(nD), ip_WBlocks(nD), AB, CD
integer(kind=iwp) :: MAB, MCD, l_Int, iD, ipV, ipW
real(kind=wp), allocatable :: FJK(:)
integer(kind=iwp), external :: LDF_nBasAux_Pair
#include "WrkSpc.fh"

! Get row and column dimension of integrals
MAB = LDF_nBasAux_Pair(AB)
MCD = LDF_nBasAux_Pair(CD)

! Return if nothing to do
if ((MAB < 1) .or. (MCD < 1)) return

! Allocate integrals (J_AB | K_CD)
l_Int = MAB*MCD
call mma_allocate(FJK,l_Int,label='LDFFJK')

! Compute integrals (J_AB | K_CD)
call LDF_ComputeIntegrals_JK_2P(AB,CD,l_Int,FJK)

! Compute contributions
do iD=1,nD
  ipV = iWork(ip_VBlocks(iD)-1+CD)
  ipW = iWork(ip_WBlocks(iD)-1+AB)
  call dGeMV_('N',MAB,MCD,Const,FJK,MAB,Work(ipV),1,One,Work(ipW),1)
end do

! Deallocate integrals
call mma_deallocate(FJK)

end subroutine LDF_Fock_CoulombOnly0_3
