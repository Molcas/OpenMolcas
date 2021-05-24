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

implicit none
real*8 Const
integer nD
integer ip_VBlocks(nD)
integer ip_WBlocks(nD)
integer AB, CD
#include "WrkSpc.fh"

integer LDF_nBasAux_Pair
external LDF_nBasAux_Pair

integer MAB, MCD
integer ip_Int, l_Int
integer iD
integer ipV, ipW

! Get row and column dimension of integrals
MAB = LDF_nBasAux_Pair(AB)
MCD = LDF_nBasAux_Pair(CD)

! Return if nothing to do
if ((MAB < 1) .or. (MCD < 1)) return

! Allocate integrals (J_AB | K_CD)
l_Int = MAB*MCD
call GetMem('LDFFJK','Allo','Real',ip_Int,l_Int)

! Compute integrals (J_AB | K_CD)
call LDF_ComputeIntegrals_JK_2P(AB,CD,l_Int,Work(ip_Int))

! Compute contributions
do iD=1,nD
  ipV = iWork(ip_VBlocks(iD)-1+CD)
  ipW = iWork(ip_WBlocks(iD)-1+AB)
  call dGeMV_('N',MAB,MCD,Const,Work(ip_Int),MAB,Work(ipV),1,1.0d0,Work(ipW),1)
end do

! Deallocate integrals
call GetMem('LDFFJK','Free','Real',ip_Int,l_Int)

end subroutine LDF_Fock_CoulombOnly0_3
