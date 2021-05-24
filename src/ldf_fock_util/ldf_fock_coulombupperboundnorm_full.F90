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
! Copyright (C) 2011, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_Fock_CoulombUpperBoundNorm_Full(PrintNorm,PackedD,nD,FactC,ip_D,UBFNorm)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute norm of the upper bound to the Coulomb Fock
!          matrix error in LDF.

implicit none
logical PrintNorm
logical PackedD
integer nD
real*8 FactC(nD)
integer ip_D(nD)
real*8 UBFNorm(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

integer ip_DBlkP, l_DBlkP
integer iD

if (nD < 1) return
if (NumberOfAtomPairs < 1) return

l_DBlkP = nD
call GetMem('CUBDBP','Allo','Inte',ip_DBlkP,l_DBlkP)
do iD=1,nD
  call LDF_AllocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
  call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,iWork(ip_DBlkP-1+iD))
  call LDF_ScaleOffDiagonalMatrixBlocks(iWork(ip_DBlkP-1+iD),2.0d0)
end do
call LDF_Fock_CoulombUpperBoundNorm(PrintNorm,nD,FactC,iWork(ip_DBlkP),UBFNorm)
do iD=1,nD
  call LDF_DeallocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
end do
call GetMem('CUBDBP','Free','Inte',ip_DBlkP,l_DBlkP)

end subroutine LDF_Fock_CoulombUpperBoundNorm_Full
