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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: PrintNorm, PackedD
integer(kind=iwp), intent(in) :: nD, ip_D(nD)
real(kind=wp), intent(in) :: FactC(nD)
real(kind=wp), intent(out) :: UBFNorm(nD)
integer(kind=iwp) :: l_DBlkP, iD
integer(kind=iwp), allocatable :: DBlkP(:)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

if (nD < 1) return
if (NumberOfAtomPairs < 1) return

l_DBlkP = nD
call mma_allocate(DBlkP,l_DBlkP,label='CUBDBP')
do iD=1,nD
  call LDF_AllocateBlockMatrix('UBD',DBlkP(iD))
  call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,DBlkP(iD))
  call LDF_ScaleOffDiagonalMatrixBlocks(DBlkP(iD),Two)
end do
call LDF_Fock_CoulombUpperBoundNorm(PrintNorm,nD,FactC,DBlkP,UBFNorm)
do iD=1,nD
  call LDF_DeallocateBlockMatrix('UBD',DBlkP(iD))
end do
call mma_deallocate(DBlkP)

end subroutine LDF_Fock_CoulombUpperBoundNorm_Full
