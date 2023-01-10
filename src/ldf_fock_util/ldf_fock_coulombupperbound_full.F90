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

subroutine LDF_Fock_CoulombUpperBound_Full(PrintNorm,Add,PackedD,PackedF,nD,FactC,ip_D,F)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute the upper bound correction to the LDF Fock
!          matrix (Coulomb only):
!
!      |Fexact(uv)-F(uv)| <= sqrt[(Delta(uv)|Delta(uv))]*U
!
!      U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: PrintNorm, Add, PackedD, PackedF
integer(kind=iwp), intent(in) :: nD, ip_D(nD)
real(kind=wp), intent(in) :: FactC(nD)
real(kind=wp), intent(inout) :: F(*)
integer(kind=iwp) :: l, iD, l_DBlkP, l_FBlkP
integer(kind=iwp), allocatable :: DBlkP(:), FBlkP(:)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

! Return if nothing to do
if (nD < 1) return
if (NumberOfAtomPairs < 1) return

! Allocate and extract density matrix blocks
l_DBlkP = nD
call mma_allocate(DBlkP,l_DBlkP,label='CUBFDBP')
do iD=1,nD
  call LDF_AllocateBlockMatrix('UBD',DBlkP(iD))
  call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,DBlkP(iD))
  call LDF_ScaleOffdiagonalMatrixBlocks(DBlkP(iD),Two)
end do

! If not Add, initialize Fock matrices
if (PackedF) then
  l = nBas_Valence*(nBas_Valence+1)/2
else
  l = nBas_Valence**2
end if
if (.not. Add) then
  F(1:nD*l) = Zero
end if

! Allocate and extract Fock matrix blocks
l_FBlkP = nD
call mma_allocate(FBlkP,l_FBlkP,label='CUBFFBP')
do iD=1,nD
  call LDF_AllocateBlockMatrix('Fck',FBlkP(iD))
  call LDF_Full2Blocked(F((iD-1)*l+1),PackedF,FBlkP(iD))
end do

! Compute upper bound and add to blocked Fock matrices
call LDF_Fock_CoulombUpperBound(PrintNorm,nD,FactC,DBlkP,FBlkP)

! Get full storage Fock matrices from blocked ones
do iD=1,nD
  call LDF_Blocked2Full(FBlkP(iD),PackedF,F((iD-1)*l+1))
end do

! Deallocations
do iD=1,nD
  call LDF_DeallocateBlockMatrix('Fck',FBlkP(iD))
end do
call mma_deallocate(FBlkP)
do iD=1,nD
  call LDF_DeallocateBlockMatrix('UBD',DBlkP(iD))
end do
call mma_deallocate(DBlkP)

end subroutine LDF_Fock_CoulombUpperBound_Full
