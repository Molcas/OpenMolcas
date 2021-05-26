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

subroutine LDF_ComputeCoulombIntermediates0(nD,ip_DBlocks,ip_VBlocks)
! Thomas Bondo Pedersen, September 2010.
!
! Purpose: Compute Coulomb intermediates
!
!      V(J) = sum_uv C(uv,J)*D(uv)
!
!      using LDF fitting coefficients.
!
! BLOCKED VERSION

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nD, ip_DBlocks(nD), ip_VBlocks(nD)
integer(kind=iwp) :: TaskListID, iAtomPair, iAtom, jAtom, l_C, iD, ipV, ipD, nuv, M
real(kind=wp), allocatable :: LDFCBlk(:)
integer(kind=iwp), external :: LDF_nBas_Atom, LDF_nBasAux_Pair
logical(kind=iwp), external :: Rsv_Tsk
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

#ifdef _MOLCAS_MPP_
! Initialize V arrays
do iD=1,nD
  call LDF_ZeroBlockVector(ip_VBlocks(iD))
end do
#endif

! Allocate array for storing coefficients
l_C = 0
do iAtomPair=1,NumberOfAtomPairs
  iAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+1)
  jAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+2)
  nuv = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  M = LDF_nBasAux_Pair(iAtomPair)
  l_C = max(l_C,nuv*M)
end do
call mma_allocate(LDFCBlk,l_C,label='LDFCBlk')

! Compute atom pair blocks of V
call Init_Tsk(TaskListID,NumberOfAtomPairs)
do while (Rsv_Tsk(TaskListID,iAtomPair))
  ! Get dimensions
  iAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+1)
  jAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+2)
  nuv = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  M = LDF_nBasAux_Pair(iAtomPair)
  ! Read coefficients for this atom pair
  call LDF_CIO_ReadC(iAtomPair,LDFCBlk,l_C)
  ! Compute V for this atom pair and for each density
  do iD=1,nD
    ipD = iWork(ip_DBlocks(iD)-1+iAtomPair)
    ipV = iWork(ip_VBlocks(iD)-1+iAtomPair)
    call dGeMV_('T',nuv,M,One,LDFCBlk,nuv,Work(ipD),1,Zero,Work(ipV),1)
  end do
end do
call Free_Tsk(TaskListID)

! Deallocation
call mma_deallocate(LDFCBlk)

#ifdef _MOLCAS_MPP_
! Get complete V on all nodes
do iD=1,nD
  call LDF_P_AddBlockVector(ip_VBlocks(iD))
end do
#endif

end subroutine LDF_ComputeCoulombIntermediates0
