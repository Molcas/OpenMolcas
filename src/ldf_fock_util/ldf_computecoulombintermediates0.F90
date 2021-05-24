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
!          V(J) = sum_uv C(uv,J)*D(uv)
!
!          using LDF fitting coefficients.
!
! BLOCKED VERSION

implicit none
integer nD
integer ip_DBlocks(nD)
integer ip_VBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

integer LDF_nBas_Atom, LDF_nBasAux_Pair
external LDF_nBas_Atom, LDF_nBasAux_Pair

logical Rsv_Tsk
external Rsv_Tsk

integer TaskListID
integer iAtomPair
integer iAtom, jAtom
integer ip_C, l_C
integer iD
integer ipV, ipD
integer nuv, M

integer i, j
integer AP_Atoms
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

#ifdef _MOLCAS_MPP_
! Initialize V arrays
do iD=1,nD
  call LDF_ZeroBlockVector(ip_VBlocks(iD))
end do
#endif

! Allocate array for storing coefficients
l_C = 0
do iAtomPair=1,NumberOfAtomPairs
  iAtom = AP_Atoms(1,iAtomPair)
  jAtom = AP_Atoms(2,iAtomPair)
  nuv = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  M = LDF_nBasAux_Pair(iAtomPair)
  l_C = max(l_C,nuv*M)
end do
call GetMem('LDFCBlk','Allo','Real',ip_C,l_C)

! Compute atom pair blocks of V
call Init_Tsk(TaskListID,NumberOfAtomPairs)
do while (Rsv_Tsk(TaskListID,iAtomPair))
  ! Get dimensions
  iAtom = AP_Atoms(1,iAtomPair)
  jAtom = AP_Atoms(2,iAtomPair)
  nuv = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  M = LDF_nBasAux_Pair(iAtomPair)
  ! Read coefficients for this atom pair
  call LDF_CIO_ReadC(iAtomPair,Work(ip_C),l_C)
  ! Compute V for this atom pair and for each density
  do iD=1,nD
    ipD = iWork(ip_DBlocks(iD)-1+iAtomPair)
    ipV = iWork(ip_VBlocks(iD)-1+iAtomPair)
    call dGeMV_('T',nuv,M,1.0d0,Work(ip_C),nuv,Work(ipD),1,0.0d0,Work(ipV),1)
  end do
end do
call Free_Tsk(TaskListID)

! Deallocation
call GetMem('LDFCBlk','Free','Real',ip_C,l_C)

#ifdef _MOLCAS_MPP_
! Get complete V on all nodes
do iD=1,nD
  call LDF_P_AddBlockVector(ip_VBlocks(iD))
end do
#endif

end subroutine LDF_ComputeCoulombIntermediates0
