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

subroutine LDF_Fock_CoulombUpperBound_Full(PrintNorm,Add,PackedD,PackedF,nD,FactC,ip_D,ip_F)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute the upper bound correction to the LDF Fock
!          matrix (Coulomb only):
!
!      |Fexact(uv)-F(uv)| <= sqrt[(Delta(uv)|Delta(uv))]*U
!
!      U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|

implicit none
logical PrintNorm
logical Add
logical PackedD
logical PackedF
integer nD
real*8 FactC(nD)
integer ip_D(nD)
integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

integer l
integer iD
integer ip_DBlkP, l_DBlkP
integer ip_FBlkP, l_FBlkP

! Return if nothing to do
if (nD < 1) return
if (NumberOfAtomPairs < 1) return

! Allocate and extract density matrix blocks
l_DBlkP = nD
call GetMem('CUBFDBP','Allo','Inte',ip_DBlkP,l_DBlkP)
do iD=1,nD
  call LDF_AllocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
  call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,iWork(ip_DBlkP-1+iD))
  call LDF_ScaleOffdiagonalMatrixBlocks(iWork(ip_DBlkP-1+iD),2.0d0)
end do

! If not Add, initialize Fock matrices
if (.not. Add) then
  if (PackedF) then
    l = nBas_Valence*(nBas_Valence+1)/2
  else
    l = nBas_Valence**2
  end if
  do iD=1,nD
    call Cho_dZero(Work(ip_F(iD)),l)
  end do
end if

! Allocate and extract Fock matrix blocks
l_FBlkP = nD
call GetMem('CUBFFBP','Allo','Inte',ip_FBlkP,l_FBlkP)
do iD=1,nD
  call LDF_AllocateBlockMatrix('Fck',iWork(ip_FBlkP-1+iD))
  call LDF_Full2Blocked(Work(ip_F(iD)),PackedF,iWork(ip_FBlkP-1+iD))
end do

! Compute upper bound and add to blocked Fock matrices
call LDF_Fock_CoulombUpperBound(PrintNorm,nD,FactC,iWork(ip_DBlkP),iWork(ip_FBlkP))

! Get full storage Fock matrices from blocked ones
do iD=1,nD
  call LDF_Blocked2Full(iWork(ip_FBlkP-1+iD),PackedF,Work(ip_F(iD)))
end do

! Deallocations
do iD=1,nD
  call LDF_DeallocateBlockMatrix('Fck',iWork(ip_FBlkP-1+iD))
end do
call GetMem('CUBFFBP','Allo','Inte',ip_FBlkP,l_FBlkP)
do iD=1,nD
  call LDF_DeallocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
end do
call GetMem('CUBFDBP','Free','Inte',ip_DBlkP,l_DBlkP)

end subroutine LDF_Fock_CoulombUpperBound_Full
