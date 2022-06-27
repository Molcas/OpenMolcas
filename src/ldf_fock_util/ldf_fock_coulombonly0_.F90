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

subroutine LDF_Fock_CoulombOnly0_(Mode,nD,FactC,ip_DBlocks,ip_VBlocks,ip_FBlocks)
! Thomas Bondo Pedersen, September 2010.
!
!          OLD CODE
!
! Purpose: Compute Coulomb contributions to the Fock matrix using
!          Coulomb intermediates V.
!
! See LDF_Fock_CoulombOnly for an outline of the algorithm.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Mode, nD, ip_DBlocks(nD), ip_VBlocks(nD), ip_FBlocks(nD)
real(kind=wp), intent(in) :: FactC(nD)
integer(kind=iwp) :: l_WBlkP, iD, TaskListID, AB, CD, l, A, B, nuv, M, ipW, ipF
integer(kind=iwp), allocatable :: WBlkP(:)
real(kind=wp), allocatable :: C_AB(:)
character(len=*), parameter :: SecNam = 'LDF_Fock_CoulombOnly0_'
logical(kind=iwp), external :: Rsv_Tsk
integer(kind=iwp), external :: LDF_nBas_Atom, LDF_nBasAux_Pair
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

! Allocate and initialize W intermediates
l_WBlkP = nD
call mma_allocate(WBlkP,l_WBlkP,label='WBlk_P')
do iD=1,nD
  call LDF_AllocateBlockVector('Win',WBlkP(iD))
  call LDF_ZeroBlockVector(WBlkP(iD))
end do

if ((Mode == 1) .or. (Mode == 3)) then
  ! Parallel loop over atom pairs AB (A>=B)
  call Init_Tsk(TaskListID,NumberOfAtomPairs)
  do while (Rsv_Tsk(TaskListID,AB))
    ! Serial loop over atom pairs CD (C>=D)
    do CD=1,NumberOfAtomPairs
      ! F(u_A v_B) = F(u_A v_B) + sum_[K_CD] (u_A v_B | K_CD)*V(K_CD)
      call LDF_Fock_CoulombOnly0_1(nD,FactC,ip_VBlocks,ip_FBlocks,AB,CD)
      ! W(J_AB) = W(J_AB) + sum_[k_C l_D] (J_AB | k_C l_D)*D(k_C l_D)
      call LDF_Fock_CoulombOnly0_2(nD,ip_DBlocks,WBlkP,AB,CD)
      if (Mode == 1) then
        ! W(J_AB) = W(J_AB) - sum_[K_CD] (J_AB | K_CD)*V(K_CD)
        call LDF_Fock_CoulombOnly0_3(-One,nD,ip_VBlocks,WBlkP,AB,CD)
      end if
    end do ! end serial loop over atom pairs CD
    ! F(u_A v_B) = F(u_A v_B) + sum_[J_AB] C(u_A v_B,J_AB)*W(J_AB)
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    M = LDF_nBasAux_Pair(AB)
    l = nuv*M
    call mma_allocate(C_AB,l,label='C_AB')
    call LDF_CIO_ReadC(AB,C_AB,l)
    do iD=1,nD
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      ipW = iWork(WBlkP(iD)-1+AB)
      call dGeMV_('N',nuv,M,FactC(iD),C_AB,nuv,Work(ipW),1,One,Work(ipF),1)
    end do
    call mma_deallocate(C_AB)
  end do ! end parallel loop over atom pairs AB
  call Free_Tsk(TaskListID)
else if (Mode == 2) then ! non-robust fitting
  ! Parallel loop over atom pairs AB (A>=B)
  call Init_Tsk(TaskListID,NumberOfAtomPairs)
  do while (Rsv_Tsk(TaskListID,AB))
    ! Serial loop over atom pairs CD (C>=D)
    do CD=1,NumberOfAtomPairs
      ! W(J_AB) = W(J_AB) + sum_[K_CD] (J_AB | K_CD)*V(K_CD)
      call LDF_Fock_CoulombOnly0_3(One,nD,ip_VBlocks,WBlkP,AB,CD)
    end do ! end serial loop over atom pairs CD
    ! F(u_A v_B) = F(u_A v_B) + sum_[J_AB] C(u_A v_B,J_AB)*W(J_AB)
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    M = LDF_nBasAux_Pair(AB)
    l = nuv*M
    call mma_allocate(C_AB,l,label='C_AB')
    call LDF_CIO_ReadC(AB,C_AB,l)
    do iD=1,nD
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      ipW = iWork(WBlkP(iD)-1+AB)
      call dGeMV_('N',nuv,M,FactC(iD),C_AB,nuv,Work(ipW),1,One,Work(ipF),1)
    end do
    call mma_deallocate(C_AB)
  end do ! end parallel loop over atom pairs AB
  call Free_Tsk(TaskListID)
else
  write(u6,'(A,A,I6)') SecNam,': unknown Mode:',Mode
  call LDF_NotImplemented()
end if

#ifdef _MOLCAS_MPP_
! Add F over nodes
do iD=1,nD
  call LDF_P_AddBlockMatrix(ip_FBlocks(iD))
end do
#endif

! Deallocate W intermediates
do iD=1,nD
  call LDF_DeallocateBlockVector('Win',WBlkP(iD))
end do
call mma_deallocate(WBlkP)

end subroutine LDF_Fock_CoulombOnly0_
