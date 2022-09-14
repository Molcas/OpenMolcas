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

subroutine LDF_Fock_CoulombOnly_(UseExactIntegralDiagonal,Timing,Mode,tau,nD,FactC,ip_DBlocks,ip_V,ip_FBlocks,CNorm,DNorm,VNorm)
! Thomas Bondo Pedersen, October 2010.
!
! Purpose: Compute Coulomb contributions to the Fock matrix using
!          Coulomb intermediates V. Integral-driven algorithm.
!
! See LDF_Fock_CoulombOnly for more details.

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: UseExactIntegralDiagonal, Timing
integer(kind=iwp), intent(in) :: Mode, nD, ip_DBlocks(nD), ip_V(nD), ip_FBlocks(nD)
real(kind=wp), intent(in) :: tau(2), FactC(nD), CNorm(4,*), DNorm(*), VNorm(*)
integer(kind=iwp) :: l_WBlkP, iD, nAtom, TaskListID, jAB, AB, A, B, CD, C, nuv, M, MA, MB, MAB, MCD, l_Int, ipD, ipV, ipF, ipW, &
                     l_C, ipC, l_tauW, ip, l, l_VNrm, l_DNrm
real(kind=wp) :: Const, tC1, tC2, tIC1, tIC2, tW1, tW2, tIW1, tIW2, tIC, tIW, tauWr, GABCD, GCDAB
integer(kind=iwp), allocatable :: WBlkP(:)
real(kind=wp), allocatable :: FckInt(:), FckCoef(:), tauW(:), VNrm(:), DNrm(:)
character(len=*), parameter :: SecNam = 'LDF_Fock_CoulombOnly_'
logical(kind=iwp), external :: Rsv_Tsk
integer(kind=iwp), external :: LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_integral_prescreening_info.fh"
#include "ldf_a2ap.fh"

! Get number of atoms
nAtom = LDF_nAtom()

! Set up prescreening info
l_VNrm = nAtom+NumberOfAtomPairs
call mma_allocate(VNrm,l_VNrm,label='VNrm')
VNrm(:) = Zero
do iD=1,nD
  ip = (nAtom+NumberOfAtomPairs)*(iD-1)
  do AB=1,nAtom+NumberOfAtomPairs
    VNrm(AB) = max(VNrm(AB),VNorm(ip+AB))
  end do
end do
l_DNrm = NumberOfAtomPairs
call mma_allocate(DNrm,l_DNrm,label='DNrm')
DNrm(:) = Zero
do iD=1,nD
  ip = NumberOfAtomPairs*(iD-1)
  do AB=1,NumberOfAtomPairs
    DNrm(AB) = max(DNrm(AB),DNorm(ip+AB))
  end do
end do
l_tauW = nAtom+NumberOfAtomPairs
call mma_allocate(tauW,l_tauW,label='tauW')
if (tau(2) > Zero) then
  call LDF_SetA2AP()
  do A=1,nAtom
    l = iWork(ip_A2AP-1+2*(A-1)+1)
    ip = iWork(ip_A2AP-1+2*(A-1)+2)-1
    tauW(A) = Zero
    do jAB=1,l
      AB = iWork(ip+jAB)
      if (iWork(ip_AP_Atoms-1+2*(AB-1)+1) == A) then
        tauW(A) = max(tauW(A),CNorm(2,AB))
      else if (iWork(ip_AP_Atoms-1+2*(AB-1)+2) == A) then
        tauW(A) = max(tauW(A),CNorm(3,AB))
      else
        call WarningMessage(2,SecNam//': logical error [A2AP]')
        call LDF_Quit(1)
      end if
    end do
    if (tauW(A) > 1.0e-16_wp) then
      tauWr = tau(2)/tauW(A)
    else
      tauWr = 1.0e99_wp
    end if
    tauW(A) = tauWr
  end do
  call LDF_UnsetA2AP()
  do AB=1,NumberOfAtomPairs
    if (iWork(ip_AP_2CFunctions-1+2*(AB-1)+1) > 0) then
      if (CNorm(4,AB) > 1.0e-16_wp) then
        tauW(nAtom+AB) = tau(2)/CNorm(4,AB)
      else
        tauW(nAtom+AB) = 1.0e99_wp
      end if
    else
      tauW(nAtom+AB) = 1.0e99_wp
    end if
  end do
else
  tauW(:) = Zero
end if

! Allocate and initialize W intermediates
l_WBlkP = nD
call mma_allocate(WBlkP,l_WBlkP,label='WBlkP')
do iD=1,nD
  call LDF_AllocateAuxBasVector('Win',WBlkP(iD))
  call LDF_ZeroAuxBasVector(WBlkP(iD))
end do

!======================
! 3-index contributions
!======================

if ((Mode == 1) .or. (Mode == 3)) then
  if (Timing) then
    tIC = Zero
    tIW = Zero
    call CWTime(tC1,tW1)
  end if
  call Init_Tsk(TaskListID,NumberOfAtomPairs)
  do while (Rsv_Tsk(TaskListID,AB))
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do C=1,nAtom
      tauWr = tauW(C)
      if ((Work(ip_IDiag_Sm-1+AB)*Work(ip_GDiag_1C_Sm-1+C)*VNorm(C) >= tau(2)) .or. &
          (Work(ip_IDiag_Sm-1+AB)*Work(ip_GDiag_1C_Sm-1+C)*DNorm(AB) >= tauWr)) then
        ! Compute integrals (u_A v_B|J_C)
        M = LDF_nBasAux_Atom(C)
        l_Int = nuv*M
        call mma_allocate(FckInt,l_Int,label='Fck3Int1')
        if (Timing) call CWTime(tIC1,tIW1)
        call LDF_Compute3IndexIntegrals_1(AB,C,tau(1),l_Int,FckInt)
        if (Timing) then
          call CWTime(tIC2,tIW2)
          tIC = tIC+(tIC2-tIC1)
          tIW = tIW+(tIW2-tIW1)
        end if
        ! Compute Fock matrix contribution
        ! F(u_A v_B)+=FactC*sum_J_C (u_A v_B|J_C)*V(J_C)
        do iD=1,nD
          ipV = iWork(ip_V(iD)-1+C)
          ipF = iWork(ip_FBlocks(iD)-1+AB)
          call dGeMV_('N',nuv,M,FactC(iD),FckInt,nuv,Work(ipV),1,One,Work(ipF),1)
        end do
        ! Compute W contribution
        ! W(J_C)+=sum_u_Av_B (u_A v_B|J_C)*D(u_A v_B)
        do iD=1,nD
          ipD = iWork(ip_DBlocks(iD)-1+AB)
          ipW = iWork(WBlkP(iD)-1+C)
          call dGeMV_('T',nuv,M,One,FckInt,nuv,Work(ipD),1,One,Work(ipW),1)
        end do
        call mma_deallocate(FckInt)
      end if
    end do
    if (LDF2) then
      ! Contributions from two-center aux functions
      do CD=1,NumberOfAtomPairs
        MCD = iWork(ip_AP_2CFunctions-1+2*(CD-1)+1)
        if (MCD > 0) then
          tauWr = tauW(nAtom+CD)
          if ((Work(ip_IDiag_Sm-1+AB)*Work(ip_GDiag_2C_Sm-1+CD)*VNorm(nAtom+CD) >= tau(2)) .or. &
              (Work(ip_IDiag_Sm-1+AB)*Work(ip_GDiag_2C_Sm-1+CD)*DNorm(AB) >= tauWr)) then
            ! Compute integrals (u_A v_B|J_CD)
            l_Int = nuv*MCD
            call mma_allocate(FckInt,l_Int,label='Fck3Int2')
            if (Timing) call CWTime(tIC1,tIW1)
            call LDF_Compute3IndexIntegrals_2(AB,CD,tau(1),l_Int,FckInt)
            if (Timing) then
              call CWTime(tIC2,tIW2)
              tIC = tIC+(tIC2-tIC1)
              tIW = tIW+(tIW2-tIW1)
            end if
            ! Compute Fock matrix contribution
            ! F(u_A v_B)+=
            ! FactC*sum_J_CD (u_A v_B|J_CD)*V(J_CD)
            do iD=1,nD
              ipV = iWork(ip_V(iD)-1+nAtom+CD)
              ipF = iWork(ip_FBlocks(iD)-1+AB)
              call dGeMV_('N',nuv,MCD,FactC(iD),FckInt,nuv,Work(ipV),1,One,Work(ipF),1)
            end do
            ! Compute W contribution
            ! W(J_CD)+=sum_u_Av_B (u_A v_B|J_CD)*D(u_A v_B)
            do iD=1,nD
              ipD = iWork(ip_DBlocks(iD)-1+AB)
              ipW = iWork(WBlkP(iD)-1+nAtom+CD)
              call dGeMV_('T',nuv,MCD,One,FckInt,nuv,Work(ipD),1,One,Work(ipW),1)
            end do
            call mma_deallocate(FckInt)
          end if
        end if
      end do
    end if
  end do
  call Free_Tsk(TaskListID)
  if (Timing) then
    call CWTime(tC2,tW2)
    write(u6,'(A,2(1X,F12.2),A)') 'Time spent on 3-index contributions:              ',tC2-tC1,tW2-tW1,' seconds'
    write(u6,'(A,2(1X,F12.2),A)') '      - of which integrals required:              ',tIC,tIW,' seconds'
  end if
end if

!======================
! 2-index contributions
!======================

if ((Mode == 1) .or. (Mode == 2)) then
  if (Timing) then
    tIC = Zero
    tIW = Zero
    call CWTime(tC1,tW1)
  end if
  if (Mode == 1) then
    Const = -One
  else
    Const = One
  end if
  call Init_Tsk(TaskListID,nAtom)
  do while (Rsv_Tsk(TaskListID,A))
    MA = LDF_nBasAux_Atom(A)
    if (MA > 0) then
      do B=1,A-1
        MB = LDF_nBasAux_Atom(B)
        if (MB > 0) then
          if ((Work(ip_GDiag_1C_Sm-1+A)*Work(ip_GDiag_1C_Sm-1+B)*VNorm(B) >= tauW(A)) .or. &
              (Work(ip_GDiag_1C_Sm-1+A)*Work(ip_GDiag_1C_Sm-1+B)*VNorm(A) >= tauW(B))) then
            ! Compute integrals (J_A|K_B)
            l_Int = MA*MB
            call mma_allocate(FckInt,l_Int,label='Fck2Int11')
            if (Timing) call CWTime(tIC1,tIW1)
            call LDF_Compute2IndexIntegrals_11(A,B,tau(1),l_Int,FckInt)
            if (Timing) then
              call CWTime(tIC2,tIW2)
              tIC = tIC+(tIC2-tIC1)
              tIW = tIW+(tIW2-tIW1)
            end if
            ! Compute W contribution
            ! W(J_A)+=Const*sum_K_B (J_A|K_B)*V(K_B)
            do iD=1,nD
              ipV = iWork(ip_V(iD)-1+B)
              ipW = iWork(WBlkP(iD)-1+A)
              call dGeMV_('N',MA,MB,Const,FckInt,MA,Work(ipV),1,One,Work(ipW),1)
            end do
            ! Compute W contribution
            ! W(K_B)+=Const*sum_J_A (J_A|K_B)*V(J_A)
            do iD=1,nD
              ipV = iWork(ip_V(iD)-1+A)
              ipW = iWork(WBlkP(iD)-1+B)
              call dGeMV_('T',MA,MB,Const,FckInt,MA,Work(ipV),1,One,Work(ipW),1)
            end do
            call mma_deallocate(FckInt)
          end if
        end if
      end do
      if (Work(ip_GDiag_1C_Sm-1+A)*Work(ip_GDiag_1C_Sm-1+A)*VNorm(A) >= tauW(A)) then
        ! Compute integrals (J_A|K_A)
        l_Int = MA**2
        call mma_allocate(FckInt,l_Int,label='Fck2Int11')
        if (Timing) call CWTime(tIC1,tIW1)
        call LDF_Compute2IndexIntegrals_11(A,A,tau(1),l_Int,FckInt)
        if (Timing) then
          call CWTime(tIC2,tIW2)
          tIC = tIC+(tIC2-tIC1)
          tIW = tIW+(tIW2-tIW1)
        end if
        ! Compute W contribution
        ! W(J_A)+=Const*sum_K_A (J_A|K_A)*V(K_A)
        do iD=1,nD
          ipV = iWork(ip_V(iD)-1+A)
          ipW = iWork(WBlkP(iD)-1+A)
          call dGeMV_('N',MA,MA,Const,FckInt,MA,Work(ipV),1,One,Work(ipW),1)
        end do
        call mma_deallocate(FckInt)
      end if
      if (LDF2) then
        ! Two-center contributions
        do CD=1,NumberOfAtomPairs
          MCD = iWork(ip_AP_2CFunctions-1+2*(CD-1)+1)
          if (MCD > 0) then
            if ((Work(ip_GDiag_1C_Sm-1+A)*Work(ip_GDiag_2C_Sm-1+CD)*VNorm(nAtom+CD) >= tauW(A)) .or. &
                (Work(ip_GDiag_1C_Sm-1+A)*Work(ip_GDiag_2C_Sm-1+CD)*VNorm(A) >= tauW(nAtom+CD))) then
              ! Compute integrals (J_A|K_CD)
              l_Int = MA*MCD
              call mma_allocate(FckInt,l_Int,label='Fck2Int12')
              if (Timing) call CWTime(tIC1,tIW1)
              call LDF_Compute2IndexIntegrals_12(A,CD,tau(1),l_Int,FckInt)
              if (Timing) then
                call CWTime(tIC2,tIW2)
                tIC = tIC+(tIC2-tIC1)
                tIW = tIW+(tIW2-tIW1)
              end if
              ! Compute W contribution
              ! W(J_A)+=Const*sum_K_CD (J_A|K_CD)*V(K_CD)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+nAtom+CD)
                ipW = iWork(WBlkP(iD)-1+A)
                call dGeMV_('N',MA,MCD,Const,FckInt,MA,Work(ipV),1,One,Work(ipW),1)
              end do
              ! Compute W contribution
              ! W(K_CD)+=Const*sum_J_A (J_A|K_CD)*V(J_A)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+A)
                ipW = iWork(WBlkP(iD)-1+nAtom+CD)
                call dGeMV_('T',MA,MCD,Const,FckInt,MA,Work(ipV),1,One,Work(ipW),1)
              end do
              call mma_deallocate(FckInt)
            end if
          end if
        end do
      end if
    end if
  end do
  call Free_Tsk(TaskListID)
  if (LDF2) then
    call Init_Tsk(TaskListID,NumberOfAtomPairs)
    do while (Rsv_Tsk(TaskListID,AB))
      MAB = iWork(ip_AP_2CFunctions-1+2*(AB-1)+1)
      if (MAB > 0) then
        do CD=1,AB-1
          MCD = iWork(ip_AP_2CFunctions-1+2*(CD-1)+1)
          if (MCD > 0) then
            GABCD = Work(ip_GDiag_2C_Sm-1+AB)*Work(ip_GDiag_2C_Sm-1+CD)*VNorm(nAtom+CD)
            GCDAB = Work(ip_GDiag_2C_Sm-1+AB)*Work(ip_GDiag_2C_Sm-1+CD)*VNorm(nAtom+AB)
            if ((GABCD >= tauW(nAtom+AB)) .or. (GCDAB >= tauW(nAtom+CD))) then
              ! Compute integrals (J_AB|K_CD)
              l_Int = MAB*MCD
              call mma_allocate(FckInt,l_Int,label='Fck2Int22')
              if (Timing) call CWTime(tIC1,tIW1)
              call LDF_Compute2IndexIntegrals_22(AB,CD,tau(1),l_Int,FckInt)
              if (Timing) then
                call CWTime(tIC2,tIW2)
                tIC = tIC+(tIC2-tIC1)
                tIW = tIW+(tIW2-tIW1)
              end if
              ! Compute W contribution
              ! W(J_AB)+=Const*sum_K_CD (J_AB|K_CD)*V(K_CD)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+nAtom+CD)
                ipW = iWork(WBlkP(iD)-1+nAtom+AB)
                call dGeMV_('N',MAB,MCD,Const,FckInt,MAB,Work(ipV),1,One,Work(ipW),1)
              end do
              ! Compute W contribution
              ! W(K_CD)+=Const*sum_J_AB (J_AB|K_CD)*V(J_AB)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+nAtom+AB)
                ipW = iWork(WBlkP(iD)-1+nAtom+CD)
                call dGeMV_('T',MAB,MCD,Const,FckInt,MAB,Work(ipV),1,One,Work(ipW),1)
              end do
              call mma_deallocate(FckInt)
            end if
          end if
        end do
        ! Compute integrals (J_AB|K_AB)
        if (Work(ip_GDiag_2C_Sm-1+AB)*Work(ip_GDiag_2C_Sm-1+AB)*VNorm(nAtom+AB) >= tauW(nAtom+AB)) then
          l_Int = MAB**2
          call mma_allocate(FckInt,l_Int,label='Fck2Int22')
          if (Timing) call CWTime(tIC1,tIW1)
          call LDF_Compute2IndexIntegrals_22(AB,AB,tau(1),l_Int,FckInt)
          if (Timing) then
            call CWTime(tIC2,tIW2)
            tIC = tIC+(tIC2-tIC1)
            tIW = tIW+(tIW2-tIW1)
          end if
          ! Compute W contribution
          ! W(J_AB)+=Const*sum_K_AB (J_AB|K_AB)*V(K_AB)
          do iD=1,nD
            ipV = iWork(ip_V(iD)-1+nAtom+AB)
            ipW = iWork(WBlkP(iD)-1+nAtom+AB)
            call dGeMV_('N',MAB,MAB,Const,FckInt,MAB,Work(ipV),1,One,Work(ipW),1)
          end do
          call mma_deallocate(FckInt)
        end if
      end if
    end do
    call Free_Tsk(TaskListID)
  end if
  if (Timing) then
    call CWTime(tC2,tW2)
    write(u6,'(A,2(1X,F12.2),A)') 'Time spent on 2-index contributions:              ',tC2-tC1,tW2-tW1,' seconds'
    write(u6,'(A,2(1X,F12.2),A)') '      - of which integrals required:              ',tIC,tIW,' seconds'
  end if
end if

! Deallocate prescreening info
call mma_deallocate(tauW)
if (l_DNrm > 0) then
  call mma_deallocate(DNrm)
end if
if (l_VNrm > 0) then
  call mma_deallocate(VNrm)
end if

#ifdef _MOLCAS_MPP_
! Add W over nodes
if ((nProcs > 1) .and. Is_Real_Par()) then
  if (Timing) call CWTime(tC1,tW1)
  do iD=1,nD
    call LDF_P_AddAuxBasVector(WBlkP(iD))
  end do
  if (Timing) then
    call CWTime(tC2,tW2)
    write(u6,'(A,2(1X,F12.2),A)') 'Parallel overhead for W intermediates:            ',tC2-tC1,tW2-tW1,' seconds'
  end if
end if
#endif

!===============
! Compute F+=C*W
!===============

if (Timing) call CWTime(tC1,tW1)
call Init_Tsk(TaskListID,NumberOfAtomPairs)
do while (Rsv_Tsk(TaskListID,AB))
  ! Read coefficients
  A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
  B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
  nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
  l_C = nuv*LDF_nBasAux_Pair_wLD(AB)
  call mma_allocate(FckCoef,l_C,label='FckCoef')
  call LDF_CIO_ReadC_wLD(AB,FckCoef,l_C)
  ! Compute F(u_A v_B) += FactC*sum_J C(u_A v_B,J)*W(J)
  ipC = 1
  MA = LDF_nBasAux_Atom(A)
  if (MA > 0) then
    do iD=1,nD
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      ipW = iWork(WBlkP(iD)-1+A)
      call dGeMV_('N',nuv,MA,FactC(iD),FckCoef(ipC),nuv,Work(ipW),1,One,Work(ipF),1)
    end do
    ipC = ipC+nuv*MA
  end if
  if (B /= A) then
    MB = LDF_nBasAux_Atom(B)
    if (MB > 0) then
      do iD=1,nD
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        ipW = iWork(WBlkP(iD)-1+B)
        call dGeMV_('N',nuv,MB,FactC(iD),FckCoef(ipC),nuv,Work(ipW),1,One,Work(ipF),1)
      end do
      ipC = ipC+nuv*MB
    end if
  end if
  MAB = iWork(ip_AP_2CFunctions-1+2*(AB-1)+1)
  if (MAB > 0) then
    do iD=1,nD
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      ipW = iWork(WBlkP(iD)-1+nAtom+AB)
      call dGeMV_('N',nuv,MAB,FactC(iD),FckCoef(ipC),nuv,Work(ipW),1,One,Work(ipF),1)
    end do
  end if
  call mma_deallocate(FckCoef)
end do
call Free_Tsk(TaskListID)
if (Timing) then
  call CWTime(tC2,tW2)
  write(u6,'(A,2(1X,F12.2),A)') 'Time spent computing C*W contribution:            ',tC2-tC1,tW2-tW1,' seconds'
end if

! Deallocate W intermediates
do iD=1,nD
  call LDF_DeallocateAuxBasVector('Win',WBlkP(iD))
end do
call mma_deallocate(WBlkP)

! Use exact integral diagonal blocks if requested.
! Computed by adding the correction
! dF(uv) = FactC * sum_kl { (uv|kl) - [uv|kl] }*D(kl)
! where |uv) and |kl) belong to the same atom pair
! and [uv|kl] is an LDF integral approximation.
if (UseExactIntegralDiagonal) then
  if (Timing) call CWTime(tC1,tW1)
  call LDF_Fock_CoulombOnly_XIDI(Mode,tau(1),nD,FactC,ip_DBlocks,ip_FBlocks)
  if (Timing) then
    call CWTime(tC2,tW2)
    write(u6,'(A,2(1X,F12.2),A)') 'Time spent computing XIDI corrections:            ',tC2-tC1,tW2-tW1,' seconds'
  end if
end if

#ifdef _MOLCAS_MPP_
if ((nProcs > 1) .and. Is_Real_Par()) then
  if (Timing) call CWTime(tC1,tW1)
  do iD=1,nD
    call LDF_P_AddBlockMatrix(ip_FBlocks(iD))
  end do
  if (Timing) then
    call CWTime(tC2,tW2)
    write(u6,'(A,2(1X,F12.2),A)') 'Parallel overhead for F blocks:                   ',tC2-tC1,tW2-tW1,' seconds'
  end if
end if
#endif

if (Timing) call xFlush(u6)

end subroutine LDF_Fock_CoulombOnly_
