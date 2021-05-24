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

subroutine LDF_Fock_CoulombOnly_(UseExactIntegralDiagonal,Timing,Mode,tau,nD,FactC,ip_DBlocks,ip_V,ip_FBlocks,ip_CNorm,ip_DNorm, &
                                 ip_VNorm)
! Thomas Bondo Pedersen, October 2010.
!
! Purpose: Compute Coulomb contributions to the Fock matrix using
!          Coulomb intermediates V. Integral-driven algorithm.
!
! See LDF_Fock_CoulombOnly for more details.

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif

implicit none
logical UseExactIntegralDiagonal
logical Timing
integer Mode
real*8 tau(2)
integer nD
real*8 FactC(nD)
integer ip_DBlocks(nD)
integer ip_V(nD)
integer ip_FBlocks(nD)
integer ip_CNorm
integer ip_DNorm
integer ip_VNorm
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_integral_prescreening_info.fh"
#include "ldf_a2ap.fh"

character*21 SecNam
parameter(SecNam='LDF_Fock_CoulombOnly_')

logical Rsv_Tsk
external Rsv_Tsk

integer LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom
integer LDF_nBasAux_Pair_wLD
external LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom
external LDF_nBasAux_Pair_wLD

integer ip_WBlkP, l_WBlkP
integer iD
integer nAtom
integer TaskListID
integer jAB, AB, A, B
integer CD, C
integer nuv, M, MA, MB, MAB, MCD
integer ip_Int, l_Int
integer ipD, ipV, ipF, ipW
integer ip_C, l_C, ipC
integer ip_tauW, l_tauW, ip_tauWA
integer ip, l
integer ip_VNrm, l_VNrm
integer ip_DNrm, l_DNrm

real*8 Const
real*8 tC1, tC2, tIC1, tIC2
real*8 tW1, tW2, tIW1, tIW2
real*8 tIC, tIW
real*8 tauW
real*8 GABCD, GCDAB

integer i, j
integer ip_W
integer AP_Atoms
integer AP_2CFunctions
integer A2AP
real*8 IAB
real*8 GA, GAB
real*8 VA, VAB
real*8 DAB
real*8 CAB_A, CAB_B, CAB_AB
ip_W(i) = iWork(ip_WBlkP-1+i)
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)
AP_2CFunctions(i,j) = iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
A2AP(i,j) = iWork(ip_A2AP-1+2*(j-1)+i)
IAB(i) = Work(ip_IDiag_Sm-1+i)
GA(i) = Work(ip_GDiag_1C_Sm-1+i)
GAB(i) = Work(ip_GDiag_2C_Sm-1+i)
VA(i) = Work(ip_VNrm-1+i)
VAB(i) = Work(ip_VNrm-1+nAtom+i)
DAB(i) = Work(ip_DNrm-1+i)
CAB_A(i) = Work(ip_CNorm+4*(i-1)+1)
CAB_B(i) = Work(ip_CNorm+4*(i-1)+2)
CAB_AB(i) = Work(ip_CNorm+4*(i-1)+3)

! Get number of atoms
nAtom = LDF_nAtom()

! Set up prescreening info
if (nD == 1) then
  l_VNrm = 0
  ip_VNrm = ip_VNorm
  l_DNrm = 0
  ip_DNrm = ip_DNorm
else
  l_VNrm = nAtom+NumberOfAtomPairs
  call GetMem('VNrm','Allo','Real',ip_VNrm,l_VNrm)
  call Cho_dZero(Work(ip_VNrm),l_VNrm)
  do iD=0,nD-1
    ip = ip_VNorm+(nAtom+NumberOfAtomPairs)*iD
    do AB=0,nAtom+NumberOfAtomPairs-1
      Work(ip_VNrm+AB) = max(Work(ip_VNrm+AB),Work(ip+AB))
    end do
  end do
  l_DNrm = NumberOfAtomPairs
  call GetMem('DNrm','Allo','Real',ip_DNrm,l_DNrm)
  call Cho_dZero(Work(ip_DNrm),l_DNrm)
  do iD=0,nD-1
    ip = ip_DNorm+NumberOfAtomPairs*iD
    do AB=0,NumberOfAtomPairs-1
      Work(ip_DNrm+AB) = max(Work(ip_DNrm+AB),Work(ip+AB))
    end do
  end do
end if
l_tauW = nAtom+NumberOfAtomPairs
call GetMem('tauW','Allo','Real',ip_tauW,l_tauW)
if (tau(2) > 0.0d0) then
  call LDF_SetA2AP()
  do A=1,nAtom
    l = A2AP(1,A)
    ip = A2AP(2,A)-1
    ip_tauWA = ip_tauW-1+A
    Work(ip_tauWA) = 0.0d0
    do jAB=1,l
      AB = iWork(ip+jAB)
      if (AP_Atoms(1,AB) == A) then
        Work(ip_tauWA) = max(Work(ip_tauWA),CAB_A(AB))
      else if (AP_Atoms(2,AB) == A) then
        Work(ip_tauWA) = max(Work(ip_tauWA),CAB_B(AB))
      else
        call WarningMessage(2,SecNam//': logical error [A2AP]')
        call LDF_Quit(1)
      end if
    end do
    if (Work(ip_tauWA) > 1.0d-16) then
      tauW = tau(2)/Work(ip_tauWA)
    else
      tauW = 1.0d99
    end if
    Work(ip_tauWA) = tauW
  end do
  call LDF_UnsetA2AP()
  do AB=1,NumberOfAtomPairs
    if (AP_2CFunctions(1,AB) > 0) then
      if (CAB_AB(AB) > 1.0d-16) then
        Work(ip_tauW-1+nAtom+AB) = tau(2)/CAB_AB(AB)
      else
        Work(ip_tauW-1+nAtom+AB) = 1.0d99
      end if
    else
      Work(ip_tauW-1+nAtom+AB) = 1.0d99
    end if
  end do
else
  call Cho_dZero(Work(ip_tauW),l_tauW)
end if

! Allocate and initialize W intermediates
l_WBlkP = nD
call GetMem('WBlkP','Allo','Inte',ip_WBlkP,l_WBlkP)
do iD=1,nD
  call LDF_AllocateAuxBasVector('Win',iWork(ip_WBlkP-1+iD))
  call LDF_ZeroAuxBasVector(ip_W(iD))
end do

!======================
! 3-index contributions
!======================

if ((Mode == 1) .or. (Mode == 3)) then
  if (Timing) then
    tIC = 0.0d0
    tIW = 0.0d0
    call CWTime(tC1,tW1)
  end if
  call Init_Tsk(TaskListID,NumberOfAtomPairs)
  do while (Rsv_Tsk(TaskListID,AB))
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do C=1,nAtom
      tauW = Work(ip_tauW-1+C)
      if ((IAB(AB)*GA(C)*VA(C) >= tau(2)) .or. (IAB(AB)*GA(C)*DAB(AB) >= tauW)) then
        ! Compute integrals (u_A v_B|J_C)
        M = LDF_nBasAux_Atom(C)
        l_Int = nuv*M
        call GetMem('Fck3Int1','Allo','Real',ip_Int,l_Int)
        if (Timing) call CWTime(tIC1,tIW1)
        call LDF_Compute3IndexIntegrals_1(AB,C,tau(1),l_Int,Work(ip_Int))
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
          call dGeMV_('N',nuv,M,FactC(iD),Work(ip_Int),nuv,Work(ipV),1,1.0d0,Work(ipF),1)
        end do
        ! Compute W contribution
        ! W(J_C)+=sum_u_Av_B (u_A v_B|J_C)*D(u_A v_B)
        do iD=1,nD
          ipD = iWork(ip_DBlocks(iD)-1+AB)
          ipW = iWork(ip_W(iD)-1+C)
          call dGeMV_('T',nuv,M,1.0d0,Work(ip_Int),nuv,Work(ipD),1,1.0d0,Work(ipW),1)
        end do
        call GetMem('Fck3Int1','Free','Real',ip_Int,l_Int)
      end if
    end do
    if (LDF2) then
      ! Contributions from two-center aux functions
      do CD=1,NumberOfAtomPairs
        MCD = AP_2CFunctions(1,CD)
        if (MCD > 0) then
          tauW = Work(ip_tauW-1+nAtom+CD)
          if ((IAB(AB)*GAB(CD)*VAB(CD) >= tau(2)) .or. (IAB(AB)*GAB(CD)*DAB(AB) >= tauW)) then
            ! Compute integrals (u_A v_B|J_CD)
            l_Int = nuv*MCD
            call GetMem('Fck3Int2','Allo','Real',ip_Int,l_Int)
            if (Timing) call CWTime(tIC1,tIW1)
            call LDF_Compute3IndexIntegrals_2(AB,CD,tau(1),l_Int,Work(ip_Int))
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
              call dGeMV_('N',nuv,MCD,FactC(iD),Work(ip_Int),nuv,Work(ipV),1,1.0d0,Work(ipF),1)
            end do
            ! Compute W contribution
            ! W(J_CD)+=sum_u_Av_B (u_A v_B|J_CD)*D(u_A v_B)
            do iD=1,nD
              ipD = iWork(ip_DBlocks(iD)-1+AB)
              ipW = iWork(ip_W(iD)-1+nAtom+CD)
              call dGeMV_('T',nuv,MCD,1.0d0,Work(ip_Int),nuv,Work(ipD),1,1.0d0,Work(ipW),1)
            end do
            call GetMem('Fck3Int2','Free','Real',ip_Int,l_Int)
          end if
        end if
      end do
    end if
  end do
  call Free_Tsk(TaskListID)
  if (Timing) then
    call CWTime(tC2,tW2)
    write(6,'(A,2(1X,F12.2),A)') 'Time spent on 3-index contributions:              ',tC2-tC1,tW2-tW1,' seconds'
    write(6,'(A,2(1X,F12.2),A)') '      - of which integrals required:              ',tIC,tIW,' seconds'
  end if
end if

!======================
! 2-index contributions
!======================

if ((Mode == 1) .or. (Mode == 2)) then
  if (Timing) then
    tIC = 0.0d0
    tIW = 0.0d0
    call CWTime(tC1,tW1)
  end if
  if (Mode == 1) then
    Const = -1.0d0
  else
    Const = 1.0d0
  end if
  call Init_Tsk(TaskListID,nAtom)
  do while (Rsv_Tsk(TaskListID,A))
    MA = LDF_nBasAux_Atom(A)
    if (MA > 0) then
      do B=1,A-1
        MB = LDF_nBasAux_Atom(B)
        if (MB > 0) then
          if ((GA(A)*GA(B)*VA(B) >= Work(ip_tauW-1+A)) .or. (GA(A)*GA(B)*VA(A) >= Work(ip_tauW-1+B))) then
            ! Compute integrals (J_A|K_B)
            l_Int = MA*MB
            call GetMem('Fck2Int11','Allo','Real',ip_Int,l_Int)
            if (Timing) call CWTime(tIC1,tIW1)
            call LDF_Compute2IndexIntegrals_11(A,B,tau(1),l_Int,Work(ip_Int))
            if (Timing) then
              call CWTime(tIC2,tIW2)
              tIC = tIC+(tIC2-tIC1)
              tIW = tIW+(tIW2-tIW1)
            end if
            ! Compute W contribution
            ! W(J_A)+=Const*sum_K_B (J_A|K_B)*V(K_B)
            do iD=1,nD
              ipV = iWork(ip_V(iD)-1+B)
              ipW = iWork(ip_W(iD)-1+A)
              call dGeMV_('N',MA,MB,Const,Work(ip_Int),MA,Work(ipV),1,1.0d0,Work(ipW),1)
            end do
            ! Compute W contribution
            ! W(K_B)+=Const*sum_J_A (J_A|K_B)*V(J_A)
            do iD=1,nD
              ipV = iWork(ip_V(iD)-1+A)
              ipW = iWork(ip_W(iD)-1+B)
              call dGeMV_('T',MA,MB,Const,Work(ip_Int),MA,Work(ipV),1,1.0d0,Work(ipW),1)
            end do
            call GetMem('Fck2Int11','Free','Real',ip_Int,l_Int)
          end if
        end if
      end do
      if (GA(A)*GA(A)*VA(A) >= Work(ip_tauW-1+A)) then
        ! Compute integrals (J_A|K_A)
        l_Int = MA**2
        call GetMem('Fck2Int11','Allo','Real',ip_Int,l_Int)
        if (Timing) call CWTime(tIC1,tIW1)
        call LDF_Compute2IndexIntegrals_11(A,A,tau(1),l_Int,Work(ip_Int))
        if (Timing) then
          call CWTime(tIC2,tIW2)
          tIC = tIC+(tIC2-tIC1)
          tIW = tIW+(tIW2-tIW1)
        end if
        ! Compute W contribution
        ! W(J_A)+=Const*sum_K_A (J_A|K_A)*V(K_A)
        do iD=1,nD
          ipV = iWork(ip_V(iD)-1+A)
          ipW = iWork(ip_W(iD)-1+A)
          call dGeMV_('N',MA,MA,Const,Work(ip_Int),MA,Work(ipV),1,1.0d0,Work(ipW),1)
        end do
        call GetMem('Fck2Int11','Free','Real',ip_Int,l_Int)
      end if
      if (LDF2) then
        ! Two-center contributions
        do CD=1,NumberOfAtomPairs
          MCD = AP_2CFunctions(1,CD)
          if (MCD > 0) then
            if ((GA(A)*GAB(CD)*VAB(CD) >= Work(ip_tauW-1+A)) .or. (GA(A)*GAB(CD)*VA(A) >= Work(ip_tauW-1+nAtom+CD))) then
              ! Compute integrals (J_A|K_CD)
              l_Int = MA*MCD
              call GetMem('Fck2Int12','Allo','Real',ip_Int,l_Int)
              if (Timing) call CWTime(tIC1,tIW1)
              call LDF_Compute2IndexIntegrals_12(A,CD,tau(1),l_Int,Work(ip_Int))
              if (Timing) then
                call CWTime(tIC2,tIW2)
                tIC = tIC+(tIC2-tIC1)
                tIW = tIW+(tIW2-tIW1)
              end if
              ! Compute W contribution
              ! W(J_A)+=Const*sum_K_CD (J_A|K_CD)*V(K_CD)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+nAtom+CD)
                ipW = iWork(ip_W(iD)-1+A)
                call dGeMV_('N',MA,MCD,Const,Work(ip_Int),MA,Work(ipV),1,1.0d0,Work(ipW),1)
              end do
              ! Compute W contribution
              ! W(K_CD)+=Const*sum_J_A (J_A|K_CD)*V(J_A)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+A)
                ipW = iWork(ip_W(iD)-1+nAtom+CD)
                call dGeMV_('T',MA,MCD,Const,Work(ip_Int),MA,Work(ipV),1,1.0d0,Work(ipW),1)
              end do
              call GetMem('Fck2Int12','Free','Real',ip_Int,l_Int)
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
      MAB = AP_2CFunctions(1,AB)
      if (MAB > 0) then
        do CD=1,AB-1
          MCD = AP_2CFunctions(1,CD)
          if (MCD > 0) then
            GABCD = GAB(AB)*GAB(CD)*VAB(CD)
            GCDAB = GAB(AB)*GAB(CD)*VAB(AB)
            if ((GABCD >= Work(ip_tauW-1+nAtom+AB)) .or. (GCDAB >= Work(ip_tauW-1+nAtom+CD))) then
              ! Compute integrals (J_AB|K_CD)
              l_Int = MAB*MCD
              call GetMem('Fck2Int22','Allo','Real',ip_Int,l_Int)
              if (Timing) call CWTime(tIC1,tIW1)
              call LDF_Compute2IndexIntegrals_22(AB,CD,tau(1),l_Int,Work(ip_Int))
              if (Timing) then
                call CWTime(tIC2,tIW2)
                tIC = tIC+(tIC2-tIC1)
                tIW = tIW+(tIW2-tIW1)
              end if
              ! Compute W contribution
              ! W(J_AB)+=Const*sum_K_CD (J_AB|K_CD)*V(K_CD)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+nAtom+CD)
                ipW = iWork(ip_W(iD)-1+nAtom+AB)
                call dGeMV_('N',MAB,MCD,Const,Work(ip_Int),MAB,Work(ipV),1,1.0d0,Work(ipW),1)
              end do
              ! Compute W contribution
              ! W(K_CD)+=Const*sum_J_AB (J_AB|K_CD)*V(J_AB)
              do iD=1,nD
                ipV = iWork(ip_V(iD)-1+nAtom+AB)
                ipW = iWork(ip_W(iD)-1+nAtom+CD)
                call dGeMV_('T',MAB,MCD,Const,Work(ip_Int),MAB,Work(ipV),1,1.0d0,Work(ipW),1)
              end do
              call GetMem('Fck2Int22','Free','Real',ip_Int,l_Int)
            end if
          end if
        end do
        ! Compute integrals (J_AB|K_AB)
        if (GAB(AB)*GAB(AB)*VAB(AB) >= Work(ip_tauW-1+nAtom+AB)) then
          l_Int = MAB**2
          call GetMem('Fck2Int22','Allo','Real',ip_Int,l_Int)
          if (Timing) call CWTime(tIC1,tIW1)
          call LDF_Compute2IndexIntegrals_22(AB,AB,tau(1),l_Int,Work(ip_Int))
          if (Timing) then
            call CWTime(tIC2,tIW2)
            tIC = tIC+(tIC2-tIC1)
            tIW = tIW+(tIW2-tIW1)
          end if
          ! Compute W contribution
          ! W(J_AB)+=Const*sum_K_AB (J_AB|K_AB)*V(K_AB)
          do iD=1,nD
            ipV = iWork(ip_V(iD)-1+nAtom+AB)
            ipW = iWork(ip_W(iD)-1+nAtom+AB)
            call dGeMV_('N',MAB,MAB,Const,Work(ip_Int),MAB,Work(ipV),1,1.0d0,Work(ipW),1)
          end do
          call GetMem('Fck2Int22','Free','Real',ip_Int,l_Int)
        end if
      end if
    end do
    call Free_Tsk(TaskListID)
  end if
  if (Timing) then
    call CWTime(tC2,tW2)
    write(6,'(A,2(1X,F12.2),A)') 'Time spent on 2-index contributions:              ',tC2-tC1,tW2-tW1,' seconds'
    write(6,'(A,2(1X,F12.2),A)') '      - of which integrals required:              ',tIC,tIW,' seconds'
  end if
end if

! Deallocate prescreening info
call GetMem('tauW','Free','Real',ip_tauW,l_tauW)
if (l_DNrm > 0) then
  call GetMem('DNrm','Free','Real',ip_DNrm,l_DNrm)
end if
if (l_VNrm > 0) then
  call GetMem('VNrm','Free','Real',ip_VNrm,l_VNrm)
end if

#ifdef _MOLCAS_MPP_
! Add W over nodes
if ((nProcs > 1) .and. Is_Real_Par()) then
  if (Timing) call CWTime(tC1,tW1)
  do iD=1,nD
    call LDF_P_AddAuxBasVector(ip_W(iD))
  end do
  if (Timing) then
    call CWTime(tC2,tW2)
    write(6,'(A,2(1X,F12.2),A)') 'Parallel overhead for W intermediates:            ',tC2-tC1,tW2-tW1,' seconds'
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
  A = AP_Atoms(1,AB)
  B = AP_Atoms(2,AB)
  nuv = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
  l_C = nuv*LDF_nBasAux_Pair_wLD(AB)
  call GetMem('FckCoef','Allo','Real',ip_C,l_C)
  call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l_C)
  ! Compute F(u_A v_B) += FactC*sum_J C(u_A v_B,J)*W(J)
  ipC = ip_C
  MA = LDF_nBasAux_Atom(A)
  if (MA > 0) then
    do iD=1,nD
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      ipW = iWork(ip_W(iD)-1+A)
      call dGeMV_('N',nuv,MA,FactC(iD),Work(ipC),nuv,Work(ipW),1,1.0d0,Work(ipF),1)
    end do
    ipC = ipC+nuv*MA
  end if
  if (B /= A) then
    MB = LDF_nBasAux_Atom(B)
    if (MB > 0) then
      do iD=1,nD
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        ipW = iWork(ip_W(iD)-1+B)
        call dGeMV_('N',nuv,MB,FactC(iD),Work(ipC),nuv,Work(ipW),1,1.0d0,Work(ipF),1)
      end do
      ipC = ipC+nuv*MB
    end if
  end if
  MAB = AP_2CFunctions(1,AB)
  if (MAB > 0) then
    do iD=1,nD
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      ipW = iWork(ip_W(iD)-1+nAtom+AB)
      call dGeMV_('N',nuv,MAB,FactC(iD),Work(ipC),nuv,Work(ipW),1,1.0d0,Work(ipF),1)
    end do
  end if
  call GetMem('FckCoef','Free','Real',ip_C,l_C)
end do
call Free_Tsk(TaskListID)
if (Timing) then
  call CWTime(tC2,tW2)
  write(6,'(A,2(1X,F12.2),A)') 'Time spent computing C*W contribution:            ',tC2-tC1,tW2-tW1,' seconds'
end if

! Deallocate W intermediates
do iD=0,nD-1
  call LDF_DeallocateAuxBasVector('Win',iWork(ip_WBlkP+iD))
end do
call GetMem('WBlkP','Free','Inte',ip_WBlkP,l_WBlkP)

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
    write(6,'(A,2(1X,F12.2),A)') 'Time spent computing XIDI corrections:            ',tC2-tC1,tW2-tW1,' seconds'
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
    write(6,'(A,2(1X,F12.2),A)') 'Parallel overhead for F blocks:                   ',tC2-tC1,tW2-tW1,' seconds'
  end if
end if
#endif

if (Timing) call xFlush(6)

end subroutine LDF_Fock_CoulombOnly_
