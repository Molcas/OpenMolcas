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
! Copyright (C) 1992,2007, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PGet0(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,n1,n2,n3,n4,MemPSO,Mem2,nMem2,iShell_A,iShell_B,iShell_C, &
                 iShell_D,nQuad,PMax)
!***********************************************************************
!                                                                      *
! Object: to act as a shell towards the manipulations of generating or *
!         accessing the 2nd order density matrix.                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '92.                                             *
!                                                                      *
!             Modified for RI Feb. 2007                                *
!***********************************************************************

use setup, only: nSOs
use pso_stuff, only: Bin, Case_2C, Case_3C, D0, DS, DSVar, DVar, G_Toc, Gamma_MRCISD, Gamma_On, lBin, lPSO, lSA, LuGamma, nDens, &
                     nGamma, nNP, nV_k, nZ_p_k, SO2CI, U_K, V_K, Z_P_K
use iSD_data, only: iSO2Sh
use Sizes_of_Seward, only: S
use RICD_Info, only: Do_RI
use Symmetry_Info, only: nIrrep
use EtWas, only: CoulFac, ExFac, nAsh, nCRED, nScr1, nScr2
use mspdft_grad, only: DoGradMSPD
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iBas, jBas, kBas, lBas, iAO(4), iAOst(4), ijkl, nPSO, n1, n2, n3, n4, MemPSO, nMem2, &
                                 iShell_A, iShell_B, iShell_C, iShell_D, nQuad
real(kind=wp), intent(out) :: PSO(ijkl,nPSO), Mem2(nMem2), PMax
integer(kind=iwp) :: ipC, ipiPam, ipMAP, ipPAM, ipS1, ipS2, kOp(4), nSA

!                                                                      *
!***********************************************************************
!                                                                      *

!write(u6,*) 'Print out in integral_util/pget0 starting'
!call RecPrt('DSO in PGet0',' ',D0,ndens,5)  ! ====== yma ======

PMax = One
nSA = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! RASSCF wavefunction

if (lPSO) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ipPam = 1
  ipiPam = ipPam+MemPSO
  ipMap = ipiPam+n1+n2+n3+n4
  ipC = ipMap+4*S%nDim
  ipS1 = ipC+nCred
  ipS2 = ipS1+2*nScr1

  if (lSA) nSA = 4
  if (DoGradMSPD) nSA = 5

  if (nIrrep == 1) then
    kOp(1) = 0
    kOp(2) = 0
    kOp(3) = 0
    kOp(4) = 0
    if (Case_2C) then
      if (Do_RI) then
        call PGet1_RI2(PSO,ijkl,nPSO,iCmp,iAO,iAOst,jBas,lBas,kOp,ExFac,CoulFac,PMax,V_K,U_K,nV_K,Z_p_k,nSA)
        !write(u6,*) 'PGet1_RI2 ===========' ! yma
      else
        ! Not modified yet
        call Abend()
      end if
    else if (Case_3C) then
      if (Do_RI) then
        call PGet1_RI3(PSO,ijkl,nPSO,iCmp,iAO,iAOst,jBas,kBas,lBas,kOp,D0,DVar,nDens,ExFac,CoulFac,PMax,V_K,U_K,nV_K,Z_p_k,nnP(0), &
                       nSA,nAsh)
      else
        ! Not modified yet
        call Abend()
      end if
    else
      call PGet3(PSO,ijkl,nPSO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,Mem2(ipPam),n1,n2,n3,n4,Mem2(ipiPam),Mem2(ipMap), &
                 S%nDim,Mem2(ipC),nCred,Mem2(ipS1),nScr1,Mem2(ipS2),nScr2,PMax)
      !yma write(u6,*) 'PGet3 ==========='
    end if
  else
    if (Case_2C) then
      if (Do_RI) then
        call PGet2_RI2(iCmp,jBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,ExFac,CoulFac,PMax,V_K,nV_K,Z_p_k,nSA,nZ_p_k)
        !yma write(u6,*) 'PGet2_RI2 ==========='
      else
        ! Not modified yet
        call Abend()
      end if
    else if (Case_3C) then
      if (Do_RI) then
        call PGet2_RI3(iCmp,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,D0,nDens,ExFac,CoulFac,PMax,V_K,nV_K,Z_p_k,nSA,nAsh)
      else
        ! Not modified yet
        call Abend()
      end if
    else
      PSO(:,:) = Zero  !yma for testing

      !write(u6,*) 'Print out in integral_util/pget0 before'
      !call RecPrt('DSO in PGet0',' ',D0,ndens,5)  ! ====== yma ======

      call PGet4(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,Mem2(ipPam),n1,n2,n3,n4,Mem2(ipiPam),Mem2(ipMap),S%nDim, &
                 Mem2(ipC),nCred,Mem2(ipS1),nScr1,Mem2(ipS2),nScr2,PMax)
      !yma write(u6,*) 'PGet4 ============'
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! SCF and DFT wavefunction

  if (Gamma_On .and. (nGamma > nMem2)) then
    write(u6,*) 'pGet0: nGamma < nMem2'
    call abend()
  end if

  if (nIrrep == 1) then
    kOp(1) = 0
    kOp(2) = 0
    kOp(3) = 0
    kOp(4) = 0
    if (Gamma_On) then
      if (Do_RI) call Abend()
      if (gamma_mrcisd) then
        call Read_Bin_Columbus(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,nQuad,Mem2,nGamma,LuGamma,Bin,lBin)
      else
        call Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,nQuad,Mem2,nGamma,LuGamma,Bin,lBin)
      end if
      call PGet1_Aces(PSO,ijkl,nPSO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,D0,DVar,DS,DSVar,nDens,Mem2,nGamma,SO2cI,nSOs,iSO2Sh, &
                      PMax)
    else
      if (Case_2C) then
        if (Do_RI) then
          call PGet1_RI2(PSO,ijkl,nPSO,iCmp,iAO,iAOst,jBas,lBas,kOp,ExFac,CoulFac,PMax,V_K,U_K,nV_K,Z_p_k,nSA)
        else
          call PGet1_CD2(PSO,ijkl,nPSO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,ExFac,CoulFac,PMax,V_K,U_K,nV_K,Z_p_k,nnP(0))
        end if
      else if (Case_3C) then
        if (Do_RI) then
          call PGet1_RI3(PSO,ijkl,nPSO,iCmp,iAO,iAOst,jBas,kBas,lBas,kOp,D0,DVar,nDens,ExFac,CoulFac,PMax,V_K,U_K,nV_K,Z_p_k, &
                         nnP(0),nSA,nAsh)
        else
          call PGet1_CD3(PSO,ijkl,nPSO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,D0,DVar,nDens,ExFac,CoulFac,PMax,V_K,U_K,nV_K)
        end if
      else
        call PGet1(PSO,ijkl,nPSO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,D0,DS,nDens,ExFac,CoulFac,PMax)
      end if
    end if
  else
    if (Gamma_On) then
      if (Do_RI) call Abend()
      if (gamma_mrcisd) then
        call Read_Bin_Columbus(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,nQuad,Mem2,nGamma,LuGamma,Bin,lBin)
      else
        call Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,G_Toc,nQuad,Mem2,nGamma,LuGamma,Bin,lBin)
      end if
      call PGet2_Aces(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,D0,DVar,DS,DSVar,nDens,Mem2,nGamma,SO2cI,nSOs,iSO2Sh,PMax)
    else
      if (Case_2C) then
        if (Do_RI) then

          call PGet2_RI2(iCmp,jBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,ExFac,CoulFac,PMax,V_K,nV_K,Z_p_k,nSA,nZ_p_k)
        else
          call PGet2_CD2(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,CoulFac,PMax,V_K,nV_K)
        end if
      else if (Case_3C) then
        if (Do_RI) then
          call PGet2_RI3(iCmp,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,D0,nDens,ExFac,CoulFac,PMax,V_K,nV_K,Z_p_k,nSA,nAsh)
        else
          call PGet2_CD3(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,D0,nDens,CoulFac,PMax,V_K,nV_K)
        end if

      else
        call PGet2(iCmp,iBas,jBas,kBas,lBas,iAO,iAOst,ijkl,PSO,nPSO,D0,DS,nDens,ExFac,CoulFac,PMax)
      end if
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('PSO in PGet0',' ',PSO,ijkl,nPSO)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PGet0
