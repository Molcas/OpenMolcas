************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992,2007, Roland Lindh                                *
************************************************************************
      SubRoutine PGet0(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                 Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                 n1,n2,n3,n4,MemPSO,ipMem2,
     &                 iShell_A,iShell_B,iShell_C,iShell_D,nQuad,
     &                 PMax)
************************************************************************
*                                                                      *
* Object: to act as a shell towards the manipulations of generating or *
*         accessing the 2nd order density matrix.                      *
*                                                                      *
* Called from: Twoel                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              PGet1                                                   *
*              PGet2                                                   *
*              PGet3                                                   *
*              PGet4                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '92.                                             *
*                                                                      *
*             Modified for RI Feb. 2007                                *
************************************************************************
      use aces_stuff
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
#include "print.fh"
#include "real.fh"
#include "shinf.fh"
#include "setup.fh"
#include "WrkSpc.fh"
#include "etwas.fh"
#include "columbus_gamma.fh"
      Real*8 PSO(ijkl,nPSO)
      Integer iShell(4), iAO(4), iCmp(4), kOp(4), iAOst(4)
      Logical Shijij
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      iRout = 248
      iPrint = nPrint(iRout)
      Call qEnter('PGet0')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
!      write(*,*)"Print out in integral_util/pget0 starting"
!      Call RecPrt('DSO in PGet0',' ',D0,ndens,5)  ! ====== yma ======

      PMax=One
      nSA=1
*
      ipPam=ipMem2
      ipiPam = ipPam + MemPSO
      ipMap = ipiPam + n1+n2+n3+n4
      ipC = ipMap + 4*nDim
      ipS1 = ipC + nCred
      ipS2 = ipS1 + 2*nScr1
*
*     Correction to get correct types in subsequent calls
*
      ipiPam = ip_of_iWork_d(Work(ipiPam))
      ipMap = ip_of_iWork_d(Work(ipMap))
*                                                                      *
************************************************************************
*                                                                      *
*       RASSCF wavefunction
*
        If (lPSO) Then
         If (lSA) nSA=4
         If (nIrrep.eq.1) Then
            kOp(1) = 0
            kOp(2) = 0
            kOp(3) = 0
            kOp(4) = 0
            If (Case_2C) Then
               If (Do_RI) Then
                  Call PGet1_RI2(PSO,ijkl,nPSO,iCmp,
     &                           iShell,iAO,iAOst,Shijij,
     &                           iBas,jBas,kBas,lBas,kOp,
     &                           ExFac,CoulFac,PMax,
     &                           Work(ip_V_K),Work(ip_U_K),nV_K,
     &                           Z_p_k,nSA)
!                  write(*,*)"PGet1_RI2 ===========" yma
               Else
*Not modified yet
                  Call Abend()
               End If
            Else If (Case_3C) Then
               If (Do_RI) Then
                  Call PGet1_RI3(PSO,ijkl,nPSO,iCmp,
     &                           iShell,iAO,iAOst,Shijij,
     &                           iBas,jBas,kBas,lBas,kOp,D0,
     &                           DS,DVar,nDens,
     &                           ExFac,CoulFac,PMax,
     &                           Work(ip_V_K),Work(ip_U_K),nV_K,
     &                           Z_p_k,nnP(0),nSA,nAsh)
               Else
*Not modified yet
                  Call Abend()
               End If
            Else
               Call PGet3(PSO,ijkl,nPSO,iCmp,
     &                    iShell,iAO,iAOst,Shijij,
     &                    iBas,jBas,kBas,lBas,kOp,D0,nDens,
     &                    Work(ipPam),n1,n2,n3,n4,
     &                    iWork(ipiPam),iWork(ipMap),nDim,
     &                    Work(ipC),nCred,Work(ipS1),nScr1,
     &                    Work(ipS2),nScr2,PMax)
!yma                  write(*,*)"PGet3 ==========="
            End If
         Else
            If (Case_2C) Then
               If (Do_RI) Then
                  Call PGet2_RI2(iCmp,iShell,
     &                           iBas,jBas,kBas,lBas,
     &                           Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                           ExFac,CoulFac,
     &                           PMax,Work(ip_V_K),Work(ip_U_K),nV_K,
     &                           Z_p_k,nSA,nZ_p_k)
!yma                  write(*,*)"PGet2_RI2 ==========="
               Else
*Not modified yet
                  Call Abend()
               EndIf
            Else If (Case_3C) Then
               If (Do_RI) Then
                  Call PGet2_RI3(iCmp,iShell,
     &                           iBas,jBas,kBas,lBas,
     &                           Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                           D0,DS,nDens,ExFac,
     &                           CoulFac,PMax,Work(ip_V_K),nV_K,
     &                           Z_p_k,nSA,nAsh)

               Else
*Not modified yet
                  Call Abend()
               EndIf
            Else
            do i=1,ijkl  !yma for testing
              do j=1,nPSO
                PSO(i,j)=0.0d0
              end do
            end do



!      write(*,*)"Print out in integral_util/pget0 before"
!      Call RecPrt('DSO in PGet0',' ',D0,ndens,5)  ! ====== yma ======

               Call PGet4(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                    Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                    D0,nDens,
     &                    Work(ipPam),n1,n2,n3,n4,
     &                    iWork(ipiPam),iWork(ipMap),nDim,
     &                    Work(ipC),nCred,Work(ipS1),nScr1,
     &                    Work(ipS2),nScr2,PMax)
!yma                  write(*,*)"PGet4 ============"
            End If
         End If
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        SCF and DFT wavefunction
*
         iComp = 1
         If (nIrrep.eq.1) Then
            kOp(1) = 0
            kOp(2) = 0
            kOp(3) = 0
            kOp(4) = 0
            If (Gamma_On) Then
               If (Do_RI) Call Abend()
               ipGamma=ipMem2
               If (gamma_mrcisd) Then
               Call Read_Bin_Columbus
     &                (iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Work(ipGamma),nGamma,
     &                       LuGamma,Bin,lBin)
               Else
               Call Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Work(ipGamma),nGamma,
     &                       LuGamma,Bin,lBin)
               Endif
               Call PGet1_Aces(PSO,ijkl,nPSO,iCmp,
     &                         iShell,iAO,iAOst,Shijij,
     &                         iBas,jBas,kBas,lBas,kOp,D0,
     &                         DVar,DS,DSVar,
     &                         nDens,Work(ipGamma),nGamma,
     &                         SO2cI,nSOs,iWork(ipSOsh),PMax)
            Else
               If (Case_2C) Then
                  If (Do_RI) Then
                     Call PGet1_RI2(PSO,ijkl,nPSO,iCmp,
     &                              iShell,iAO,iAOst,Shijij,
     &                              iBas,jBas,kBas,lBas,kOp,
     &                              ExFac,CoulFac,PMax,
     &                              Work(ip_V_K),Work(ip_U_K),nV_K,
     &                              Z_p_k, nSA)
                  Else
                     Call PGet1_CD2(PSO,ijkl,nPSO,iCmp,
     &                              iShell,iAO,iAOst,Shijij,
     &                              iBas,jBas,kBas,lBas,kOp,
     &                              ExFac,CoulFac,PMax,
     &                              Work(ip_V_K),Work(ip_U_K),nV_K,
     &                              Z_p_k,nnP(0))
                  End If
               Else If (Case_3C) Then
                  If (Do_RI) Then
                     Call PGet1_RI3(PSO,ijkl,nPSO,iCmp,
     &                              iShell,iAO,iAOst,Shijij,
     &                              iBas,jBas,kBas,lBas,kOp,D0,
     &                              DS,DVar,nDens,
     &                              ExFac,CoulFac,PMax,
     &                              Work(ip_V_K),Work(ip_U_K),nV_K,
     &                              Z_p_k,nnP(0),nSA,nAsh)
                  Else
                     Call PGet1_CD3(PSO,ijkl,nPSO,iCmp,
     &                              iShell,iAO,iAOst,Shijij,
     &                              iBas,jBas,kBas,lBas,kOp,D0,
     &                              DS,DVar,nDens,
     &                              ExFac,CoulFac,PMax,
     &                              Work(ip_V_K),Work(ip_U_K),nV_K)
                  End If
               Else
                  Call PGet1(PSO,ijkl,nPSO,iCmp,
     &                       iShell,iAO,iAOst,Shijij,
     &                       iBas,jBas,kBas,lBas,kOp,D0,
     &                       DS,nDens,ExFac,CoulFac,PMax)
               End If
            End If
         Else
            If (Gamma_On) Then
               If (Do_RI) Call Abend()
               ipGamma=ipMem2
               If (gamma_mrcisd) Then
               Call Read_Bin_Columbus
     &                (iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Work(ipGamma),nGamma,
     &                       LuGamma,Bin,lBin)
               else
               Call Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Work(ipGamma),nGamma,
     &                       LuGamma,Bin,lBin)
               endif
               Call PGet2_Aces(iCmp,iShell,
     &                         iBas,jBas,kBas,lBas,
     &                         Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                         D0,DVar,DS,
     &                         DSVar,nDens, Work(ipGamma),
     &                         nGamma,
     &                         SO2cI,nSOs,iWork(ipSOsh),PMax)
            Else
               If (Case_2C) Then
                  If (Do_RI) Then

                     Call PGet2_RI2(iCmp,iShell,
     &                              iBas,jBas,kBas,lBas,
     &                              Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                              ExFac,CoulFac,
     &                              PMax,Work(ip_V_K),Work(ip_U_K),nV_K,
     &                              Z_p_k, nSA,nZ_p_k)
                  Else
                     Call PGet2_CD2(iCmp,iShell,
     &                              iBas,jBas,kBas,lBas,
     &                              Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                              ExFac,CoulFac,
     &                              PMax,Work(ip_V_K),nV_K)
                  End If
               Else If (Case_3C) Then
                  If (Do_RI) Then
                     Call PGet2_RI3(iCmp,iShell,
     &                              iBas,jBas,kBas,lBas,
     &                              Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                              D0,DS,nDens,ExFac,
     &                              CoulFac,PMax,Work(ip_V_K),nV_K,
     &                              Z_p_k,nSA,nAsh)
                  Else
                     Call PGet2_CD3(iCmp,iShell,
     &                              iBas,jBas,kBas,lBas,
     &                              Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                              D0,DS,nDens,ExFac,
     &                              CoulFac,PMax,Work(ip_V_K),nV_K)
                  End If
*
               Else
                  Call PGet2(iCmp,iShell,
     &                       iBas,jBas,kBas,lBas,
     &                       Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                       D0,DS,nDens,ExFac,CoulFac,
     &                       PMax)
               End If
            End If
         End If
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      If (iPrint.ge.99) Call RecPrt('PSO in PGet0',' ',
     &                               PSO,ijkl,nPSO)
      Call qExit('PGet0')
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
