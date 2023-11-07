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
      SubRoutine PGet0(iCmp,iBas,jBas,kBas,lBas,
     &                 Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                 n1,n2,n3,n4,MemPSO,Mem2,nMem2,
     &                 iShell_A,iShell_B,iShell_C,iShell_D,nQuad,PMax)
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
      use setup
      use pso_stuff, only: lPSO, lSA, Case_2C, Case_3C, Gamma_On,
     &                     nGamma, Gamma_MRCISD, Bin, D0, DS, DVar,
     &                     G_Toc, lBin, LuGamma, nDens, nNP, nV_k,
     &                     nZ_p_k, SO2CI, U_K, V_K, Z_P_K, DSVar
      use iSD_data, only: iSO2Sh
      use Sizes_of_Seward, only: S
      use RICD_Info, only: Do_RI
      use Symmetry_Info, only: nIrrep
      use Constants, only: One
      use EtWas, only: nCRED, nScr1, nScr2, CoulFac, ExFac, nAsh
      use mspdft_grad, only: DoGradMSPD
      Implicit None

      Integer iBas, jBas, kBas, lBas, ijkl, nPSO, n1, n2, n3, n4,
     &        MemPSO, nMem2, iShell_A, iShell_B, iShell_C, iShell_D,
     &        nQuad
      Real*8 PMax
      Real*8 PSO(ijkl,nPSO), Mem2(nMem2)
      Integer iAO(4), iCmp(4), kOp(4), iAOst(4)
      Logical Shijij

      Integer nSA, ipPAM, ipiPam, ipC, ipS1, ipS2, i, j, ipMAP
!                                                                      *
!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
!
!     write(*,*)"Print out in integral_util/pget0 starting"
!     Call RecPrt('DSO in PGet0',' ',D0,ndens,5)  ! ====== yma ======

      PMax=One
      nSA=1
!                                                                      *
!***********************************************************************
!                                                                      *
!     RASSCF wavefunction
!
      If (lPSO) Then
!                                                                      *
!***********************************************************************
!                                                                      *
!
         ipPam=1
         ipiPam = ipPam + MemPSO
         ipMap = ipiPam + n1+n2+n3+n4
         ipC = ipMap + 4*S%nDim
         ipS1 = ipC + nCred
         ipS2 = ipS1 + 2*nScr1
!
         If (lSA) nSA=4
         If (DoGradMSPD) nSA=5
!
         If (nIrrep.eq.1) Then
            kOp(1) = 0
            kOp(2) = 0
            kOp(3) = 0
            kOp(4) = 0
            If (Case_2C) Then
               If (Do_RI) Then
                  Call PGet1_RI2(PSO,ijkl,nPSO,iCmp,
     &                           iAO,iAOst,
     &                           jBas,lBas,kOp,
     &                           ExFac,CoulFac,PMax,
     &                           V_K,U_K,nV_K,
     &                           Z_p_k,nSA)
!                  write(*,*)"PGet1_RI2 ===========" yma
               Else
!Not modified yet
                  Call Abend()
               End If
            Else If (Case_3C) Then
               If (Do_RI) Then
                  Call PGet1_RI3(PSO,ijkl,nPSO,iCmp,
     &                           iAO,iAOst,
     &                           jBas,kBas,lBas,kOp,D0,
     &                           DVar,nDens,
     &                           ExFac,CoulFac,PMax,
     &                           V_K,U_K,nV_K,
     &                           Z_p_k,nnP(0),nSA,nAsh)
               Else
!Not modified yet
                  Call Abend()
               End If
            Else
               Call PGet3(PSO,ijkl,nPSO,iCmp,
     &                    iAO,iAOst,Shijij,
     &                    iBas,jBas,kBas,lBas,kOp,
     &                    Mem2(ipPam),n1,n2,n3,n4,
     &                    Mem2(ipiPam),Mem2(ipMap),S%nDim,
     &                    Mem2(ipC),nCred,Mem2(ipS1),nScr1,
     &                    Mem2(ipS2),nScr2,PMax)
!yma                  write(*,*)"PGet3 ==========="
            End If
         Else
            If (Case_2C) Then
               If (Do_RI) Then
                  Call PGet2_RI2(iCmp,
     &                           jBas,lBas,
     &                           iAO, iAOst, ijkl, PSO, nPSO,
     &                           ExFac,CoulFac,
     &                           PMax,V_K,nV_K,
     &                           Z_p_k,nSA,nZ_p_k)
!yma                  write(*,*)"PGet2_RI2 ==========="
               Else
!Not modified yet
                  Call Abend()
               EndIf
            Else If (Case_3C) Then
               If (Do_RI) Then
                  Call PGet2_RI3(iCmp,
     &                           jBas,kBas,lBas,
     &                           iAO, iAOst, ijkl, PSO, nPSO,
     &                           D0,nDens,ExFac,
     &                           CoulFac,PMax,V_K,nV_K,
     &                           Z_p_k,nSA,nAsh)

               Else
!Not modified yet
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

               Call PGet4(iCmp,iBas,jBas,kBas,lBas,
     &                    Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                    Mem2(ipPam),n1,n2,n3,n4,
     &                    Mem2(ipiPam),Mem2(ipMap),S%nDim,
     &                    Mem2(ipC),nCred,Mem2(ipS1),nScr1,
     &                    Mem2(ipS2),nScr2,PMax)
!yma                  write(*,*)"PGet4 ============"
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
      Else
!                                                                      *
!***********************************************************************
!                                                                      *
!        SCF and DFT wavefunction
!
         If (Gamma_On .and. nGamma.gt.nMem2) Then
            Write (6,*) 'pGet0: nGamma.lt.nMem2'
            Call abend()
         End If
!
         If (nIrrep.eq.1) Then
            kOp(1) = 0
            kOp(2) = 0
            kOp(3) = 0
            kOp(4) = 0
            If (Gamma_On) Then
               If (Do_RI) Call Abend()
               If (gamma_mrcisd) Then
               Call Read_Bin_Columbus
     &                (iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Mem2,nGamma,
     &                       LuGamma,Bin,lBin)
               Else
               Call Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Mem2,nGamma,
     &                       LuGamma,Bin,lBin)
               Endif
               Call PGet1_Aces(PSO,ijkl,nPSO,iCmp,
     &                         iAO,iAOst,Shijij,
     &                         iBas,jBas,kBas,lBas,kOp,D0,
     &                         DVar,DS,DSVar,
     &                         nDens,Mem2,nGamma,
     &                         SO2cI,nSOs,iSO2Sh,PMax)
            Else
               If (Case_2C) Then
                  If (Do_RI) Then
                     Call PGet1_RI2(PSO,ijkl,nPSO,iCmp,
     &                              iAO,iAOst,
     &                              jBas,lBas,kOp,
     &                              ExFac,CoulFac,PMax,
     &                              V_K,U_K,nV_K,
     &                              Z_p_k, nSA)
                  Else
                     Call PGet1_CD2(PSO,ijkl,nPSO,iCmp,
     &                              iAO,iAOst,
     &                              iBas,jBas,kBas,lBas,kOp,
     &                              ExFac,CoulFac,PMax,
     &                              V_K,U_K,nV_K,
     &                              Z_p_k,nnP(0))
                  End If
               Else If (Case_3C) Then
                  If (Do_RI) Then
                     Call PGet1_RI3(PSO,ijkl,nPSO,iCmp,
     &                              iAO,iAOst,
     &                              jBas,kBas,lBas,kOp,D0,
     &                              DVar,nDens,
     &                              ExFac,CoulFac,PMax,
     &                              V_K,U_K,nV_K,
     &                              Z_p_k,nnP(0),nSA,nAsh)
                  Else
                     Call PGet1_CD3(PSO,ijkl,nPSO,iCmp,
     &                              iAO,iAOst,
     &                              iBas,jBas,kBas,lBas,kOp,D0,
     &                              DVar,nDens,
     &                              ExFac,CoulFac,PMax,
     &                              V_K,U_K,nV_K)
                  End If
               Else
                  Call PGet1(PSO,ijkl,nPSO,iCmp,
     &                       iAO,iAOst,Shijij,
     &                       iBas,jBas,kBas,lBas,kOp,D0,
     &                       DS,nDens,ExFac,CoulFac,PMax)
               End If
            End If
         Else
            If (Gamma_On) Then
               If (Do_RI) Call Abend()
               If (gamma_mrcisd) Then
               Call Read_Bin_Columbus
     &                (iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Mem2,nGamma,
     &                       LuGamma,Bin,lBin)
               else
               Call Read_Bin(iShell_A,iShell_B,iShell_C,iShell_D,
     &                       G_Toc,nQuad,
     &                       Mem2,nGamma,
     &                       LuGamma,Bin,lBin)
               endif
               Call PGet2_Aces(iCmp,
     &                         iBas,jBas,kBas,lBas,
     &                         Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                         D0,DVar,DS,
     &                         DSVar,nDens,Mem2,nGamma,
     &                         SO2cI,nSOs,iSO2Sh,PMax)
            Else
               If (Case_2C) Then
                  If (Do_RI) Then

                     Call PGet2_RI2(iCmp,
     &                              jBas,lBas,
     &                              iAO, iAOst, ijkl, PSO, nPSO,
     &                              ExFac,CoulFac,
     &                              PMax,V_K,nV_K,
     &                              Z_p_k, nSA,nZ_p_k)
                  Else
                     Call PGet2_CD2(iCmp,
     &                              iBas,jBas,kBas,lBas,
     &                              iAO, iAOst, ijkl, PSO, nPSO,
     &                              CoulFac,
     &                              PMax,V_K,nV_K)
                  End If
               Else If (Case_3C) Then
                  If (Do_RI) Then
                     Call PGet2_RI3(iCmp,
     &                              jBas,kBas,lBas,
     &                              iAO, iAOst, ijkl, PSO, nPSO,
     &                              D0,nDens,ExFac,
     &                              CoulFac,PMax,V_K,nV_K,
     &                              Z_p_k,nSA,nAsh)
                  Else
                     Call PGet2_CD3(iCmp,
     &                              iBas,jBas,kBas,lBas,
     &                              iAO, iAOst, ijkl, PSO, nPSO,
     &                              D0,nDens,
     &                              CoulFac,PMax,V_K,nV_K)
                  End If
!
               Else
                  Call PGet2(iCmp,
     &                       iBas,jBas,kBas,lBas,
     &                       Shijij, iAO, iAOst, ijkl, PSO, nPSO,
     &                       D0,DS,nDens,ExFac,CoulFac,
     &                       PMax)
               End If
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
      End If
!
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('PSO in PGet0',' ',PSO,ijkl,nPSO)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End SubRoutine PGet0
