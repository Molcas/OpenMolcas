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
      SubRoutine PGet1_RI3(PAO,ijkl,nPAO,iCmp,
     &                 iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                 DSO,DSSO,DSO_Var,nDSO,ExFac,CoulFac,PMax,V_K,
     &                 U_K,mV_k,ZpK,nnP1,nSA,nAct)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density.                        *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified for 3-center RI gradients, March 2007           *
************************************************************************
      use Basis_Info, only: nBas
      use SOAO_Info, only: iAOtSO
      use pso_stuff, only: lPSO, lsa, ipAorb, Thpkl
      use ExTerm, only: CijK, CilK, BklK, BMP2, iMP2prpt, LuBVector
      use ExTerm, only: Ymnij, ipYmnij, nYmnij, CMOi
#ifdef _DEBUGPRINT_
      use ExTerm, only: iOff_Ymnij
#endif
      use ExTerm, only: Yij
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "exterm.fh"
#include "WrkSpc.fh"
      Real*8 PAO(ijkl,nPAO), DSO(nDSO,nSA), DSSO(nDSO), V_k(mV_k,nSA),
     &       U_k(mV_k), DSO_Var(nDSO),ZpK(nnP1,mV_K,*)
      Integer iAO(4), kOp(4), iAOst(4), iCmp(4)
      Integer nj(4), jSkip(4), NumOrb(4), nAct(0:7)
      Logical Shijij,Found

      Real*8, Pointer :: Xli(:)=>Null(), Xki(:)=>Null()
      Type V2
        Real*8, Pointer:: A2(:,:)=>Null()
      End Type V2
      Type (V2):: Xli2(2), Xki2(2)
      Type (V2):: Xli3(2), Xki3(2)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      kYmnij(l,iDen)=Ymnij(ipYmnij(iDen)-1+l)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      iComp = 1
      Call PrMtrx('DSO     ',[iD0Lbl],iComp,1,D0)
      Write (6,*)
      Write (6,*) 'Distribution of Ymnij'
      iSym=1
      If (nYmnij(iSym,1).gt.0) Then
         Write (6,*) 'iSym=',iSym
         Do i = iOff_Ymnij(iSym,1)+1, iOff_Ymnij(iSym,1)+nYmnij(iSym,1)
            Write (6,*) 'kYmnij=',kYmnij(i,1)
         End Do
       End If
       Write (6,*) 'jbas,kbas,lbas=',jBas,kBas,lBas
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     DeSymP will compensate for degeneracy due to permutational
*     symmetry. We will have to compensate for that here!
*
      Call CWTime(Cpu1,Wall1)
*
      iOff1 = nBas(0)
      Fac = One / Four
      PMax=Zero
      iPAO = 0
*
      jSym = 1
      kSym = 1
      lSym = iEor(jSym-1,kSym-1)+1
      NumOrb(1) = nChOrb(kSym-1,1)
*
      Call Qpg_iScalar('SCF mode',Found)
      If (Found) Then
         Call Get_iScalar('SCF mode',iUHF) ! either 0 or 1
      Else
         iUHF=0
      EndIf
*
*     Test if we have any exchange contribution of significance
*
*
      ExFac_=ExFac
      If (ExFac.ne.0) Then
*
*        Pick up the number of MOs which passed the threshold test.
*
         nj2=0
         Do iSO=1,nKdens
           jSkip(iSO)=0
           nj(iSO)=nYmnij(jSym,iSO)
           NumOrb(iSO) = nChOrb(kSym-1,iSO)
*
*        If all included skip presceening.
*
!          trick for skipping unnecessary overhead
           If (-nj(iSO).eq.NumOrb(iSO)) Then
              jSkip(iSO)=1
              nj(iSO)=NumOrb(iSO)
           EndIf
*
*        If all excluded process only for Coulombic contributions.
*
           nj2=nj2+nj(iSO)
         End Do
         If ((nj2.eq.0).and.(.not.lPSO))  ExFac=Zero
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (ExFac.ne.Zero .and. NumOrb(1).gt.0 .and. iMP2prpt.ne.2
     &    .and. .not. lPSO  .and. iUHF.eq.0 ) Then
*                                                                      *
************************************************************************
*                                                                      *
*        HF and Hybrid DFT
*
*        number of functions in the kS and lS shell
*
         nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)
*
         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)
*
*        Pointers to the full list of the X_mu,i elements.
*
         lda = SIZE(CMOi(1)%SB(1)%A2,1)
         ik = 1 + lda*(kSO-1)
         il = 1 + lda*(lSO-1)
         Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
         Xli(1:) => CMOi(1)%SB(1)%A1(il:)
*
*        Collect the X_mu,i which survived the prescreening.
*        Replace the pointers above, i.e. Xki, Xli.
*
         If (nj(1).le.NumOrb(1) .and. jSkip(1).eq.0) Then
*
*           Note that the X_mu,i are stored as X_i,mu!
*
            imo=1
            Do k=1,nj(1)
               kmo=kYmnij(k,1) ! CD-MO index
*
*              Pick up X_mu,i for all mu's that belong to shell k
*
               call dcopy_(nKBas,Xki(kmo),NumOrb(1),
     &                           Yij(imo,1,1),nj(1))
*
*              Pick up X_mu,i for all mu's that belong to shell l
*
               call dcopy_(nLBas,Xli(kmo),NumOrb(1),
     &                           Yij(imo,2,1),nj(1))
*
               imo=imo+1
            End Do
*           Reset pointers!
            Xki(1:nj(1)*nKBas) => Yij(1:nj(1)*nKBas,1,1)
            Xli(1:nj(1)*nLBas) => Yij(1:nj(1)*nLBas,2,1)
         ElseIf (nj(1).gt.NumOrb(1)) Then
            Call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
            Call Abend()
         EndIf

         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSO_off = jSO - iOff1
*
*           Read a block of C_kl^J
*
            lCVec = nIJR(kSym,lSym,1)*jBas ! Block size
            iAdr = nIJR(kSym,lSym,1)*(jSO_off-1) + iAdrCVec(jSym,kSym,1)
            Call dDaFile(LuCVector(jSym,1),2,CijK,lCVec,iAdr)
*
*           Extract only those C_kl^Js for which we deem k and l to
*           belong to the shell-pair and to be of significance.
*
            If (nj(1).le.NumOrb(1) .and. jSkip(1).eq.0) Then
               ij=1
               Do j=1,nj(1)
                 jmo=kYmnij(j,1)
                 Do i=1,nj(1)
                   imo=kYmnij(i,1)
                   jC=imo+NumOrb(1)*(jmo-1)
                   call dcopy_(jBas,CijK(jC),NumOrb(1)**2,
     &                              CilK(ij),nj(1)**2)
                   ij=ij+1
                 End Do
               End Do
               n2j=nj(1)**2*jBas
               CijK(1:n2j)=CilK(1:n2j)
            End If
*
*           Transform according to Eq. 16 (step 4) and generate B_kl^J
*
*** ----    E(jK,m) = Sum_i C(i,jK)' * X(i,m)
*
            Call dGEMM_('T','N',nj(1)*jBas,nKBas,nj(1),
     &                   1.0d0,CijK,nj(1),
     &                         Xki,nj(1),
     &                   0.0d0,CilK,nj(1)*jBas)
*
*** ----    B(Km,n) = Sum_j E(j,Km)' * X(j,n)
*
            Call dGEMM_('T','N',jBas*nKBas,nLBas,nj(1),
     &                   1.0d0,CilK,nj(1),
     &                         Xli,nj(1),
     &                   0.0d0,BklK,jBas*nKBas)
*
            Do i3 = 1, iCmp(3)
               kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                  iPAO = iPAO + 1
                  nijkl = 0
*
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl
                     Do kAOk = 0, kBas-1
                        kSOk = kSO + kAOk
*
                        indexB = (kAOk + (i3-1)*kBas)*jBas
     &                         + (lAOl + (i4-1)*lBas)*nKBas*jBas
                        Indk=Max(kSOk,lSOl)
                        Indl=kSOk+lSOl-Indk
                        Indkl=(Indk-1)*Indk/2+Indl
*
                        Do jAOj = 0, jBas-1
                           jSOj = jSO + jAOj - iOff1
                           nijkl = nijkl + 1
                           indexB = indexB + 1
*
*-----------------------Coulomb contribution: V_k(j)*D(kl)
*
                           temp=CoulFac*V_k(jSOj,1)*DSO(Indkl,1)
*
*-----------------------Exchange contribution: B(K,m,n)
*
                           temp = temp - ExFac*Half*BklK(indexB)
*
                           PMax=Max(PMax,Abs(temp))
                           PAO(nijkl,iPAO) =  Fac * temp
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Xki=>Null()
         Xli=>Null()
*                                                                      *
************************************************************************
*                                                                      *
      Else If (ExFac.ne.Zero .and. NumOrb(1).gt.0 .and. iMP2prpt.ne.2
     &    .and. .not. lPSO  .and. iUHF.eq.1 ) Then
*                                                                      *
************************************************************************
*                                                                      *
*        UHF and Hybrid UDFT
*
*        number of functions in the kS and lS shell
*
         nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)
*
         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)
*
*        Pointers to the full list of the X_mu,i elements.
*
         Do iSO=1,2
           If (nIJR(kSym,lSym,iSO).ne.0) Then
*
             Xki2(iSO)%A2(1:,1:) => CMOi(iSO)%SB(1)%A2(1:,kSO:)
             Xli2(iSO)%A2(1:,1:) => CMOi(iSO)%SB(1)%A2(1:,lSO:)
*
*        Collect the X_mu,i which survived the prescreening.
*        Replace the pointers above, i.e. Xki, Xli.
*
             If (nj(iSO).le.NumOrb(iSO) .and. jSkip(iSO).eq.0) Then
*
*           Note that the X_mu,i are stored as X_i,mu!
*
               imo=1
               Do k=1,nj(iSO)
                 kmo=kYmnij(k,iSO) ! CD-MO index
*
*              Pick up X_mu,i for all mu's that belong to shell k
*
                 call dcopy_(nKBas,Xki2(iSO)%A2(kmo,1),NumOrb(iSO),
     &                             Yij(imo,1,iSO),nj(iSO))
*
*              Pick up X_mu,i for all mu's that belong to shell l
*
                 call dcopy_(nLBas,Xli2(iSO)%A2(kmo,1),NumOrb(iSO),
     &                             Yij(imo,2,iSO),nj(iSO))
*
                 imo=imo+1
               End Do
*           Reset pointers!
               Xki2(iSO)%A2(1:nj(iSO),1:nKBas) =>
     &                    Yij(1:nj(iSO)*nKBas,1,iSO)
               Xli2(iSO)%A2(1:nj(iSO),1:nLBas) =>
     &                    Yij(1:nj(iSO)*nLBas,2,iSO)
             ElseIf (nj(iSO).gt.NumOrb(iSO)) Then
               Call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
               Call Abend()
             EndIf
           EndIf
         End Do

         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSO_off = jSO - iOff1
*
            Factor=0.0d0
*
            Do iSO=1,2
               If ((nIJR(kSym,lSym,iSO).ne.0).and.(nj(iSO).ne.0)) Then
*
*           Read a block of C_kl^J
*
                 lCVec = nIJR(kSym,lSym,iSO)*jBas ! Block size
                 iAdr = nIJR(kSym,lSym,iSO)*(jSO_off-1) +
     &                  iAdrCVec(jSym,kSym,iSO)
                 Call dDaFile(LuCVector(jSym,iSO),2,CijK,lCVec,iAdr)
*
*           Extract only those C_kl^Js for which we deem k and l to
*           belong to the shell-pair and to be of significance.
*
                 If (nj(iSO).le.NumOrb(iSO) .and. jSkip(iSO).eq.0) Then
                    ij=1
                    Do j=1,nj(iSO)
                      jmo=kYmnij(j,iSO)
                      Do i=1,nj(iSO)
                        imo=kYmnij(i,iSO)
                        jC=imo+NumOrb(iSO)*(jmo-1)
                        call dcopy_(jBas,CijK(jC),NumOrb(iSO)**2,
     &                                   Cilk(ij),nj(iSO)**2)
                        ij=ij+1
                      End Do
                    End Do
                    n2j=nj(iSO)**2*jBas
                    CijK(1:n2j)=CilK(1:n2j)
                 End If
*
*           Transform according to Eq. 16 (step 4) and generate B_kl^J
*
*** ----    E(jK,m) = Sum_i C(i,jK)' * X(i,m)
*
                 Call dGEMM_('T','N',nj(iSO)*jBas,nKBas,nj(iSO),
     &                        1.0d0,CijK,nj(iSO),
     &                              Xki2(iSO)%A2,nj(iSO),
     &                        0.0d0,CilK,nj(iSO)*jBas)
*
*** ----    B(Km,n) = Sum_j E(j,Km)' * X(j,n)
*
                 Call dGEMM_('T','N',jBas*nKBas,nLBas,nj(iSO),
     &                        1.0d0,CilK,nj(iSO),
     &                              Xli2(iSO)%A2,nj(iSO),
     &                        Factor,BklK,jBas*nKBas)
                 Factor=1.0d0
              EndIf
            End Do
*
            Do i3 = 1, iCmp(3)
               kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                  iPAO = iPAO + 1
                  nijkl = 0
*
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl
                     Do kAOk = 0, kBas-1
                        kSOk = kSO + kAOk
*
                        indexB = (kAOk + (i3-1)*kBas)*jBas
     &                         + (lAOl + (i4-1)*lBas)*nKBas*jBas
                        Indk=Max(kSOk,lSOl)
                        Indl=kSOk+lSOl-Indk
                        Indkl=(Indk-1)*Indk/2+Indl
*
                        Do jAOj = 0, jBas-1
                           jSOj = jSO + jAOj - iOff1
                           nijkl = nijkl + 1
                           indexB = indexB + 1
*
*-----------------------Coulomb contribution: V_k(j)*D(kl)
*
                           temp=CoulFac*V_k(jSOj,1)*DSO(Indkl,1)
*
*-----------------------Exchange contribution: B(K,m,n)
*
                           temp = temp - ExFac*BklK(indexB)
*
                           PMax=Max(PMax,Abs(temp))
                           PAO(nijkl,iPAO) =  Fac * temp
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Do iSO=1,2
             Xki2(iSO)%A2 => Null()
             Xli2(iSO)%A2 => Null()
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      Else If (ExFac.ne.Zero .and. NumOrb(1).gt.0 .and. iMP2prpt.ne.2
     &    .and. lPSO .and. .not. LSA ) Then
*                                                                      *
************************************************************************
*                                                                      *
*        CASSCF case
*
*        number of functions in the kS and lS shell
*
         nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)
*
         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)
*
*        Pointers to the full list of the X_mu,i elements.
*
         lda = SIZE(CMOi(1)%SB(1)%A2,1)
         ik  = 1 + lda*(kSO-1)
         il  = 1 + lda*(lSO-1)
         Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
         Xli(1:) => CMOi(1)%SB(1)%A1(il:)
*
*        Collect the X_mu,i which survived the prescreening.
*        Replace the pointers above, i.e. Xki, Xli.
*
         If (nj(1).le.NumOrb(1) .and. jSkip(1).eq.0.and.nj(1).ne.0) Then
*
*           Note that the X_mu,i are stored as X_i,mu!
*
            imo=1
            Do k=1,nj(1)
               kmo=kYmnij(k,1) ! CD-MO index
*
*              Pick up X_mu,i for all mu's that belong to shell k
*
               call dcopy_(nKBas,Xki(kmo),NumOrb(1),
     &                           Yij(imo,1,1),nj(1))
*
*              Pick up X_mu,i for all mu's that belong to shell l
*
               call dcopy_(nLBas,Xli(kmo),NumOrb(1),
     &                           Yij(imo,2,1),nj(1))
*
               imo=imo+1
            End Do
*           Reset pointers!
            Xki(1:nj(1)*nKBas) => Yij(1:nj(1)*nKBas,1,1)
            Xli(1:nj(1)*nLBas) => Yij(1:nj(1)*nLBas,2,1)
         ElseIf (nj(1).gt.NumOrb(1)) Then
            Call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
            Call Abend()
         EndIf

         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSO_off = jSO - iOff1
*
*           Read a block of C_kl^J
*
            lCVec = nIJR(kSym,lSym,1)*jBas ! Block size
            iAdr = nIJR(kSym,lSym,1)*(jSO_off-1) + iAdrCVec(jSym,kSym,1)
            Call dDaFile(LuCVector(jSym,1),2,CijK,lCVec,iAdr)
*
*           Extract only those C_kl^Js for which we deem k and l to
*           belong to the shell-pair and to be of significance.
*
            If (nj(1).ne.0) Then
              If (nj(1).le.NumOrb(1) .and. jSkip(1).eq.0) Then
                 ij=1
                 Do j=1,nj(1)
                   jmo=kYmnij(j,1)
                   Do i=1,nj(1)
                     imo=kYmnij(i,1)
                     jC=imo+NumOrb(1)*(jmo-1)
                     call dcopy_(jBas,CijK(jC),NumOrb(1)**2,CilK(ij),
     &                          nj(1)**2)
                     ij=ij+1
                   End Do
                 End Do
                 n2j=nj(1)**2*jBas
                 CijK(1:n2j)=CilK(1:n2j)
              End If
*
*             Transform according to Eq. 16 (step 4) and generate B_kl^J
*
*** ----      E(jK,m) = Sum_i C(i,jK)' * X(i,m)
*
              Call dGEMM_('T','N',nj(1)*jBas,nKBas,nj(1),
     &                     1.0d0,CijK,nj(1),
     &                           Xki,nj(1),
     &                     0.0d0,CilK,nj(1)*jBas)
*
*** ----      B(Km,n) = Sum_j E(j,Km)' * X(j,n)
*
              Call dGEMM_('T','N',jBas*nKBas,nLBas,nj(1),
     &                     1.0d0,CilK,nj(1),
     &                           Xli,nj(1),
     &                     0.0d0,BklK,jBas*nKBas)
*
            Else
               BklK(1:jBas*nKBas*nLBas)=Zero
            EndIf
*
**          Active term
*
            Do jAOj=0,jBas-1
              jSOj = jSO_off + jAOj
              Do i4 = 1, iCmp(4)
                lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
                Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl-1
                  lp=ipAOrb(0,1)+lSOl*nAct(lSym-1)
                  Do kAct=1,nAct(kSym-1)
                    tmp=ddot_(kact,Zpk(kAct*(kAct-1)/2+1,jSOj,1),1,
     &                       Work(lp),1)
*
*
                    Do lAct=kAct+1,nAct(lSym-1)
                      tmp=tmp+Zpk(lAct*(lAct-1)/2+kAct,jSOj,1)*
     &                        Work(lp+lAct-1)
                    End Do
                    CilK(kAct)=tmp
                  End Do
*
                  Do i3 = 1, iCmp(3)
                    kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                    kp=ipAOrb(0,1)+(kSO-1)*nAct(kSym-1)
                    iThpkl= jAOj+ (i3-1)*kBas*jBas
     &                         + (lAOl + (i4-1)*lBas)*nKBas*jBas+1
                    Call dGeMV_('T',nAct(kSym-1),kBas,1.0d0,
     &                         Work(kp),
     &                         nAct(kSym-1),Cilk,1,0.0d0,
     &                         Thpkl(iThpkl),jBas)
                  End Do
                End Do
              End Do
            End Do
*
            Do i3 = 1, iCmp(3)
               kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                  iPAO = iPAO + 1
                  nijkl = 0
*
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl
                     Do kAOk = 0, kBas-1
                        kSOk = kSO + kAOk
*
                        iThpkl=(kAOk + (i3-1)*kBas)*jBas
     &                         + (lAOl + (i4-1)*lBas)*nKBas*jBas
                        indexB=iThpkl
                        Indk=Max(kSOk,lSOl)
                        Indl=kSOk+lSOl-Indk
                        Indkl=(Indk-1)*Indk/2+Indl
*
                        Do jAOj = 0, jBas-1
                           jSOj = jSO + jAOj - iOff1
                           nijkl = nijkl + 1
                           indexB = indexB + 1
                           iThpkl = iThpkl + 1
*
*-----------------------Coulomb contribution: V_k(j)*D(kl)
*
                           temp=CoulFac*V_k(jSOj,1)*DSO(Indkl,1)
*
*-----------------------Exchange contribution: B(K,m,n)
*
                           temp = temp - ExFac*Half*Bklk(indexB)
*
*-----------------------Active space contribution: Sum_p Z(p,K)*Th(p,m,n)
*
                           temp=temp+Thpkl(iThpkl)
*
                           PMax=Max(PMax,Abs(temp))
                           PAO(nijkl,iPAO) =  Fac * temp
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Xki=>Null()
         Xli=>Null()
*                                                                      *
************************************************************************
*                                                                      *
      Else If (ExFac.ne.Zero .and. iMP2prpt.ne.2
     &   .and. lPSO .and. lSA ) Then
*                                                                      *
************************************************************************
*                                                                      *
*        SA-CASSCF case
*
*        number of functions in the kS and lS shell
*
         nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)
*
         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)
*
*        Pointers to the full list of the X_mu,i elements.
*
         Do iSO=1,2
           If (nIJR(kSym,lSym,iSO).ne.0) Then
             iMOleft=iSO
             iMOright=iSO+2
*
             Xki2(iSO)%A2(1:,1:) => CMOi(iSO+2)%SB(1)%A2(1:,kSO:)
             Xki3(iSO)%A2(1:,1:) => CMOi(iSO  )%SB(1)%A2(1:,kSO:)
             Xli2(iSO)%A2(1:,1:) => CMOi(iSO  )%SB(1)%A2(1:,lSO:)
             Xli3(iSO)%A2(1:,1:) => CMOi(iSO+2)%SB(1)%A2(1:,lSO:)
*
*            Collect the X_mu,i which survived the prescreening.
*            Replace the pointers above, i.e. Xki, Xli.
*
             If ((nj(iMOright).le.NumOrb(iMOright))
     &           .and.(jSkip(iMOright).eq.0)) Then

*               Note that the X_mu,i are stored as X_i,mu!
*
                imo=1
                Do k=1,nj(iMOright)
                   kmo=kYmnij(k,iMOright) ! CD-MO index
*
*                  Pick up X_mu,i for all mu's that belong to shell k
*
                   call dcopy_(nKBas,
     &                         Xki2(iSO)%A2(kmo,1),NumOrb(iMOright),
     &                         Yij(imo,1,iMOright),nj(iMOright))

                   call dcopy_(nLBas,
     &                         Xli3(iSO)%A2(kmo,1),NumOrb(iMOright),
     &                         Yij(imo,2,iMOright),nj(iMOright))

                   imo=imo+1
                End Do
*               Reset pointers!
                nk = nj(iMOright)
                Xki2(iSO)%A2(1:nk,1:nKBas) => Yij(1:nk*nKBas,1,iMOright)
                Xli3(iSO)%A2(1:nk,1:nLBas) => Yij(1:nk*nLBas,2,iMOright)
             ElseIf (nj(iMOright).gt.NumOrb(iMOright)) Then
                Call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
                Call Abend()
             End If
*
             If ((nj(iMOleft).le.NumOrb(iMOleft))
     &            .and.(jSkip(iMOleft).eq.0)) Then
*
                imo=1
                Do k=1,nj(iMOleft)
                   kmo=kYmnij(k,iMOleft) ! CD-MO index
*
*                  Pick up X_mu,i for all mu's that belong to shell l
*
                   call dcopy_(nLBas,
     &                         Xli2(iSO)%A2(kmo,1),NumOrb(iMOleft),
     &                         Yij(imo,2,iMOleft),nj(iMOleft))

                   call dcopy_(nKBas,
     &                         Xki3(iSO)%A2(kmo,1),NumOrb(iMOleft),
     &                         Yij(imo,1,iMOleft),nj(iMOleft))
                   imo=imo+1
                End Do
                nk = nj(iMOleft)
                Xli2(iSO)%A2(1:nk,1:nLBas) => Yij(1:nk*nLBas,2,iMOleft)
                Xki3(iSO)%A2(1:nk,1:nKBas) => Yij(1:nk*nKBas,1,iMOleft)
             ElseIf (nj(iMOleft).gt.NumOrb(iMOleft)) Then
                Call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
                Call Abend()
             EndIf
           EndIf
         End Do

         Do i2 = 1, iCmp(2)
*
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSO_off = jSO - iOff1
*
            Factor=0.0d0
*
            Do iSO=1,2
              iMOleft=iSO
              iMOright=iSO+2
              If ((nIJR(kSym,lSym,iSO).ne.0).and.(nj(iMOright).ne.0)
     &           .and.(nj(iMOleft).ne.0))  Then
*
*           Read a block of C_kl^J
*
                lCVec = nIJR(kSym,lSym,iSO)*jBas ! Block size
                iAdr = nIJR(kSym,lSym,iSO)*(jSO_off-1) +
     &                 iAdrCVec(jSym,kSym,iSO)
                Call dDaFile(LuCVector(jSym,iSO),2,CijK,
     &               lCVec,iAdr)
*
*           Extract only those C_kl^Js for which we deem k and l to
*           belong to the shell-pair and to be of significance.
*
*MGD skipped jSkip() since not used and complicated in this case
                 If (nj(iMOright).le.NumOrb(iMOright).or.
     &               nj(iMOleft ).le.NumOrb(iMOleft )) Then
                   ij=1
                   Do j=1,nj(iMOleft)
                     jmo=kYmnij(j,iMOleft)
                     Do i=1,nj(iMOright)
                       imo=kYmnij(i,iMOright)
                       jC=imo+NumOrb(iMOright)*(jmo-1)
                       call dcopy_(jBas,CijK(jC),NumOrb(iMOright)*
     &                             NumOrb(iMOleft),CilK(ij),
     &                             nj(iMOright)*nj(iMOleft))
                       ij=ij+1
                     End Do
                   End Do
                   n2j=nj(iMOright)*nj(iMOleft)*jBas
                   CijK(1:n2j)=CilK(1:n2j)
                 End If
*
*           Transform according to Eq. 16 (step 4) and generate B_kl^J
*
*** ----    E(jK,m) = Sum_i C(i,jK)' * X(i,m)
*
                Call dGEMM_('T','N',nj(iMOleft)*jBas,nKBas,nj(iMOright),
     &                     1.0d0,CijK,nj(iMOright),
     &                           Xki2(iSO)%A2,nj(iMOright),
     &                     0.0d0,CilK,nj(iMOleft)*jBas)
*
*** ----    B(Km,n) = Sum_j E(j,Km)' * X(j,n)
*
                Call dGEMM_('T','N',jBas*nKBas,nLBas,nj(iMOleft),
     &                     1.0d0,CilK,nj(iMOleft),
     &                           Xli2(iSO)%A2,nj(iMOleft),
     &                     Factor,BklK,jBas*nKBas)
                Factor=1.0d0
*
** Add transpose
*
*Transpose Cijk->Cjik
                Do ijBas=1,jBas
                  nnk =nj(iMOleft)*nj(iMOright)*(ijBas-1)

                  Do ileft=1,nj(iMOleft)
                    njk = nj(iMOright)*(ileft-1) + nnk

                    Do iright=1,nj(iMOright)
                      nik= nj(iMOleft)*(iright-1)+ nnk

                      ijk = iright + njk
                      jik = ileft  + nik
*
                      CilK(jik)=CijK(ijk)
                    End Do
                  End Do
                End Do
*
*** ----    E(iK,m) = Sum_j C(j,iK)' * X(j,m)
*
                Call dGEMM_('T','N',nj(iMOright)*jBas,nKBas,nj(iMOleft),
     &                     1.0d0,CilK,nj(iMOleft),
     &                           Xki3(iSO)%A2,nj(iMOleft),
     &                     0.0d0,CijK,nj(iMOright)*jBas)
*
*** ----    B(Km,n) = Sum_j E(i,Km)' * X(i,n)
*
                Call dGEMM_('T','N',jBas*nKBas,nLBas,nj(iMOright),
     &                     1.0d0,CijK,nj(iMOright),
     &                           Xli3(iSO)%A2,nj(iMOright),
     &                     Factor,BklK,jBas*nKBas)
              EndIf
            End Do
*
**          Active term
*
            Call dzero(Thpkl,jBas*nKBas*nLBas)
            Do iVec=1,4
              iMO1=1
              iMO2=1
              iVec_=iVec
              fact=1.0d0
              If (iVec.eq.2) iMO2=2
              If (iVec.eq.3) fact=2.0d0
              If (iVec.eq.4) Then
                iMO1=2
                iVec_=2
              EndIf

              Do jAOj=0,jBas-1
                jSOj = jSO_off + jAOj
                Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
                  Do lAOl = 0, lBas-1
                    lSOl = lSO + lAOl-1
                    lp=ipAOrb(0,iMO1)+lSOl*nAct(lSym-1)
                    Do kAct=1,nAct(kSym-1)
                     tmp=ddot_(kact,Zpk(kAct*(kAct-1)/2+1,jSOj,iVec_),1,
     &                         Work(lp),1)
                      Do lAct=kAct+1,nAct(lSym-1)
                        tmp=tmp+Zpk(lAct*(lAct-1)/2+kAct,jSOj,iVec_)*
     &                    Work(lp+lAct-1)
                      End Do
                      CilK(kAct)=tmp
                    End Do
*
                    Do i3 = 1, iCmp(3)
                      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                      kp=ipAOrb(0,iMO2)+(kSO-1)*nAct(kSym-1)
                      iThpkl= jAOj+ (i3-1)*kBas*jBas
     &                           + (lAOl + (i4-1)*lBas)*nKBas*jBas+1
                      Call dGeMV_('T',nAct(kSym-1),kBas,fact,
     &                           Work(kp),
     &                           nAct(kSym-1),Cilk,1,1.0d0,
     &                           Thpkl(iThpkl),jBas)
                    End Do
                  End Do
                End Do
              End Do
            End Do
*
            Do i3 = 1, iCmp(3)
               kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                  iPAO = iPAO + 1
                  nijkl = 0
*
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl
                     Do kAOk = 0, kBas-1
                        kSOk = kSO + kAOk
*
                        iThpkl=(kAOk + (i3-1)*kBas)*jBas
     &                         + (lAOl + (i4-1)*lBas)*nKBas*jBas
                        indexB=iThpkl
                        Indk=Max(kSOk,lSOl)
                        Indl=kSOk+lSOl-Indk
                        Indkl=(Indk-1)*Indk/2+Indl
*
                        Do jAOj = 0, jBas-1
                           jSOj = jSO + jAOj - iOff1
                           nijkl = nijkl + 1
                           indexB = indexB + 1
                           iThpkl = iThpkl + 1
*
*-----------------------SA-CASSCF Coulomb contribution
*
                           temp=CoulFac*(V_k(jSOj,1)*DSO(Indkl,2)+
     &                                   V_k(jSOj,2)*DSO(Indkl,1)+
     &                                   V_k(jSOj,3)*DSO(Indkl,4)+
     &                                   V_k(jSOj,4)*DSO(Indkl,3)+
     &                                   V_k(jSOj,1)*DSO(Indkl,5)+
     &                                   V_k(jSOj,5)*DSO(Indkl,1))
*
*-----------------------Exchange contribution: B(K,m,n)
*
*
                           temp = temp - Factor*ExFac*Half*BklK(indexB)
*
*-----------------------Active space contribution: Sum_p Z(p,K)*Th(p,m,n)
*
                           temp=temp+Thpkl(iThpkl)
*
                           PMax=Max(PMax,Abs(temp))
                           PAO(nijkl,iPAO) =  Fac * temp
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Do iSO=1,2
            Xki2(iSO)%A2 => Null()
            Xki3(iSO)%A2 => Null()
            Xli2(iSO)%A2 => Null()
            Xli3(iSO)%A2 => Null()
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      Else If (ExFac.ne.Zero.and.NumOrb(1).gt.0.and.iMP2prpt.eq.2) Then
*                                                                      *
************************************************************************
*                                                                      *
*        MP2 case
*
         nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)

         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

         lda = SIZE(CMOi(1)%SB(1)%A2,1)
         ik  = 1 + lda*(kSO-1)
         il  = 1 + lda*(lSO-1)
         Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
         Xli(1:) => CMOi(1)%SB(1)%A1(il:)

         If (nj(1).le.NumOrb(1) .and. jSkip(1).eq.0) Then
            imo=1
            Do k=1,nj(1)
               kmo=kYmnij(k,1)
               call dcopy_(nKBas,Xki(kmo),NumOrb(1),
     &                           Yij(imo,1,1),nj(1))
               call dcopy_(nLBas,Xli(kmo),NumOrb(1),
     &                           Yij(imo,2,1),nj(1))
               imo=imo+1
            End Do
            Xki(1:nj(1)*nKBas) => Yij(1:nj(1)*nKBas,1,1)
            Xli(1:nj(1)*nLBas) => Yij(1:nj(1)*nLBas,2,1)
         ElseIf (nj(1).gt.NumOrb(1)) Then
            Call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
            Call Abend()
         EndIf

         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSO_off = jSO - iOff1
*
            lCVec = nIJR(kSym,lSym,1)*jBas
            iAdr = nIJR(kSym,lSym,1)*(jSO_off-1) + iAdrCVec(jSym,kSym,1)
            Call dDaFile(LuCVector(jSym,1),2,CijK,lCVec,iAdr)
*
            If (nj(1).le.NumOrb(1) .and. jSkip(1).eq.0) Then
               ij=1
               Do j=1,nj(1)
                 jmo=kYmnij(j,1)
                 Do i=1,nj(1)
                   imo=kYmnij(i,1)
                   jC=imo+NumOrb(1)*(jmo-1)
                   call dcopy_(jBas,CijK(jC),NumOrb(1)**2,CilK(ij),
     &                        nj(1)**2)
                   ij=ij+1
                 End Do
               End Do
               n2j=nj(1)**2*jBas
               CijK(1:n2j)=CilK(1:n2j)
            EndIf
*
*** ---- C(jK,m) = sum_i C(i,jK)' * X(i,m)
*
            Call dGEMM_('T','N',nj(1)*jBas,nKBas,nj(1),
     &                   1.0d0,CijK,nj(1),
     &                         Xki,nj(1),
     &                   0.0d0,CilK,nj(1)*jBas)
*
*** ---- B(Km,n) = sum_j C(j,Km)' * X(j,n)
*
            Call dGEMM_('T','N',jBas*nKBas,nLBas,nj(1),
     &                   1.0d0,CilK,nj(1),
     &                         Xli,nj(1),
     &                   0.0d0,BklK,jBas*nKBas)
*
****
            lBVec = nBas(0)*nBas(0)*jBas
            Do i = 1,2
               iAdr = 1 + nBas(0)*nBas(0)*(jSO_off-1)
               Call dDaFile(LuBVector(i),2,Bmp2(:,i),lBVec,iAdr)
            End Do
*
            Do i3 = 1, iCmp(3)
               kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                  iPAO = iPAO + 1
                  nijkl = 0
*
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl
                     Do kAOk = 0, kBas-1
                        kSOk = kSO + kAOk
*
                        indexB = (kAOk + (i3-1)*kBas)*jBas
     &                         + (lAOl + (i4-1)*lBas)*nKBas*jBas
                        Indk=Max(kSOk,lSOl)
                        Indl=kSOk+lSOl-Indk
                        Indkl=(Indk-1)*Indk/2+Indl
*
                        Do jAOj = 0, jBas-1
                           jSOj = jSO + jAOj - iOff1
                           nijkl = nijkl + 1
                           indexB = indexB + 1
*
*-----------------------Coulomb contribution: V_k(j)*D(kl)
*
                           temp = CoulFac*(V_k(jSOj,1)*DSO(Indkl,1)
     &                       + U_k(jSOj)*DSO(Indkl,1)
     &                       + V_k(jSOj,1)*(DSO_Var(Indkl)-DSO(indkl,1))
     &                       + Compute_B(irc,kSOk,lSOl,jAOj,iOff1,2))
*
*-----------------------Exchange contribution: B(K,m,n)
*
                           temp = temp - ExFac*Half*(BklK(indexB)
     &                        + Compute_B(irc,kSOk,lSOl,jAOj,iOff1,1))
*                          temp = temp - ExFac*Half*(
*    &                        + Compute_B(irc,kSOk,lSOl,jAOj,iOff1,1))

                           PMax=Max(PMax,Abs(temp))
                           PAO(nijkl,iPAO) =  Fac * temp
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Xki=>Null()
         Xli=>Null()
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        Pure DFT or case when no exhange
*
         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            Do i3 = 1, iCmp(3)
               kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i4 = 1, iCmp(4)
                  lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                  iPAO = iPAO + 1
                  nijkl = 0
*
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl
                     Do kAOk = 0, kBas-1
                        kSOk = kSO + kAOk
*
                        Indk=Max(kSOk,lSOl)
                        Indl=kSOk+lSOl-Indk
                        Indkl=(Indk-1)*Indk/2+Indl
*
                        Do jAOj = 0, jBas-1
                           jSOj = jSO + jAOj - iOff1
                           nijkl = nijkl + 1
*
*-----------------------Coulomb contribution: V_k(j)*D(kl)
*
                           temp=CoulFac*V_k(jSOj,1)*DSO(Indkl,1)
*                          temp=0.0D0
*
                           PMax=Max(PMax,Abs(temp))
                           PAO(nijkl,iPAO) =  Fac * temp
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*     Reset ExFac always.
*
      ExFac=ExFac_
*
      If (iPAO.ne.nPAO) Then
         Write (6,*) ' Error in PGet1_RI3!'
         Call Abend
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt(' In PGet1_RI3:PAO ',' ',PAO,ijkl,nPAO)
      Do i = 1, ijkl
         Write (6,*) DDot_(nPAO,PAO(i,1),ijkl,PAO(i,1),ijkl)
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tbvec(1) = tbvec(1) + Cpu
      tbvec(2) = tbvec(2) + Wall
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Shijij)
         Call Unused_integer(iBas)
         Call Unused_real_array(DSSO)
      End If
      End
