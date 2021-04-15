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
*               2009, Francesco Aquilante                              *
************************************************************************
      SubRoutine PGet1_RI2(PAO,ijkl,nPAO,iCmp,iAO,iAOst,
     &                     Shijij,iBas,jBas,kBas,lBas,kOp,ExFac,
     &                     CoulFac,PMax,V_K,U_K,mV_K,Z_p_K,nSA)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density.                        *
*                                                                      *
*          (Only for use with C1 point group symmetry)                 *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified for RI-DFT, March 2007                          *
*                                                                      *
*             Modified for RI-HF/CAS, Dec 2009 (F. Aquilante)          *
************************************************************************
      use Basis_Info, only: nBas
      use SOAO_Info, only: iAOtSO
      use pso_stuff, only: nnP, lPSO, lsa, DMdiag, nPos
      use ExTerm, only: CijK, AMP2
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
      Real*8 PAO(ijkl,nPAO), V_K(mV_K,nSA), U_K(mV_K),
     &       Z_p_K(nnP(0),mV_K,*)
      Integer iAO(4), kOp(4), iAOst(4), iCmp(4)
      Logical Shijij,Found

      Real*8, Pointer:: V2(:)=>Null()
      Real*8, Pointer:: CiKj(:,:)=>Null()
      Real*8, Pointer:: CiKl(:)=>Null()
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Do i=1,nSA
         Call RecPrt('PGet1_RI2: V_k',' ',V_k(1,i),1,mV_k)
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     DeSymP will treat up to eight fold degeneracy due to permutational
*     symmetry of shell quadruplets. We will have to compensate for that
*     here since we only have shell doublets.
*
*
      Call CWTime(Cpu1,Wall1)
*
      If (Min(lBas,jBas) .eq.0) Return
*
      Fac = One / Four
      PMax=Zero
      iPAO=0
      iOffA=nBas(0)
*
      Call Qpg_iScalar('SCF mode',Found)
      If (Found) Then
         Call Get_iScalar('SCF mode',iUHF) ! either 0 or 1
      Else
         iUHF=0
      EndIf
*
*                                                                      *
************************************************************************
*                                                                      *
*     Pure DFT
*
      If (ExFac.eq.Zero) Then
*                                                                      *
************************************************************************
*                                                                      *
*
*        Pure DFT
*
         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
*
            Do i4 = 1, iCmp(4)
*
               lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
               iPAO = iPAO + 1
               nijkl = 0
*
               Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl - iOffA
                  Do jAOj = 0, jBas-1
*
                     jSOj = jSO + jAOj - iOffA
                     nijkl = nijkl + 1
*
*----- Coulomb contribution
                     temp=CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
*                    temp=Zero
*
                     PMax=Max(PMax,Abs(temp))
                     PAO(nijkl,iPAO) = Fac * temp
*
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Hybrid DFT and HF
*
      Else If(iMP2prpt .ne. 2 .and. .not. lPSO .and. iUHF.eq.0) Then
*                                                                      *
************************************************************************
*                                                                      *
         iSO=1
*
         jSym = 1
         kSym = jSym
         iSym = 1
         lSym = iEor(jSym-1,iSym-1)+1
*
         nik = nIJ1(iSym,kSym,iSO)

         n = nik*jBas
         iS = 1
         iE = n
         CiKj(1:n,1:1) => CijK(iS:iE)
         n = nik*lBas
         iS = iE + 1
         iE = iE + n
         CiKl(1:n) => CijK(iS:iE)
*
         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
*
*           Pick up the MO transformed fitting coefficients, C_ik^J
            jSOj = jSO - iOffA
            iAdrJ = nik*(jSOj-1) + iAdrCVec(jSym,iSym,1)
            Call dDaFile(LuCVector(jSym,1),2,CikJ(:,1),nik*jBas,iAdrJ)

            Do i4 = 1, iCmp(4)
               lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
               iPAO = iPAO + 1
               nijkl = 0

               If (lSO.ne.jSO) Then
                  lSOl = lSO - iOffA
                  iAdrL = nik*(lSOl-1) + iAdrCVec(jSym,iSym,1)
                  Call dDaFile(LuCVector(jSym,1),2,CiKl,nik*lBas,iAdrL)

                  V2(1:)=>CiKl(1:)
               Else
                  V2(1:)=>CiKj(1:,1)
               EndIf

               Call FZero(Work(ip_A),jBas*lBas)
               CALL DGEMM_('T','N',jBas,lBas,nik,
     &                    1.0d0,CiKj,nik,
     &                          V2,nik,
     &                    0.0d0,Work(ip_A),jBas)

               Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl - iOffA
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj - iOffA
                     nijkl = nijkl + 1

                     temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
                     temp = temp - ExFac*Work(ip_A+nijkl-1)
*
                     PMax = Max(PMax,Abs(temp))
                     PAO(nijkl,iPAO) = Fac*temp
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Hybrid UDFT and UHF
*
      Else If(iMP2prpt .ne. 2 .and. .not. lPSO .and. iUHF.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
         jSym = 1
         kSym = jSym
         iSym = 1
         lSym = iEor(jSym-1,iSym-1)+1
         nik1= nIJ1(iSym,kSym,1)
         nik2= nIJ1(iSym,kSym,2)
         nik = Max(nik1,nik2)
*
         n = nik * jBas
         iS = 1
         iE = n * 2
         CiKj(1:n,1:2) => CijK(iS:iE)
*
         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSOj = jSO - iOffA
*
*           Pick up the MO transformed fitting coefficients, C_ik^J
            If (nik1.ne.0) Then
              iAdrJ = nik1*(jSOj-1) + iAdrCVec(jSym,iSym,1)
              Call dDaFile(LuCVector(jSym,1),2,CiKj(:,1),nik1*jBas,
     &                     iAdrJ)
            EndIf
*
            If (nik2.ne.0) Then
              iAdrJ = nik2*(jSOj-1) + iAdrCVec(jSym,iSym,2)
              Call dDaFile(LuCVector(jSym,2),2,CikJ(:,2),nik2*jBas,
     &                     iAdrJ)
            EndIf
*
            Do i4 = 1, iCmp(4)
               lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
               iPAO = iPAO + 1
               nijkl = 0

               Factor=0.0d0
               Call FZero(Work(ip_A),jBas*lBas)
               Do iSO=1,nKVec
                 nik = nIJ1(iSym,kSym,iSO)

                 CiKl(1:nik*lBas) => CijK(iE+1:iE+nik*lBas)

                 If (nik==0) Cycle

                 If (lSO.ne.jSO) Then
                    lSOl = lSO - iOffA
                    iAdrL = nik*(lSOl-1) + iAdrCVec(jSym,iSym,iSO)
                    Call dDaFile(LuCVector(jSym,iSO),2,CiKl,nik*lBas,
     &                           iAdrL)
                    V2(1:) => CiKl(1:)
                 Else
                    V2(1:) => CiKj(1:,iSO)
                 EndIf
*
                 CALL DGEMM_('T','N',jBas,lBas,nik,
     &                       1.0d0,CikJ(:,iSO),nik,
     &                             V2,nik,
     &                       Factor,Work(ip_A),jBas)
                 Factor=1.0d0
               End Do

               Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl - iOffA
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj - iOffA
                     nijkl = nijkl + 1

                     temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
                     temp = temp - 2.0d0*ExFac*Work(ip_A+nijkl-1)
*
                     PMax = Max(PMax,Abs(temp))
                     PAO(nijkl,iPAO) = Fac*temp
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*     CASSCF
*
      Else If(iMP2prpt .ne. 2 .and. lPSO .and. .not.LSA) Then
*                                                                      *
************************************************************************
*                                                                      *
         iSO=1
*
         jSym = 1
         kSym = jSym
         iSym = 1
         lSym = iEor(jSym-1,iSym-1)+1
*
         nik = nIJ1(iSym,kSym,iSO)
         iS = 1
         iE = nik*jBas
         CiKj(1:nik*jBas,1:1) => CijK(iS:iE)
         iS = iE + 1
         iE = iE + nik*lBas
         CiKl(1:nik*lBas) => CijK(iS:iE)
*
         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
*
*           Pick up the MO transformed fitting coefficients, C_ik^J
            jSOj = jSO - iOffA
            iAdrJ = nik*(jSOj-1) + iAdrCVec(jSym,iSym,1)
            Call dDaFile(LuCVector(jSym,1),2,CikJ(:,1),nik*jBas,iAdrJ)

            Do i4 = 1, iCmp(4)
               lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
               iPAO = iPAO + 1
               nijkl = 0

               If (lSO.ne.jSO) Then
                  lSOl = lSO - iOffA
                  iAdrL = nik*(lSOl-1) + iAdrCVec(jSym,iSym,1)
                  Call dDaFile(LuCVector(jSym,1),2,CiKl,nik*lBas,iAdrL)

                  V2(1:) = CiKl(1:)
               Else
                  V2(1:) = CiKj(1:,1)
               EndIf

               Call FZero(Work(ip_A),jBas*lBas)
               CALL DGEMM_('T','N',jBas,lBas,nik,
     &                    1.0d0,CiKj(:,1),nik,
     &                          V2,nik,
     &                    0.0d0,Work(ip_A),jBas)

               Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl - iOffA
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj - iOffA
                     nijkl = nijkl + 1

                     temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
                     temp = temp - ExFac*Work(ip_A+nijkl-1)
*
*----- Active space contribution
                     temp2=0.0d0
                     Do jp=1,nnP(0)
                       temp2 = temp2 +
     &                         sign(1.0d0,DMdiag(jp,1))*
     &                         Z_p_K(jp,jSOj,1)*Z_p_K(jp,lSOl,1)
                     End Do
                     temp=temp+temp2
*
                     PMax = Max(PMax,Abs(temp))
                     PAO(nijkl,iPAO) = Fac*temp
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*     SA-CASSCF
*
      Else If( iMP2prpt .ne. 2 .and. lPSO .and. lSA ) Then
*                                                                      *
************************************************************************
*                                                                      *
         jSym = 1
         kSym = jSym
         iSym = 1
         lSym = iEor(jSym-1,iSym-1)+1
*
         nik1= nIJ1(iSym,kSym,1)
         nik2= nIJ1(iSym,kSym,2)
         nik = Max(nik1,nik2)

         iS = 1
         iE = 2 * nik * jBas
         CiKj(1:nik*jBas,1:2) => CijK(iS:iE)
*
         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
            jSOj = jSO - iOffA
*
*           Pick up the MO transformed fitting coefficients, C_ik^J
            If (nik1.ne.0) Then
              iAdrJ = nik1*(jSOj-1) + iAdrCVec(jSym,iSym,1)
              Call dDaFile(LuCVector(jSym,1),2,CikJ(:,1),nik1*jBas,
     &                     iAdrJ)
            EndIf
*
            If (nik2.ne.0) Then
              iAdrJ = nik2*(jSOj-1) + iAdrCVec(jSym,iSym,2)
*
              Call dDaFile(LuCVector(jSym,2),2,CikJ(:,2),nik2*jBas,
     &                     iAdrJ)
            EndIf
*
            Do i4 = 1, iCmp(4)
               lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
               iPAO = iPAO + 1
               nijkl = 0

               Factor=0.0d0
               Call FZero(Work(ip_A),jBas*lBas)

               Do iSO=1,nKVec
                 nik = nIJ1(iSym,kSym,iSO)

                 CiKl(1:nik*lBas) => CijK(iE+1:iE+nik*lBas)

                 If (nik==0) Cycle

                 If (lSO.ne.jSO) Then
                    lSOl = lSO - iOffA
                    iAdrL = nik*(lSOl-1) + iAdrCVec(jSym,iSym,iSO)
                    Call dDaFile(LuCVector(jSym,iSO),2,CikL,nik*lBas,
     &                           iAdrL)
                    V2(1:) => CiKl(1:)
                 Else
                    V2(1:) => CiKj(1:,iSO)
                 EndIf
*
** Here one should keep track of negative eigenvalues of the densities
*
                 iSO2=iSO+2
*
                 Do l=1,lBas
                   Do k=1,jBas

                     tmp=0.0d0

                     Do i=1,nChOrb(0,iSO)
                        Do j=1,nChOrb(0,iSO2)

                           jik = j + nChOrb(0,iSO2)*(i-1) + nik*(k-1)
                           jil = j + nChOrb(0,iSO2)*(i-1) + nik*(l-1)
                           If (j<=npos(0,iSO)) Then
                              tmp = tmp + CiKj(jik,iSO)* V2(jil)
                           Else
                              tmp = tmp - CiKj(jik,iSO)* V2(jil)
                           End If
                        End Do
                     End Do

                     kl = (ip_A-1)  + k + jBas*(l-1)
                     Work(kl)= Factor*Work(kl)+tmp

                   End Do
                 End Do
                 Factor=1.0d0

               End Do

               Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl - iOffA
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj - iOffA
                     nijkl = nijkl + 1

                     temp=CoulFac*(V_K(lSOl,1)*V_K(jSOj,2)+
     &                             V_K(lSOl,2)*V_K(jSOj,1)+
     &                             V_K(lSOl,3)*V_K(jSOj,4)+
     &                             V_K(lSOl,4)*V_K(jSOj,3)+
     &                             V_K(lSOl,1)*V_K(jSOj,5)+
     &                             V_K(lSOl,5)*V_K(jSOj,1))
                     temp = temp - ExFac*Work(ip_A+nijkl-1)
*
*----- Active space contribution
                     temp2=0.0d0
                     Do jp=1,nnP(0)
                       temp2 = temp2 +
     &                           sign(1.0d0,DMdiag(jp,1))*
     &                         Z_p_K(jp,jSOj,1)*Z_p_K(jp,lSOl,1)+
     &                           sign(2.0d0,DMdiag(jp,2))*
     &                         (Z_p_K(jp,jSOj,2)*Z_p_K(jp,lSOl,3)+
     &                          Z_p_K(jp,jSOj,3)*Z_p_K(jp,lSOl,2))
                     End Do
                     temp=temp+temp2
*
                     PMax = Max(PMax,Abs(temp))
                     PAO(nijkl,iPAO) = Fac*temp
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*     MP2
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
         iSO=1
*
         jSym = 1
         kSym = jSym
         iSym = 1
         lSym = iEor(jSym-1,iSym-1)+1

         nik = nIJ1(iSym,lSym,iSO)

         iS = 1
         iE = nik * Max(jBas,lBas) * 2
         CiKj(1:nik * Max(jBas,lBas),1:2) => CijK(iS:iE)

         Do i2 = 1, iCmp(2)
            jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

            jSOj = jSO - iOffA
            iAdrJ = nik*(jSOj-1) + iAdrCVec(jSym,iSym,1)
            Call dDaFile(LuCVector(jSym,1),2,CiKj(:,1),nik*jBas,iAdrJ)

            Do i4 = 1, iCmp(4)
               lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
               iPAO = iPAO + 1
               nijkl = 0
*
               If (lSO.ne.jSO) Then
                  lSOl = lSO - iOffA
                  iAdrL = nik*(lSOl-1) + iAdrCVec(jSym,iSym,1)
                  Call dDaFile(LuCVector(jSym,1),2,CiKj(:,2),nik*lBas,
     &                         iAdrL)

                  V2(1:) => CiKj(1:,2)
               Else
                  V2(1:) => CiKj(1:,1)
               EndIf

               Call FZero(Work(ip_A),jBas*lBas)
               CALL DGEMM_('T','N',jBas,lBas,nik,
     &                    1.0d0,CiKj(:,1),nik,
     &                          V2,nik,
     &                    0.0d0,Work(ip_A),jBas)
*
               Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl - iOffA
*
*                 While the I/O here has been moved outside the
*                 inner loop this needs to be reconsidered and
*                 improved such that it can be moved out yet
*                 another loop (or more.)
*
                  lTot = jBas
                  iAdrA = nAuxVe*(lSOl-1) + (jSO - iOffA)
                  Call dDaFile(LuAVector(1),2,AMP2(:,1),lTot,iAdrA)
                  iAdrA = nAuxVe*(lSOl-1) + (jSO - iOffA)
                  Call dDaFile(LuAVector(2),2,AMP2(:,2),lTot,iAdrA)
*
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj - iOffA
                     nijkl = nijkl + 1
*
                     temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
     &                    + CoulFac*V_K(jSOj,1)*U_K(lSOl)
     &                    + CoulFac*U_K(jSOj)*V_K(lSOl,1)
     &                    - ExFac*Work(ip_A+nijkl-1)
*
                     tempJ_mp2=AMP2(1+jAOj,2)
                     temp = temp + tempJ_mp2*CoulFac
*
                     tempK_mp2=AMP2(1+jAOj,1)
                     temp = temp - ExFac*half*tempK_mp2
*
                     PMax = Max(PMax,Abs(temp))
                     PAO(nijkl,iPAO) = Fac*temp
                  End Do
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If

      CiKj => Null()
      CiKl => Null()
      V2   => Null()
*                                                                      *
************************************************************************
*                                                                      *
      If (iPAO.ne.nPAO) Then
         Write (6,*) ' Error in PGet1_RI2!'
         Call Abend
      End If
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In PGet1_RI2:PAO ',' ',PAO,ijkl,nPAO)
#endif
      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tavec(1) = tavec(1) + Cpu
      tavec(2) = tavec(2) + Wall
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Shijij)
         Call Unused_integer(iBas)
         Call Unused_integer(kBas)
      End If
      End
