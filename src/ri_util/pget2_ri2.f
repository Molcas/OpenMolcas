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
      SubRoutine PGet2_RI2(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                  Shijij, iAO, iAOst, nijkl,PSO,nPSO,
     &                  ExFac,CoulFac,PMax,V_K,U_K,mV_K,Z_p_K,nSA)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density matrix.                 *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
* Called from: PGet0                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified to RI-DFT, March 2007                           *
*                                                                      *
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "lundio.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
      Real*8 PSO(nijkl,nPSO), V_K(mV_K,nSA),Z_p_K(nZ_p_k,*)
      Integer iCmp(4), iShell(4), iAO(4), iAOst(4)
      Logical Shijij, Found
*     Local Array
      Integer jSym(0:7), lSym(0:7)
      Integer iTwoj(0:7),CumnnP(0:7),CumnnP2(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUG_
#ifdef _DEBUG_
      iRout = 39
      iPrint = nPrint(iRout)
      iPrint=99
      Call qEnter('PGet2_RI2')
      Call RecPrt('V_K',' ',V_K,1,mV_K)
#endif

      Call CWTime(Cpu1,Wall1)
*                                                                      *
************************************************************************
*                                                                      *
      lOper=1
      PMax=Zero
      iSO=1
      ip_CikJ = ip_CijK
      nMaxBas = Max(jBas,lBas)
*
      Call FZero(PSO,nijkl*nPSO)
*
      If (lPSO) Then
        CumnnP(0)=0
        CumnnP2(0)=0
        Do i=1,nIrrep-1
          nB = nBas_Aux(i-1)
          If (i.eq.1) nB = nB-1
          CumnnP(i)=CumnnP(i-1)+nnP(i-1)
          CumnnP2(i)=CumnnP2(i-1)+nnP(i-1)*nB
        End Do
      End If
*
      Call Qpg_iScalar('SCF mode',Found)
      If (Found) Then
         Call Get_iScalar('SCF mode',iUHF) ! either 0 or 1
      Else
         iUHF=0
      EndIf

*                                                                      *
************************************************************************
*                                                                      *
      Fac = One/Four
      MemSO2 = 0
*                                                                      *
************************************************************************
*                                                                      *
*     Pure DFT
*
      If (ExFac.eq.Zero) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell(2))+i2),
     &         iTwoj(j)).ne.0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndS(iShell(4))+i4),
     &             iTwoj(j)).ne.0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2.ne.j4) Go To 410
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
                  If (j2.ne.0) Go To 410
*
                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
 410              Continue
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Hybrid DFT and HF
*
      Else If (iMP2prpt .ne. 2 .and. .Not.lPSO .and. iUHF.eq.0 ) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell(2))+i2),
     &         iTwoj(j)).ne.0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndS(iShell(4))+i4),
     &             iTwoj(j)).ne.0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2.ne.j4) Go To 510
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  Call FZero(Work(ip_A),jBas*lBas)
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik.eq.0) Go To 520
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,Work(ip_CikJ),
     &                            nik*jBas,iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        ip_CikL = ip_CikJ + nik*nMaxBas
                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,
     &                               Work(ip_CikL),nik*lBas,iAdrL)
                        ip_V2 = ip_CikL
                     Else
                        ip_V2 = ip_CikJ
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,Work(ip_CikJ),nik,
     &                           Work(ip_V2),nik,
     &                           1.0D0,Work(ip_A),jBas)
*
 520                 Continue
*
                  End Do

                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp
     &                       - ExFac*Work(ip_A+mijkl-1)
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
 510              Continue
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Hybrid UDFT and UHF
*
      Else If (iMP2prpt .ne. 2 .and. .Not.lPSO .and. iUHF.eq.1 ) Then
*
         Write (6,*) 'Pget2_RI2: UDFT/UHF not implemented yet.'
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
*     CASSCF
*
      Else If (iMP2prpt .ne. 2 .and. lPSO .and. .Not. LSA) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell(2))+i2),
     &         iTwoj(j)).ne.0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndS(iShell(4))+i4),
     &             iTwoj(j)).ne.0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2.ne.j4) Go To 610
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  Call FZero(Work(ip_A),jBas*lBas)
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik.eq.0) Go To 620
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,Work(ip_CikJ),
     &                            nik*jBas,iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        ip_CikL = ip_CikJ + nik*nMaxBas
                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,
     &                               Work(ip_CikL),nik*lBas,iAdrL)
                        ip_V2 = ip_CikL
                     Else
                        ip_V2 = ip_CikJ
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,Work(ip_CikJ),nik,
     &                           Work(ip_V2),nik,
     &                           1.0D0,Work(ip_A),jBas)
*
 620                 Continue
*
                  End Do
*
                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp
     &                       - ExFac*Work(ip_A+mijkl-1)
*
                        temp2=0.0d0
                        jpSOj=CumnnP2(j2)+(jSOj-1)*nnP(j2)
                        jpSOl=CumnnP2(j2)+(lSOl-1)*nnP(j2)
                        Do jp=1,nnP(j2)
                          temp2=temp2+sign(1.0d0,
     &                          Work(ipDMdiag+CumnnP(j2)+jp-1))*
     &                          Z_p_K(jpSOj+jp,1)*Z_p_K(jpSOl+jp,1)
                        End Do
                        temp=temp+temp2
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
 610              Continue
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     SA-CASSCF
*
      Else If ( iMP2prpt .ne. 2 .and. lPSO .and. lSA ) Then
*                                                                      *
************************************************************************
*                                                                      *
      Write (6,*) 'Pget2_ri2: SA-CASSCF not implemented yet'
      Call Abend()
*
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell(2))+i2),
     &         iTwoj(j)).ne.0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndS(iShell(4))+i4),
     &             iTwoj(j)).ne.0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2.ne.j4) Go To 710
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  Call FZero(Work(ip_A),jBas*lBas)
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik.eq.0) Go To 720
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,Work(ip_CikJ),
     &                            nik*jBas,iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        ip_CikL = ip_CikJ + nik*nMaxBas
                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,
     &                               Work(ip_CikL),nik*lBas,iAdrL)
                        ip_V2 = ip_CikL
                     Else
                        ip_V2 = ip_CikJ
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,Work(ip_CikJ),nik,
     &                           Work(ip_V2),nik,
     &                           1.0D0,Work(ip_A),jBas)
*
 720                 Continue
*
                  End Do

                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=CoulFac*(V_K(lSOl,1)*V_K(jSOj,2)+
     &                                   V_K(lSOl,2)*V_K(jSOj,1)+
     &                                   V_K(lSOl,3)*V_K(jSOj,4)+
     &                                   V_K(lSOl,4)*V_K(jSOj,3))
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp
     &                       - ExFac*Work(ip_A+mijkl-1)
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
 710              Continue
               End Do
            End Do
*
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
      Write (6,*) 'Pget2_ri2: MP2 not implemented yet'
      Call Abend()
*
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell(2))+i2),
     &         iTwoj(j)).ne.0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndS(iShell(4))+i4),
     &             iTwoj(j)).ne.0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2.ne.j4) Go To 810
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  Call FZero(Work(ip_A),jBas*lBas)
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik.eq.0) Go To 820
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,Work(ip_CikJ),
     &                            nik*jBas,iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        ip_CikL = ip_CikJ + nik*nMaxBas
                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,
     &                               Work(ip_CikL),nik*lBas,iAdrL)
                        ip_V2 = ip_CikL
                     Else
                        ip_V2 = ip_CikJ
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,Work(ip_CikJ),nik,
     &                           Work(ip_V2),nik,
     &                           1.0D0,Work(ip_A),jBas)
*
 820                 Continue
*
                  End Do
*
                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp
     &                       - ExFac*Work(ip_A+mijkl-1)
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
 810              Continue
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (nPSO.ne.MemSO2) Then
        Write (6,*) ' PGet2: nPSO.ne.MemSO2'
        Write (6,*) nPSO, MemSO2
        Call Abend()
      End If
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' In PGet2_RI2:PSO ',' ',PSO,nijkl,nPSO)
      End If
      Call GetMem(' Exit PGet2_RI2','CHECK','REAL',iDum,iDum)
      Call qExit('PGet2_RI2')
#endif
*
      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tavec(1) = tavec(1) + Cpu
      tavec(2) = tavec(2) + Wall
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_integer(kBas)
         Call Unused_logical(Shijij)
         Call Unused_real(U_K)
      End If
      End
