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
      SubRoutine PGet1_CD3(PAO,ijkl,nPAO,iCmp,iShell,
     &                 iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                 DSO,DSSO,DSO_Var,nDSO,ExFac,CoulFac,PMax,V_k,
     &                 U_k,mV_k)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density.                        *
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
*             Modified for Cholesky 1-center gradients May 2007 by     *
*             R. Lindh                                                 *
*                                                                      *
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "chomp2g_alaska.fh"
#include "exterm.fh"
#include "WrkSpc.fh"
      Real*8 PAO(ijkl,nPAO), DSO(nDSO), DSSO(nDSO), V_k(mV_k),
     &       U_k(mV_k), DSO_Var(nDSO)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iCmp(4)
      Logical Shijij, skip
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 39
      iPrint = nPrint(iRout)
*define _DEBUG_
#ifdef _DEBUG_
      Call qEnter('PGet1_CD3')
      iPrint=99
      If (iPrint.ge.99) Then
         iComp = 1
         Call PrMtrx('DSO     ',[iD0Lbl],iComp,[ipD0],Work)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in PAO in a
*     sequential way.
*
C     Fac = One / Four

      Call CWTime(Cpu1,Wall1)

      Fac = One / Two
      PMax=Zero
      skip = .false.
      iPAO = 0
      jSym = 1
      kSym = 1
      lSym = 1
      NumOrb = nChOrb(kSym-1,1)

      If(ExFac .ne. Zero .and. NumOrb .gt.0 .and. iMP2prpt .ne. 2) Then

         nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)

         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         index2k= NumOrb*(kSO-1)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)
         index2l= NumOrb*(lSO-1)

         Do i1 = 1, iCmp(1)
            iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
            Do iAOi = 0, iBas-1
               iSOi = iSO + iAOi
               Do i2 = 1, iCmp(2)
                  jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj
                     Indi=Max(iSOi,jSOj)
                     Indj=iSOi+jSOj-Indi
                     If(Indi .eq. Indj) Then
                        Fac_ij = 1.0d0
                     Else
                        Fac_ij = 0.5d0
                     End If
                     Indij=(Indi-1)*Indi/2+Indj
                     ijVec=mn2K(Indij,1)

                     If(ijVec.ne.0) Then
                        iAdr = nIJR(kSym,lSym,1)*(ijVec-1) +
     &                       iAdrCVec(jSym,kSym,1)
                        Call dDaFile(LuCVector(jSym,1),2,Work(ip_CijK),
     &                       nIJR(kSym,lSym,1),iAdr)

                        Call dGEMM_('T','N',NumOrb,nKBas,NumOrb,
     &                             1.0d0,Work(ip_CijK),NumOrb,
     &                             Work(ip_CMOi(1)+index2k),NumOrb,
     &                             0.0d0,Work(ip_CilK),Max(1,NumOrb))

                        Call dGEMM_('T','N',nKBas,nLBas,NumOrb,
     &                             1.0d0,Work(ip_CilK),NumOrb,
     &                             Work(ip_CMOi(1)+index2l),NumOrb,
     &                             0.0d0,Work(ip_BklK),Max(1,nKBas))
                     End If

                     Do i3 = 1, iCmp(3)
                        kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                        Do i4 = 1, iCmp(4)
                           lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

                           iPAO = i4 + (i3-1)*iCmp(4)
     &                               + (i2-1)*iCmp(4)*iCmp(3)
     &                               + (i1-1)*iCmp(4)*iCmp(3)*iCmp(2)

                           Do lAOl = 0, lBas-1
                              lSOl = lSO + lAOl
                              Do kAOk = 0, kBas-1
                                 kSOk = kSO + kAOk
                                 indexB = ip_BklK +
     &                                    (kAOk + (i3-1)*kBas)
     &                                  + (lAOl + (i4-1)*lBas)*nKBas
                                 nijkl = iAOi + jAOj*iBas
     &                                 + kAOk*iBas*jBas
     &                                 + lAOl*iBas*jBas*kBas + 1

                                 Indk=Max(kSOk,lSOl)
                                 Indl=kSOk+lSOl-Indk
                                 Indkl=(Indk-1)*Indk/2+Indl
                                 temp=V_k(Indij)*DSO(Indkl)*coulfac
                                 If(ijVec .ne. 0) Then
                                    tempK = Work(indexB)
                                 Else
                                    tempK = 0.0d0
                                 End If

                                 temp = temp - tempK*ExFac*Half*fac_ij
                                 PMax=Max(PMax,Abs(Temp))
                                 PAO(nijkl,iPAO) = Fac * temp
*
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do

      Else If(iMP2prpt .eq. 2 .and. NumOrb .gt. 0) Then

                  nKBas = kBas*iCmp(3)
         nLBas = lBas*iCmp(4)

         kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
         index2k= NumOrb*(kSO-1)
         lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)
         index2l= NumOrb*(lSO-1)

         Do i1 = 1, iCmp(1)
            iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
            Do iAOi = 0, iBas-1
               iSOi = iSO + iAOi
               Do i2 = 1, iCmp(2)
                  jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                  Do jAOj = 0, jBas-1
                     jSOj = jSO + jAOj
                     Indi=Max(iSOi,jSOj)
                     Indj=iSOi+jSOj-Indi
                     If(Indi .eq. Indj) Then
                        Fac_ij = 1.0d0
                     Else
                        Fac_ij = 0.5d0
                     End If
                     Indij=(Indi-1)*Indi/2+Indj
                     ijVec=mn2K(Indij,1)
                     If(ijVec.ne.0) Then
                        iAdr = nIJR(kSym,lSym,1)*(ijVec-1) +
     &                       iAdrCVec(jSym,kSym,1)
                        Call dDaFile(LuCVector(jSym,1),2,Work(ip_CijK),
     &                       nIJR(kSym,lSym,1),iAdr)

                        Call dGEMM_('T','N',NumOrb,nKBas,NumOrb,
     &                             1.0d0,Work(ip_CijK),NumOrb,
     &                             Work(ip_CMOi(1)+index2k),NumOrb,
     &                             0.0d0,Work(ip_CilK),Max(1,NumOrb))

                        Call dGEMM_('T','N',nKBas,nLBas,NumOrb,
     &                             1.0d0,Work(ip_CilK),NumOrb,
     &                             Work(ip_CMOi(1)+index2l),NumOrb,
     &                             0.0d0,Work(ip_BklK),Max(1,nKBas))
                        lBVec = nBas(0)*nBas(0)
                        Do i = 1,2
                           iAdr = 1 + nBas(0)*nBas(0)*(ijVec-1)
                           Call dDaFile(LuBVector(i),2,Work(ip_B_mp2(i))
     &                                  , lBVec,iAdr)
                        End Do

                     End If
                     Do i3 = 1, iCmp(3)
                        kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                        Do i4 = 1, iCmp(4)
                           lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

                           iPAO = i4 +  (i3-1)*iCmp(4)
     &                               + (i2-1)*iCmp(4)*iCmp(3)
     &                               + (i1-1)*iCmp(4)*iCmp(3)*iCmp(2)

                           Do lAOl = 0, lBas-1
                              lSOl = lSO + lAOl
                              Do kAOk = 0, kBas-1
                                 kSOk = kSO + kAOk
                                 indexB = ip_BklK +
     &                                    (kAOk + (i3-1)*kBas)
     &                                  + (lAOl + (i4-1)*lBas)*nKBas
                                 nijkl = iAOi + jAOj*iBas
     &                                 + kAOk*iBas*jBas
     &                                 + lAOl*iBas*jBas*kBas + 1

                                 Indk=Max(kSOk,lSOl)
                                 Indl=kSOk+lSOl-Indk
                                 Indkl=(Indk-1)*Indk/2+Indl
                                 temp=V_k(Indij)*DSO(Indkl)*coulfac

                                 If(ijVec.ne.0) Then
                                    tempK = Work(indexB)
                                 Else
                                    tempK = 0.0d0
                                 End If
                                 temp = temp
     &                                + U_k(indij)*DSO(indkl)*CoulFac
                                 temp = temp + V_k(indij)*
     &                               (DSO_Var(indkl)-DSO(indkl))*CoulFac
                                 if(ijVec.ne.0) Then
                                    tempJ = Compute_B_4(irc,kSOk,
     &                                   lSOl,0,nBas(0),2)
                                    temp = temp + tempJ*CoulFac*
     &                                   fac_ij

                                    tempK = tempK +
     &                                   Compute_B_4(irc,kSOk,lSOl,
     &                                   0,nBas(0),1)
                                 End If
                                 temp = temp - tempK*ExFac*Half*fac_ij
                                 PMax=Max(PMax,Abs(Temp))
                                 PAO(nijkl,iPAO) = Fac * temp
*
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do


      Else
         Do i1 = 1, iCmp(1)
            Do i2 = 1, iCmp(2)
               Do i3 = 1, iCmp(3)
                  Do i4 = 1, iCmp(4)
                     iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
                     jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                     kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                     lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                     iPAO = iPAO + 1
*
                     nijkl = 0
                     Do lAOl = 0, lBas-1
                        lSOl = lSO + lAOl
                        Do kAOk = 0, kBas-1
                           kSOk = kSO + kAOk
                           Indk=Max(kSOk,lSOl)
                           Indl=kSOk+lSOl-Indk
                           Indkl=(Indk-1)*Indk/2+Indl
                           Do jAOj = 0, jBas-1
                              jSOj = jSO + jAOj
                              Do iAOi = 0, iBas-1
                                 iSOi = iSO + iAOi
                                 nijkl = nijkl + 1
*
*---------------------------V_k(ij)*D(kl)
*
                                 Indi=Max(iSOi,jSOj)
                                 Indj=iSOi+jSOj-Indi
                                 Indij=(Indi-1)*Indi/2+Indj
*
                                 temp=V_k(Indij)*DSO(Indkl)*coulfac
*
                                 PMax=Max(PMax,Abs(Temp))
                                 PAO(nijkl,iPAO) = Fac * temp
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If
      If (iPAO.ne.nPAO) Then
         Write (6,*) ' Error in PGet1_CD3!'
         Call Abend
      End If
*
#ifdef _DEBUG_
      Call RecPrt(' In PGet1_CD3:PAO ',' ',PAO,ijkl,nPAO)
      Call GetMem(' Exit PGet1_CD3','CHECK','REAL',iDum,iDum)
      Call qExit('PGet1_CD3')
#endif

      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tbvec(1) = tbvec(1) + Cpu
      tbvec(2) = tbvec(2) + Wall

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
         Call Unused_real_array(DSSO)
      End If
      End
