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
* Copyright (C) 1992,2000, Roland Lindh                                *
************************************************************************
      SubRoutine PGet1_Aces(PAO,ijkl,nPAO,iCmp,iShell,
     &                      iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                      DSO,DSO_Var,DSSO,DSSO_Var,nDSO,
     &                      Gamma,nGamma,iSO2cI,nSOs,
     &                      iSO2Sh,PMax)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density.                        *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*          DSO: HF 1st order density                                   *
*          DSO_Var: 1st order density of correlated wf.                *
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
*     Modified to Aces 2 by RL, July 2000, Gainesville, FL, USA        *
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
************ columbus interface ****************************************
#include "columbus_gamma.fh"
      parameter (exfac=1d0)
      Real*8 PAO(ijkl,nPAO), DSO(nDSO),  DSO_Var(nDSO),
     &       Gamma(nGamma), DSSO(nDSO), DSSO_Var(nDSO)
      Integer iSO2cI(2,nSOs), iSO2Sh(nSOs)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iCmp(4)
      Logical Shijij
*                                                                      *
************************************************************************
*                                                                      *
*     Statement Function
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 39
      iPrint = nPrint(iRout)
#ifdef _DEBUG_
      Call qEnter('PGet1   ')
      If (iPrint.ge.99) Then
         iComp = 1
         Call PrMtrx('DSO     ',[iD0Lbl],iComp,[ipD0],Work)
         Call PrMtrx('DSO_Var ',[iD0Lbl],iComp,[ipDVar],Work)
         Write (6,*) ' nBases..=',iBas,jBas,kBas,lBas
         Write (6,*) 'iSO2Sh=',iSO2Sh
         Write (6,*) 'iSO2cI(1)',(iSO2cI(1,i),i=1,nSOs)
         Write (6,*) 'iSO2cI(2)',(iSO2cI(2,i),i=1,nSOs)
         Call RecPrt('PGet1: Gamma',' ',Gamma,1,nGamma)
      End If
#endif
*
*     Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in PAO in a
*     sequential way.
*
      PMax=Zero
      iPAO=0
      t14 = Quart * exfac
      Do 100 i1 = 1, iCmp(1)
         Do 200 i2 = 1, iCmp(2)
            Do 300 i3 = 1, iCmp(3)
               Do 400 i4 = 1, iCmp(4)
*
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                iPAO = iPAO + 1
                nijkl = 0
                Do 120 lAOl = 0, lBas-1
                   lSOl = lSO + lAOl
                   iShell_D=iSO2Sh(lSOl)
                   Index_D =iSO2cI(1,lSOl)
                   nDim_D  =iSO2cI(2,lSOl)
                   Do 220 kAOk = 0, kBas-1
                      kSOk = kSO + kAOk
                      iShell_C=iSO2Sh(kSOk)
                      Index_C =iSO2cI(1,kSOk)
                      nDim_C  =iSO2cI(2,kSOk)
                      nDim_CD=nDim_C*nDim_D
                      iShell_CD=iTri(iShell_C,iShell_D)
                      If (iShell_C.gt.iShell_D) Then
                         Index_CD=(Index_D-1)*nDim_C + Index_C
                      Else If (iShell_C.eq.iShell_D) Then
                         Index_CD=iTri(Index_C,Index_D)
                      Else
                         Index_CD=(Index_C-1)*nDim_D + Index_D
                      End If
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         iShell_B=iSO2Sh(jSOj)
                         Index_B =iSO2cI(1,jSOj)
                         nDim_B  =iSO2cI(2,jSOj)
                         Do 420 iAOi = 0, iBas-1
                            iSOi = iSO + iAOi
                            iShell_A=iSO2Sh(iSOi)
                            Index_A =iSO2cI(1,iSOi)
                            nDim_A  =iSO2cI(2,iSOi)
                            nDim_AB=nDim_A*nDim_B
                            iShell_AB=iTri(iShell_A,iShell_B)
                            If (iShell_A.gt.iShell_B) Then
                               Index_AB=(Index_B-1)*nDim_A + Index_A
                            Else If (iShell_A.eq.iShell_B) Then
                               Index_AB=iTri(Index_A,Index_B)
                            Else
                               Index_AB=(Index_A-1)*nDim_B + Index_B
                            End If
                            If (iShell_AB.gt.iShell_CD) Then
                               Index_ABCD=(Index_CD-1)*nDim_AB+Index_AB
                            Else If (iShell_AB.eq.iShell_CD) Then
                               Index_ABCD=iTri(Index_AB,Index_CD)
                            Else
                               Index_ABCD=(Index_AB-1)*nDim_CD+Index_CD
                            End If
                            nijkl = nijkl + 1

************ columbus interface ****************************************
*do not reconstruct the two-particle density from the one-particle
*density or partial two-particle densities but simply read them from
*file
                            if (gamma_mrcisd) goto 95
*
*---------------------------D(ij)*D(kl)
*
                            Indi=Max(iSOi,jSOj)
                            Indj=iSOi+jSOj-Indi
                            Indk=Max(kSOk,lSOl)
                            Indl=kSOk+lSOl-Indk
                            Indij=(Indi-1)*Indi/2+Indj
                            Indkl=(Indk-1)*Indk/2+Indl
                            temp= DSO(Indij)*DSO(Indkl)
     &                          +(DSO_Var(Indij)-DSO(Indij))*DSO(Indkl)
     &                          +DSO(Indij)*(DSO_Var(Indkl)-DSO(Indkl))
*
*--------------------------- -0.25*D(ik)*D(jl)
*
                            Indi=Max(iSOi,kSOk)
                            Indk=iSOi+kSOk-Indi
                            Indj=Max(jSOj,lSOl)
                            Indl=jSOj+lSOl-Indj
                            Indik=(Indi-1)*Indi/2+Indk
                            Indjl=(Indj-1)*Indj/2+Indl
                            temp=temp - t14*(
     &                           DSO(Indik)*DSO(Indjl)
     &                        +(DSO_Var(Indik)-DSO(Indik))*DSO(Indjl)
     &                        +DSO(Indik)*(DSO_Var(Indjl)-DSO(Indjl))
     &                        +DSSO(Indik)*DSSO(Indjl)
     &                        +(DSSO_Var(Indik)-DSSO(Indik))*DSSO(Indjl)
     &                        +DSSO(Indik)*(DSSO_Var(Indjl)-DSSO(Indjl))
     &                                      )
*
*--------------------------- -0.25*D(il)*D(jk)
*
                            Indi=Max(iSOi,lSOl)
                            Indl=iSOi+lSOl-Indi
                            Indj=Max(jSOj,kSOk)
                            Indk=jSOj+kSOk-Indj
                            Indil=(Indi-1)*Indi/2+Indl
                            Indjk=(Indj-1)*Indj/2+Indk
                            temp=temp - t14*(
     &                           DSO(Indil)*DSO(Indjk)
     &                        +(DSO_Var(Indil)-DSO(Indil))*DSO(Indjk)
     &                        +DSO(Indil)*(DSO_Var(Indjk)-DSO(Indjk))
     &                        +DSSO(Indil)*DSSO(Indjk)
     &                        +(DSSO_Var(Indil)-DSSO(Indil))*DSSO(Indjk)
     &                        +DSSO(Indil)*(DSSO_Var(Indjk)-DSSO(Indjk))
     &                                      )
*
                            temp = temp + Four*Gamma(Index_ABCD)
 95                         If(gamma_mrcisd) Then
                               temp = Gamma(Index_ABCD)
                            End If
*
                            PMax=Max(PMax,Abs(temp))
                            PAO(nijkl,iPAO) = temp
*
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
      If (iPAO.ne.nPAO) Then
         Call WarningMessage(2,' Error in PGet1_Aces!')
         Call Abend()
      End If
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' In PGet1:PAO ',' ',PAO,ijkl,nPAO)
         Do 3333 i = 1, ijkl
            Write (6,*) DDot_(nPAO,PAO(i,1),ijkl,
     &                            PAO(i,1),ijkl)
 3333    Continue
      End If
      Call GetMem(' Exit PGet1','CHECK','REAL',iDum,iDum)
      Call qExit('PGet1')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
      End If
      End
