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
      SubRoutine PGet2_Aces(iCmp,iBas,jBas,kBas,lBas,
     &                      Shijij, iAO, iAOst, nijkl,PSO,nPSO,
     &                      DSO,DSO_Var,DSSO,DSSO_Var,nDSO,
     &                      Gamma,nGamma,iSO2cI,nSOs,
     &                      iSO2Sh,PMax)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density matrix.                 *
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
      use SOAO_Info, only: iAOtSO, iOffSO
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "lundio.fh"
#include "print.fh"
#include "WrkSpc.fh"
************ columbus interface ****************************************
#include "columbus_gamma.fh"
      parameter (exfac=1d0)
      Real*8 PSO(nijkl,nPSO), DSO(nDSO),  DSO_Var(nDSO),
     &       Gamma(nGamma),  DSSO(nDSO), DSSO_Var(nDSO)
      Integer iSO2cI(2,nSOs), iSO2Sh(nSOs)
      Integer iCmp(4), iAO(4), iAOst(4)
      Logical Shijij
*     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
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
      Call qEnter('PGet2')
      If (iPrint.ge.99) Then
         Write (6,*) 'nSOs=',nSOs
         Write (6,*) 'iSO2Sh=',iSO2Sh
         iComp = 1
         Call PrMtrx(' In PGet2:DSO ',[iD0Lbl],iComp,1,D0)
         Call PrMtrx(' In PGet2:DSO_Var ',[iD0Lbl],iComp,1,DVar)
      End If
#endif
      lOper = 1
      t14 = Quart * ExFac
      PMax=Zero
*
*-----Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      MemSO2 = 0
      Do 100 i1 = 1, iCmp(1)
         niSym = 0
         Do 101 j = 0, nIrrep-1
            If (iAOtSO(iAO(1)+i1,j)>0) Then
               iSym(niSym) = j
               niSym = niSym + 1
            End if
101      Continue
         Do 200 i2 = 1, iCmp(2)
            njSym = 0
            Do 201 j = 0, nIrrep-1
               If (iAOtSO(iAO(2)+i2,j)>0) Then
                  jSym(njSym) = j
                  njSym = njSym + 1
               End If
201         Continue
            Do 300 i3 = 1, iCmp(3)
               nkSym = 0
               Do 301 j = 0, nIrrep-1
                  If (iAOtSO(iAO(3)+i3,j)>0) Then
                     kSym(nkSym) = j
                     nkSym = nkSym + 1
                  End If
301            Continue
               Do 400 i4 = 1, iCmp(4)
                  nlSym = 0
                  Do 401 j = 0, nIrrep-1
                     If (iAOtSO(iAO(4)+i4,j)>0) Then
                        lSym(nlSym) = j
                        nlSym = nlSym + 1
                     End If
401               Continue
*
*------Loop over irreps which are spanned by the basis function.
*
       Do 110 is = 0, niSym-1
          j1 = iSym(is)
*
          Do 210 js = 0, njSym-1
             j2 = jSym(js)
             j12 = iEor(j1,j2)
*
             Do 310 ks = 0, nkSym-1
                j3 = kSym(ks)
                j123 = iEor(j12,j3)
                Do 410 ls = 0, nlSym-1
                   j4 = lSym(ls)
                   If (j123.ne.j4) Go To 410
*
                MemSO2 = MemSO2 + 1
*
*               Unfold the way the eight indices have been reordered.
                iSO_r = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO_r = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO_r = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO_r = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
                iSO_a = iSO_r+iOffSO(j1)
                jSO_a = jSO_r+iOffSO(j2)
                kSO_a = kSO_r+iOffSO(j3)
                lSO_a = lSO_r+iOffSO(j4)
*
                mijkl = 0
                Do 120 lAOl = 0, lBas-1
                   lSOl   = lSO_r + lAOl
                   lSOl_a = lSO_a + lAOl
                   iShell_D=iSO2Sh(lSOl_a)
                   Index_D =iSO2cI(1,lSOl_a)
                   nDim_D  =iSO2cI(2,lSOl_a)
                   Do 220 kAOk = 0, kBas-1
                      kSOk   = kSO_r + kAOk
                      kSOk_a = kSO_a + kAOk
                      iShell_C=iSO2Sh(kSOk_a)
                      Index_C =iSO2cI(1,kSOk_a)
                      nDim_C  =iSO2cI(2,kSOk_a)
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
                         jSOj   = jSO_r + jAOj
                         jSOj_a = jSO_a + jAOj
                         iShell_B=iSO2Sh(jSOj_a)
                         Index_B =iSO2cI(1,jSOj_a)
                         nDim_B  =iSO2cI(2,jSOj_a)
                         Do 420 iAOi = 0, iBas-1
                            iSOi   = iSO_r + iAOi
                            iSOi_a = iSO_a + iAOi
                            iShell_A=iSO2Sh(iSOi_a)
                            Index_A =iSO2cI(1,iSOi_a)
                            nDim_A  =iSO2cI(2,iSOi_a)
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
                            mijkl = mijkl + 1
*
************ columbus interface ****************************************
*do not reconstruct the two-particle density from the one-particle
*density or partial two-particle densities but simply read them from
*file

                            if (gamma_mrcisd) goto 95
*---------------------------Contribution D(ij)*D(kl) to P(ijkl)
                            If (j1.eq.j2) Then
*------------------------------j3.eq.j4 also
                               Indi=Max(iSOi,jSOj)
                               Indj=iSOi+jSOj-Indi
                               Indk=Max(kSOk,lSOl)
                               Indl=kSOk+lSOl-Indk
                               iPntij=iPntSO(j1,j2,lOper,nbas)
                               iPntkl=iPntSO(j3,j4,lOper,nbas)
                               Indij=iPntij+(Indi-1)*Indi/2+Indj
                               Indkl=iPntkl+(Indk-1)*Indk/2+Indl
                               temp=DSO(Indij)*DSO(Indkl)
     &                           +(DSO_Var(Indij)-DSO(Indij))*DSO(Indkl)
     &                           +DSO(Indij)*(DSO_Var(Indkl)-DSO(Indkl))
                            Else
                               temp = Zero
                            End If
*
*---------------------------Contribution -1/4*D(ik)*D(jl) to P(ijkl)
                            If (j1.eq.j3) Then
*------------------------------j2.eq.j4 also
                               Indi=Max(iSOi,kSOk)
                               Indk=iSOi+kSOk-Indi
                               Indj=Max(jSOj,lSOl)
                               Indl=jSOj+lSOl-Indj
                               iPntik=iPntSO(j1,j3,lOper,nbas)
                               iPntjl=iPntSO(j2,j4,lOper,nbas)
                               Indik=iPntik+(Indi-1)*Indi/2+Indk
                               Indjl=iPntjl+(Indj-1)*Indj/2+Indl
                               temp=temp-t14*(
     &                              DSO(Indik)*DSO(Indjl)
     &                        +(DSO_Var(Indik)-DSO(Indik))*DSO(Indjl)
     &                        +DSO(Indik)*(DSO_Var(Indjl)-DSO(Indjl))
     &                        +DSSO(Indik)*DSSO(Indjl)
     &                        +(DSSO_Var(Indik)-DSSO(Indik))*DSSO(Indjl)
     &                        +DSSO(Indik)*(DSSO_Var(Indjl)-DSSO(Indjl))
     &                                       )
                            End If
*
*---------------------------Contribution -1/4*D(il)*D(jk) to P(ijkl)
                            If (j1.eq.j4) Then
*------------------------------j2.eq.j3 also
                               Indi=Max(iSOi,lSOl)
                               Indl=iSOi+lSOl-Indi
                               Indj=Max(jSOj,kSOk)
                               Indk=jSOj+kSOk-Indj
                               iPntil=iPntSO(j1,j4,lOper,nbas)
                               iPntjk=iPntSO(j2,j3,lOper,nbas)
                               Indil=iPntil+(Indi-1)*Indi/2+Indl
                               Indjk=iPntjk+(Indj-1)*Indj/2+Indk
                               temp=temp-t14*(
     &                              DSO(Indil)*DSO(Indjk)
     &                        +(DSO_Var(Indil)-DSO(Indil))*DSO(Indjk)
     &                        +DSO(Indil)*(DSO_Var(Indjk)-DSO(Indjk))
     &                        +DSSO(Indil)*DSSO(Indjk)
     &                        +(DSSO_Var(Indil)-DSSO(Indil))*DSSO(Indjk)
     &                        +DSSO(Indil)*(DSSO_Var(Indjk)-DSSO(Indjk))
     &                                       )
                            End If
*
*                           Write (*,*) 'iSO:',
*    &                                  iSOi_a,jSOj_a,kSOk_a,lSOl_a
*                           Write (*,*) 'iShell:',iShell_A,iShell_B,
*    &                                            iShell_C,iShell_D
*                           Write (*,*) 'nDim:',nDim_A,nDim_B,
*    &                                          nDim_C,nDim_D
*                           Write (*,*)  temp , Gamma(Index_ABCD),
*    &                                   Index_ABCD
                            temp = temp + Four*Gamma(Index_ABCD)
 95                         If(gamma_mrcisd) Then
                               temp=Gamma(Index_ABCD)
                            End If
*
                            PMax=Max(PMax,Abs(Temp))
                            PSO(mijkl,MemSO2) =  temp
*
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
*
 410            Continue
 310         Continue
 210      Continue
 110   Continue
*
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
      If (nPSO.ne.MemSO2) Then
         Call WarningMessage(2,'PGet2_Aces: nPSO.ne.MemSO2')
         Write (6,*) nPSO, MemSO2
         Call Abend()
      End If
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' In PGet2:PSO ',' ',PSO,nijkl,nPSO)
      End If
      Call GetMem(' Exit PGet2','CHECK','REAL',iDum,iDum)
      Call qExit('PGet2')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Shijij)
      End
