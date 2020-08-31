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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine PGet1(PAO,ijkl,nPAO,iCmp,
     &                 iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                 DSO,DSSO,nDSO,ExFac,CoulFac,PMax)
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
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8 PAO(ijkl,nPAO), DSO(nDSO), DSSO(nDSO)
      Integer iAO(4), kOp(4), iAOst(4), iCmp(4)
      Logical Shijij
*
      iRout = 39
      iPrint = nPrint(iRout)
#ifdef _DEBUG_
      Call qEnter('PGet1   ')
      If (iPrint.ge.99) Then
         iComp = 1
         Call PrMtrx('DSO     ',[iD0Lbl],iComp,1,D0)
         Write (6,*) ' nBases..=',iBas,jBas,kBas,lBas
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
      t14 = Quart * ExFac
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
                   Do 220 kAOk = 0, kBas-1
                      kSOk = kSO + kAOk
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         Do 420 iAOi = 0, iBas-1
                            iSOi = iSO + iAOi
                            nijkl = nijkl + 1

*
*---------------------------D(ij)*D(kl)
*
                            Indi=Max(iSOi,jSOj)
                            Indj=iSOi+jSOj-Indi
                            Indk=Max(kSOk,lSOl)
                            Indl=kSOk+lSOl-Indk
                            Indij=(Indi-1)*Indi/2+Indj
                            Indkl=(Indk-1)*Indk/2+Indl
                            temp=DSO(Indij)*DSO(Indkl)*coulfac
*
*--------------------------- -0.25*D(ik)*D(jl)
*
                            Indi=Max(iSOi,kSOk)
                            Indk=iSOi+kSOk-Indi
                            Indj=Max(jSOj,lSOl)
                            Indl=jSOj+lSOl-Indj
                            Indik=(Indi-1)*Indi/2+Indk
                            Indjl=(Indj-1)*Indj/2+Indl
                            temp=temp - t14* (
     &                           DSO(Indik) *DSO(Indjl)
     &                          +DSSO(Indik)*DSSO(Indjl) )
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
     &                           DSO(Indil) *DSO(Indjk)
     &                          +DSSO(Indil)*DSSO(Indjk) )
*
                            PMax=Max(PMax,Abs(Temp))
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
         Call WarningMessage(2,' Error in PGet1!')
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
         Call Unused_logical(Shijij)
      End If
      End
