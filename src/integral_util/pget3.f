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
      SubRoutine PGet3(PAO,ijkl,nPAO,iCmp,iShell,
     &                 iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                 DAO,nDAO,
     &                 PAOPam,n1,n2,n3,n4,iPam,MapPam,mDim,
     &                 Cred,nCred,Scr1,nScr1,Scr2,nScr2,PMax)
************************************************************************
*  Object: to assemble the index list of the batch of the 2nd order    *
*          density matrix.                                             *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
* Called from: PGet0                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*             Modified from PGet1, June '92                            *
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8 PAO(ijkl,nPAO), PAOPam(n1,n2,n3,n4), DAO(nDAO),
     &       Cred(nCred), Scr1(nScr1,2), Scr2(nScr2)
      Integer iShell(4), iAO(4), kOp(4),
     &          iAOst(4), nPam(4), iPam(n1+n2+n3+n4), iiBas(4),
     &          MapPam(4,mDim), iCmp(4)
      Logical Shijij
*
      iRout = 39
      iPrint = nPrint(iRout)
*     Call qEnter('PGet3   ')
      If (iPrint.ge.99) Then
         iComp = 1
         Write (6,*) ' nBases..=',iBas,jBas,kBas,lBas
      End If
*
*-----Prepare some data for Pam
*
      iiBas(1) = iBas
      iiBas(2) = jBas
      iiBas(3) = kBas
      iiBas(4) = lBas
      nPSOPam=n1*n2*n3*n4
*
*-----Set up table with SO indices in iPam and a table
*     with the number of basis functions in nPam.
*
      Call ICopy(4,[0],0,nPam,1)
      in1 = 0
      Do 9 jPam = 1, 4
         in2 = 0
         Do 11 i1 = 1, iCmp(jPam)
            iSO = iAOtSO(iAO(jPam)+i1,0)
     &          + iAOst(jPam)
            nPam(jPam) = nPam(jPam) + iiBas(jPam)
            Do 12 iAOi = 0, iiBas(jPam)-1
               iSOi = iSO + iAOi
               in2 = in2 + 1
               iPam(in1+in2) = iSOi
               MapPam(jPam,iSOi) = in2
 12         Continue
 11      Continue
         in1 = in1 + in2
 9    Continue
*
*     Get the scrambled 2nd order density matrix
*
      If (LSA) Then
      Call PTrans_sa(Work(ipCMo),nPam,iPam,n1+n2+n3+n4,
     &            DAO,PAOPam,nPSOPam,Work(ipG1),nG1,Work(ipG2),nG2,
     &            Cred,nCred/2,Scr1(1,1),nScr1,Scr2,nScr2,Scr1(1,2),
     &            nScr1)
      Else
       Call PTrans(Work(ipCMo),nPam,iPam,n1+n2+n3+n4,
     &            DAO,PAOPam,nPSOPam,Work(ipG1),nG1,Work(ipG2),nG2,
     &            Cred,nCred,Scr1,nScr1,Scr2,nScr2)
      End If
*
*     Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in PAO in a
*     sequential way.
*
      PMax=Zero
      iPAO=0
      Do 100 i1 = 1, iCmp(1)
         Do 200 i2 = 1, iCmp(2)
            Do 300 i3 = 1, iCmp(3)
               Do 400 i4 = 1, iCmp(4)
*
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)+iOffSO(kOp(1))
                jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)+iOffSO(kOp(2))
                kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)+iOffSO(kOp(3))
                lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)+iOffSO(kOp(4))
*
                iPAO = iPAO + 1
                nijkl = 0
                Do 120 lAOl = 0, lBas-1
                   lSOl = lSO + lAOl
                   k4 = MapPam(4,lSOl)
                   Do 220 kAOk = 0, kBas-1
                      kSOk = kSO + kAOk
                      k3 = MapPam(3,kSOk)
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         k2 = MapPam(2,jSOj)
                         Do 420 iAOi = 0, iBas-1
                            iSOi = iSO + iAOi
                            k1 = MapPam(1,iSOi)
                            nijkl = nijkl + 1
*
                            PMax=Max(PMax,Abs(PAOPam(k1,k2,k3,k4)))
                            PAO(nijkl,iPAO) = PAOPam(k1,k2,k3,k4)
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
        Call WarningMessage(2,' Error in PGet3!')
        Call Abend()
      End If
*
*     Call GetMem(' Exit PGet3','CHECK','REAL',iDum,iDum)
*     Call qExit('PGet3')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
      End If
      End
