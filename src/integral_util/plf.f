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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine PLF(AOInt,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &               iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp)
************************************************************************
*  Object: to sift and index the Petite List Format integrals.         *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
* Called from: Twoel                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             May '90                                                  *
************************************************************************
      use SOAO_Info, only: iAOtSO
      use LundIO
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4)
      Logical Shijij, Shij, Shkl, Qijij, iShij, iShkl, Qij, Qkl,
     &        iQij, iQkl
*
      iRout = 109
      iPrint = nPrint(iRout)
      iShij = iShell(1).eq.iShell(2)
      iShkl = iShell(3).eq.iShell(4)
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      Do 100 i1 = 1, iCmp
         jCmpMx = jCmp
         If (Shij) jCmpMx = i1
         Do 200 i2 = 1, jCmpMx
            Qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) Then
               i12 = jCmp*(i1-1) + i2
            Else
               i12 = iCmp*(i2-1) + i1
            End If
            Do 300 i3 = 1, kCmp
               lCmpMx = lCmp
               If (Shkl) lCmpMx = i3
               Do 400 i4 = 1, lCmpMx
                  Qkl = i3.eq.i4
                  If (iShell(4).gt.iShell(3)) Then
                     i34 = lCmp*(i3-1) + i4
                  Else
                     i34 = kCmp*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) Go To 400
                  Qijij = Shijij .and. i12.eq.i34
                  iQij = iShij .and. i1.eq.i2 .and. kOp(1).eq.kOp(2)
                  iQkl = iShkl .and. i3.eq.i4 .and. kOp(3).eq.kOp(4)
*
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                nkl = 0
                fac = one
                Do 120 lAOl = 0, lBas-1
                   lSOl = lSO + lAOl
                   Do 220 kAOk = 0, kBas-1
                      kSOk = kSO + kAOk
                      nkl = nkl + 1
                      If (iQkl .and. kSOk.lt.lSOl) Go To 220
                      If (kSOk.lt.lSOl) Then
                          kSOkk = lSOl
                          lSOll = kSOk
                      Else
                          kSOkk = kSOk
                          lSOll = lSOl
                      End If
                      nijkl = (nkl-1)*iBas*jBas
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         Do 420 iAOi = 0, iBas-1
                            nijkl = nijkl + 1
                            iSOi = iSO + iAOi
                            If (iSOi.lt.jSOj .and. iQij) Go To 420
                            If (iSOi.lt.jSOj) Then
                               iSOii = jSOj
                               jSOjj = iSOi
                            Else
                               iSOii = iSOi
                               jSOjj = jSOj
                            End If
                            iSOij = iSOii*(iSOii-1)/2 + jSOjj
                            iSOkl = kSOkk*(kSOkk-1)/2 + lSOll
                            If (iSOij.lt.iSOkl .and. Qijij) Go To 420
*
                            If (Abs(AOInt(nijkl,i1,i2,i3,i4)).gt.
     &                          ThrInt) Then
                               If (iSOij.lt.iSOkl) Then
                                  iiSOii = kSOkk
                                  jjSOjj = lSOll
                                  kkSOkk = iSOii
                                  llSOll = jSOjj
                               Else
                                  iiSOii = iSOii
                                  jjSOjj = jSOjj
                                  kkSOkk = kSOkk
                                  llSOll = lSOll
                               End If
                               IntTot = IntTot + 1
*
                               nUt=nUt + 1
                               Buf%Buf(nUt) = AOInt(nijkl,i1,i2,i3,i4)
                               Buf%iBuf(nUt) = llSOll + kkSOkk*2**8 +
     &                                     jjSOjj*2**16
                               Buf%iBuf(nUt)=iOr( iShft(iiSOii,24),
     &                                            Buf%iBuf(nUt) )
                               If (nUt.eq.nBuf-1) Then
                                  Call iDafile(Lu_28,1,Buf,lBuf,iDisk)
                                  nUt=0
                               End If
*                              XInt=AOInt(nijkl,i1,i2,i3,i4)
*                              Sum = Sum + XInt
*                              SumAbs = SumAbs+Abs(XInt)
*                              SumSq= SumSq + XInt**2
                            End If
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
      Return
      End
