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
      SubRoutine IndSft(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                  Shijij,iAO,iAOst,ijkl,SOInt,nSOInt)
************************************************************************
*  Object: to sift and index the SO integrals.                         *
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
*             April '90                                                *
************************************************************************
      use SOAO_Info, only: iAOtSO, iOffSO
      use LundIO
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8 SOInt(ijkl,nSOInt)
      Integer iCmp(4), iShell(4), iAO(4), iAOst(4)
      Logical Shijij, Shij, Shkl, Qijij, Qij, Qkl,
     &        iQij, iQkl, Wijij
*     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
*
      iRout = 39
      iPrint = nPrint(iRout)
      k12=0
      k34=0
      If (iPrint.ge.99)
     &   Call RecPrt(' In IndSft:SOInt ',' ',SOInt,ijkl,nSOInt)
      MemSO2 = 0
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      Do 100 i1 = 1, iCmp(1)
         Do 101 j = 0, nIrrep-1
            ix = 0
            If (iAOtSO(iAO(1)+i1,j)>0) ix = 2**j
            iSym(j) = ix
101      Continue
         jCmpMx = iCmp(2)
         If (Shij) jCmpMx = i1
         Do 200 i2 = 1, jCmpMx
            Do 201 j = 0, nIrrep-1
               ix = 0
               If (iAOtSO(iAO(2)+i2,j)>0) ix = 2**j
               jSym(j) = ix
201         Continue
            Qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) Then
               i12 = iCmp(2)*(i1-1) + i2
            Else
               i12 = iCmp(1)*(i2-1) + i1
            End If
            Do 300 i3 = 1, iCmp(3)
               Do 301 j = 0, nIrrep-1
                  ix = 0
                  If (iAOtSO(iAO(3)+i3,j)>0) ix = 2**j
                  kSym(j) = ix
301            Continue
               lCmpMx = iCmp(4)
               If (Shkl) lCmpMx = i3
               Do 400 i4 = 1, lCmpMx
                  Do 401 j = 0, nIrrep-1
                     ix = 0
                     If (iAOtSO(iAO(4)+i4,j)>0) ix = 2**j
                     lSym(j) = ix
401               Continue
                  Qkl = i3.eq.i4
                  If (iShell(4).gt.iShell(3)) Then
                     i34 = iCmp(4)*(i3-1) + i4
                  Else
                     i34 = iCmp(3)*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) Go To 400
                  Qijij = Shijij .and. i12.eq.i34
                  iQij = Shij .and. i1.eq.i2
                  iQkl = Shkl .and. i3.eq.i4
*
*      Loop over irreps which are spanned by the basis function.
*      Again, the loop structure is restricted to ensure unique
*      integrals.
*
       Do 110 j1 = 0, nIrrep-1
          If (iSym(j1).eq.0) Go To 110
*
          j2Max = nIrrep-1
          If (Shij .and. Qij) j2Max = j1
          Do 210 j2 = 0, j2Max
             If (jSym(j2).eq.0) Go To 210
             j12 = iEor(j1,j2)
             If (Qijij) Then
                If (Shij .and. Qij) Then
                    k12 = j1*(j1+1)/2 + j2+1
                Else If (Shij) Then
                    k12 = nIrrep*j1 + j2+1
                Else If (iShell(1).gt.iShell(2)) Then
                    k12 = nIrrep*j1 + j2+1
                Else
                    k12 = nIrrep*j2 + j1+1
                End If
             End If
*
             Do 310 j3 = 0, nIrrep-1
                If (kSym(j3).eq.0) Go To 310
                j4 = iEor(j12,j3)
                If (lSym(j4).eq.0) Go To 310
                If (Shkl .and. Qkl .and. j4.gt.j3) Go To 310
                If (Qijij) Then
                   If (Shkl .and. Qkl) Then
                      k34 = j3*(j3+1)/2 + j4+1
                   Else If (Shkl) Then
                      k34 = nIrrep*j3 + j4+1
                   Else If (iShell(3).gt.iShell(4)) Then
                      k34 = nIrrep*j3 + j4+1
                   Else
                      k34 = nIrrep*j4 + j3+1
                   End If
                   If (Qijij .and. k34.gt.k12) Go To 310
                End If
                Wijij = Qijij .and. k12.eq.k34
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)
*
                MemSO2 = MemSO2 + 1
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
                      nij = 0
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         Do 420 iAOi = 0, iBas-1
                            nij = nij + 1
                            nijkl = (nkl-1)*iBas*jBas + nij
                            If (Abs(SOInt(nijkl,MemSO2)).gt.ThrInt) Then
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
                               If (iSOij.lt.iSOkl .and. Wijij)
     &                              Go To 420
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
                               Buf%nUt=Buf%nUt + 1
                               Buf%Buf(Buf%nUt) = SOInt(nijkl,MemSO2)
                               Buf%iBuf(Buf%nUt) = llSOll + kkSOkk*2**8
     &                                           +  jjSOjj*2**16
                               Buf%iBuf(Buf%nUt)=iOr( iShft(iiSOii,24),
     &                                            Buf%iBuf(Buf%nUt) )
                               If (Buf%nUt.eq.nBuf-1) Then
                                  Call dDafile(Lu_28,1,Buf%Buf,lBuf,
     &                                         iDisk)
                                  Buf%nUt=0
                               End If
                            End If
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
*
 310         Continue
 210      Continue
 110   Continue
*
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
      Return
      End
