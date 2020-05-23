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
      SubRoutine PGet4(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                 Shijij, iAO, iAOst, ijkl,PSO,nPSO,DSO,nDSO,
     &                 PSOPam,n1,n2,n3,n4,iPam,MapPam,mDim,
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
*             Modified from PGet2, October '92.                        *
************************************************************************
      use pso_stuff
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "lundio.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8 PSO(ijkl,nPSO), PSOPam(n1,n2,n3,n4), DSO(nDSO),
     &       Cred(nCred), Scr1(nScr1,2), Scr2(nScr2)
      Integer nPam(4,0:7), iPam(n1+n2+n3+n4), iiBas(4),
     &          iCmp(4), iShell(4), iAO(4),
     &          iAOst(4), MapPam(4,mDim)
      Logical Shijij
*     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*
      iRout = 39
      iPrint = nPrint(iRout)
*     Call qEnter('PGet4')
      lOper = 1
*
*     Prepare some data for Pam
*
      iiBas(1) = iBas
      iiBas(2) = jBas
      iiBas(3) = kBas
      iiBas(4) = lBas
      nPSOPam=n1*n2*n3*n4
*
*-----Set up table with SO indices in iPam and a table
*     with number of basis functions in each irrep in nPam.
*     Observe that the SO index is only within a given irrep.
*
      Call ICopy(4*8,[0],0,nPam,1)
      in1 = 0
      Do 9 jPam = 1, 4
         in2 = 0
         Do 10 j = 0, nIrrep-1
            Do 11 i1 = 1, iCmp(jPam)
               If (iAnd(IrrCmp(IndS(iShell(jPam))+i1),
     &             iTwoj(j)).ne.0) Then
                   iSO = iAOtSO(iAO(jPam)+i1,j)
     &                 + iAOst(jPam)
                   nPam(jPam,j) = nPam(jPam,j) + iiBas(jPam)
                   Do 12 iAOi = 0, iiBas(jPam)-1
                      iSOi = iSO + iAOi
                      in2 = in2 + 1
                      iPam(in1+in2) = iSOi
                      MapPam(jPam,iSOi+iOffSO(j)) = in2
 12                Continue
               End If
 11         Continue
 10      Continue
         in1 = in1 + in2
 9    Continue
*
*     Get the scrambled 2nd order density matrix
*
      If (LSA) Then
!      write(*,*)"This or ??? in pget4"  !yma

!      do i=1,nG1
!        write(*,*)i,"V-ipG1",Work(ipG1+i-1)
!      end do
!      write(*,*)
!      do i=1,nG2
!        write(*,*)i,"V-ipG2",Work(ipG2+i-1)
!      end do

      Call PTrans_sa(CMO(1,1),nPam,iPam,n1+n2+n3+n4,
     &            DSO,PSOPam,nPSOPam,G1,nG1,G2,nG2,
     &            Cred,nCred/2,Scr1(1,1),nScr1,Scr2,nScr2,Scr1(1,2),
     &            nScr1)
      Else
      Call PTrans(CMO(1,1),nPam,iPam,n1+n2+n3+n4,
     &            DSO,PSOPam,nPSOPam,G1,nG1,G2,nG2,
     &            Cred,nCred,Scr1,nScr1,Scr2,nScr2)
      End If
*
*-----Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      PMax=Zero
      MemSO2 = 0
      Do 100 i1 = 1, iCmp(1)
         niSym = 0
         Do 101 j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndS(iShell(1))+i1),
     &          iTwoj(j)).ne.0) Then
               iSym(niSym) = j
               niSym = niSym + 1
            End if
101      Continue
         Do 200 i2 = 1, iCmp(2)
            njSym = 0
            Do 201 j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndS(iShell(2))+i2),
     &             iTwoj(j)).ne.0) Then
                  jSym(njSym) = j
                  njSym = njSym + 1
               End If
201         Continue
            Do 300 i3 = 1, iCmp(3)
               nkSym = 0
               Do 301 j = 0, nIrrep-1
                  If (iAnd(IrrCmp(IndS(iShell(3))+i3),
     &                iTwoj(j)).ne.0) Then
                     kSym(nkSym) = j
                     nkSym = nkSym + 1
                  End If
301            Continue
               Do 400 i4 = 1, iCmp(4)
                  nlSym = 0
                  Do 401 j = 0, nIrrep-1
                     If (iAnd(IrrCmp(IndS(iShell(4))+i4),
     &                   iTwoj(j)).ne.0) Then
                        lSym(nlSym) = j
                        nlSym = nlSym + 1
                     End If
401               Continue
*
*------Loop over irreps which are spanned by the basis function.
*      Again, the loop structure is restricted to ensure unique
*      integrals.
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
                   If (j123.ne.j4) Go To 411
*
                MemSO2 = MemSO2 + 1
*
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)
*
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
                            nijkl = nijkl + 1
                            k1 = MapPam(1,iSOi)
*
*---------------------------Pick up the contribution.
*
                            PMax=Max(PMax,Abs(PSOPam(k1,k2,k3,k4)))
                            PSO(nijkl,MemSO2) = PSOPam(k1,k2,k3,k4)
*
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
*
 411               Continue
 410            Continue
 310         Continue
 210      Continue
 110   Continue
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
      If (nPSO.ne.MemSO2) Then
         Call WarningMessage(2,'PGet4: nPSO.ne.MemSO2')
         Call Abend()
      End If
*
*     Call GetMem(' Exit PGet4','CHECK','REAL',iDum,iDum)
*     Call qExit('PGet4')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Shijij)
      End
