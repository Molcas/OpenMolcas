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
      SubRoutine IndSft_RI_3(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                       Shijij, iAO, iAOst, ijkl,SOint,nSOint,
     &                       iSOSym,nSOs,
     &                       TInt,nTInt,iOff,iShlSO,nBasSh,
     &                       iSOShl,nSO,nShell,nSym,iSSOff)
************************************************************************
*  object: to sift and index the SO integrals.                         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          april '90                                                   *
*                                                                      *
************************************************************************
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "WrkSpc.fh"
*
      Real*8 SOint(ijkl,nSOint), TInt(nTInt)
      Integer iSOShl(nSO), iShlSO(nSO), nBasSh(0:nSym-1,nShell),
     &        iSSOff(0:nIrrep-1,0:nIrrep-1)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), iSOSym(2,nSOs), jOffSO(0:7)
      Logical Shijij, Shkl, qkl
*     local array
      Integer jSym(0:7), kSym(0:7), lSym(0:7), iOff(3,0:7)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      jOffSO(0)=0
      Do iIrrep = 1, nIrrep-1
         jOffSO(iIrrep)=jOffSO(iIrrep-1)+nBas(iIrrep-1)
      End Do
      memSO2 = 0
*
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      Shkl = iShell(3).eq.iShell(4)
      If (iShell(4).gt.iShell(3)) Then
         Call WarningMessage(2,'Error in IndSft_RI_3')
         Write (6,*) 'iShell(4).gt.iShell(3)'
         Call Abend()
      End If
*
      i1=1
      j1=0
      Do i2 = 1, iCmp(2)
         Do 201 j = 0, nIrrep-1
            ix = 0
            If (iAOtSO(iAO(2)+i2,j)>0) ix = 2**j
            jSym(j) = ix
201      Continue
         If (iShell(2).gt.iShell(1)) then
            i12 = iCmp(2)*(i1-1) + i2
         else
            i12 = iCmp(1)*(i2-1) + i1
         End If
         Do i3 = 1, iCmp(3)
            Do 301 j = 0, nIrrep-1
               ix = 0
               If (iAOtSO(iAO(3)+i3,j)>0) ix = 2**j
               kSym(j) = ix
301         Continue
            lCmpMx = iCmp(4)
            If (Shkl) lCmpMx = i3
            Do 400 i4 = 1, lCmpMx
               Do 401 j = 0, nIrrep-1
                  ix = 0
                  If (iAOtSO(iAO(4)+i4,j)>0) ix = 2**j
                  lSym(j) = ix
401            Continue
               qkl = i3.eq.i4
               If (iShell(4).gt.iShell(3)) then
                  i34 = iCmp(4)*(i3-1) + i4
               else
                  i34 = iCmp(3)*(i4-1) + i3
               End If
*
*      loop over Irreps which are spanned by the basis function.
*      again, the loop structure is restricted to ensure unique
*      integrals.
*
       Do 210 j2 = 0, nIrrep-1
          If (jSym(j2).eq.0) go to 210
          j12 = iEor(j1,j2)
*
          Do 310 j3 = 0, nIrrep-1
             If (kSym(j3).eq.0) go to 310
             j4 = iEor(j12,j3)
             If (lSym(j4).eq.0) go to 310
             If (Shkl .and. qkl .and. j4.gt.j3) go to 310
*
             memSO2 = memSO2 + 1
             If ( (nSkip(j2+1)+
     &             nSkip(j3+1)+nSkip(j4+1) ).ne.0 ) GoTo 310
*                                                                      *
************************************************************************
*                                                                      *
*            Number of auxiliary basis functions in this symmetry block.
             mm       = iOff(1,j12)
             If (mm.eq.0) Go  To 310
*            Offset to J12 symmetry block.
             iOff_J12 = iOff(2,j12)
*            Effective number of valence basis products in this symmetry
*            block.
             n3C      = iOff(3,j12)
             If (n3C.eq.0) Go To 310
*            Offset to the symmetry block of this shell pair.
             iOff_L   = iSSOff(j3,j4)
*                                                                      *
************************************************************************
*                                                                      *
*            Compute index within the irrep. Keep the indexation
*            of the two basis set sets apart.
*
             jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)-nBas(j2)
             kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+jOffSO(j3)
             lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+jOffSO(j4)
*
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1
                   iD = iShlSO(lSOl)
                   iShD = iSOShl(lSOl)
                   nD = nBasSh(j4,iShD)
                   Do kSOk = kSO, kSO+kBas-1
                      iC = iShlSO(kSOk)
                      iShC = iSOShl(kSOk)
                      nC = nBasSh(j3,iShC)
*
                      If (iShC.eq.iShD) Then
                         If (j12.eq.0) Then
                            kl=iTri(iC,iD)
                         Else If (j3.gt.j4) Then
                            kl = (iC-1)*nD + iD
                         Else
                            kl = (iD-1)*nC + iC
                         End If
                      Else
                         If (iShC.ge.iShD) Then
                            kl = (iD-1)*nC + iC
                         Else
                            kl = (iC-1)*nD + iD
                         End If
                      End If
*
                      Do jSOj = jSO, jSO+jBas-1
                         iAux= jSOj
                         nijkl = nijkl + 1
                         AInt=SOint(nijkl,memSO2)
*                                                                      *
************************************************************************
*                                                                      *
                         If (j12.eq.0) Then
                            If (kSOk.ge.lSOl) Then
                               kl_B=(iAux-1)*n3C+kl + iOff_L
                               TInt(kl_B)=AInt
                            End If
                         Else
                            kl_B=(iAux-1)*n3C+kl + iOff_L
                            TInt(kl_B)=AInt
                         End If
*                                                                      *
************************************************************************
*                                                                      *
                      End Do
                   End Do
                End Do
*
310       Continue
210    Continue
*
400         Continue
         End Do
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_logical(Shijij)
         Call Unused_integer_array(iSOSym)
      End If
      End
