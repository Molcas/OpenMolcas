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
      Subroutine PLF_RI_3(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                    iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                    TInt,nTInt,iOff,iShlSO,nBasSh,iSOShl,
     &                    nSO,nShell,nSym,iSSOff)
************************************************************************
*                                                                      *
*  object: to sift and index the petite list format integrals.         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          May '90                                                     *
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
      Real*8 AOint(ijkl,jCmp,kCmp,lCmp), TInt(nTInt)
      Integer iSOShl(nSO), iShlSO(nSO), nBasSh(0:nSym-1,nShell)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iSOs(4), iOff(3)
      Logical Shijij, Shkl
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      Shkl = iShell(3).eq.iShell(4)
      iOff1 = nBas(0)
      mm = iOff(1)
      n3C= iOff(3)
      If (iShell(4).gt.iShell(3)) Then
         Write (6,*) 'iShell(4).gt.iShell(3)'
         Call Abend()
      End If
*
      Do i2 = 1, jCmp
         iSOs(2)=iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
         Do i3 = 1, kCmp
            iSOs(3)=iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
            lCmp_Max=lCmp
            If (Shkl) lCmp_Max=i3
            Do i4 = 1, lCmp_Max
               iSOs(4)=iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*                                                                      *
************************************************************************
*                                                                      *
               If (Shkl.and.i3.eq.i4) Then
*                                                                      *
************************************************************************
*                                                                      *
                  nijkl = 0
                  Do lSOl = iSOs(4), iSOs(4)+lBas-1
                     iD = iShlSO(lSOl)                  ! Relative index
                     Do kSOk = iSOs(3), iSOs(3)+kBas-1
                        iC = iShlSO(kSOk)
                        iShC=iSOShl(kSOk)
                        nC = nBasSh(0,iShC)
*
                        kl = iTri(iC,iD) + iSSOff
*
                        Do jSOj = iSOs(2), iSOs(2)+jBas-1
                           nijkl = nijkl + 1
                           If (lSOl.gt.kSOk) Go To 99
                           iAux = jSOj - iOff1
                           AInt=AOint(nijkl,i2,i3,i4)
*
                           kl_B = (iAux-1)*n3C + kl
                           TInt(kl_B) = AInt
*
 99                        Continue
*
                        End Do
                     End Do
                  End Do
*                                                                      *
************************************************************************
*                                                                      *
               Else
*                                                                      *
************************************************************************
*                                                                      *
                  nijkl = 0
                  Do lSOl = iSOs(4), iSOs(4)+lBas-1
                     iD = iShlSO(lSOl)
                     Do kSOk = iSOs(3), iSOs(3)+kBas-1
                        iC = iShlSO(kSOk)
                        iShC=iSOShl(kSOk)
                        nC = nBasSh(0,iShC)
*
                        If (Shkl) Then
                           kl = iTri(iC,iD)
                        Else
                           kl = (iD-1)*nC + iC
                        End If
                        kl = kl + iSSOff
*
                        Do jSOj = iSOs(2), iSOs(2)+jBas-1
                           iAux = jSOj - iOff1
                           nijkl = nijkl + 1
                           AInt=AOint(nijkl,i2,i3,i4)
*
                           kl_B = (iAux-1)*n3C + kl
                           TInt(kl_B) = AInt
*
                        End Do
                     End Do
                  End Do
*                                                                      *
************************************************************************
*                                                                      *
               End If
*                                                                      *
************************************************************************
*                                                                      *
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iCmp)
         Call Unused_logical(Shijij)
         Call Unused_integer(iBas)
      End If
      End
