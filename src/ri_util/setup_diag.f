************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Setup_Diag(nBas,nIrrep,nShell,iShij,nij_Shell,
     &                      ip_iSSOff,iSOInf,MaxAO,ip_ShlSO,
     &                      ip_SOShl,ip_nBasSh)
      use iSD_data
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "setup.fh"
#include "nsd.fh"
      Integer nBas(0:nIrrep-1), iShij(2,nij_Shell), iSOInf(3,4*MaxAO)
*                                                                      *
************************************************************************
*                                                                      *
      nSO = 0
      Do iIrrep = 0, nIrrep-1
         nSO = nSO + nBas(iIrrep)
      End Do
*
      Call GetMem('SOShl','Allo','Integer',ip_SOShl,nSO)
      Call GetMem('ShlSO','Allo','Integer',ip_ShlSO,nSO)
      Call GetMem('nBasSh','Allo','Integer',ip_nBasSh,nShell*nIrrep)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate SOShl, SO -> shell index
*
      Do iSO = 1, nSO
         iCnttp=iSOInf(1,iSO)
         iCnt  =iSOInf(2,iSO)
         iAng  =iSOInf(3,iSO)
*
*        Find the Shell from which this basis function is derived.
*
         Do iSkal = 1, nShell
            jCnttp=iSD(13,iSkal)
            jCnt  =iSD(14,iSkal)
            jAng  =iSD( 1,iSkal)
            If (jCnttp.eq.iCnttp .and.
     &          jCnt  .eq.iCnt   .and.
     &          jAng  .eq.iAng         ) Then
               iWork(ip_SOShl-1+iSO)=iSkal
               Go To 99
            End If
         End Do
 99      Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      nSSOff = nIrrep**2 * nij_Shell
      Call GetMem('iSSOff','Allo','Integer',ip_iSSOff,nSSOff)
*                                                                      *
************************************************************************
*                                                                      *
      Call Setup_Diag_(iWork(ip_SOShl),nSO,iWork(ip_ShlSO),
     &                 iWork(ip_nBasSh),nShell,nIrrep,nBas,
     &                 iWork(ip_iSSOff),nij_Shell,iShij)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine Setup_Diag_(iSOShl,nSO,iShlSO,nBasSh,nShell,nIrrep,
     &                       nBas,iSSOff,nij_Shell,iShij)
      Implicit Real*8 (a-h,o-z)
      Integer nBasSh(0:nIrrep-1,nShell), iSOShl(nSO), iShlSO(nSO),
     &        nBas(0:nIrrep-1), iShij(2,nij_Shell), iTtmp(0:7),
     &        iSSOff(0:nIrrep-1,0:nIrrep-1,nij_Shell)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate index array for relative index within the shell and irrep
*
      iSO = 0
      Do iIrrep = 0, nIrrep-1
         Do iShell = 1, nShell
*
            iSO_Shl = 0
            Do iBas = iSO+1 , iSO+nBas(iIrrep)
               If (iSOShl(iBas).eq.iShell) Then
                  iSO_Shl = iSO_Shl + 1
*                                                                      *
************************************************************************
*                                                                      *
*                 Save the relative index within the shell and irrep
*                 of a given absolute SO index.
*
                  iShlSO(iBas) = iSO_Shl
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Save the total number of basis functions a specific shell
*           has in a given irrep.
*
            nBasSh(iIrrep,iShell)=iSO_Shl
*                                                                      *
************************************************************************
*                                                                      *
         End Do
         iSO = iSO + nBas(iIrrep)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      iAcc = 0
      Do ijShell = 1, nij_Shell
         iShell = iShij(1,ijShell)
         jShell = iShij(2,ijShell)
         Call ICopy(nIrrep,0,0,iTtmp,1)
*                                                                      *
************************************************************************
*                                                                      *
         If (iShell.gt.jShell) Then   ! iShell > jShell
*
            Do jIrrep = 0, nIrrep-1
               nB = nBasSh(jIrrep,jShell)
               Do iIrrep = 0, nIrrep-1
                  nA = nBasSh(iIrrep,iShell)
*
                  ijIrrep = iEor(iIrrep,jIrrep)
                  iSSOff(iIrrep,jIrrep,ijShell)=iTtmp(ijIrrep)
                  nab = na*nb
                  iTtmp(ijIrrep)=iTtmp(ijIrrep)+nab
*
               End Do
            End Do
*
         Else                        ! iShell = jShell
*
            Do iIrrep = 0, nIrrep-1
               nA = nBasSh(iIrrep,iShell)
               Do jIrrep = 0, iIrrep
                  nB = nBasSh(jIrrep,jShell)
*
                  ijIrrep = iEor(iIrrep,jIrrep)
                  iSSOff(iIrrep,jIrrep,ijShell)=iTtmp(ijIrrep)
                  iSSOff(jIrrep,iIrrep,ijShell)=iTtmp(ijIrrep)
                  nab = na*nb
                  If (iIrrep.eq.jIrrep) nab = na*(na+1)/2
                  iTtmp(ijIrrep) = iTtmp(ijIrrep) + nab
               End Do
            End Do
*
         End If
*
*        Now update the index to be the total offset within a slice
*        for a fixed shell-pair
*
         Do ijIrrep = 0, nIrrep-1
            Do iIrrep = 0, nIrrep-1
               jIrrep = iEor(ijIrrep,iIrrep)
               iSSOff(iIrrep,jIrrep,ijShell)=
     &                            iSSOff(iIrrep,jIrrep,ijShell)+iAcc
            End Do
            iAcc = iAcc + iTtmp(ijIrrep)
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do   ! ijShell = 1, nij_Shell
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
