!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Setup_Aux(nIrrep,nBas,nShell,nShell_Aux,nSO,           &
     &                     TMax,CutOff,nij_Shell,                       &
     &                     nBas_Aux,nChV,iTOffs)
      use iSD_data
      use SOAO_Info, only: iSOInf
      use j12, only: ShlSO, SOShl, nBasSh, iSSOff, iShij
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "setup.fh"
#include "nsd.fh"
      Integer nBas(0:nIrrep-1), nBas_Aux(0:nIrrep-1),                   &
     &        nChV(0:nIrrep-1), iTOffs(3,0:nIrrep-1)
      Real*8  TMax(nShell,nShell)
!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*) 'Setup_Aux:nIrrep:   ',nIrrep
      Write (6,*) 'Setup_Aux:nBas:     ',nBas
      Write (6,*) 'Setup_Aux:nBas_Aux: ',nBas_Aux
      Write (6,*) 'Setup_Aux:nChV:     ',nChV
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      nSO = 0
      nSO_Aux = 0
      Do iIrrep = 0, nIrrep-1
         nSO = nSO + nBas(iIrrep)
         nSO_Aux = nSO_Aux + nBas_Aux(iIrrep)
      End Do
!
      Call mma_allocate(SOShl,nSO+nSO_Aux,Label='SOShl')
      Call mma_allocate(ShlSO,nSO+nSO_Aux,Label='ShlSO')
      Call mma_allocate(nBasSh,[0,nIrrep-1],                            &
     &                         [1,nShell+nShell_Aux],Label='nBasSh')
!                                                                      *
!***********************************************************************
!                                                                      *
      Do iSO = 1, nSO+nSO_Aux
         iCnttp=iSOInf(1,iSO)
         iCnt  =iSOInf(2,iSO)
         iAng  =iSOInf(3,iSO)
!        Write (*,*) 'iCnttp,iCnt,iAng=',iCnttp,iCnt,iAng
!
!        Find the Shell from which this basis function is derived.
!
         Do iSkal = 1, nShell+nShell_Aux
            jCnttp=iSD(13,iSkal)
            jCnt  =iSD(14,iSkal)
            jAng  =iSD( 1,iSkal)
            If (jCnttp.eq.iCnttp .and.                                  &
     &          jCnt  .eq.iCnt   .and.                                  &
     &          jAng  .eq.iAng         ) Then
               SOShl(iSO)=iSkal
!              Write (*,*) 'Found in shell=',iSkal
               Go To 99
            End If
         End Do
 99      Continue
      End Do
!     Call iVcPrt('SOShl',' ',SOShl,nSO+nSO_Aux)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute the number of effective shell pairs.
!
      TMax_ij=Zero
      Do iSkal = 1, nShell
         Do jSkal = 1, iSkal
            TMax_ij=Max(TMax_ij,TMax(iSkal,jSkal))
         End Do
      End Do
!
      nij_Shell=0
      Do iSkal = 1, nShell
         Do jSkal = 1, iSkal
            If (TMax(iSkal,jSkal)*TMax_ij.ge.CutOff) Then
               nij_Shell = nij_Shell + 1
            End If
         End Do
      End Do
      Call mma_allocate(iShij,2,nij_Shell,Label='iShij')
!
      ij_Shell = 0
      Do iSkal = 1, nShell
         Do jSkal = 1, iSkal
            If (TMax(iSkal,jSkal)*TMax_ij.ge.CutOff) Then
               ij_Shell = ij_Shell + 1
               iShij(1,ij_Shell)=iSkal
               iShij(2,ij_Shell)=jSkal
#ifdef _DEBUGPRINT_
               Write (6,*) 'ij_Shell,iSkal,jSkal=',                     &
     &                      ij_Shell,iSkal,jSkal
#endif
            End If
         End Do
      End Do
!                                                                      *
!***********************************************************************
!                                                                      *
      Call mma_allocate(iSSOff,[0,nIrrep-1],[0,nIrrep-1],               &
     &                  [1,nij_Shell],Label='iSSOff')
!                                                                      *
!***********************************************************************
!                                                                      *
!
      Call Setup_Aux_(SOShl,nSO+nSO_Aux,ShlSO,                          &
     &                        nBasSh,nShell+nShell_Aux,nIrrep,nBas,     &
     &                        iSSOff,nij_Shell,iShij,                   &
     &                        nBas_Aux,nChV,iTOffs)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*) 'nij_Shell=',nij_Shell
      Write (6,*)
      Do ij_Shell = 1, nij_Shell
         Write (6,*) iShij(1,ij_Shell),iShij(2,ij_Shell)
      End Do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End
