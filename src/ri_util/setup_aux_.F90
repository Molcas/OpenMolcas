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
      Subroutine Setup_Aux_(iSOShl,nSO,iShlSO,nBasSh,nShell,            &
     &                              nIrrep,nBas,iSSOff,nij_Shell,iShij, &
     &                              nBas_Aux,nChV,iTOffs)
      Implicit Real*8 (a-h,o-z)
      Integer iSOShl(nSO), iShlSO(nSO), nBasSh(0:nIrrep-1,nShell),      &
     &        nBas(0:nIrrep-1), nBas_Aux(0:nIrrep-1), nChV(0:nIrrep-1), &
     &        iSSOff(0:nIrrep-1,0:nIrrep-1,nij_Shell),                  &
     &        iShij(2,nij_Shell), iTOffs(3,0:nIrrep-1), iTtmp(0:7)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Generate index array for relative index within the shell and irrep
!
!     Call iVcPrt('iSOShl',' ',iSOShl,nSO)
      iSO = 0
      Do iIrrep = 0, nIrrep-1
         Do iShell = 1, nShell
!
            iSO_Shl = 0
            Do iBas = iSO+1 , iSO+nBas(iIrrep)
               If (iSOShl(iBas).eq.iShell) Then
                  iSO_Shl = iSO_Shl + 1
!                                                                      *
!***********************************************************************
!                                                                      *
!                 Save the relative index within the shell and irrep
!                 of a given absolute SO index.
!
                  iShlSO(iBas) = iSO_Shl
               End If
            End Do
!                                                                      *
!***********************************************************************
!                                                                      *
!           Save the total number of basis functions a specific shell
!           has in a given irrep.
!
            nBasSh(iIrrep,iShell)=iSO_Shl
!                                                                      *
!***********************************************************************
!                                                                      *
         End Do
         iSO = iSO + nBas(iIrrep)
      End Do
!                                                                      *
!***********************************************************************
!
!     Initialize
!
      Call ICopy(3*nIrrep,[0],0,iTOffs,1)
      Call iCopy(nij_Shell*nIrrep**2,[0],0,iSSOff,1)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute offsets within the symmetry block for a fixed shell pair.
!
!     Note that for each pair of valence shells all the products which
!     are of the same irrep are consecutive.
!
!
      Do ijShell = 1, nij_Shell
         iShell = iShij(1,ijShell)
         jShell = iShij(2,ijShell)
         Call ICopy(nIrrep,[0],0,iTtmp,1)
!                                                                      *
!***********************************************************************
!                                                                      *
         If (iShell.gt.jShell) Then   ! iShell > jShell
!
            Do jIrrep = 0, nIrrep-1
               nB = nBasSh(jIrrep,jShell)
               Do iIrrep = 0, nIrrep-1
                  nA = nBasSh(iIrrep,iShell)
!
                  ijIrrep = iEor(iIrrep,jIrrep)
                  iSSOff(iIrrep,jIrrep,ijShell)=iTtmp(ijIrrep)
                  nab = na*nb
                  iTtmp(ijIrrep)=iTtmp(ijIrrep)+nab
!
               End Do
            End Do
!
         Else                        ! iShell = jShell
!
            Do iIrrep = 0, nIrrep-1
               nA = nBasSh(iIrrep,iShell)
               Do jIrrep = 0, iIrrep
                  nB = nBasSh(jIrrep,jShell)
!
                  ijIrrep = iEor(iIrrep,jIrrep)
                  iSSOff(iIrrep,jIrrep,ijShell)=iTtmp(ijIrrep)
                  iSSOff(jIrrep,iIrrep,ijShell)=iTtmp(ijIrrep)
                  nab = na*nb
                  If (iIrrep.eq.jIrrep) nab = na*(na+1)/2
                  iTtmp(ijIrrep) = iTtmp(ijIrrep) + nab
               End Do
            End Do
!
         End If
!
         Do iIrrep = 0, nIrrep-1
            iTOffs(3,iIrrep) = iTOffs(3,iIrrep) + iTtmp(iIrrep)
         End Do
!
!        Now update the index to be the total offset within a slice
!        for a fixed shell-pair
!
         iAcc = 0
         Do ijIrrep = 0, nIrrep-1
            Do iIrrep = 0, nIrrep-1
               jIrrep = iEor(ijIrrep,iIrrep)
               iSSOff(iIrrep,jIrrep,ijShell)=                           &
     &                            iSSOff(iIrrep,jIrrep,ijShell)+iAcc
            End Do
            nI = nBas_Aux(ijIrrep)
            If (ijIrrep.eq.0) nI=nI-1
            iAcc = iAcc + nI*iTtmp(ijIrrep)
         End Do
!                                                                      *
!***********************************************************************
!                                                                      *
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*)  'iSSOff'
      Write (6,*)
      Do ijShell = 1, nij_Shell
         iShell = iShij(1,ijShell)
         jShell = iShij(2,ijShell)
         Write (6,*)
         Write (6,*) 'iShell,jShell=',iShell,jShell
         Write (6,*)
         Do i = 0, nIrrep-1
            Write (6,'(8I4)') (iSSOff(i,j,ijShell),j=0,nIrrep-1)
         End Do
      End Do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Set up pointers for the J12 matrix and compute total size of the
!     3-center integrals.
!
      iOff_V12=0
      Do iIrrep = 0, nIrrep-1
         iTOffs(1,iIrrep) = nChV(iIrrep) ! # of vectors
         nAux = nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nAux = nAux - 1
         iTOffs(2,iIrrep) = iOff_V12
         iOff_V12 = iOff_V12 + nAux**2
      End Do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) ' iSO, iShlSO(iSO), relative index in irrep'
      Do jSO = 1, iSO
         Write (6,*) jSO, iShlSO(jSO)
      End Do
!
      Write (6,*)
      Write (6,*) ' iShell: number of basis functions in each irrep'
      Do iShell = 1, nShell
         Write (6,*) iShell,':',                                        &
     &         (nBasSh(iIrrep,iShell),iIrrep=0,nIrrep-1)
      End Do
      Write (6,*)
      Write (6,*)  'iSSOff'
      Write (6,*)
      Do ijShell = 1, nij_Shell
         iShell = iShij(1,ijShell)
         jShell = iShij(2,ijShell)
         Write (6,*)
         Write (6,*) 'iShell,jShell=',iShell,jShell
         Write (6,*)
         Do i = 0, nIrrep-1
            Write (6,'(8I4)') (iSSOff(i,j,ijShell),j=0,nIrrep-1)
         End Do
      End Do
      Write (6,*)
      Write (6,*)  'iTOffs'
      Write (6,*)
      Do i = 0, nIrrep-1
         Write (6,*) (iTOffs(j,i),j=1,3)
      End Do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End
