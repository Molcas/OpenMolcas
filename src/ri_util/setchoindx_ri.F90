!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
      SubRoutine SetChoIndx_RI(iiBstRSh,nnBstRSh,IndRed,IndRsh,iRS2F,   &
     &                         I_nSym,I_nnShl,I_mmBstRT,I_3,I_2,        &
     &                         iShij,nShij)
      use ChoArr, only: iSP2F, iBasSh, nBasSh, nBstSh
      Implicit Real*8 (a-h,o-z)
      Integer iiBstRSh(I_nSym,I_nnShl,I_3), nnBstRSh(I_nSym,I_nnShl,I_3)
      Integer IndRed(I_mmBstRT,I_3), IndRsh(I_mmBstRT)
      Integer iRS2F(I_2,I_mmBstRT), iShij(2,nShij)
#include "choorb.fh"
#include "cholesky.fh"

      Integer  Cho_iSAOSh
      External Cho_iSAOSh

      Integer iRS(8)

      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

!     nnBstRSh(iSym,iSh_ij,1) = #elements in compound sym. iSym of
!                               shell-pair ab in 1st reduced set.
!     IndRSh(jRS): shell-pair to which element jRS of first reduced set
!                  belongs.
!     IndRed(jRS,1): address (without symmetry) in shell-pair of element
!                    jRS of first reduced set.
!     ------------------------------------------------------------------

      Call iCopy(nSym*nnShl,[0],0,nnBstRSh(1,1,1),1)
      Call iCopy(nSym,iiBstR(1,1),1,iRS,1)
      Do iSh_ij= 1,nShij
         iShla=iShij(1,iSh_ij)
         iShlb=iShij(2,iSh_ij)
         iShlab=iTri(iShla,iShlb)
!        Write (*,*) 'iSh_ij,iShlab,iShla,iShlb=',
!    &               iSh_ij,iShlab,iShla,iShlb
         If (iShlab .ne. iSP2F(iSh_ij)) Then
            Call SysAbendMsg('SetChoIndx_RI','SP2F setup error',' ')
         End If
!
         If (iShla.gt.iShlb) Then
!
!        code for shell a > shell b
!
            Do iSymb = 1,nSym
               Do ibb = 1,nBasSh(iSymb,iShlb)
                  ib = iBasSh(iSymb,iShlb) + ibb
                  Do iSyma = 1,nSym
                     iSym = MulD2h(iSyma,iSymb)
                     Do iaa = 1,nBasSh(iSyma,iShla)
                        ia = iBasSh(iSyma,iShla) + iaa
                        iab = nBstSh(iShla)*(ib-1) + ia
                        nnBstRSh(iSym,iSh_ij,1) =                       &
     &                     nnBstRSh(iSym,iSh_ij,1) + 1
                        iRS(iSym) = iRS(iSym) + 1
                        IndRSh(iRS(iSym)) = iShlab
                        IndRed(iRS(iSym),1) = iab
                     End Do
                  End Do
               End Do
            End Do
!
         Else
!
!            code for shell a = shell b follows
!
            Do ia = 1,nBstSh(iShla)
               iSyma = Cho_iSAOSh(ia,iShla)
               Do ib = 1,ia
                  iab = iTri(ia,ib)
                  iSymb = Cho_iSAOSh(ib,iShlb)
                  iSym = MulD2h(iSyma,iSymb)
                  nnBstRSh(iSym,iSh_ij,1) = nnBstRSh(iSym,iSh_ij,1) + 1
                  iRS(iSym) = iRS(iSym) + 1
                  IndRSh(iRS(iSym)) = iShlab
                  IndRed(iRS(iSym),1) = iab
               End Do
            End Do
!
         End If
      End Do   ! iSh_ij

!     Check.
!     ------

      nErr = 0
      Do iSym = 1,nSym
         iCount = nnBstRSh(iSym,1,1)
         Do iSh_ij = 2,nnShl
            iCount = iCount + nnBstRSh(iSym,iSh_ij,1)
         End Do
         If (iCount .ne. nnBstR(iSym,1)) Then
            nErr = nErr + 1
         End If
      End Do
      If (nErr .ne. 0) Then
         Call SysAbendMsg('SetChoIndx_RI','Setup error',                &
     &                    'iCount vs. nnBstR')
      End If
      Do iSym = 1,nSym
         If ((iRS(iSym)-iiBstR(iSym,1)) .ne. nnBstR(iSym,1)) Then
            nErr = nErr + 1
         End If
      End Do
      If (nErr .ne. 0) Then
         Call SysAbendMsg('SetChoIndx_RI','Setup error','ShP RS1 count')
      End If

!     iiBstRSh(iSym,iSh_ij,1) = offset to elements in compound sym. iSym
!                               of shell-pair ab in 1st reduced set.
!     ------------------------------------------------------------------

      Do iSym = 1,nSym
         iiBstRSh(iSym,1,1) = 0
         Do iSh_ij = 2,nnShl
            iiBstRSh(iSym,iSh_ij,1) = iiBstRSh(iSym,iSh_ij-1,1)         &
     &                              + nnBstRSh(iSym,iSh_ij-1,1)
         End Do
      End Do

!     Check.
!     ------

      nErr = 0
      Do iSym = 1,nSym
         Do iSh_ij = 1,nnShl
            jRS1 = iiBstR(iSym,1) + iiBstRSh(iSym,iSh_ij,1) + 1
            jRS2 = jRS1 + nnBstRSh(iSym,iSh_ij,1) - 1
            Do jRS = jRS1,jRS2
               If (IndRSh(jRS) .ne. iSP2F(iSh_ij)) Then
                  nErr = nErr + 1
               End If
            End Do
         End Do
      End Do
      If (nErr .ne. 0) Then
         Call SysAbendMsg('SetChoIndx_RI','Setup error','IndRSh')
      End If

!     Copy index arrays to "locations" 2 and 3.
!     Note: IndRed here returns the index in 1st reduced set.
!     -------------------------------------------------------

      Do i = 2,3
         Do jRS = 1,nnBstRT(1)
            IndRed(jRS,i) = jRS
         End Do
         Call iCopy(nSym*nnShl,iiBstRSh(1,1,1),1,iiBstRSh(1,1,i),1)
         Call iCopy(nSym*nnShl,nnBstRSh(1,1,1),1,nnBstRSh(1,1,i),1)
      End Do

      Call Cho_RStoF(iRS2F,2,nnBstRT(1),1)

      Return
      End
