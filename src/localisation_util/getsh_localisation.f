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
      SubRoutine GetSh_Localisation(X,nBas,m,XSh,nShell,iSO2Sh,iOpt,
     &                              Norm)
      Implicit Real*8 (a-h,o-z)
      Real*8  X(nBas,m), XSh(nShell,*)
      Integer iSO2Sh(*)
      Character*3 Norm  ! 'MAX' or 'FRO'

      Character*3 myNorm

      If (nBas.lt.1 .or. nShell.lt.1) Return

      myNorm = Norm
      Call UpCase(myNorm)

      If (iOpt .eq. 1) Then
         Call dCopy_(nShell*m,[0.0d0],0,XSh,1)
         If (myNorm .eq. 'MAX') Then
            Do j = 1,m
               Do i = 1,nBas
                  iShell = iSO2Sh(i)
                  XSh(iShell,j) = max(abs(X(i,j)),XSh(iShell,j))
               End Do
            End Do
         Else If (myNorm .eq. 'FRO') Then
            Do j = 1,m
               Do i = 1,nBas
                  iShell = iSO2Sh(i)
                  XSh(iShell,j) = XSh(iShell,j) + X(i,j)**2
               End Do
               Do iShell = 1,nShell
                  XSh(iShell,j) = sqrt(XSh(iShell,j))
               End Do
            End Do
         End If
      Else
         If (m .ne. nBas) Then
            Call SysAbendMsg('GetSh_Localisation','Fatal error',
     &                       'm != nBas')
         End If
         Call dCopy_(nShell*nShell,[0.0d0],0,XSh,1)
         If (myNorm .eq. 'MAX') Then
            Do j = 1,nBas
               jShell = iSO2Sh(j)
               Do i = 1,nBas
                  iShell = iSO2Sh(i)
                  XSh(iShell,jShell) =
     &               max(abs(X(i,j)),XSh(iShell,jShell))
               End Do
            End Do
         Else If (myNorm .eq. 'FRO') Then
            Do j = 1,nBas
               jShell = iSO2Sh(j)
               Do i = 1,nBas
                  iShell = iSO2Sh(i)
                  XSh(iShell,jShell) = XSh(iShell,jShell) + X(i,j)**2
               End Do
            End Do
            Do jShell = 1,nShell
               Do iShell = 1,nShell
                  Xsh(iShell,jShell) = sqrt(Xsh(iShell,jShell))
               End Do
            End Do
         End If
      End If

      End
