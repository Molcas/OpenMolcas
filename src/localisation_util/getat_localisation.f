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
      SubRoutine GetAt_Localisation(X,nBas,m,XAt,nAtoms,iOpt,
     &                              nBas_per_Atom,nBas_Start,
     &                              Norm)
      Implicit Real*8 (a-h,o-z)
      Real*8  X(nBas,m), XAt(nAtoms,*)
      Integer nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
      Character*3 Norm  ! 'MAX' or 'FRO'

      Character*3 myNorm

      If (nBas.lt.1 .or. nAtoms.lt.1) Return

      myNorm = Norm
      Call UpCase(myNorm)

      If (iOpt .eq. 1) Then
         Call dCopy_(nAtoms*m,0.0d0,0,XAt,1)
         If (myNorm .eq. 'MAX') Then
            Do j = 1,m
               Do iAt = 1,nAtoms
                  i1 = nBas_Start(iAt)
                  i2 = i1 + nBas_per_Atom(iAt) - 1
                  Do i = i1,i2
                     XAt(iAt,j) = max(abs(X(i,j)),XAt(iAt,j))
                  End Do
               End Do
            End Do
         Else If (myNorm .eq. 'FRO') Then
            Do j = 1,m
               Do iAt = 1,nAtoms
                  i1 = nBas_Start(iAt)
                  i2 = i1 + nBas_per_Atom(iAt) - 1
                  Do i = i1,i2
                     XAt(iAt,j) = XAt(iAt,j) + X(i,j)**2
                  End Do
                  XAt(iAt,j) = sqrt(XAt(iAt,j))
               End Do
            End Do
         End If
      Else
         If (m .ne. nBas) Then
            Call SysAbendMsg('GetAt_Localisation','Fatal error',
     &                       'm != nBas')
         End If
         Call dCopy_(nAtoms*nAtoms,0.0d0,0,XAt,1)
         If (myNorm .eq. 'MAX') Then
            Do jAt = 1,nAtoms
               j1 = nBas_Start(jAt)
               j2 = j1 + nBas_per_Atom(jAt) - 1
               Do j = j1,j2
                  Do iAt = 1,nAtoms
                     i1 = nBas_Start(iAt)
                     i2 = i1 + nBas_per_Atom(iAt) - 1
                     Do i = i1,i2
                        XAt(iAt,jAt) = max(abs(X(i,j)),XAt(iAt,jAt))
                     End Do
                  End Do
               End Do
            End Do
         Else If (myNorm .eq. 'FRO') Then
            Do jAt = 1,nAtoms
               j1 = nBas_Start(jAt)
               j2 = j1 + nBas_per_Atom(jAt) - 1
               Do j = j1,j2
                  Do iAt = 1,nAtoms
                     i1 = nBas_Start(iAt)
                     i2 = i1 + nBas_per_Atom(iAt) - 1
                     Do i = i1,i2
                        XAt(iAt,jAt) = XAt(iAt,jAt) + X(i,j)**2
                     End Do
                  End Do
               End Do
            End Do
            Do jAt = 1,nAtoms
               Do iAt = 1,nAtoms
                  XAt(iAt,jAt) = sqrt(XAt(iAt,jAt))
               End Do
            End Do
         End If
      End If

      End
