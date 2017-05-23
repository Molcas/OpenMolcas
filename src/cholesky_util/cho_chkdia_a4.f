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
      SubRoutine Cho_ChkDia_A4(Diag,Dmax,iSym,nNeg,nNegT,nConv,xM,yM,zM)
C
C     Purpose: check for negative diagonals.
C              Dmax is the max. diagonal (global), used for screening.
C              On exit,
C              nNeg = #negative zeroed diagonals.
C              nNegT = #negative diagonals
C              nConv = #diagonal elements < ThrCom
C              xM = max. element in Diag
C              yM = min. element in Diag
C              zM = max. abs. element in Diag
C
      Implicit None
      Real*8  Diag(*), Dmax
      Integer iSym, nNeg, nNegT, nConv
      Real*8  xM, yM, zM
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*13 SecNam
      Parameter (SecNam = 'Cho_ChkDia_A4')

      Real*8  Tst
      Integer i, j, j1, j2

      Integer IndRed
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)

      nNeg = 0
      nNegT = 0
      nConv = 0
      If (nnBstR(iSym,2) .gt. 0) Then
         xM = -9.9d9
         yM = 9.9d9
      Else
         xM = 0.0d0
         yM = 0.0d0
      End If

      j1 = iiBstR(iSym,2) + 1
      j2 = j1 + nnBstR(iSym,2) - 1
      Do j = j1,j2
         i = IndRed(j,2) ! addr in 1st reduced set
         xM = max(xM,Diag(i))
         yM = min(yM,Diag(i))
         If (Diag(i) .lt. 0.0d0) Then
            nNegT = nNegT + 1
            If (Diag(i) .lt. ThrNeg) Then
               nNeg = nNeg + 1
               If (Diag(i) .lt. TooNeg) Then
                  Write(Lupri,'(A,A,I12,1X,1P,D16.8)')
     &            SecNam,': diagonal too negative: ',i,Diag(i)
                  Write(Lupri,'(A,A)')
     &            SecNam,': shutting down Cholesky decomposition!'
                  Call Cho_Quit('Diagonal too negative in '//SecNam,104)
               End If
               If (Diag(i) .lt. WarNeg) Then
                  Write(Lupri,'(A,A,I12,1X,1P,D16.8,A)')
     &            SecNam,': Negative diagonal: ',i,Diag(i),
     &            ' (zeroed)'
               End If
               Diag(i) = 0.0d0
            End If
         End If
      End Do

      zM = max(abs(xM),abs(yM))

      Do j = j1,j2
         i = IndRed(j,2)
         Tst = sqrt(abs(Diag(i))*Dmax)*Damp(2)
         If (Tst .le. ThrCom) Then
            nConv = nConv + 1
            If (ScDiag) Diag(i) = 0.0d0
         End If
      End Do

      End
