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
      SubRoutine Cho_SimRI_Z1CDia(Diag,Thr,Indx)
C
C     Purpose: Zero 1-center diagonals smaller then Thr.
C              Diag is the diagonal (1st reduced set).
C              On exit, Indx(i)=1 if diagonal i was zeroed, else
C              Indx(i)=0 (thus, Indx must have same dimension as Diag).
C
      Implicit None
      Real*8  Diag(*)
      Real*8  Thr
      Integer Indx(*)
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer iShlAB, iShlA, iShlB
      Integer iAB1, iAB2, iAB
      Integer n
      Real*8  zmx

      Real*8 Zero
      Parameter (Zero = 0.0d0)

      Integer Inf_SimRI
      Parameter (Inf_SimRI = 0)

      Integer iSP2F, iAtomShl, iiBstRSh, nnBstRSh
      Integer i, j, k
      iSP2F(i)=iWork(ip_iSP2F-1+i)
      iAtomShl(i)=iWork(ip_iAtomShl-1+i)
      iiBstRSh(i,j,k)=iWork(ip_iiBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)
      nnBstRSh(i,j,k)=iWork(ip_nnBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)

      Call Cho_iZero(Indx,nnBstR(1,1))

      zmx = 0.0d0
      n = 0
      Do iShlAB = 1,nnShl
         Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.True.)
         If (iAtomShl(iShlA) .eq. iAtomShl(iShlB)) Then
            iAB1 = iiBstR(1,1)
     &           + iiBstRSh(1,iShlAB,1) + 1
            iAB2 = iAB1 + nnBstRSh(1,iShlAB,1) - 1
            Do iAB = iAB1,iAB2
               If (Diag(iAB) .lt. Thr) Then
                  zmx = max(zmx,Diag(iAB))
                  n = n + 1
                  Indx(iAB) = 1
                  Diag(iAB) = 0.0d0
               End If
            End Do
         End If
      End Do

      If (iPrint .gt. Inf_SimRI) Then
         Write(LuPri,'(/,A,I7,A,1P,D10.2,A)')
     &   'Simulating RI:',n,' 1-center diagonals < ',Thr,
     &   ' have been zeroed'
         If (n .gt. 0) Then
            Write(LuPri,'(A,1P,D15.7)')
     &      'Largest zeroed diagonal: ',zmx
         End If
      End If

      End
