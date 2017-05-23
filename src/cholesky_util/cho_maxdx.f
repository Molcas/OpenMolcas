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
      SubRoutine Cho_MaxDX(Diag,Dmax)
C
C     Purpose: get max. diagonal elements in each sym. block,
C              qualified diagonals excluded.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 Diag(*), Dmax(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      iQuAB(i,j)=iWork(ip_iQuAB-1+MaxQual*(j-1)+i)
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)

      MxQ = nQual(1)
      Do jSym = 2,nSym
         MxQ = max(MxQ,nQual(jSym))
      End Do
      Call GetMem('DXQ','Allo','Real',ipExQ,MxQ)

      Do jSym=1,nSym

         Dmax(jSym) = 0.0d0
         If (nQual(jSym) .lt. 1) goto 10  ! next symm

         Do iQ=1,nQual(jSym)
            iab=IndRed(iQuAB(iQ,jSym),2) ! addr in 1st red set
            Work(ipExQ+iQ-1) = Diag(iab)
            Diag(iab) = 0.0d0
         End Do

         jRab1 = iiBstr(jSym,2) + 1
         jRab2 = jRab1 + nnBstR(jSym,2) - 1
         Do jRab=jRab1,jRab2
            iRab = IndRed(jRab,2) ! addr in 1st red set
            Dmax(jSym) = Max(Dmax(jSym),Diag(iRab))
         End Do

C --- Restore the qualified
         Do iQ=1,nQual(jSym)
            iab=IndRed(iQuAB(iQ,jSym),2) ! addr in 1st red set
            Diag(iab) = Work(ipExQ+iQ-1)
         End Do

10       Continue

      End Do

      Call GetMem('DXQ','Free','Real',ipExQ,MxQ)

      End
