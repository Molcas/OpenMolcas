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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine GetU_ER(U,R,n)
C
C     Thomas Bondo Pedersen, November 2005.
C
C     Purpose: compute U = R*[R^T*R]^(-1/2).
C
C     (used by ER orbital localisation - hence the _ER)
C
      Implicit None
      Integer n
      Real*8 U(n,n), R(n,n)
#include "WrkSpc.fh"

      Integer nn,n2
      Integer ipRTR, lRTR
      Integer ipSqrt, lSqrt, ipISqrt, lISqrt
      Integer ipScr, lScr
      Integer iTask

      If (n .lt. 1) Return

C     Allocations.
C     ------------

      nn = n*(n+1)/2
      n2 = n**2

      lRTR = n2
      lSqrt = n2
      lISqrt = n2
      lScr = 2*n2 + nn
      Call GetMem('RTR','Allo','Real',ipRTR,lRTR)
      Call GetMem('Sqrt','Allo','Real',ipSqrt,lSqrt)
      Call GetMem('ISqrt','Allo','Real',ipISqrt,lISqrt)
      Call GetMem('Scr','Allo','Real',ipScr,lScr)

C     Compute R^T*R.
C     --------------

      Call DGEMM_('T','N',n,n,n,1.0d0,R,n,R,n,0.0d0,Work(ipRTR),n)

C     Compute inverse square root of R^T*R.
C     -------------------------------------

      iTask = 2 ! compute sqrt as well as inverse sqrt
      Call SqrtMt(Work(ipRTR),n,iTask,Work(ipSqrt),Work(ipISqrt),
     &            Work(ipScr))

C     Compute U.
C     ----------

      Call DGEMM_('N','N',n,n,n,1.0d0,R,n,Work(ipISqrt),n,0.0d0,U,n)

C     De-allocations.
C     ---------------

      Call GetMem('Scr','Free','Real',ipScr,lScr)
      Call GetMem('ISqrt','Free','Real',ipISqrt,lISqrt)
      Call GetMem('Sqrt','Free','Real',ipSqrt,lSqrt)
      Call GetMem('RTR','Free','Real',ipRTR,lRTR)

      End
