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
      SubRoutine Chk_Unitary(irc,U,n,Thr)
C
C     Thomas Bondo Pedersen, November 2005.
C
C     Purpose: check that U is unitary.
C
      Implicit None
      Integer irc, n
      Real*8  U(n,n), Thr
#include "WrkSpc.fh"

      Integer n2, ipUTU, lUTU, i, ip0
      Real*8  RMS, x2

      real*8 ddot_
      external ddot_

      If (n .lt. 1) Then
         irc = 0
         Return
      End If

      n2 = n**2
      lUTU = n2
      Call GetMem('UTU','Allo','Real',ipUTU,lUTU)

      Call dCopy_(n2,[0.0d0],0,Work(ipUTU),1)
      ip0 = ipUTU - 1
      Do i = 1,n
         Work(ip0+n*(i-1)+i) = 1.0d0
      End Do
      Call DGEMM_('T','N',n,n,n,-1.0d0,U,n,U,n,1.0d0,Work(ipUTU),n)

      x2 = dble(n2)
      RMS = sqrt(dDot_(n2,Work(ipUTU),1,Work(ipUTU),1)/x2)
      If (RMS .gt. Thr) Then
         irc = 1
      Else
         irc = 0
      End If

      Call GetMem('UTU','Free','Real',ipUTU,lUTU)

      End
