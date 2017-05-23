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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine Anasize_Localisation(Den,CMO,XMO,nShell,nOrb,iSym)
C
C     Author: T.B. Pedersen
C
C     Purpose: sparsity analysis of shell-based matrices.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 Den(nShell,nShell), CMO(nShell,nOrb), XMO(nShell,nOrb)
#include "WrkSpc.fh"

      Character*17 XHead
      Character*20 CHead
      Character*36 DHead

C     Return if nothing to do.
C     ------------------------

      If (nShell .lt. 0) Return

C     Set up bins.
C     ------------

      lBin = 9
      Call GetMem('Bin','Allo','Real',ipBin,lBin)
      StpSiz = 1.0d-1
      Work(ipBin) = 1.0d0
      ip0 = ipBin - 1
      Do iBin = 2,lBin
         Work(ip0+iBin) = Work(ip0+iBin-1)*StpSiz
      End Do

C     Density.
C     --------

      lDLT = nShell*(nShell+1)/2
      Call GetMem('LTDen','Allo','Real',ipDLT,lDLT)
      Call Sq2Tri(Den,Work(ipDLT),nShell)
      Write(DHead,'(A34,I2)') 'Histogram of density matrix , sym.',iSym
      Call Cho_Head(DHead,'=',80,6)
      Call Cho_Anasize(Work(ipDLT),lDlt,Work(ipBin),lBin,6)
      Call GetMem('LTDen','Free','Real',ipDLT,lDLT)

      If (nOrb .lt. 1) Go To 1 ! return after de-allocation

C     Original MOs.
C     -------------

      Write(CHead,'(A18,I2)') 'Original MOs, sym.',iSym
      Call Cho_Head(CHead,'=',80,6)
      Do i = 1,nOrb
         Write(6,'(/,2X,A,I5)') 'Original MO no.',i
         Call Cho_Anasize(CMO(1,i),nShell,Work(ipBin),lBin,6)
      End Do

C     Local MOs.
C     ----------

      Write(XHead,'(A15,I2)') 'Local MOs, sym.',iSym
      Call Cho_Head(XHead,'=',80,6)
      Do i = 1,nOrb
         Write(6,'(/,2X,A,I5)') 'Local MO no.',i
         Call Cho_Anasize(XMO(1,i),nShell,Work(ipBin),lBin,6)
      End Do

C     De-allocate and return.
C     -----------------------

    1 Call GetMem('Bin','Free','Real',ipBin,lBin)

      End
