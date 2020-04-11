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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SubRoutine Wavelet_Transform(irc,ipCMO,nSym,nBas,nFro,nOrb2Loc,
     &                                 inv,Silent,xNrm)
C
C     Author: F. Aquilante
C
C     Purpose: wavelet transform of the MO basis (inv=0)
C              "       backtransform (inv=1)
C
      Implicit Real*8 (a-h,o-z)
      Integer irc, ipCMO, nSym, nBas(nSym), nFro(nSym), nOrb2Loc(nSym)
      Integer inv
      Logical Silent
      Real*8  xNrm

#include "WrkSpc.fh"

      Character*17 SecNam
      Parameter (SecNam = 'Wavelet_Transform')

      Integer  Log2
      External Log2

      real*8 ddot_
      external ddot_
*
*
      irc = 0
      xNrm = 0.0d0
      If (.not.Silent) Then
         If (inv.eq.0) Write(6,'(/,1X,A)')'Wavelet transform of the MOs'
         If (inv.eq.1) Write(6,'(/,1X,A)')'Inverse wavelet transform'//
     &                                    ' of the MOs'
         Write(6,'(1X,A,8(1X,I6))')
     &   'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
         Write(6,'(1X,A,8(1X,I6))')
     &   'Orbitals to transform:',(nOrb2Loc(iSym),iSym=1,nSym)
      End If

      If (inv.eq.1) go to 1000 ! Inverse wavelet transform

      njOrb = Log2(nOrb2Loc(1))
      l_Scr = nBas(1)*(2**njOrb-1)
      Do iSym = 2,nSym
         njOrb = Log2(nOrb2Loc(iSym))
         l_Scr = max(l_Scr,nBas(iSym)*(2**njOrb-1))
      End Do
      Call GetMem('Scratch','Allo','Real',ipScr,l_Scr)
      kOffC = ipCMO
      Do iSym = 1,nSym
         If (nOrb2Loc(iSym) .gt. 0) Then
            kOff1 = kOffC + nBas(iSym)*nFro(iSym)
            kOff2 = kOff1
            njOrb = Log2(nOrb2Loc(iSym))
            Do While (njOrb .ge. 1)
              Call FWT_Haar(nBas(iSym),njOrb,Work(ipScr),Work(kOff2))
              njOrb = 2**njOrb
              kOff2 = kOff2 + nBas(iSym)*njOrb
              njOrb = Log2(nOrb2Loc(iSym)-njOrb)
            End Do
            xNrm = xNrm + dDot_(nBas(iSym)*nOrb2Loc(iSym),Work(kOff1),1,
     &                                                   Work(kOff1),1)
            If (irc .ne. 0) Then
               irc  = 1
               xNrm = -9.9d9
               Return
            End If
         End If
         kOffC = kOffC + nBas(iSym)**2
      End Do
      xNrm = sqrt(xNrm)
      Call GetMem('Scratch','Free','Real',ipScr,l_Scr)
      Return
*
1000  Continue
      njOrb = Log2(nOrb2Loc(1))
      l_Scr = nBas(1)*2**njOrb
      Do iSym = 2,nSym
         njOrb = Log2(nOrb2Loc(iSym))
         l_Scr = max(l_Scr,nBas(iSym)*2**njOrb)
      End Do
      Call GetMem('Scratch','Allo','Real',iScr,l_Scr)
      kOffC = ipCMO
      Do iSym = 1,nSym
         If (nOrb2Loc(iSym) .gt. 0) Then
            kOff1 = kOffC + nBas(iSym)*nFro(iSym)
            kOff2 = kOff1
            njOrb = Log2(nOrb2Loc(iSym))
            Do While (njOrb .ge. 1)
              Call Inv_FWT_Haar(nBas(iSym),njOrb,Work(iScr),Work(kOff2))
              njOrb = 2**njOrb
              kOff2 = kOff2 + nBas(iSym)*njOrb
              njOrb = Log2(nOrb2Loc(iSym)-njOrb)
            End Do
            xNrm = xNrm + dDot_(nBas(iSym)*nOrb2Loc(iSym),Work(kOff1),1,
     &                                                   Work(kOff1),1)
            If (irc .ne. 0) Then
               irc  = 1
               xNrm = -9.9d9
               Return
            End If
         End If
         kOffC = kOffC + nBas(iSym)**2
      End Do
      xNrm = sqrt(xNrm)
      Call GetMem('Scratch','Free','Real',iScr,l_Scr)
*
      End
*                                                                      *
************************************************************************
*                                                                      *
      Integer Function Log2(n)

      Implicit none
      Integer n, m

      m=n
      Log2=0
      Do while (m .gt. 1)
         m=m/2
         Log2=Log2+1
      End Do

      Return
      End
