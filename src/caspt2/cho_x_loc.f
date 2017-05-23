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
      Subroutine Cho_x_Loc(irc,Thrs,nSym,nBas,nFro,nIsh,
     &                                        nAsh,nSsh,CMO)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym), nAsh(nSym)
      Integer nIsh(nSym), nSsh(nSym)
      Real*8  Thrs, CMO(*)
#include "WrkSpc.fh"

      irc=0
      l_Dens = 0
      Do iSym = 1,nSym
         l_Dens = max(l_Dens,nBas(iSym)**2)
      End Do
      Call GetMem('Density','Allo','Real',ip_Dens,l_Dens)
      kOffC = 0
      Do iSym = 1,nSym
         If (nIsh(iSym) .gt. 0) Then
            kOff1 = 1 + kOffC + nBas(iSym)*nFro(iSym)
            Call GetDens_Localisation(Work(ip_Dens),CMO(kOff1),
     &                                nBas(iSym),nIsh(iSym))
            Call FZero(CMO(kOff1),nBas(iSym)*nIsh(iSym))
            Call ChoLoc(irc,Work(ip_Dens),CMO(kOff1),Thrs,yNrm,
     &                  nBas(iSym),nIsh(iSym))
            If (irc .ne. 0) Then
               Call GetMem('Density','Free','Real',ip_Dens,l_Dens)
               irc  = 1
               Return
            End If
         End If
         If (nSsh(iSym) .gt. 0) Then
            kOff1= 1+kOffC+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
            Call GetDens_Localisation(Work(ip_Dens),CMO(kOff1),
     &                                nBas(iSym),nSsh(iSym))
            Call FZero(CMO(kOff1),nBas(iSym)*nSsh(iSym))
            Call ChoLoc(irc,Work(ip_Dens),CMO(kOff1),Thrs,yNrm,
     &                  nBas(iSym),nSsh(iSym))
            If (irc .ne. 0) Then
               Call GetMem('Density','Free','Real',ip_Dens,l_Dens)
               irc  = 1
               Return
            End If
         End If
         kOffC = kOffC + nBas(iSym)**2
      End Do

      Call GetMem('Density','Free','Real',ip_Dens,l_Dens)
      Return
      End
