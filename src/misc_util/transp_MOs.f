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
      SubRoutine Transp_MOs(CMO1,CMO2,nSym,nFro,nIsh,nAsh,
     &                                     nSsh,nBas)
C
C        CMO1(alpha,p) -> CMO2(p,alpha)
C
C        At the same time, exclude Fro and Del orbitals
C
#include "implicit.fh"
      Integer nSym, nFro(nSym), nIsh(nSym), nAsh(nSym)
      Integer nSsh(nSym), nBas(nSym)
      Real*8  CMO1(*), CMO2(*)

      iCount = 0
      lCount = 0
      Do iSym = 1,nSym

         jCount = iCount + nBas(iSym)*nFro(iSym)

         nOrb=nIsh(iSym)+nAsh(iSym)+nSsh(iSym)

         Do i = 1,nOrb
            kOff1 = jCount + nBas(iSym)*(i-1) + 1
            kOff2 = lCount + i
            Call dCopy_(nBas(iSym),CMO1(kOff1),1,CMO2(kOff2),nOrb)
         End Do

         lCount = lCount + nOrb*nBas(iSym)
         iCount = iCount + nBas(iSym)**2

      End Do

      End
