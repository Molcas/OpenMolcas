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
      SUBROUTINE Compute_Shanks(E1,E2,EOrb,lthEOrb,
     &                          nBas,nFro,nOcc,nSym,E0,Shanks1_E)

      Implicit Real*8 (a-h,o-z)
      Real*8 E1, E2, EOrb(lthEOrb), E0, Shanks1_E
      Integer nSym, nBas(nSym), nFro(nSym), nOcc(nSym)

      E0=0.0d0
      ioff=0
      Do iSym=1,nSym
         nOrb=nFro(iSym)+nOcc(iSym)
         Do iorb=1,nOrb
            jorb=ioff+iorb
            E0=E0+EOrb(jorb)
         End Do
         ioff=ioff+nBas(iSym)
      End Do
      E0=2.0d0*E0
*
      Call Peek_dScalar('PotNuc',PotNuc)
      E0 = E0 + PotNuc

*   Shanks formula

      Shanks1_E = (E2*E0 - E1**2) / (E2 - 2.0d0*E1 + E0)

      Return
      End
