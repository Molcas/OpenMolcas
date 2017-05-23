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
      Subroutine Set_iOff(nSym,nA,nB,jStart,iOff)
      Implicit Real*8(a-h,o-z)
      Integer nSym, jStart, nA(nSym), nB(nSym), iOff(nSym)

      iOff(1)=jStart
      Do i=2,nSym
         iOff(i)=iOff(i-1)+nA(i-1)*nB(i-1)
      End Do

      Return
      End
