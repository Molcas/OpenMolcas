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

      SUBROUTINE set_nnA(nSym,nAorb,nnA)

      Implicit Real*8 (a-h,o-z)
      Integer nSym,nAorb(8)
      Integer nnA(8,8)

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
**************************************************

      do j=1,nSym
         do i=j,nSym

            kSym     = MulD2h(i,j)

            nnA(i,j) = nAorb(i)*nAorb(j)
     &               + Min(0,kSym-2)*nAorb(i)*(nAorb(i)-1)/2

            nnA(j,i) = nnA(i,j)

         end do
      end do

      Return
      END
