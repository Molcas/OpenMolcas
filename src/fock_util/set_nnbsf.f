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

      SUBROUTINE set_nnBSF(nSym,Nbas,nnBSF,n2BSF)

      Implicit Real*8 (a-h,o-z)
      Integer nSym,nBas(8)
      Integer nnBSF(8,8),n2BSF(8,8)



      do j=1,nSym
         do i=j,nSym

            kSym = iEOR(i-1,j-1) + 1

            nnBSF(i,j) = nBas(i)*nBas(j)
     &                 + Min(0,kSym-2)*nBas(i)*(nBas(i)-1)/2

            nnBSF(j,i) = nnBSF(i,j)

            n2BSF(i,j) = nBas(i)*nBas(j)

            n2BSF(j,i) = n2BSF(i,j)

         end do
      end do

      Return
      END
