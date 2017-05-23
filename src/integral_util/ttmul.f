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
      Subroutine TTMul(A,B,C,nRowA,nColA,nRowB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 A(nRowA,nColA), B(nRowB,nRowA),
     &       C(nColA,nRowB)
*
      nCache_=(64/8)*1024
      mCache=(nCache_*3)/4 - nRowA*nColA
      Incj=mCache/(nRowA+nColA)
*
*-----Sectioning of long index
*
      Do jj = 1, nRowB, Incj
         njVec=Min(Incj,nRowB-jj+1)
*
         Do i = 1, nColA
*-----------Set target to zero
            Do j = jj, jj+njVec-1
               C(i,j) = Zero
            End Do
            Do k = 1, nRowA
               If (A(k,i).ne.Zero) Then
                  Do j = jj, jj+njVec-1
                     C(i,j) = C(i,j) + A(k,i)*B(j,k)
                  End Do
               End If
            End Do
         End Do
*
      End Do    ! End of sectioning
*
      Return
      End
