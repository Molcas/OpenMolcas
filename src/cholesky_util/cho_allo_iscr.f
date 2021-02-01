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
      SubRoutine Cho_Allo_iScr(DoDummy)
C
C     Purpose: allocate iScr array for reading and reordering vectors.
C              If (DoDummy): make dummy (length 1) allocation.
C
      use ChoArr, only: iScr
      Implicit None
      Logical DoDummy
#include "cholesky.fh"
#include "stdalloc.fh"

      Integer iSym, l_iScr

      If (DoDummy) Then
         l_iScr = 1
      Else
         l_iScr = nnBstR(1,1)
         Do iSym = 2,nSym
            l_iScr = max(l_iScr,nnBstR(iSym,1))
         End Do
      End If
      Call mma_allocate(iScr,l_iScr,Label='iScr')

      End
