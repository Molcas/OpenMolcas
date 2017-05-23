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
      Subroutine Log2Int(A,B,N)
      Logical A(N)
      Integer B(N)
*
      Do i = 1, N
         B(i)=0
         If (A(i)) B(i)=1
      End Do
*
      Return
      End
      Subroutine Int2Log(A,B,N)
      Integer A(N)
      Logical B(N)
*
      Do i = 1, N
         B(i)= A(i).eq.1
      End Do
*
      Return
      End
