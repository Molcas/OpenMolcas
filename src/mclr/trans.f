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
      SubRoutine Trans(A,n,m,B)
*
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 A(n,m), B(m,n)
*
      Do i=1, m
         Do j=1, n
            B(i,j) = A(j,i)
         End Do
      End Do
*
      Return
      End
