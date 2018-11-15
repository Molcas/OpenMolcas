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
      Subroutine pa_sort (a, n)
      Implicit None
      Integer n, a(n), temp, k, l
      Do k = 1, n - 1
        Do l = k+1, n
         If (a(k).gt.a(l)) Then
         temp = a(k)
         a(k) = a(l)
         a(l) = temp
         End If
        End Do
      End Do
      Return
      End

