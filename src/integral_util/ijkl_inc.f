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
      Subroutine ijkl_inc(i,j,k,l)
      Integer i,j,k,l
      l = l + 1
      If ( i.eq.k.and.l.le.j) Return
      If ( i.ne.k.and.l.le.k) Return
      l = 1
      k = k +1
      If (k.le.i) Return
      k = 1
      j = j + 1
      If (j.le.i) Return
      j=1
      i=i+1
      Return
      End
