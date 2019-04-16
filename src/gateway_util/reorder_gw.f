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
       Subroutine Reorder_GW(A,B,k,l,n,m)
       Implicit Real*8 (a-h,o-z)
       Real*8 A(k,l,n,m), B(k,n,l,m)
*
       Do ik = 1, k
          Do il = 1, l
             Do in = 1, n
                Do im = 1, m
*
                    B(ik,in,il,im) = A(ik,il,in,im)
*
                End Do
             End Do
          End Do
       End Do
*
       Return
       End
