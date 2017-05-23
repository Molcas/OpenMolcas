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
      Subroutine Filler(N,M,A)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Dimension A(N,M)
      k=0
      Do i=1,N
        Do j=1,M
          k=k+1
          A(i,j)=1.000d0*j+0.100d0*i+0.001d0*k
        EndDo
      EndDo
      Return
      End
