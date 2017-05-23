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
C
C----------------------------------------------------------------------|
C
      subroutine dmxma(n,transa,transb,a,b,c,alpha)
C
C Square real matrices multiplication
C
      implicit none
      integer n
      Real*8 a(*),b(*),c(*),alpha
      character*1 transa,transb
      call dgemm_(transa,transb,n,n,n,alpha,a(1),n,b(1),n,0.d0,c(1),n)
      return
      end
