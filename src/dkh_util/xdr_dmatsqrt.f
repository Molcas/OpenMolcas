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
      subroutine XDR_dmatsqrt(a,n)
C
C Compute the inverse square root of a real symmetric matrix : A**(-1/2)
C
      implicit none
#include "WrkSpc.fh"
      Real*8 a(*),dia
      integer n,iCr,iW,iTmp,i,j,k,info
      call getmem('tmp ','ALLOC','REAL',iTmp,8*n)
      call getmem('Cr  ','ALLOC','REAL',iCr ,n*n+4)
      call getmem('Eig ','ALLOC','REAL',iW  ,n+4)
      do i = 0,n*n-1
        Work(iCr+i) = a(i+1)
      end do
      call dsyev_('V','L',n,Work(iCr),n,Work(iW),Work(iTmp),8*n,info)
      do i = 0,n-1
        dia = 1.d0/sqrt(sqrt(Work(iW+i)))
        do j = 0,n-1
          k = i*n + j
          Work(iCr+k) = Work(iCr+k) * dia
        end do
      end do
      call dmxma(n,'N','T',Work(iCr),Work(iCr),a,1.d0)
      call getmem('tmp ','FREE','REAL',iTmp,8*n)
      call getmem('Cr  ','FREE','REAL',iCr ,n*n+4)
      call getmem('Eig ','FREE','REAL',iW  ,n+4)
      return
      end
