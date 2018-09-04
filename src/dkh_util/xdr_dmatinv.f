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
      subroutine XDR_dmatinv(a,n)
C
C Invert a real square matrix
C
      implicit none
#include "WrkSpc.fh"
      Real*8 a(*)
      integer n,ipiv,iTMp,info1,info2
#ifdef _MOLCAS_MPP_
      logical :: isCloneQ
      call check_parallel_data(a,n*n,isCloneQ,"C")
#endif
      call getmem('ipiv','ALLOC','INTE',ipiv,n+4)
      call getmem('tmp ','ALLOC','REAL',iTmp,n+4)
      call dgetrf_(n,n,a(1),n,iWork(ipiv),info1)
      call dgetri_(n  ,a(1),n,iWork(ipiv),Work(iTmp),n,info2)
      call getmem('ipiv','FREE','INTE',ipiv,n+4)
      call getmem('tmp ','FREE','REAL',iTmp,n+4)
#ifdef _MOLCAS_MPP_
      if(isCloneQ)then
         call check_parallel_data(a,n*n,isCloneQ,"S")
      end if
#endif
      return
      end
