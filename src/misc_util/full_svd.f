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
      subroutine full_svd(m,n,amat,umat,vmat,svals)
      implicit real*8 (a-h,o-z)
      dimension amat(m,*)
      dimension umat(m,*)
      dimension vmat(n,*)
      dimension svals(*)
#include "WrkSpc.fh"
      dimension wrk1_lapack(1000)

*  Note that dgesvd returns V**T, not V.
*      nm=min(n,m)
      write(6,*)' In full_svd. Calling dgesvd:'
      call xflush(6)
      call dgesvd_('A','A',m,n,amat,m,svals,umat,m,vmat,n,
     &                                     wrk1_lapack,-1,info)
      write(6,*)' full_svd back from dgesvd'
      call xflush(6)
      lwork=int(wrk1_lapack(1))
      write(6,*)' lwork:',lwork
      call xflush(6)
      call getmem('lapckwrk','allo','real',ipwork,lwork)
      write(6,*)' Calling dgesvd again:'
      call xflush(6)
      call dgesvd_('A','A',m,n,amat,m,svals,umat,m,vmat,n,
     &                                 work(ipwork),lwork,info)
      write(6,*)' full_svd back from dgesvd'
      call xflush(6)
      call getmem('lapckwrk','free','real',ipwork,lwork)

      return
      end
