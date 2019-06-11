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
      dimension wrk1_lapack(1)

*  Note that dgesvd returns V**T, not V.
c     write(6,*)' In full_svd. Calling dgesvd:'
      call dgesvd_('A','A',m,n,amat,m,svals,umat,m,vmat,n,
     &                                     wrk1_lapack,-1,info)
c     write(6,*)' full_svd back from dgesvd'
      lwork=int(wrk1_lapack(1))
c     write(6,*)' lwork:',lwork
      call getmem('lapckwrk','allo','real',ipwork,lwork)
c     write(6,*)' Calling dgesvd again:'
      call dgesvd_('A','A',m,n,amat,m,svals,umat,m,vmat,n,
     &                                 work(ipwork),lwork,info)
c     write(6,*)' full_svd back from dgesvd'
      call getmem('lapckwrk','free','real',ipwork,lwork)

      return
      end
