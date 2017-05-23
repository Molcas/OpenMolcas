************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine mxdiag_cvb(a,eigval,n)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n,n),eigval(n)

      itmp = mstackr_cvb(n*3)
cvv MKL blas dsyev doesn't work properly. use this simple workaround!!
      call dsyev_vv('V','L',n,a,n,eigval,w(itmp),n*3,ierr)
      call mfreer_cvb(itmp)
      if(ierr.ne.0)then
        write(6,*)' Fatal error in mxdiag, ierr :',ierr
        call abend_cvb()
      endif
      return
      end
