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
      subroutine ddrhs_cvb(vec,ndim,ioffs)
      implicit real*8(a-h,o-z)
#include "malloc_cvb.fh"
#include "direct_cvb.fh"
      dimension vec(ndim)

      nvrhs=nvrhs+1
      if(nvrhs.gt.mxrhs)then
        write(6,*)' Too many RHS vectors in Davidson!',nvrhs,mxrhs
        call abend_cvb()
      endif
      if(ndim+ioffs.gt.nparm)then
        write(6,*)' Illegal call to DDRHS :',ndim,ioffs,nparm
        call abend_cvb()
      endif
      call fzero(w(idd(ivrhs)+(nvrhs-1)*nparm),ioffs)
      call fmove_cvb(vec,w(ioffs+idd(ivrhs)+(nvrhs-1)*nparm),ndim)
      call fzero(w(ndim+ioffs+idd(ivrhs)+(nvrhs-1)*nparm),
     >  nparm-ioffs-ndim)
      return
      end
