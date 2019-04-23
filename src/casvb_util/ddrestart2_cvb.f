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
      subroutine ddrestart2_cvb(c,axc,vec,
     >  hp,solp,
     >  maxdav,n,
     >  nvguess1,nvrestart1,
     >  eigval,eigvec)
      implicit real*8 (a-h,o-z)
#include "direct_cvb.fh"
      dimension c(n,maxdav),axc(n,maxdav),vec(n)
      dimension hp(maxdav,maxdav)
      dimension solp(maxdav)
      dimension eigval(maxdav),eigvec(maxdav,maxdav)

      call fmove_cvb(hp,eigvec,maxdav*maxdav)
      call mxdiag_cvb(eigvec,eigval,maxdav)

      nvrestart1=0
      nvguess1=0
      call mxatb_cvb(c,solp,n,maxdav,1,vec)
      if(ifollow.le.2)then
c  (Put lower-lying solutions in AxC :)
        do 200 ir=1,nroot-1
        if(ifollow.eq.1)then
          ir_use=maxdav-ir+1
        else
          ir_use=ir
        endif
200     call mxatb_cvb(c,eigvec(1,ir_use),n,maxdav,1,axc(1,ir+1))
        call fmove_cvb(axc(1,2),c(1,2),n*(nroot-1))
      endif
      call fmove_cvb(vec,c(1,1),n)
      return
      end
