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
      subroutine mxgendiag_cvb(a,s,eigval,n)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n,n),s(n,n),eigval(n)

      info=0
      lwrk=-1
      call dsygv_(1,'V','U',n,a,n,s,n,eigval,wrk,lwrk,info)
      lwrk=nint(wrk)
      i1 = mstackr_cvb(lwrk)
      call dsygv_(1,'V','U',n,a,n,s,n,eigval,w(i1),lwrk,info)
      call mfreer_cvb(i1)
      if(info.ne.0)then
        write(6,*)' Error in generalized diagonalization!'
        write(6,*)' Dsygv exited with code:',info
        call abend_cvb()
      endif
      return
      end
