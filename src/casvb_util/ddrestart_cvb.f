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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine ddrestart_cvb(c,axc,vec,
     >  hp,solp,
     >  maxdav,n,
     >  nvguess1,nvrestart1)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension c(n,maxdav),axc(n,maxdav),vec(n)
      dimension hp(maxdav,maxdav)
      dimension solp(maxdav)

      i1 = mstackr_cvb(maxdav)
      i2 = mstackr_cvb(maxdav*maxdav)

      call ddrestart2_cvb(c,axc,vec,
     >  hp,solp,
     >  maxdav,n,
     >  nvguess1,nvrestart1,
     >  w(i1),w(i2))

      call mfreer_cvb(i1)
      return
      end
