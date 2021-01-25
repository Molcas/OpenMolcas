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
      subroutine seth_cvb(iarr,n)
      implicit real*8(a-h,o-z)
#include "malloc_cvb.fh"
#include "seth_cvb.fh"
      dimension iarr(n)

      call wrbis_cvb(iarr,n,ncnt)
      return
      end
      subroutine geth_cvb(iarr,n)
      implicit real*8(a-h,o-z)
#include "malloc_cvb.fh"
#include "seth_cvb.fh"
      dimension iarr(n)
      if(icnt.lt.ncnt)then
        call rdbis_cvb(iarr,n,icnt)
      else
        call izero(iarr,n)
      endif
      return
      end
