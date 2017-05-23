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
      subroutine cvprt_cvb(a,l)
      implicit real*8(a-h,o-z)
#include "formats_cvb.fh"
      character*20 a
      logical l
      character*16 a1
      save huge
      data huge/1d20/

      if(l)then
        write(6,'(2a)')a,'     Converged.'
      else
        write(6,'(2a)')a,' Not converged.'
      endif
      return
      entry cvprt2_cvb(a1,f1,f2,ic)
      if(abs(f2).ne.huge)then
        if(ic.eq.1.and.f1.lt.f2)then
          write(6,formcvp)a1,f1,'     smaller than',f2
        elseif(ic.eq.1)then
          write(6,formcvp)a1,f1,' not smaller than',f2
        elseif(ic.eq.2.and.f1.gt.f2)then
          write(6,formcvp)a1,f1,'     greater than',f2
        elseif(ic.eq.2)then
          write(6,formcvp)a1,f1,' not greater than',f2
        endif
      endif
      return
      end
