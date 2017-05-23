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
      subroutine change0_cvb()
#include "rls_cvb.fh"

      call touch_cvb('MEM0')
      release(1)=.false.
      return
      end
      function chpcmp_cvb(itst)
      implicit real*8 (a-h,o-z)
      logical chpcmp_cvb
#include "lstprm_cvb.fh"

      iprm=iprm+1
      if(iprm.gt.mxprm)then
        write(6,*)' Dimensioning error in CHPCMP!',iprm,mxprm
        call abend_cvb()
      endif
      chpcmp_cvb=(lstprm(iprm).ne.itst)
      lstprm(iprm)=itst
      return
      end
      function lchpcmp_cvb(ltst)
      implicit real*8 (a-h,o-z)
      logical lchpcmp_cvb,chpcmp_cvb
      logical ltst

      if(ltst)then
        itst=1
      else
        itst=0
      endif
      lchpcmp_cvb=chpcmp_cvb(itst)
      return
      end
