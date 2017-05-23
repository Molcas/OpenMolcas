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
      subroutine mkcifree2_cvb(cvb,ifxstr,tconstr)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension cvb(nvb),ifxstr(nfxvb),tconstr(nvb,nvb)

      if(strucopt)then
        nfrvb=nvb
      else
        nfrvb=0
      endif
      nfr=nfrorb+nfrvb
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(cvb)
        call Unused_integer_array(ifxstr)
        call Unused_real_array(tconstr)
      end if
      end
