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
      subroutine mksymcvb2_cvb(cvb,tconstr,cvbdet)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension cvb(nvb)
      dimension tconstr(nvb*nvb)
      dimension cvbdet(ndetvb)
      save thresh
      data thresh/1.d-15/

c  Constraints on struc coeffs - either symmetry or deleted
      if(iconstruc.gt.0)then
        if(ip(1).ge.0)
     >  write(6,'(/,2a)')' Imposing constraints on ',
     >    'the structure coefficients.'
        call symtrizcvb_cvb(cvb)
        psnrm=ddot_(nvb,cvb,1,cvb,1)
        if(psnrm.lt.thresh)then
          write(6,*)' Fatal error - structure coefficients',
     >      ' null after symmetrization!'
          call abend_cvb()
        endif
        if(ip(1).ge.0)then
          write(6,'(/,a)')' Constrained structure coefficients :'
          write(6,'(a)')' ------------------------------------'
          call vecprint_cvb(cvb,nvb)
        endif
      endif
      call str2vbc_cvb(cvb,cvbdet)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_real_array(tconstr)
      end
