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
      subroutine mktrnspn2_cvb(cvb,cvbdet)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "spinb_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb)

      if(ip(1).ge.1)
     >  write(6,'(/,4a)')' Changing spin basis : ',
     >  spinb(kbasiscvb)(1:len_trim_cvb(spinb(kbasiscvb))),' --> ',
     >  spinb(kbasis)(1:len_trim_cvb(spinb(kbasis)))
      call str2vbc_cvb(cvb,cvbdet)
      kbasiscvb=kbasis
      nvb=nvb_cvb(kbasiscvb)
      call vb2strc_cvb(cvbdet,cvb)
      return
      end
