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
      subroutine asonc2_cvb(c,axc,sxc,nvec,
     >   civbh,civbs,orbs,gjorb,gjorb2,gjorb3,cvbdet)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension c(nvb,nvec),axc(nvb,nvec),sxc(nvb,nvec)
      dimension civbh(ndet),civbs(ndet)
      dimension orbs(norb,norb),gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvbdet(ndetvb)

      do 100 ivec=1,nvec
      call str2vbf_cvb(c(1,ivec),cvbdet)
      call vb2cif_cvb(cvbdet,civbs)
      call vb2cif_cvb(cvbdet,civbh)
      call makecivbhs_cvb(civbh,civbs,orbs,gjorb,gjorb2,gjorb3)
      call ci2vbg_cvb(civbh,cvbdet)
      call vb2strg_cvb(cvbdet,axc(1,ivec))
      call ci2vbg_cvb(civbs,cvbdet)
100   call vb2strg_cvb(cvbdet,sxc(1,ivec))
      return
      end
