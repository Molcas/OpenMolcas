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
      subroutine makecivecp_cvb(civec,civecp,orbs)
c  Construct CIVECP :
      implicit real*8 (a-h,o-z)
c ... Content of CI vectors ...
      logical, external :: tstcnt_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb)
      dimension civec(ndet),civecp(ndet)

      if(tstcnt_cvb(civecp,3))return

      iowrk  = mstackr_cvb(norb*norb)
      igjorb = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))
      call transp_cvb(orbs,w(iowrk),norb,norb)
      call gaussj_cvb(w(iowrk),w(igjorb))
      if(memplenty)then
        call getci_cvb(civec)
        call cicopy_cvb(civec,civecp)
      else
        call cird_cvb(civecp,61001.2d0)
      endif
      call applyt_cvb(civecp,w(igjorb))
      call mfreer_cvb(iowrk)

      call setcnt_cvb(civecp,3)
      return
      end
