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
      subroutine makecivb_cvb(civec,civb,cvbdet,orbs,cvb,ic)
c  Construct CIVB and CVBDET :
c  IC=0 : CIVB will contain full set of structures (if PROJCAS).
c  IC=1 : CIVB will contain only VB structures.
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)
      dimension civec(ndet),civb(ndet)
      dimension cvbdet(ndetvb)

      if(tstcnt_cvb(civb,3-ic))return

      if(.not.projcas)then
        call str2vbc_cvb(cvb,cvbdet)
        call vb2cic_cvb(cvbdet,civb)
      elseif(projcas)then
        iorbinv = mstackr_cvb(norb*norb)
        igjorb  = mstackr_cvb(norb*norb+ihlf_cvb(norb+2*norb*norb))

        if(memplenty)then
          call getci_cvb(civec)
          call cicopy_cvb(civec,civb)
        else
          call cird_cvb(civb,61001.2d0)
        endif
        call fmove(orbs,w(iorbinv),norb*norb)
        call mxinv_cvb(w(iorbinv),norb)
        call gaussj_cvb(w(iorbinv),w(igjorb))
        call applyt_cvb(civb,w(igjorb))
        call ci2vbc_cvb(civb,cvbdet)
        call vb2strc_cvb(cvbdet,cvb)
        if(ic.eq.1)call vb2cic_cvb(cvbdet,civb)

        call mfreer_cvb(iorbinv)
      endif
      call setcnt_cvb(civb,3-ic)
      return
      end
