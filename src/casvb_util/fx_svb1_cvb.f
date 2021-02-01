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
      subroutine fx_svb1_cvb(fx,fast,orbstry,cvbtry,
     >   civec,civecp,civbs,civb,
     >   gjorb,gjorb2,gjorb3,cvbdet)
      implicit real*8 (a-h,o-z)
      logical fast
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"
#include "fx_cvb.fh"
      dimension orbstry(norb,norb),cvbtry(nvb)
      dimension civec(ndet),civecp(ndet),civbs(ndet),civb(ndet)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvbdet(ndetvb)

      call str2vbc_cvb(cvbtry,cvbdet)
      if(fast)then
        call makecivb_cvb(civec,civb,cvbdet,orbstry,cvbtry,1)
        call gaussj_cvb(orbstry,gjorb)
        call applyt_cvb(civb,gjorb)
        call proj_cvb(civb)
        call cinorm_cvb(civb,ovraa_try)
        if(memplenty)then
          call cidot_cvb(civec,civb,ovrab_try)
        else
          call cird_cvb(civecp,61001.2d0)
          call cidot_cvb(civecp,civb,ovrab_try)
        endif
      else
        call makecivb_cvb(civec,civb,cvbdet,orbstry,cvbtry,0)
        call makecivecp_cvb(civec,civecp,orbstry)
        call makecivbs_cvb(civbs,orbstry,gjorb,gjorb2,gjorb3,cvbdet)

        call pvbdot_cvb(civb,civbs,ovraa_try)
        call pvbdot_cvb(civb,civecp,ovrab_try)
      endif
      svb=ovrab_try/sqrt(ovraa_try)
      fx=svb
      if(fast.and.ip(3).ge.2)
     >  write(6,formE)' Svb :      ',svb
      return
      end
