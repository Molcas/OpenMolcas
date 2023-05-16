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
      subroutine input_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"

#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"

      ibase = mstacki_cvb(0)

      mxpair=mxorb_cvb*(mxorb_cvb+1)/2
      mxdimrel=mxpair*(3+mxops)
      iorbrel = mstacki_cvb(mxpair*(3+mxops))

      ifxorb = mstacki_cvb(mxorb_cvb)

      iorts = mstacki_cvb(2*mxpair)
      irots = mstacki_cvb(2*mxpair)

      izeta = mstacki_cvb(mxsyme)

      iorbs = mstackrz_cvb(mxaobf*mxorb_cvb)
      irdorbs = mstackiz_cvb(mxorb_cvb)

      call input2_cvb(
     >  iwork(iorbrel),mxdimrel,iwork(ifxorb),
     >  iwork(iorts),iwork(irots),iwork(izeta),work(iorbs),
     >  iwork(irdorbs))
      call mfreei_cvb(ibase)
      return
      end
