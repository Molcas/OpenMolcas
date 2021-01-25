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

#include "malloc_cvb.fh"

      ibase = mstacki_cvb(0)

      mxpair=mxorb*(mxorb+1)/2
      mxdimrel=mxpair*(3+mxops)
      iorbrel = mstacki_cvb(mxpair*(3+mxops))

      ifxorb = mstacki_cvb(mxorb)

      iorts = mstacki_cvb(2*mxpair)
      irots = mstacki_cvb(2*mxpair)

      izeta = mstacki_cvb(mxsyme)

      iorbs = mstackrz_cvb(mxaobf*mxorb)
      irdorbs = mstackiz_cvb(mxorb)

      call input2_cvb(
     >  iw(iorbrel),mxdimrel,iw(ifxorb),
     >  iw(iorts),iw(irots),iw(izeta),w(iorbs),iw(irdorbs))
      call mfreei_cvb(ibase)
      return
      end
