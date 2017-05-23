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
      subroutine applythmes_cvb(civbh,orbs,gjorb,gjorb2,gjorb3)
c  Apply T(O) (H - E) T(O) to CIVBH :
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "fx_cvb.fh"
      dimension civbh(ndet)
      dimension orbs(norb,norb),gjorb(*),gjorb2(*),gjorb3(*)

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

      call applyt_cvb(civbh,gjorb)
      call proj_cvb(civbh)
      call applyhpcx_cvb(civbh,-ww/ovraa)
      call applyt_cvb(civbh,gjorb2)
      return
      end
