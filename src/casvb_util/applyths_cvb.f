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
      subroutine applyths_cvb(civbh,civbs,orbs,gjorb,gjorb2,gjorb3)
c  Apply T(O) H T(O) to CIVBH & T(s) to CIVBS :
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension civbh(ndet),civbs(ndet)
      dimension orbs(norb,norb),gjorb(*),gjorb2(*),gjorb3(*)

      call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

      if(.not.proj)then
        call cicopy_cvb(civbs,civbh)
        call applyth_cvb(civbh,orbs,gjorb,gjorb2,gjorb3)
        call applyt_cvb(civbs,gjorb3)
      else
        call applyt_cvb(civbs,gjorb)
        call proj_cvb(civbs)
        call cicopy_cvb(civbs,civbh)
        call applyh_cvb(civbh)
        call applyt_cvb(civbs,gjorb2)
        call applyt_cvb(civbh,gjorb2)
      endif
      return
      end
