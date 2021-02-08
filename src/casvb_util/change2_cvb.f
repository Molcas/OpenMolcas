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
      subroutine change2_cvb()
      implicit real*8 (a-h,o-z)
      logical changed
c ... Change of dimensioning variables ...
      logical, external :: chpcmp_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "rls_cvb.fh"

      changed=.false.
c  VB wavefunction
c  (MXNVB should be upper bound on KBASIS & KBASISCVB ) :
      nvb_alloc=max(nvb_cvb(kbasiscvb),nvb_cvb(kbasis),mxnvb)
      if(chpcmp_cvb(norb))changed=.true.
      if(chpcmp_cvb(nvb_alloc))changed=.true.
      if(changed)call touch_cvb('MEM2')
      return

      entry chop2_cvb()
      if(release(2))call mfreer_cvb(lv(1))
      release(2)=.true.
      release(3)=.false.

c  Note zeroing of ORBS and CVB :
      lv(1) = mstackrz_cvb(norb*norb)
c  (MXNVB should be upper bound on KBASIS & KBASISCVB ) :
      nvb_alloc=max(nvb_cvb(kbasiscvb),nvb_cvb(kbasis),mxnvb)
      lv(2)= mstackrz_cvb(nvb_alloc)
      return
      end
