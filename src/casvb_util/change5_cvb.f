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
      subroutine change5_cvb()
      implicit real*8 (a-h,o-z)
      logical changed,construc
c ... Change of dimensioning variables ...
      logical, external :: chpcmp_cvb,lchpcmp_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "rls_cvb.fh"

c  Dimensioning for symmetry handling :
      changed=.false.
      if(chpcmp_cvb(nsyme))changed=.true.
      if(chpcmp_cvb(ndimrel))changed=.true.
      if(chpcmp_cvb(norbrel))changed=.true.
      if(chpcmp_cvb(nvb))changed=.true.
      if(chpcmp_cvb(nzrvb))changed=.true.
      if(chpcmp_cvb(nort))changed=.true.
      if(chpcmp_cvb(ndrot))changed=.true.

      orbfr_is_unit=(ndimrel.eq.0.and.nfxorb.eq.0
     >  .and.nort.eq.0.and.(.not.plc_const))
c  Set ORBFR_IS_UNIT if optimization method is 'NONE' :
      if(imethod.eq.11)orbfr_is_unit=.true.
      if(lchpcmp_cvb(orbfr_is_unit))changed=.true.
      nfxvbr=nfxvb
      if(lfxvb.eq.1)nfxvbr=nvb-nfxvb
      nzrvbr=nzrvb
      if(lzrvb.eq.1)nzrvbr=nvb-nzrvb
c  NPR arrays -- depend on nature of opt procedure
      construc=((nzrvbr.gt.0).or.
     >          (nfxvbr.gt.0.and.nfxvbr.lt.nvb).or.
     >          (nzeta.gt.0))
      if(construc)then
        if(nvb.gt.20.or..not.strucopt)then
          iconstruc=1
        else
          iconstruc=2
        endif
      else
        iconstruc=0
      endif
      if(chpcmp_cvb(iconstruc))changed=.true.
      if(changed)call touch_cvb('MEM5')
      return

      entry chop5_cvb()
      if(release(5))call mfreer_cvb(ls(1))
      release(5)=.true.
      release(6)=.false.

      ls(1)= mstackr_cvb(nsyme*norb*norb)
      ls(2)= mstacki_cvb(ndimrel)
      ls(3)= mstacki_cvb(norb)
      ls(4)= mstackr_cvb(norb*norb*norb)
      ls(5)= mstacki_cvb(2*(norb-1))
      ls(6)= mstackr_cvb(min(norb-1,norbrel)*norb*norb)
      ls(8)= mstacki_cvb(norb)
      ls(9)= mstacki_cvb(nvb)
      ls(10)= mstacki_cvb(nzrvb)
      ls(11)= mstacki_cvb(2*nort)
      ls(12)= mstacki_cvb(2*ndrot)
      ls(13)= mstacki_cvb(nsyme)
      if(.not.orbfr_is_unit)then
        ls(14)= mstackr_cvb(nprorb*nprorb)
      else
        ls(14)= mstackr_cvb(0)
      endif
      if(iconstruc.eq.2)then
        ls(15)= mstackr_cvb(nvb*nvb)
      else
        ls(15)= mstackr_cvb(0)
      endif
      ls(16)= mstacki_cvb(norb*nzeta)
      return
      end
