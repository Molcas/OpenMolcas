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
      subroutine change3_cvb()
      implicit real*8 (a-h,o-z)
      logical changed
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "rls_cvb.fh"
#include "malloc_cvb.fh"

      changed=.false.
c Spin functions coefficients (BIKCOF) + inverse (AIKCOF)
c  Get KBASISCVB if we don't know it already (eqv. to
c  GETGUESS later):
c  Need to reserve enough space for both KBASIS and KBASISCVB
c  --> figure out which one needs most :
      if((kbasis.gt.2.and.kbasis.ne.6).or.
     >   (kbasiscvb.gt.2.and.kbasiscvb.ne.6))then
        kmost=3
      elseif(kbasis.le.2.or.kbasiscvb.le.2)then
        kmost=1
      else
        kmost=6
      endif
      if(chpcmp_cvb(kmost))changed=.true.
      if(changed)call touch_cvb('MEM3')
      return

      entry chop3_cvb()
      if(release(3))call mfreer_cvb(lb(1))
      release(3)=.true.
      release(4)=.false.

      call icomb_cvb(nel,nbet,iretval1)
      call icomb_cvb(nel,nbet-1,iretval2)
      mxfns=iretval1-iretval2
      if(kbasis.eq.5)call icomb_cvb(nel,nalf,mxfns)
      call icomb_cvb(nel,nalf,mxdetvb)
      if((kbasis.gt.2.and.kbasis.ne.6).or.
     >   (kbasiscvb.gt.2.and.kbasiscvb.ne.6))then
        kmost=3
      elseif(kbasis.le.2.or.kbasiscvb.le.2)then
        kmost=1
      else
        kmost=6
      endif
      call bspset_cvb(kmost,1,need)
      if(kmost.eq.3)then
        lb(1) = mstackr_cvb(1+need)
        lb(2) = mstackr_cvb(1+need)
      elseif(kmost.eq.1)then
        lb(1) = mstackr_cvb(1+need)
        lb(2) = lb(1)
      else
        lb(1) = mstackr_cvb(1)
        lb(2) = lb(1)
      endif
c  Flag AIKCOF/BIKCOF as unset :
      w(lb(1))=zero
      w(lb(2))=zero

      lb(3) = mstacki_cvb((nel+1)*(nel+1)*(nel+1))
      lb(4) = mstacki_cvb((nel+1)*(nel+1))
      lb(5) = mstacki_cvb((nel+1)*(nel+1))
      lb(6) = mstacki_cvb((nel+1)*(nel+1))
      call bspset_cvb(kbasiscvb,2,need)
      return
      end
