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
      subroutine change4_cvb()
      implicit real*8 (a-h,o-z)
      logical changed
      logical ndres_ok
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "casinfo_cvb.fh"
#include "rls_cvb.fh"
#include "malloc_cvb.fh"
      save ndres_ok

      changed=.false.
c CI vectors
      call icomb_cvb(norb,nalf,nda)
      call icomb_cvb(norb,nbet,ndb)
      ndet = nda*ndb
      ndres=3+ndet

      memplenty=(mavailr_cvb().gt.9*ndet)

      if(chpcmp_cvb(ndres))changed=.true.
      ndres_ok=(.not.changed).and.(up2date_cvb('MEM3'))

      if(changed)call touch_cvb('RDCAS')
      if(chpcmp_cvb(nint(strtcas*1d1)))call touch_cvb('RDCAS')
      call chpcmp2_cvb(icrit,icritold)
      call chpcmp2_cvb(ifinish,ifin_old)
      if(.not.(icritold.eq.1.and.ifin_old.eq.0))call touch_cvb('RDCAS')

      if(ifinish.eq.1.or.ifinish.eq.2)then
c  Logical variables for wavefunction analysis :
c  Always evaluate svb/evb when possible
c    (a) to get exact values (wfn is updated by optim in last iteration)
c    (b) to generate ESYM (variational calculations)
        lcalcsvb=ifcasci_cvb()
        lcalcevb=ifhamil_cvb()
        lcalccivbs=.true.
        lciweights=(npcf.gt.-2.and.((.not.variat).or.endvar)
     >    .and.iciweights.gt.0)
      endif
      if(imethod.eq.11.and.ifinish.eq.0)then
        nv=2
        icase=2
      elseif(imethod.ne.4.and.ifinish.eq.0)then
        if(memplenty)then
          nv=8
          icase=1
        else
          if(icrit.eq.2.and.imethod.ne.6)then
c  No need for CITMP :
            nv=3
            icase=5
          else
            nv=5
            icase=2
          endif
        endif
      elseif(imethod.eq.4.and.ifinish.eq.0)then
        nv=2
        icase=2
      elseif(imethod.eq.12.and.ifinish.eq.0)then
        nv=8
        icase=7
      elseif(ifinish.eq.1.or.ifinish.eq.2)then
        nv=5
        if(.not.lciweights)then
          nv=nv-1
          if((.not.lcalcevb).or.lcalccivbs)nv=nv-1
          if((.not.lcalcsvb).or.lcalccivbs)nv=nv-1
        endif
        icase=3
      else
        nv=3
        icase=4
      endif
      if(chpcmp_cvb(nv))changed=.true.
      if(chpcmp_cvb(icase))changed=.true.
      if(changed)call touch_cvb('MEM4')
      return

      entry chop4_cvb()
      if(release(4))call mfreer_cvb(lc(1)-2)
      release(4)=.true.
      release(5)=.false.


c  CIVECP and CIVBH share memory --> LC(2)
      do iv=1,nv
      lc(iv) = mstackr_cvb(ndres)
      enddo
        do iv=1,nv
        w(lc(iv))=zero
        w(lc(iv)+1)=zero
        enddo
      do iv=1,nv
      lc(iv)=lc(iv)+2
      enddo
c Fix to put in "objects" :
      do iv=1,nv
      call creatci_cvb(iv,w(lc(iv)),lc(iv)+1,nint(w(lc(iv)-2)),
     >  w(lc(iv)-1))
      if(.not.ndres_ok)call setcnt_cvb(w(lc(iv)),0)
      enddo
c-- ins
      if((ifinish.eq.1.or.ifinish.eq.2).and..not.lciweights)then
        if(((.not.lcalcevb).or.lcalccivbs).and.
     >     ((.not.lcalcsvb).or.lcalccivbs))then
          lc(3)=lc(1)
          lc(4)=lc(1)
        elseif((.not.lcalcevb).or.lcalccivbs)then
          lc(4)=lc(3)
        elseif((.not.lcalcsvb).or.lcalccivbs)then
          lc(4)=lc(3)
          lc(3)=lc(1)
        endif
      endif
      if(.not.memplenty)lc(nv+1)=lc(1)
      if(nv.eq.3.or.nv.eq.4.or.nv.eq.5)then
        lc(6)=lc(2)
        lc(7)=lc(3)
        lc(8)=lc(4)
      endif
      return
      end
