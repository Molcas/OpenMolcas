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
      subroutine defs_cvb()
      implicit real*8 (a-h,o-z)
      parameter (iunset=-1357924680,unset=-1357924680d0)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
      logical ifploc
      common /ifploc_complc/ifploc
      dimension ifxorb(mxorb)

c  Default settings :
      strtvb=zero
      savvb=zero
      savvbci=zero
      kbasis=1
      mxiter=iunset
      icrit=iunset
      imethod=iunset
      isaddle=iunset
      initial=-1
      opposite=.false.
      projcas=.false.
      projsym=.false.
      npcf=iunset
      ishstruc=iunset
c  +1=CHIRGWIN +2=LOWDIN +4=INVERSE
      ivbweights=iunset
      iciweights=iunset
      opposite=.false.
      projcas=.false.
      projsym=.false.
      sij=.false.
      anyslater=.false.
      service=.false.
      do 50 i=1,10
      ip(i)=1
50    continue

      call tunedefs_cvb()
      if(variat)then
        ploc=ifploc
      else
        ploc=.false.
      endif

c  Counters
      norbrel=0
      ndimrel=0
      nfxorb=0
      nfxvb=0
      nzrvb=0
      nort=0
      lfxvb=0
      lzrvb=0
      return

      entry defs2_cvb(ifxorb)
      if(icrit.eq.iunset)then
        if((.not.variat).and.imethod.ne.6)then
          icrit=1
        else
          icrit=2
        endif
      endif
c Set ifxorb
      nfxorb=0
      do i=1,norb
      if(ifxorb(i).eq.1)nfxorb=nfxorb+1
      enddo
      call izero(ifxorb(norb+1),mxorb-norb)

c Set STRUCOPT :
      if(projcas)then
        strucopt=.false.
      elseif(imethod.eq.11)then
        strucopt=.false.
      elseif(nvb.eq.1)then
        strucopt=.false.
      else
        nfxvbr=nfxvb
        if(lfxvb.eq.1)nfxvbr=nvb-nfxvb
        nzrvbr=nzrvb
        if(lzrvb.eq.1)nzrvbr=nvb-nzrvb
        if(nfxvbr.ge.nvb)then
          strucopt=.false.
        elseif(nfxvbr+nzrvbr.ge.nvb)then
          write(6,*)' Should check!'
          call abend_cvb()
        else
          strucopt=.true.
        endif
      endif
      if(isaddle.eq.iunset)isaddle=0
c  If unset, set default for IMETHOD :
      if(imethod.eq.iunset)then
        if(isaddle.eq.0)then
          imethod=10
        else
          imethod=7
        endif
        if(nfxorb.eq.norb.and.strucopt.and.nfrag.le.1)imethod=4
      elseif(isaddle.ne.0.and.imethod.eq.1)then
        call abend_cvb()
      endif
c  If unset, set default for MXITER :
      if(mxiter.eq.iunset)then
        if(imethod.ne.4)then
          mxiter=50
        else
          mxiter=200
        endif
      endif
      if(ishstruc.eq.iunset)ishstruc=0
      if(npcf.eq.iunset)npcf=-2
      if(ivbweights.eq.iunset)ivbweights=-1
      if(iciweights.eq.iunset)iciweights=0
      call tunedefs2_cvb(imethod,.false.)

      return
      end
