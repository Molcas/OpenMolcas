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
      subroutine cimol2vb_cvb(vec,civec)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "casinfo_cvb.fh"
      dimension vec(*)
      dimension ncix(mxirrep),civec(*)

      iwr=0
      icivec=nint(civec(1))

      if(iwr.eq.0)call fzero(w(iaddr_ci(icivec)),ndet)

      icioffs=0
      do 1000 istsym_d=1,nstsym_d
      isyml=istsy_d(istsym_d)
      if(isymv(isyml).ne.1)goto 1000
      call getnci_cvb(ncix,istnel_d(istsym_d),istms2_d(istsym_d),
     >  istsy_d(istsym_d))
      nci=ncix(1)
      lcim=mstackr_cvb(nci)
      if(iwr.eq.0)then
        do 1100 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call fmove(vec(1+icioffs),w(lcim),nci)
          icioffs=icioffs+nci
          fac=sqrt(weight_d(istate,istsym_d))
          call mol2vbma_cvb(w(iaddr_ci(icivec)),w(lcim),isyml,fac)
        endif
1100    continue
      elseif(iwr.eq.1)then
        do 1200 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call vb2mol_cvb(w(iaddr_ci(icivec)),w(lcim),isyml)
          call fmove(w(lcim),vec(1+icioffs),nci)
          icioffs=icioffs+nci
        endif
1200    continue
      endif
      call mfreer_cvb(lcim)
1000  continue
      return
      end
      subroutine civb2mol_cvb(civec,vec)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "casinfo_cvb.fh"
      dimension vec(*),civec(*)
      dimension ncix(mxirrep)
      iwr=1
      icivec=nint(civec(1))

      if(iwr.eq.0)call fzero(w(iaddr_ci(icivec)),ndet)

      icioffs=0
      do 1000 istsym_d=1,nstsym_d
      isyml=istsy_d(istsym_d)
      if(isymv(isyml).ne.1)goto 1000
      call getnci_cvb(ncix,istnel_d(istsym_d),istms2_d(istsym_d),
     >  istsy_d(istsym_d))
      nci=ncix(1)
      lcim=mstackr_cvb(nci)
      if(iwr.eq.0)then
        do 1100 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call fmove(vec(1+icioffs),w(lcim),nci)
          icioffs=icioffs+nci
          fac=sqrt(weight_d(istate,istsym_d))
          call mol2vbma_cvb(w(iaddr_ci(icivec)),w(lcim),isyml,fac)
        endif
1100    continue
      elseif(iwr.eq.1)then
        do 1200 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call vb2mol_cvb(w(iaddr_ci(icivec)),w(lcim),isyml)
          call fmove(w(lcim),vec(1+icioffs),nci)
          icioffs=icioffs+nci
        endif
1200    continue
      endif
      call mfreer_cvb(lcim)
1000  continue
      return
      end
