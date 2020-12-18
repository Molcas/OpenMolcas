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
      subroutine getci_cvb(civec)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "casinfo_cvb.fh"
#include "io_cvb.fh"
      dimension ncix(mxirrep),civec(*)

      icivec=nint(civec(1))
      if(igetcnt2_cvb(icivec).eq.1)return
      if(.not.ifcasci_cvb())return
      call setcnt2_cvb(icivec,1)
      iwr=0

      if(iform_ci(icivec).ne.0)then
        write(6,*)' Unsupported format in GETCI :',iform_ci(icivec)
        call abend_cvb()
      endif

      if(iwr.eq.0)then
        if(ip(1).ge.1)then
          write(6,'(a)')' '
          call prtfid_cvb(' Restoring CI vector from ',strtci)
        endif
        call fzero(w(iaddr_ci(icivec)),ndet)
      elseif(iwr.eq.1)then
        if(ip(5).ge.1.and.valid_cvb(savvbci))then
          write(6,'(a)')' '
          call prtfid_cvb(' Saving VB CI vector to ',savvbci)
        endif
      endif

      do 1000 istsym_d=1,nstsym_d
      isyml=istsy_d(istsym_d)
      call getnci_cvb(ncix,istnel_d(istsym_d),istms2_d(istsym_d),
     >  istsy_d(istsym_d))
      nci=ncix(1)
      lcim=mstackr_cvb(nci)
      if(iwr.eq.0)then
        do 1100 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call mkfn_cvb(strtci,ibf)
          call rdcivec_cvb(w(lcim),filename(ibf),.true.)
          fac=sqrt(weight_d(istate,istsym_d))
          call mol2vbma_cvb(w(iaddr_ci(icivec)),w(lcim),isyml,fac)
        endif
1100    continue
      elseif(iwr.eq.1)then
        do 1200 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call vb2mol_cvb(w(iaddr_ci(icivec)),w(lcim),isyml)
          cnrm=one/dnrm2_(nci,w(lcim),1)
          call dscal_(nci,cnrm,w(lcim),1)
          call mkfn_cvb(savvbci,ibf)
          call wrcivec_cvb(w(lcim),filename(ibf),.not.variat)
        endif
1200    continue
      endif
      call mfreer_cvb(lcim)
1000  continue
      return
      end
      subroutine putci_cvb(civec)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "casinfo_cvb.fh"
#include "io_cvb.fh"
      dimension ncix(mxirrep),civec(*)

      icivec=nint(civec(1))
      iwr=1

      if(iform_ci(icivec).ne.0)then
        write(6,*)' Unsupported format in GETCI :',iform_ci(icivec)
        call abend_cvb()
      endif

      if(iwr.eq.0)then
        if(ip(1).ge.1)then
          write(6,'(a)')' '
          call prtfid_cvb(' Restoring CI vector from ',strtci)
        endif
        call fzero(w(iaddr_ci(icivec)),ndet)
      elseif(iwr.eq.1)then
        if(ip(5).ge.1.and.valid_cvb(savvbci))then
          write(6,'(a)')' '
          call prtfid_cvb(' Saving VB CI vector to ',savvbci)
        endif
      endif

      do 1000 istsym_d=1,nstsym_d
      isyml=istsy_d(istsym_d)
      call getnci_cvb(ncix,istnel_d(istsym_d),istms2_d(istsym_d),
     >  istsy_d(istsym_d))
      nci=ncix(1)
      lcim=mstackr_cvb(nci)
      if(iwr.eq.0)then
        do 1100 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call mkfn_cvb(strtci,ibf)
          call rdcivec_cvb(w(lcim),filename(ibf),.true.)
          fac=sqrt(weight_d(istate,istsym_d))
          call mol2vbma_cvb(w(iaddr_ci(icivec)),w(lcim),isyml,fac)
        endif
1100    continue
      elseif(iwr.eq.1)then
        do 1200 istate=1,nstats_d(istsym_d)
        if(abs(weight_d(istate,istsym_d)).gt.1.d-20)then
          call vb2mol_cvb(w(iaddr_ci(icivec)),w(lcim),isyml)
          cnrm=one/dnrm2_(nci,w(lcim),1)
          call dscal_(nci,cnrm,w(lcim),1)
          call mkfn_cvb(savvbci,ibf)
          call wrcivec_cvb(w(lcim),filename(ibf),.not.variat)
        endif
1200    continue
      endif
      call mfreer_cvb(lcim)
1000  continue
      return
      end
