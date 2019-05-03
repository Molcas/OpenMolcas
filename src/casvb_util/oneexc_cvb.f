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
      subroutine oneexc_cvb(cfrom,cto,vij,diag,iPvb)
      implicit real*8 (a-h,o-z)
      logical diag
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension vij(*),cfrom(*),cto(*)

      idens=0
      icfrom=nint(cfrom(1))
      icto=nint(cto(1))

      if(iform_ci(icfrom).ne.0)then
        write(6,*)' Unsupported format in ONEEXC/ONEDENS :',
     >    iform_ci(icfrom)
        call abend_cvb()
      elseif(iform_ci(icto).ne.0)then
        write(6,*)' Unsupported format in ONEEXC/ONEDENS :',
     >    iform_ci(icto)
        call abend_cvb()
      endif

      call oneexc2_cvb(w(iaddr_ci(icfrom)),w(iaddr_ci(icto)),vij,
     >  iw(ll(1)),iw(ll(2)),iw(ll(5)),iw(ll(6)),w(ll(9)),w(ll(10)),
     >  iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >  nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3),diag,idens,
     >  iPvb)

c  If projcas and iPvb=0 we asume the normal density/1-ex. is required:
      if(projcas.and.iPvb.ne.0)then
        if(diag)then
          nvij=norb*norb
        else
          nvij=norb*(norb-1)
        endif
        ivij2=mstackr_cvb(nvij)
        if(idens.eq.0)then
          call fmove_cvb(vij,w(ivij2),nvij)
          call dscal_(nvij,-1d0,w(ivij2),1)
        else
          call fzero(w(ivij2),nvij)
        endif
        call oneexc2_cvb(w(iaddr_ci(icfrom)),w(iaddr_ci(icto)),w(ivij2),
     >    iw(ll(1)),iw(ll(2)),iw(ll(5)),iw(ll(6)),w(ll(9)),w(ll(10)),
     >    iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >    nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3),diag,idens,
     >    3-iPvb)
        if(idens.eq.1)call daxpy_(nvij,-1d0,w(ivij2),1,vij,1)
        call mfreer_cvb(ivij2)
      endif
      return
      end
      subroutine onedens_cvb(cfrom,cto,vij,diag,iPvb)
      implicit real*8 (a-h,o-z)
      logical diag
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension vij(*),cfrom(*),cto(*)
      idens=1
      icfrom=nint(cfrom(1))
      icto=nint(cto(1))

      if(iform_ci(icfrom).ne.0)then
        write(6,*)' Unsupported format in ONEEXC/ONEDENS :',
     >    iform_ci(icfrom)
        call abend_cvb()
      elseif(iform_ci(icto).ne.0)then
        write(6,*)' Unsupported format in ONEEXC/ONEDENS :',
     >    iform_ci(icto)
        call abend_cvb()
      endif

      call oneexc2_cvb(w(iaddr_ci(icfrom)),w(iaddr_ci(icto)),vij,
     >  iw(ll(1)),iw(ll(2)),iw(ll(5)),iw(ll(6)),w(ll(9)),w(ll(10)),
     >  iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >  nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3),diag,idens,
     >  iPvb)

c  If projcas and iPvb=0 we asume the normal density/1-ex. is required:
      if(projcas.and.iPvb.ne.0)then
        if(diag)then
          nvij=norb*norb
        else
          nvij=norb*(norb-1)
        endif
        ivij2=mstackr_cvb(nvij)
        if(idens.eq.0)then
          call fmove_cvb(vij,w(ivij2),nvij)
          call dscal_(nvij,-1d0,w(ivij2),1)
        else
          call fzero(w(ivij2),nvij)
        endif
        call oneexc2_cvb(w(iaddr_ci(icfrom)),w(iaddr_ci(icto)),w(ivij2),
     >    iw(ll(1)),iw(ll(2)),iw(ll(5)),iw(ll(6)),w(ll(9)),w(ll(10)),
     >    iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >    nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3),diag,idens,
     >    3-iPvb)
        if(idens.eq.1)call daxpy_(nvij,-1d0,w(ivij2),1,vij,1)
        call mfreer_cvb(ivij2)
      endif
      return
      end
