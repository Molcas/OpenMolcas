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
      subroutine dpci2vb_cvb(civec,cvbdet,dvbdet,ic1,ret,ic)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "WrkSpc.fh"
      dimension civec(nda,ndb),cvbdet(ndetvb),dvbdet(ndetvb),dum(1)

      k1 = mstacki_cvb(nfrag)
      k2 = mstacki_cvb(nfrag)
      k3 = mstacki_cvb(nfrag+1)
      k4 = mstacki_cvb(nfrag+1)
      mxstack=100
      k5 = mstacki_cvb(mxstack)
      k6 = mstackr_cvb(nfrag+1)
      k7 = mstacki_cvb(nfrag)
      k8 = mstacki_cvb(nfrag)
      k9 = mstacki_cvb(nfrag)
      k10= mstacki_cvb(nfrag)
      call dpci2vb2_cvb(civec,cvbdet,dvbdet,dum,ic1,ret,ic,
     >  nda,ndb,ndetvb,
     >  nfrag,nda_fr(1,1),ndb_fr(1,1),
     >  iwork(ll(7)),iwork(ll(8)),iwork(k1),iwork(k2),iwork(k3),
     >  iwork(k4),
     >  iwork(k5),mxstack,work(k6),iwork(k7),
     >  iwork(ll(20)),iwork(ll(21)),
     >  iwork(k8),iwork(k9),iwork(k10),
     >  ndetvb_fr,ndavb)
      call mfreei_cvb(k1)
      return
      end
      subroutine ci2ordr_cvb(civec,cvbdet,evbdet)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "WrkSpc.fh"
      dimension civec(nda,ndb),cvbdet(ndetvb)
      dimension evbdet(*)
      icivec=nint(civec(1,1))
      if(nfrag.le.1)then
        call fzero(evbdet,ndetvb)
        return
      endif
      k1 = mstacki_cvb(nfrag)
      k2 = mstacki_cvb(nfrag)
      k3 = mstacki_cvb(nfrag+1)
      k4 = mstacki_cvb(nfrag+1)
      mxstack=100
      k5 = mstacki_cvb(mxstack)
      k6 = mstackr_cvb(nfrag+1)
      k7 = mstacki_cvb(nfrag)
      k8 = mstacki_cvb(nfrag)
      k9 = mstacki_cvb(nfrag)
      k10= mstacki_cvb(nfrag)
      call dpci2vb2_cvb(work(iaddr_ci(icivec)),cvbdet,work(lv(5)),
     >  evbdet,
     >  0,dum,5,
     >  nda,ndb,ndetvb,
     >  nfrag,nda_fr(1,1),ndb_fr(1,1),
     >  iwork(ll(7)),iwork(ll(8)),iwork(k1),iwork(k2),iwork(k3),
     >  iwork(k4),
     >  iwork(k5),mxstack,work(k6),iwork(k7),
     >  iwork(ll(20)),iwork(ll(21)),
     >  iwork(k8),iwork(k9),iwork(k10),
     >  ndetvb_fr,ndavb)
      call mfreei_cvb(k1)
      return
      end
