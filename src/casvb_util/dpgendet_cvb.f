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
      subroutine dpgendet_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "WrkSpc.fh"

      ibase=mstacki_cvb(0)

      iapr_add=ll(20)
      ixapr_add=ll(21)
      ibpr_add=ll(22)
      ixbpr_add=ll(23)
      iconfs_add=ll(15)
      idetvb_add=ll(17)
      do 100 ifrag=1,nfrag
      nalf_l=nalf_fr(1,ifrag)
      nbet_l=nbet_fr(1,ifrag)
      call icomb_cvb(norb,nalf_l,nda_l)
      call icomb_cvb(norb,nbet_l,ndb_l)

      nda_fr(1,ifrag)=nda_l
      nda_fr(1,ifrag)=nda_l

      call vbgendet_cvb(
     >  iwork(iapr_add),iwork(ixapr_add),iwork(ibpr_add),
     >  iwork(ixbpr_add),
     >  iwork(iconfs_add),iwork(idetvb_add),
     >  nconf_fr(ifrag),nconfion_fr(0,ifrag),
     >  nda_l,ndb_l,ndetvb_fr(ifrag),nel_fr(ifrag),
     >  noe,nalf_l,nbet_l,norb)
      iapr_add=iapr_add+ndetvb_fr(ifrag)
      ibpr_add=ibpr_add+ndetvb_fr(ifrag)
      ixapr_add=ixapr_add+nda_fr(1,ifrag)+1
      ixbpr_add=ixbpr_add+ndb_fr(1,ifrag)+1
      idetvb_add=idetvb_add+ndetvb_fr(ifrag)
      iconfs_add=iconfs_add+noe*nconf_fr(ifrag)
100   continue

      do 200 ifrag=1,nfrag
      iastr_fr(ifrag)=mstacki_cvb(nalf_fr(1,ifrag)*nda_fr(1,ifrag))
      ibstr_fr(ifrag)=mstacki_cvb(nbet_fr(1,ifrag)*ndb_fr(1,ifrag))
      call stringen_cvb(nel_fr(ifrag),nalf_fr(1,ifrag),
     >  iwork(iastr_fr(ifrag)),iwork(ibstr_fr(ifrag)),nda_fr(1,ifrag))
200   continue

      call izero(iwork(ll(7)),naprodvb)
      call izero(iwork(ll(8)),nbprodvb)
      k1 = mstacki_cvb((norb+1)*(nalf+1))
      k2 = mstacki_cvb(nfrag)
      k3 = mstacki_cvb(nfrag)
      k4 = mstacki_cvb(nfrag)
      k5 = mstacki_cvb(nfrag+1)
      k6 = mstacki_cvb(norb*nfrag)
      mxstack=100
      k7 = mstacki_cvb(mxstack)
      call detsort2_cvb(
     >  iwork(k1),norb,nalf,nfrag,nda_fr(1,1),nalf_fr(1,1),
     >  iwork(k2),iwork(ll(7)),iwork(k3),iwork(k4),iwork(k5),
     >  iwork(1),iastr_fr,iwork(k6),iwork(k7),mxstack)
      call mfreei_cvb(k1)
      k1 = mstacki_cvb((norb+1)*(nbet+1))
      k2 = mstacki_cvb(nfrag)
      k3 = mstacki_cvb(nfrag)
      k4 = mstacki_cvb(nfrag)
      k5 = mstacki_cvb(nfrag+1)
      k6 = mstacki_cvb(norb*nfrag)
      mxstack=100
      k7 = mstacki_cvb(mxstack)
      call detsort2_cvb(
     >  iwork(k1),norb,nbet,nfrag,ndb_fr(1,1),nbet_fr(1,1),
     >  iwork(k2),iwork(ll(8)),iwork(k3),iwork(k4),iwork(k5),
     >  iwork(1),ibstr_fr,iwork(k6),iwork(k7),mxstack)
      call mfreei_cvb(k1)

      call mfreei_cvb(ibase)
      call setiaprtot_cvb()
      return
      end
