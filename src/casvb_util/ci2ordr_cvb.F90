!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine ci2ordr_cvb(civec,cvbdet,evbdet)

use casvb_global, only: nda_fr, ndb_fr, ndetvb_fr, nfrag
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: civec(nda,ndb), cvbdet(ndetvb), evbdet(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: icivec, k1, k10, k2, k3, k4, k5, k6, k7, k8, k9, mxstack
real(kind=wp) :: dum
integer(kind=iwp), external :: mstacki_cvb, mstackr_cvb

icivec = nint(civec(1,1))
if (nfrag <= 1) then
  call fzero(evbdet,ndetvb)
  return
end if
k1 = mstacki_cvb(nfrag)
k2 = mstacki_cvb(nfrag)
k3 = mstacki_cvb(nfrag+1)
k4 = mstacki_cvb(nfrag+1)
mxstack = 100
k5 = mstacki_cvb(mxstack)
k6 = mstackr_cvb(nfrag+1)
k7 = mstacki_cvb(nfrag)
k8 = mstacki_cvb(nfrag)
k9 = mstacki_cvb(nfrag)
k10 = mstacki_cvb(nfrag)
call dpci2vb2_cvb(work(iaddr_ci(icivec)),cvbdet,work(lv(5)),evbdet,0,dum,5,nda,ndb,ndetvb,nfrag,nda_fr(1,1),ndb_fr(1,1), &
                  iwork(ll(7)),iwork(ll(8)),iwork(k1),iwork(k2),iwork(k3),iwork(k4),iwork(k5),mxstack,work(k6),iwork(k7), &
                  iwork(ll(20)),iwork(ll(21)),iwork(k8),iwork(k9),iwork(k10),ndetvb_fr,ndavb)
call mfreei_cvb(k1)

return

end subroutine ci2ordr_cvb
