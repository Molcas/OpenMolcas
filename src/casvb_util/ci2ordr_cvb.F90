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
integer(kind=iwp) :: icivec, mxstack
real(kind=wp) :: dum

icivec = nint(civec(1,1))
if (nfrag <= 1) then
  call fzero(evbdet,ndetvb)
  return
end if
mxstack = 100
call dpci2vb2_cvb(work(iaddr_ci(icivec)),cvbdet,work(lv(5)),evbdet,0,dum,5,nda,ndb,ndetvb,nfrag,nda_fr(1,1),ndb_fr(1,1), &
                  iwork(ll(7)),iwork(ll(8)),mxstack,iwork(ll(20)),iwork(ll(21)),ndetvb_fr,ndavb)

return

end subroutine ci2ordr_cvb
