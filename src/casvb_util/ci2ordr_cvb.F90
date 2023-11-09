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

use casvb_global, only: ndet, ndetvb, nfrag, vbdet
use Constants, only: Zero
use Definitions, only: wp

#include "intent.fh"

implicit none
real(kind=wp), intent(_IN_) :: civec(0:ndet), cvbdet(ndetvb)
real(kind=wp), intent(_OUT_) :: evbdet(*)
real(kind=wp) :: dum

if (nfrag <= 1) then
  evbdet(1:ndetvb) = Zero
  return
end if
call dpci2vb2_cvb(civec(1:),cvbdet,vbdet,evbdet,0,dum,5)

return

end subroutine ci2ordr_cvb
