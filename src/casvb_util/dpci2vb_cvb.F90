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

subroutine dpci2vb_cvb(civec,cvbdet,dvbdet,ic1,ret,ic)

use casvb_global, only: nda, ndb, ndetvb
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: civec(nda,ndb), cvbdet(ndetvb)
real(kind=wp), intent(in) :: dvbdet(ndetvb)
integer(kind=iwp), intent(in) :: ic1, ic
real(kind=wp), intent(_OUT_) :: ret
real(kind=wp) :: dum(1)

call dpci2vb2_cvb(civec,cvbdet,dvbdet,dum,ic1,ret,ic)

return

end subroutine dpci2vb_cvb
