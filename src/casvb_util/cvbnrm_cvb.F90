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

subroutine cvbnrm_cvb(cvb)

use casvb_global, only: nfrag, nvb_fr
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cvb(nvb)
integer(kind=iwp) :: ifrag, nvbadd
real(kind=wp), external :: dnrm2_

if (nfrag <= 1) then
  call dscal_(nvb,One/dnrm2_(nvb,cvb,1),cvb,1)
else
  nvbadd = 1
  do ifrag=1,nfrag
    call dscal_(nvb_fr(ifrag),One/dnrm2_(nvb_fr(ifrag),cvb(nvbadd),1),cvb(nvbadd),1)
    nvbadd = nvbadd+nvb_fr(ifrag)
  end do
end if

return

end subroutine cvbnrm_cvb
