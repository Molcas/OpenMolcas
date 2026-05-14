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

subroutine chgsgn_cvb(fx)

use casvb_global, only: cvb, ndetvb_fr, nfrag, nvb, nvb_fr, vbdet
use Definitions, only: wp

implicit none
real(kind=wp), intent(out) :: fx

if (nfrag <= 1) then
  cvb(1:nvb) = -cvb(1:nvb)
  vbdet(:) = -vbdet(:)
else
  cvb(1:nvb_fr(1)) = -cvb(1:nvb_fr(1))
  vbdet(1:ndetvb_fr(1)) = -vbdet(1:ndetvb_fr(1))
end if
call touch_cvb('CVB')
call fx_cvb(fx,.false.)

return

end subroutine chgsgn_cvb
