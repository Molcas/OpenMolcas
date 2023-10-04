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

use casvb_global, only: cvb, nfrag, nvb_fr, ndetvb_fr, vbdet
use Constants, only: One
use Definitions, only: wp

implicit none
real(kind=wp) :: fx
#include "main_cvb.fh"

if (nfrag <= 1) then
  call dscal_(nvb,-One,cvb,1)
  call dscal_(ndetvb,-One,vbdet,1)
else
  call dscal_(nvb_fr(1),-One,cvb,1)
  call dscal_(ndetvb_fr(1),-One,vbdet,1)
end if
call touch_cvb('CVB')
call fx_cvb(fx,.false.)

return

end subroutine chgsgn_cvb
