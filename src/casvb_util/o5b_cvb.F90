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

subroutine o5b_cvb( &
#                  define _CALLING_
#                  include "optb_interface.fh"
                  )

use casvb_global, only: hh, maxize, odx, ograd, scalesmall
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "optb_interface.fh"
integer(kind=iwp) :: ipu
real(kind=wp), external :: dnrm2_

#include "macros.fh"
unused_var(grdnrm)

call fmove_cvb(ograd,odx,nparm)
if (.not. maxize) call dscal_(nparm,-One,odx,1)
dxnrm = dnrm2_(nparm,odx,1)
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
if ((dxnrm > hh) .or. scalesmall(ipu)) then
  call dscal_(nparm,hh/dxnrm,odx,1)
  dxnrm = hh
end if

return

end subroutine o5b_cvb
