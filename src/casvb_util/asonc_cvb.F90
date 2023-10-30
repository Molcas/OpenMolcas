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

subroutine asonc_cvb( &
#                    define _CALLING_
#                    include "ddasonc_interface.fh"
                    )
! Applies H and S on c vector(s).

use casvb_global, only: civb1, civb2, cvbdet, orbs
use Definitions, only: wp, iwp

implicit none
#include "ddasonc_interface.fh"
integer(kind=iwp) :: ivec

do ivec=1,nvec
  call str2vbc_cvb(c(:,ivec),cvbdet)
  call vb2cif_cvb(cvbdet,civb2)
  call vb2cif_cvb(cvbdet,civb1)
  call makecivbhs_cvb(civb1,civb2,orbs)
  call ci2vbg_cvb(civb1,cvbdet)
  call vb2strg_cvb(cvbdet,axc(:,ivec))
  call ci2vbg_cvb(civb2,cvbdet)
  call vb2strg_cvb(cvbdet,sxc(:,ivec))
end do

return

end subroutine asonc_cvb
