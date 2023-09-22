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

subroutine o10a_cvb(nparm1)

use casvb_global, only: have_solved_it, ix, n_div, nparm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nparm1
#include "WrkSpc.fh"
integer(kind=iwp) :: ixp
real(kind=wp) :: cnrm1, cnrm2
integer(kind=iwp), external :: mstackr_cvb
real(kind=wp), external :: dnrm2_

call ddnewopt_cvb()
have_solved_it = .false.

ixp = mstackr_cvb(nparm)
call fmove_cvb(work(ix(2)),work(ixp),nparm)
call ddproj_cvb(work(ixp),nparm)
cnrm1 = dnrm2_(n_div,work(ixp),1)
cnrm2 = dnrm2_(nparm-n_div,work(n_div+ixp),1)
if (cnrm1 > cnrm2) then
  call ddguess_cvb(work(ixp),n_div,0)
  if (cnrm2 > 1.0e-8_wp) call ddguess_cvb(work(n_div+ixp),nparm-n_div,n_div)
else
  call ddguess_cvb(work(n_div+ixp),nparm-n_div,n_div)
  if (cnrm1 > 1.0e-8_wp) call ddguess_cvb(work(ixp),n_div,0)
end if
call ddrhs_cvb(work(ixp),nparm,0)
call mfreer_cvb(ixp)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nparm1)

end subroutine o10a_cvb
