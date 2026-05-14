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

subroutine o10a_cvb( &
#                   define _CALLING_
#                   include "opta_interface.fh"
                   )

use casvb_global, only: have_solved_it, n_div, nparm, nvguess, nvrestart, nvrhs, ograd
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "opta_interface.fh"
real(kind=wp) :: cnrm1, cnrm2
real(kind=wp), allocatable :: xp(:)
real(kind=wp), external :: dnrm2_

#include "macros.fh"
unused_var(nparam)

nvrestart = 0
nvguess = 0
nvrhs = 0
have_solved_it = .false.

call mma_allocate(xp,nparm,label='xp')
xp(:) = ograd(:)
call ddproj_cvb(xp,nparm)
cnrm1 = dnrm2_(n_div,xp(1:n_div),1)
cnrm2 = dnrm2_(nparm-n_div,xp(n_div+1:),1)
if (cnrm1 > cnrm2) then
  call ddguess_cvb(xp,n_div,0)
  if (cnrm2 > 1.0e-8_wp) call ddguess_cvb(xp(n_div+1:),nparm-n_div,n_div)
else
  call ddguess_cvb(xp(n_div+1:),nparm-n_div,n_div)
  if (cnrm1 > 1.0e-8_wp) call ddguess_cvb(xp(1:n_div),n_div,0)
end if
call ddrhs_cvb(xp,nparm,0)
call mma_deallocate(xp)

return

end subroutine o10a_cvb
