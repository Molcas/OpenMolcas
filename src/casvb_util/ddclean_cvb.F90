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

subroutine ddclean_cvb()

use casvb_global, only: ap, axc, c, res, rhs, rhsp, solp, solp_res, sxc
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(c,safe='*')
call mma_deallocate(axc,safe='*')
call mma_deallocate(sxc,safe='*')
call mma_deallocate(res,safe='*')
call mma_deallocate(rhs,safe='*')
call mma_deallocate(ap,safe='*')
call mma_deallocate(rhsp,safe='*')
call mma_deallocate(solp,safe='*')
call mma_deallocate(solp_res,safe='*')

return

end subroutine ddclean_cvb
