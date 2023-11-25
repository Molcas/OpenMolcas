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

if (allocated(c)) call mma_deallocate(c)
if (allocated(axc)) call mma_deallocate(axc)
if (allocated(sxc)) call mma_deallocate(sxc)
if (allocated(res)) call mma_deallocate(res)
if (allocated(rhs)) call mma_deallocate(rhs)
if (allocated(ap)) call mma_deallocate(ap)
if (allocated(rhsp)) call mma_deallocate(rhsp)
if (allocated(solp)) call mma_deallocate(solp)
if (allocated(solp_res)) call mma_deallocate(solp_res)

return

end subroutine ddclean_cvb
