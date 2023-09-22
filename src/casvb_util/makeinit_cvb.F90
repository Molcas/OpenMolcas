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
!*************************************
!** Routines to emulate unix "make" **
!*************************************

!IFG trivial
subroutine makeinit_cvb()

use casvb_global, only: ioffs, iprint, joffs, mustdeclare, ndep_ij, ndep_ji, nobj

implicit none

nobj = 0
ndep_ij = 0
ndep_ji = 0
ioffs(1) = 0
joffs(1) = 0
mustdeclare = .false.
iprint = 0

return

end subroutine makeinit_cvb
