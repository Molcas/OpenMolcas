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

subroutine hello_cvb()

use casvb_global, only: variat
use Definitions, only: u6

implicit none

if (variat) write(u6,'(a)') ' '
write(u6,10)
if (.not. variat) call date1_cvb()

return
10 format(/,'     CASVB (Valence bond MCSCF)   Authors: T. Thorsteinsson and D. L. Cooper  (1996-2006)',/)

end subroutine hello_cvb
