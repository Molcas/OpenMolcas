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
!****************************
!** Time and date routines **
!****************************

!IFG trivial
function tim0_cvb()

use Definitions, only: wp

implicit none
real(kind=wp) :: tim0_cvb
real(kind=wp) :: cpu, cpusince, wall, wallsince

call timing(cpu,cpusince,wall,wallsince)
tim0_cvb = cpu

return

end function tim0_cvb
