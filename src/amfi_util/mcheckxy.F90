!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function mcheckxy(m1,m2,m3,m4)
!bs makes a check, if there is an interaction inbetween cartesian functions
!bs with m-values m1-m4

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mcheckxy
integer(kind=iwp), intent(in) :: m1, m2, m3, m4
integer(kind=iwp) :: int12a, int12b, int34a, int34b

mcheckxy = 1
int12a = m1+m2
int12b = -m1+m2
int34a = m3+m4
int34b = -m3+m4
!bs lots of checks
if (abs(int12a+int34a) == 1) return
if (abs(int12a-int34a) == 1) return
if (abs(int12b+int34b) == 1) return
if (abs(int12b-int34b) == 1) return
if (abs(int12a+int34b) == 1) return
if (abs(int12a-int34b) == 1) return
if (abs(int12b+int34a) == 1) return
if (abs(int12b-int34a) == 1) return
mcheckxy = 0

return

end function mcheckxy
