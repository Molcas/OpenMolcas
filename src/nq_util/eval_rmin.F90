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

function Eval_RMin(Alpha,m,R_H)

use Constants, only: One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Eval_RMin
real(kind=wp), intent(in) :: Alpha, R_H
integer(kind=iwp), intent(in) :: m
real(kind=wp) :: D_m, ln_x, R_Min

!write(u6,*) 'Alpha,m,R_H=',Alpha,m,R_H
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute r_1 as a function of m and R_H
!
! Eq(25) R. Lindh, P.-A. Malmqvist, L. Gagliardi,
! TCA, 106:178-187 (2001)

D_m = -Four
if (m == 4) D_m = -2.3_wp
if (m == 2) D_m = -One
if (m == 0) D_m = 1.9_wp
if (m == -2) D_m = 9.1_wp

!x = alpha*(r_1)**2

ln_x = (Two/(real(m,kind=wp)+Three))*(D_m-log(One/R_H))

R_Min = sqrt(exp(ln_x)/Alpha)

Eval_RMin = R_Min
!write(u6,*) 'Eval_RMin=',Eval_RMin
!                                                                      *
!***********************************************************************
!                                                                      *

return

end function Eval_RMin
