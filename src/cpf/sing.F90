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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine SING(IWHY)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IWHY

goto(1,2,3),IWHY
1 write(u6,11)
call XFLUSH(u6)
goto 10
2 write(u6,12)
call XFLUSH(u6)
goto 10
3 write(u6,13)
call XFLUSH(u6)

10 return

11 format(' MATRIX WITH ZERO ROW IN DECOMPOSE.')
12 format(' SINGULAR MATRIX IN DECOMPOSE.ZERO DIVIDE IN SOLVE.')
13 format(' NO CONVERGENCE IN IMPROVE.MATRIX IS NEARLY SINGULAR.')

end subroutine SING
