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

subroutine VAM(A,LA,B,LB,C,LC,D,LD,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: LA, LB, LC, LD, N
real(kind=wp) :: A(*), B(*), C(*), D(*)
integer(kind=iwp) :: I

do I=0,N-1
  D(1+I*LD) = (A(1+I*LA)+B(1+I*LB))*C(1+I*LC)
end do

return

end subroutine VAM
