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

subroutine VSADD(A,LA,S,C,LC,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: LA, LC, N
real(kind=wp) :: A(*), S, C(*)
integer(kind=iwp) :: I

do I=0,N-1
  C(1+LC*I) = A(1+LA*I)+S
end do

return

end subroutine VSADD
