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

subroutine VNEG_CPF(A,IA,B,IB,IAB)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IA, IB, IAB
real(kind=wp) :: A(*), B(*)
integer(kind=iwp) :: I

do I=0,IAB-1
  B(1+I*IB) = -A(1+I*IA)
end do

return

end subroutine VNEG_CPF
