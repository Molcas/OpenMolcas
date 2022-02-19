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

subroutine SQUARM_CPF(A,B,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: N
real(kind=wp) :: A(*), B(N,N)
integer(kind=iwp) :: I, IIN

IIN = 1
do I=1,N
  call DCOPY_(I,A(IIN),1,B(I,1),N)
  call VNEG_CPF(A(IIN),1,B(1,I),1,I)
  IIN = IIN+I
end do

return

end subroutine SQUARM_CPF
