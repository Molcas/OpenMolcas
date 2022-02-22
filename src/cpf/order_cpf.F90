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

subroutine ORDER_CPF(C,D,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: N
real(kind=wp) :: C(*), D(*)
integer(kind=iwp) :: I, I1, IIN, IOUT, J, K, N1
real(kind=wp) :: CT, DT

if (N == 1) return
N1 = N-1
do I=1,N1
  I1 = I+1
  do J=I1,N
    if (D(I) <= D(J)) cycle
    DT = D(I)
    D(I) = D(J)
    D(J) = DT
    IIN = (I-1)*N
    IOUT = (J-1)*N
    do K=1,N
      CT = C(IIN+K)
      C(IIN+K) = C(IOUT+K)
      C(IOUT+K) = CT
    end do
  end do
end do

return

end subroutine ORDER_CPF
