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

implicit real*8(A-H,O-Z)
dimension C(*), D(*)

if (N == 1) return
N1 = N-1
do I=1,N1
  I1 = I+1
  do J=I1,N
    if (D(I) <= D(J)) GO TO 20
    DT = D(I)
    D(I) = D(J)
    D(J) = DT
    IN = (I-1)*N
    IOUT = (J-1)*N
    do K=1,N
      CT = C(IN+K)
      C(IN+K) = C(IOUT+K)
      C(IOUT+K) = CT
    end do
20  continue
  end do
end do

return

end subroutine ORDER_CPF
