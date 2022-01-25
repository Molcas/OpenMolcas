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

subroutine CPLAB(A,B,L,M,N,IA,IB,C,IC,IER)

implicit real*8(a-h,o-z)
!ulf integer L,M,N,IA,IB,IC,IER
!ulf real*8 A(IA,M),B(IB,N),C(IC,N)
!ulf real*8 TEMP
dimension A(IA,M), B(IB,N), C(IC,N)

if ((IA >= L) .and. (IB >= M) .and. (IC >= L)) GO TO 5
IER = 129
GO TO 9000
5 IER = 0
do I=1,L
  do J=1,N
    TEMP = 0.0d0
    do K=1,M
      TEMP = A(I,K)*B(K,J)+TEMP
    end do
    C(I,J) = C(I,J)+TEMP
  end do
end do
GO TO 9005
9000 continue
9005 return

end subroutine CPLAB
