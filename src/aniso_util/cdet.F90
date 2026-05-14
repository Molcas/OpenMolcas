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

function CDET(MA,N,A)
!=================================================
! MA is the MAximal dimension
! N is the dimension N<MA
! A is the Complex matrix
!=================================================

use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, iwp

implicit none
complex(kind=wp) :: CDET
integer(kind=iwp), intent(in) :: MA, N
complex(kind=wp), intent(inout) :: A(MA,MA)
integer(kind=iwp) :: I, J, K, K1, L
real(kind=wp) :: P, Q
complex(kind=wp) :: CP, CQ

I = 0
CDET = cZero

do K=1,N
  P = Zero

  do J=K,N
    Q = abs(A(J,K))
    if (Q > P) then
      P = Q
      I = J
    end if
  end do

  CP = cOne/A(I,K)

  if (I /= K) then
    CDET = -CDET
    do L=K,N
      CQ = A(I,L)
      A(I,L) = A(K,L)
      A(K,L) = CQ
    end do
  end if

  CDET = CDET*A(K,K)

  if (K < N) then
    K1 = K+1
    do I=K1,N
      CQ = A(I,K)*CP
      A(I,K1:N) = A(I,K1:N)-CQ*A(K,K1:N)
    end do
  end if

end do ! K

return

end function CDET
