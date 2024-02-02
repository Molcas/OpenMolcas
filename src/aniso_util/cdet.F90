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

complex*16 function CDET(MA,N,A)
!=================================================
! MA is the MAximal dimension
! N is the dimension N<MA
! A is the Complex matrix
!=================================================

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: N, MA
complex(kind=8), intent(inout) :: A(MA,MA)
integer :: I, J, K, L, K1
real(kind=8) :: P, Q
complex(kind=8) :: CP, CQ

I = 0
CDET = (0.0_wp,0.0_wp)

do K=1,N
  P = 0.0_wp

  do J=K,N
    Q = abs(A(J,K))
    if (Q > P) then
      P = Q
      I = J
    end if
  end do

  CP = (1.0_wp,0.0_wp)/A(I,K)

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
      do L=K1,N
        A(I,L) = A(I,L)-CQ*A(K,L)
      end do
    end do
  end if

end do ! K

return

end function CDET
