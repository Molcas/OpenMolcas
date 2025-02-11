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
! Copyright (C) 1988, Jeppe Olsen                                      *
!***********************************************************************

subroutine LULU(A,L,U,NDIM)
! LU DECOMPOSITION OF MATRIX A
!
!     A = L * U
!
! WHERE L IS A LOWER TRIANGULAR MATRIX WITH A
! UNIT DIAGONAL AND U IS AN UPPER DIAGONAL
!
! L AND U ARE STORED AS ONE DIMENSIONAL ARRAYS
!
!   L(I,J) = L(I*(I-1)/2 + J ) ( I >= J )
!
!   U(I,J) = U(J*(J-1)/2 + I ) ( J >= I )
!
! THIS ADDRESSING SCHEMES SUPPORTS VECTORIZATION OVER COLUMNS
! FOR L AND  OVER ROWS FOR U .
!
! NO PIVOTING IS DONE HERE, SO THE SCHEME GOES :
!
!     LOOP OVER R=1, NDIM
!       LOOP OVER J = R, NDIM
!         U(R,J) = A(R,J) - SUM(K=1,R-1) L(R,K) * U(K,J)
!       END OF LOOP OVER J
!
!       LOOP OVER I = R+1, NDIM
!         L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R) ) /U(R,R)
!       END OF LOOP OVER I
!     END OF LOOP OVER R
!
! JEPPE OLSEN, OCTOBER 1988

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NDIM
real(kind=wp) :: A(NDIM,NDIM), L(*), U(*)
integer(kind=iwp) :: I, J, NTEST, R
real(kind=wp) :: XFACI
real(kind=wp), external :: dDot_

do R=1,NDIM

  do J=R,NDIM
    U(J*(J-1)/2+R) = A(R,J)-dDot_(R-1,L(R*(R-1)/2+1),1,U(J*(J-1)/2+1),1)
  end do

  XFACI = One/U(R*(R+1)/2)
  L(R*(R+1)/2) = One
  do I=R+1,NDIM
    L(I*(I-1)/2+R) = (A(I,R)-dDot_(R-1,L(I*(I-1)/2+1),1,U(R*(R-1)/2+1),1))*XFACI
  end do

end do

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' L MATRIX'
  call PRSYM(L,NDIM)
  write(u6,*) ' U MATRIX ( TRANSPOSED )'
  call PRSYM(U,NDIM)
end if

end subroutine LULU
