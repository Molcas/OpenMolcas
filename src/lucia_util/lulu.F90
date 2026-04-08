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

!#define _DEBUGPRINT_
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
!   L(I,J) = L(nTri_Elem(I-1) + J) ( I >= J )
!
!   U(I,J) = U(nTri_Elem(J-1) + I) ( J >= I )
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
!         L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R)) /U(R,R)
!       END OF LOOP OVER I
!     END OF LOOP OVER R
!
! JEPPE OLSEN, OCTOBER 1988

use Index_Functions, only: nTri_Elem
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NDIM
real(kind=wp), intent(in) :: A(NDIM,NDIM)
real(kind=wp), intent(inout) :: L(nTri_Elem(NDIM)), U(nTri_Elem(NDIM))
integer(kind=iwp) :: I, J, R
real(kind=wp) :: XFACI
real(kind=wp), external :: dDot_

do R=1,NDIM

  do J=R,NDIM
    U(nTri_Elem(J-1)+R) = A(R,J)-dDot_(R-1,L(nTri_Elem(R-1)+1),1,U(nTri_Elem(J-1)+1),1)
  end do

  XFACI = One/U(nTri_Elem(R))
  L(nTri_Elem(R)) = One
  do I=R+1,NDIM
    L(nTri_Elem(I-1)+R) = (A(I,R)-dDot_(R-1,L(nTri_Elem(I-1)+1),1,U(nTri_Elem(R-1)+1),1))*XFACI
  end do

end do

#ifdef _DEBUGPRINT_
write(u6,*) ' L MATRIX'
call PRSYM(L,NDIM)
write(u6,*) ' U MATRIX (TRANSPOSED)'
call PRSYM(U,NDIM)
#endif

end subroutine LULU
