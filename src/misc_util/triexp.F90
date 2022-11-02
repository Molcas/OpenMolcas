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
! Copyright (C) 1989, Markus P. Fuelscher                              *
!***********************************************************************

subroutine TRIEXP(A,B,NDIM)
! AUTHOR:        M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, 1989
! MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
! EXPAND A REAL, SYMMETRIC MATRIX GIVEN AS THE LOWER TRIANGLE
! INTO FULL STORAGE MODE.
!
! CALLING PARAMETERS:
! A   : INPUT MATRIX (LOWER TRINAGULAR STORAGE MODE)
! B   : OUTPUT MATRIX (FULL STORAGE MODE)
! NDIM: DIMENSION OF MATRIX A AND B

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: A(*), B(*)
integer(kind=iwp) :: NDIM
integer(kind=iwp) :: I, IMODE, J, K, L1, L2, L3, L4
integer(kind=iwp), external :: ip_of_Work

IMODE = 1
! TODO: Actually, this is not allowed by the Fortran standard
if (ip_of_Work(A(1)) == ip_of_Work(B(1))) IMODE = 2

if (IMODE == 1) then
  K = 0
  do I=1,NDIM
    do J=1,I
      K = K+1
      L1 = J+(I-1)*NDIM
      B(L1) = A(K)
      L2 = I+(J-1)*NDIM
      B(L2) = A(K)
    end do
  end do
end if

if (IMODE == 2) then
  K = 1+nTri_Elem(NDIM)
  do I=NDIM,1,-1
    do J=I,1,-1
      K = K-1
      L1 = J+(I-1)*NDIM
      L2 = I+(J-1)*NDIM
      L3 = max(L2,L1)
      B(L3) = A(K)
    end do
  end do
  do I=1,NDIM
    do J=1,I
      L1 = J+(I-1)*NDIM
      L2 = I+(J-1)*NDIM
      L3 = max(L2,L1)
      L4 = min(L2,L1)
      B(L4) = B(L3)
    end do
  end do
end if

return

end subroutine TRIEXP
