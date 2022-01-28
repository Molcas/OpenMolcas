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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine TrSmr(A,B,C,N,H,W)
! TRANSFORM SYMMETRIC MATRIX A BY UNITARY TRANSFORMATION
! IN B. RESULT IS IN C

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(N*(N+1)/2), B(N,N)
real(kind=wp), intent(out) :: C(N*(N+1)/2), H(N,N), W(N,N)

call Square(A,W,n,1,n)
call DGEMM_('T','N',n,n,n,One,B,n,W,n,Zero,H,n)
call dGemm_tri('N','N',n,n,n,One,H,n,B,n,Zero,C,n)

return

end subroutine TrSmr
