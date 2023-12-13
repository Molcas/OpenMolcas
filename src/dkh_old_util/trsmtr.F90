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

subroutine TrSmtr(A,B,C,FACTOR,N,H,W)
! B*A*BT

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(N*(N+1)/2), B(N,N), FACTOR
real(kind=wp), intent(inout) :: C(N*(N+1)/2)
real(kind=wp), intent(out) :: H(N,N), W(N,N)

call Square(A,W,n,1,n)
call DGEMM_('N','N',n,n,n,One,B,n,W,n,Zero,H,n)
call dGemm_tri('N','T',n,n,n,One,H,n,B,n,Factor,C,n)

return

end subroutine TrSmtr
