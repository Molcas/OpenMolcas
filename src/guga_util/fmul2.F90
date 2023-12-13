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

subroutine FMUL2(A,B,C,NROW,NCOL,N)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NROW, NCOL, N
real(kind=wp), intent(in) :: A(NROW,N), B(NCOL,N)
real(kind=wp), intent(out) :: C(NROW,NCOL)
integer(kind=iwp) :: J, J1, K
real(kind=wp) :: CJ(1000), FAC

if (NROW > size(CJ)) then
  write(u6,*)
  write(u6,*) ' *** Error in Subroutine FMUL2 ***'
  write(u6,*) ' row dimension exceeds local buffer size'
  write(u6,*)
  call Abend()
end if

do J=1,NCOL
  CJ(1:NROW) = Zero
  if (J /= NCOL) then
    J1 = J+1
    do K=1,N
      FAC = B(J,K)
      if (FAC /= Zero) CJ(J1:NROW) = CJ(J1:NROW)+FAC*A(J1:NROW,K)
    end do
  end if
  C(:,J) = CJ(1:NROW)
end do

return

end subroutine FMUL2
