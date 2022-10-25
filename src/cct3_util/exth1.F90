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

subroutine exth1(a,b,dimp,dimq,p,nfact)
! this routine extracts A(p,q) -> B_p(q)
!
! a     - matrix a (Input)
! b     - matrix b (Output)
! dimp  - dimension of p (Input)
! dimq  - dimension of q (Input)
! p     - value of index p (Input)
! nfact - sign (+-1,0) (Input)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, p, nfact
real(kind=wp), intent(in) :: a(dimp,dimq)
real(kind=wp), intent(out) :: b(dimq)

if (nfact == 1) then
  b(:) = a(p,:)
else if (nfact == -1) then
  b(:) = -a(p,:)
else if (nfact == 0) then
  b(:) = Zero
end if

return

end subroutine exth1
