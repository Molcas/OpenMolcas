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

subroutine Cho_P_ChkInt(xInt,Diag,iSym,nErr,Tol,Report)
!
! Purpose: check diagonals in qualified integral columns.

use Cholesky, only: Cho_Real_Par
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: xInt(*), Diag(*), Tol
integer(kind=iwp), intent(in) :: iSym
integer(kind=iwp), intent(out) :: nErr
logical(kind=iwp) :: Report

if (Cho_Real_Par) then
  call Cho_P_QualSwp()
  call Cho_ChkInt(xInt,Diag,iSym,nErr,Tol,Report)
  call Cho_P_QualSwp()
else
  call Cho_ChkInt(xInt,Diag,iSym,nErr,Tol,Report)
end if

end subroutine Cho_P_ChkInt
