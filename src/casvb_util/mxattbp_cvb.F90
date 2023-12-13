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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mxattbp_cvb(a,b,n1,n2,n3,c)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n1, n2, n3
real(kind=wp), intent(in) :: a(n2,n1), b(n2,n3)
real(kind=wp), intent(inout) :: c(n1,n3)

call DGEMM_('T','N',n1,n3,n2,One,a,n2,b,n2,One,c,n1)

return

end subroutine mxattbp_cvb
