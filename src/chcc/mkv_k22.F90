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

subroutine MkV_K22(W1,W2,dim_)
! this routine does
! W1(p) = -W2(p)
!
! N.B. toto je iba plytke copy s opacnym znamienkom,
!      da sa aj s blasmi spravit
! N.B. Kvajt odflaknute

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dim_
real(kind=wp), intent(out) :: W1(dim_)
real(kind=wp), intent(in) :: W2(dim_)

W1(:) = -W2(:)

return

end subroutine MkV_K22
