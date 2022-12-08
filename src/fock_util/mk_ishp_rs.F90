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
! Copyright (C) Francesco Aquilante                                    *
!               2021, Roland Lindh                                     *
!***********************************************************************

subroutine Mk_iShp_rs(iShp_rs,nShell)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nShell
integer(kind=iwp), intent(out) :: iShp_rs(nShell*(nShell+1)/2)
integer(kind=iwp) :: iaSh, ibSh, iShp
integer(kind=iwp), external :: Cho_F2SP

! *** Mapping shell pairs from the full to the reduced set

do iaSh=1,nShell
  do ibSh=1,iaSh
    iShp = iaSh*(iaSh-1)/2+ibSh
    iShp_rs(iShp) = Cho_F2SP(iShp)
  end do
end do

end subroutine Mk_iShp_rs
