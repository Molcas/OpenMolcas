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
!***********************************************************************

subroutine INIT_NumCV(NumCV,nSym)

use Cholesky, only: NumCho
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym
integer(kind=iwp), intent(out) :: NumCV(nSym)

NumCV(:) = NumCho(1:nSym)

return

end subroutine INIT_NumCV
