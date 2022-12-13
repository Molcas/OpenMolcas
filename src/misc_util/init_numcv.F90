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

subroutine INIT_NumCV(NumCV,nSymm)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSymm
integer(kind=iwp), intent(out) :: NumCV(nSymm)
#include "cholesky.fh"

NumCV(:) = NumCho(1:nSymm)

return

end subroutine INIT_NumCV
