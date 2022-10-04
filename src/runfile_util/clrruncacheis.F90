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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************

subroutine ClrRunCacheIS()

implicit none
#include "pg_is_info.fh"
integer i

do i=1,num_IS_init
  i_IS_inmem(i) = 0
  IS_init(i) = 0
  iLbl_IS_inmem(i) = ' '
end do
num_IS_init = 0

end subroutine ClrRunCacheIS
