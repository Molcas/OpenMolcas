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

subroutine ClrRunCacheDS()

use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "pg_ds_info.fh"
integer(kind=iwp) :: i

do i=1,num_DS_init
  i_DS_inmem(i) = Zero
  DS_init(i) = 0
  iLbl_DS_inmem(i) = ' '
end do
num_DS_init = 0

end subroutine ClrRunCacheDS
