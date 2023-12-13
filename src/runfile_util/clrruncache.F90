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

subroutine ClrRunCache()

use RunFile_data, only: DS_cache, IS_cache, num_DS_init, num_IS_init
use Constants, only: Zero

implicit none

DS_cache(1:num_DS_init)%val = Zero
DS_cache(1:num_DS_init)%lab = ''
num_DS_init = 0

IS_cache(1:num_DS_init)%val = 0
IS_cache(1:num_DS_init)%lab = ''
num_IS_init = 0

return

end subroutine ClrRunCache
