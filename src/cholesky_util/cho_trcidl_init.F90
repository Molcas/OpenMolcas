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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_TrcIdl_Init()
!
! Thomas Bondo Pedersen, May 2010.
!
! Allocate and init array for tracing idle processors

use Para_Info, only: nProcs
use Cholesky, only: Cho_Real_Par, Idle
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: l_Idle

if (Cho_Real_Par) then
  l_Idle = nProcs
else
  l_Idle = 1
end if
call mma_allocate(Idle,l_Idle,Label='Idle')

Idle(:) = 0

end subroutine Cho_TrcIdl_Init
