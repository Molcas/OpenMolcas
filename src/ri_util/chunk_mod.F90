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

module Chunk_Mod

use Definitions, only: wp
#ifdef _MOLCAS_MPP_
use Definitions, only: iwp
#endif

implicit none
private

real(kind=wp), allocatable :: Chunk(:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: ip_Chunk = 0
integer(kind=iwp), allocatable :: iMap(:)
#endif

public :: Chunk
#ifdef _MOLCAS_MPP_
public :: iMap, ip_Chunk
#endif

end module Chunk_Mod
