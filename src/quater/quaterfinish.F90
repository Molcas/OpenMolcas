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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************
!  quaterFinish
!
!> @brief
!>   Clean the quater environment
!> @author Y. Carissan
!>
!> @details
!> Release the memory used by the quater program.
!***********************************************************************

subroutine quaterFinish()

use Quater_globals, only: ngeoms, list
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ig

do ig=1,ngeoms+2
  call mma_deallocate(list(ig)%geo)
  call mma_deallocate(list(ig)%geolbl)
end do

return

end subroutine quaterfinish
