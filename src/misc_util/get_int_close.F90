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

subroutine Get_Int_Close()

use GetInt_mod, only: LuCVec, Vec2
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i

! Close files.
do i=1,2
  if (LuCVec(i) /= -1) then
    call DACLOS(LuCVec(i))
    LuCVec(i) = -1
  end if
end do

if (allocated(Vec2)) call mma_deallocate(Vec2)

end subroutine Get_Int_Close
