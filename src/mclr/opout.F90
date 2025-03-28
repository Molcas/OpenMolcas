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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine opout(ii)
! opout will release the memory area of vector ii without updating the disk

use ipPage
use stdalloc, only: mma_deallocate
use Definitions, only: u6

implicit integer(a-h,o-z)

if (ii > Max_CI_Vectors) then
  write(u6,*) 'opout: ii > Max_CI_Vectors'
  write(u6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
  call Abend()
end if

if (.not. diskbased) return

if ((Status(ii) == In_Memory) .and. (ii > 0)) then
  Status(ii) = On_Disk
  call mma_deallocate(W(ii)%Vec)
end if

return

end subroutine opout
