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

integer function ipout(ii)
! ipout will page out vector ii to disk and free the memory area

use ipPage
use stdalloc, only: mma_deallocate

implicit integer(a-h,o-z)

ipout = 0
if (.not. diskbased) return

if ((Status(ii) == In_Memory) .and. (ii > 0)) then
  idisk = ida(ii)
  nn = n(ii)
  call dDafile(Lu_ip,write,W(ii)%Vec,nn,idisk)
  Status(ii) = On_Disk
  call mma_deallocate(W(ii)%Vec)
else
  ipout = -1
end if

return

end function ipout
