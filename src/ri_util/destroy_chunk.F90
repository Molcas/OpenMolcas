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

subroutine Destroy_Chunk()

use RI_glob, only: Chunk
#ifdef _MOLCAS_MPP_
use RI_glob, only: iMap, ip_Chunk
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_deallocate

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Destroy(ip_Chunk)
  call mma_deallocate(iMap)
else
  call mma_deallocate(Chunk)
end if
#else
call mma_deallocate(Chunk)
#endif

return

end subroutine Destroy_Chunk
