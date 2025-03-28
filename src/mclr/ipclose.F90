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

subroutine ipclose(ia)
! Object: release all vectors above and including the vector indexed ia.

use ipPage
use stdalloc, only: mma_deallocate
use Definitions, only: u6

real*8 rdum(1)

if (ia > Max_CI_Vectors) then
  write(u6,*) 'ipclose: ia > Max_CI_Vectors'
  write(u6,*) 'ia,Max_CI_Vectors=',ia,Max_CI_Vectors
  call Abend()
end if

! Update iDisk_Addr_End

iDisk_Addr_End = 0
if (ia < 0) then

  n_CI_Vectors = 0

else

  n_CI_Vectors = ia-1
  if (DiskBased) then
    do ii=1,ia-1
      if (Status(ii) /= Null_Vector) call dDafile(Lu_ip,dWrite,rdum(1),n(ii),iDisk_Addr_End)
    end do
  end if

end if

! Release memory and flag as a null vector

do ii=max(ia,0),Max_CI_Vectors
  if (Status(ii) == In_Memory) then
    call mma_deallocate(W(ii)%Vec)
    ida(ii) = -1
    n(ii) = 0
    Status(ii) = Null_Vector
  end if
end do

if (diskbased .and. (ia < 0)) then
  call DACLOS(Lu_ip)
  DiskBased = .false.
end if

end subroutine ipclose
