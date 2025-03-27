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

integer function ipget(nn)
! Get the index of a vector with the length nn.
! Memory or disk space is allocated.

use ipPage
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: u6

implicit integer(a-h,o-z)
character*4 Label

! Take the next memory slot.

n_CI_Vectors = n_CI_Vectors+1
ipget = n_CI_Vectors

if (n_CI_Vectors > Max_CI_Vectors) then
  write(u6,*) 'Number of CI vectors higher than Max_CI_Vectors'
  write(u6,*) 'Max_CI_Vectors=',Max_CI_Vectors
  call Abend()
end if

ida(ipget) = iDisk_Addr_End
n(ipget) = nn

! Allocate memory for vector if of non-zero  length

if (nn > 0) then
  write(Label,'(I3.3)') n_CI_Vectors
  call mma_allocate(W(ipget)%Vec,nn,Label='ipget'//Label)
  Status(ipget) = In_Memory
  W(ipget)%Vec(:) = Zero
else
  !Status(ipget) = Null_Vector

  ! The calling code doesn't have the logic to handle the
  ! case that W(i)%Vec is not allocated. Hence, we have
  ! to make a dummy allocation to make sure that the compiler
  ! doesn't puke.
  n(ipget) = 1
  write(Label,'(I3.3)') n_CI_Vectors
  call mma_allocate(W(ipget)%Vec,1,Label='ipget'//Label)
  Status(ipget) = In_Memory
  W(ipget)%Vec(:) = Zero
end if

! If diskbased mode put vector on disc and release memory

if (DiskBased) then
  if (Status(ipget) /= Null_Vector) then
    call dDafile(Lu_ip,write,W(ipget)%Vec,nn,iDisk_Addr_End)
    Status(ipget) = On_Disk
    call mma_deallocate(W(ipget)%Vec)
  end if
end if

return

end function ipget
