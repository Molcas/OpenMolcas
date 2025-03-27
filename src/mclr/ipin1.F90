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

integer function ipin1(ii,nn)
! Object: return pointer to vector ii with a length of nn and
!         make the vector available in memory as W(ii)%Vec

use ipPage
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: u6

implicit integer(a-h,o-z)
real*8, allocatable :: Tmp(:)

if (ii > Max_CI_Vectors) then
  write(u6,*) 'ipin1: ii > Max_CI_Vectors'
  write(u6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
  call Abend()
end if

if (Status(ii) == In_Memory) then

  ! ii is in memory

  ! If the size of the vector is larger than was originally set
  ! resize the reservation and copy the content

  if (nn > n(ii)) then
    call mma_allocate(Tmp,nn,Label='Tmp')
    Tmp(:) = Zero
    Tmp(1:n(ii)) = W(ii)%Vec(:)
    call mma_deallocate(W(ii)%Vec)
    call mma_allocate(W(ii)%Vec,nn,Label='ipin1')
    W(ii)%Vec(:) = Tmp(:)
    call mma_deallocate(Tmp)
    n(ii) = nn
  end if

  ip1 = ii

else if (Status(ii) == On_Disk) then

  ! ii is on disk

  call mma_allocate(W(ii)%Vec,max(n(ii),nn),Label='ipin1')
  W(ii)%Vec(:) = Zero

  nnn = min(n(ii),nn)

  ! pick up from disk

  idisk = ida(ii)
  call dDafile(Lu_ip,read,W(ii)%Vec,nnn,idisk)
  Status(ii) = In_Memory

  ip1 = ii

else if (Status(ii) == Null_Vector) then

  ip1 = -1

else

  ip1 = -1
  write(u6,*)
  write(u6,*) 'ipIn1: illegal Status(ii)'
  write(u6,*) 'ii=',ii
  write(u6,*)
  call Abend()

end if

ipin1 = ip1

return

end function ipin1
