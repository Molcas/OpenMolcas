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

subroutine ipnout(iii)
! Object: write all vectors in memory on disk but vector iii

use ipPage
use stdalloc, only: mma_deallocate
use Definitions, only: u6

implicit integer(a-h,o-z)

if (iii > Max_CI_Vectors) then
  write(u6,*) 'ipout: iii > Max_CI_Vectors'
  write(u6,*) 'iii,Max_CI_Vectors=',iii,Max_CI_Vectors
  call Abend()
end if

if (.not. DiskBased) return

do ii=1,Max_CI_Vectors

  if ((Status(ii) == In_Memory) .and. (ii /= iii)) then
    idisk = ida(ii)
    nn = n(ii)
    call dDafile(Lu_ip,write,W(ii)%Vec,nn,idisk)
    Status(ii) = On_Disk
    call mma_deallocate(W(ii)%Vec)
  end if

end do

end subroutine ipnout
