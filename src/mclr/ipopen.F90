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

logical function ipopen(nconf,page)
! Initiate the whole lot.

use ipPage
use stdalloc, only: mma_maxDBLE

implicit real*8(a-h,o-z)
logical page

! Ask how much memory is available

call mma_maxDBLE(nMax)
nmax = nmax/2

if (Page) then

  ! Initiate for disk based storage.

  if (.not. DiskBased) then
    Lu_ip = 21
    Lu_ip = IsFreeUnit(Lu_ip)
    call Daname(Lu_ip,'TEMPCIV')
    DiskBased = .true.
  end if

  ! n  : Length of CI-vector
  ! ida: disk address

  n(0:Max_CI_Vectors) = 0
  ida(0:Max_CI_Vectors) = -1
  Status(0:Max_CI_Vectors) = Null_Vector

  ! iDisk_Addr_End: next free disk address
  ! n_CI_Vectors : number of CI-vectors

  iDisk_Addr_End = 0
  n_CI_Vectors = 0

else

  if (DiskBased) then
    call ipTerm()
    DiskBased = .false.
  end if

end if

ipopen = DiskBased

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nconf)

end function ipopen
