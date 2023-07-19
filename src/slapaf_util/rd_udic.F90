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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************

subroutine Rd_UDIC(iInt,nFix,nRowH)
!***********************************************************************
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Slapaf_Info, only: iRow
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iInt, nFix, nRowH
integer(kind=iwp) :: iLines, jLines, kLines, lLines, Lu_UDIC, Skip
logical(kind=iwp) :: Found
character(len=120) :: Temp
character(len=16) :: filnam
integer(kind=iwp), external :: IsFreeUnit

Lu_UDIC = IsFreeUnit(91)
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
!open(Lu_UDIC,File=filnam,Form='FORMATTED',Status='OLD')
rewind(Lu_UDIC)

! Find begining of definitions of internal coordinates

Found = .false.
do iLines=1,iRow
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:4) == 'VARY') then
    Found = .true.
    exit
  end if
end do
if (.not. Found) then
  call WarningMessage(2,' No internal coordinates are defined!')
  call Quit_OnUserError()
end if

Skip = 2
iInt = 0
nFix = 0
nRowH = 0 ! Number of Rows of Hessian Numerically estimated
do jLines=iLines+1,iRow
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:3) == 'FIX') then
    Skip = 0
    exit
  end if
  if (Temp(1:4) == 'ROWH') then
    kLines = jLines
    Skip = 1
    exit
  end if
  ! Do not count line if continuation character
  if (index(Temp,'&') == 0) iInt = iInt+1
end do

if (Skip < 2) then
  if (Skip < 1) then
    do kLines=jLines+1,iRow
      read(Lu_UDIC,'(A)') Temp
      call UpCase(Temp)
      if (Temp(1:4) == 'ROWH') exit
      ! Do not count line if continuation character
      if (index(Temp,'&') == 0) nFix = nFix+1
    end do
  end if
  do lLines=kLines+1,iRow
    read(Lu_UDIC,'(A)') Temp
    call UpCase(Temp)
    ! Do not count line if continuation character
    if (index(Temp,'&') == 0) nRowH = nRowH+1
  end do
end if

close(Lu_UDIC)

return

end subroutine Rd_UDIC
