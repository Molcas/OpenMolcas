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

use Slapaf_Parameters, only: iRow

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
character*120 Temp
character*16 filnam

Lu_UDIC = 91
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
!open(Lu_UDIC,File=filnam,Form='FORMATTED',Status='OLD')
rewind(Lu_UDIC)

! Find begining of definitions of internal coordinates

do iLines=1,iRow
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:4) == 'VARY') Go To 100
end do
call WarningMessage(2,' No internal coordinates are defined!')
call Quit_OnUserError()

100 continue
iInt = 0
nFix = 0
nRowH = 0 ! Number of Rows of Hessian Numerically estimated
do jLines=iLines+1,iRow
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:3) == 'FIX') Go To 200
  if (Temp(1:4) == 'ROWH') then
    kLines = jLines
    Go To 300
  end if
  ! Do not count line if continuation character
  if (index(Temp,'&') == 0) iInt = iInt+1
end do
Go To 400

200 continue
do kLines=jLines+1,iRow
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:4) == 'ROWH') Go To 300
  ! Do not count line if continuation character
  if (index(Temp,'&') == 0) nFix = nFix+1
end do
300 do lLines=kLines+1,iRow
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  ! Do not count line if continuation character
  if (index(Temp,'&') == 0) nRowH = nRowH+1
end do
400 continue

close(Lu_UDIC)

return

end subroutine Rd_UDIC
