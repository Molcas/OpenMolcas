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
! Copyright (C) Giovanni Ghigo                                         *
!***********************************************************************

subroutine Rd_UDIC_RowH(nInter,nRowH,mRowH)
!***********************************************************************
!                                                                      *
! Object: Reading the Internal Coordinates required for Numerical      *
!         estimation of single rows and columns of Hessian             *
! Called from: PrePro when nRowH > 0                                   *
! Author: Giovanni Ghigo, University of Torino, Italy                  *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
character*8 Labels(nInter)
character*8 cLbl
character*120 Temp
character*16 filnam
integer mRowH(nRowH)

Lu = 6
Lu_UDIC = 91
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
rewind(Lu_UDIC)
mRowH(:) = 0

! Find begining of definitions of internal coordinates

10 read(Lu_UDIC,'(A)') Temp
call UpCase(Temp)
if (Temp(1:4) /= 'VARY') Go To 10
do iLines=1,nInter
20 read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  if (Temp(1:3) == 'FIX') Go To 20
  cLbl = '        '
  do j=1,120
    if (Temp(j:j) == ' ') goto 30
    cLbl(j:j) = Temp(j:j)
  end do
30 Labels(iLines) = cLbl
35 if (index(Temp,'&') == 0) goto 37
  read(Lu_UDIC,'(A)') Temp
  goto 35
37 continue
end do

read(Lu_UDIC,'(A)') Temp ! Skip ROWH
do iRowH=1,nRowH
  read(Lu_UDIC,'(A)') Temp
  call UpCase(Temp)
  cLbl(1:8) = Temp(1:8)
  do kLines=1,nInter
    if (cLbl == Labels(kLines)) then
      mRowH(iRowH) = kLines
      goto 40
    end if
  end do
  call WarningMessage(2,'Error in rd_udic')
  write(Lu,*) '**********************************************'
  write(Lu,*) ' ERROR: Undefined internal ROWH coordinate in '
  write(Lu,*) ' ',Temp(1:60)
  write(Lu,*) '**********************************************'
  call Quit_OnUserError()
40 continue
end do
close(Lu_UDIC)

return

end subroutine Rd_UDIC_RowH
